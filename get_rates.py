
import gzip
import argparse
import os

import pandas

from denovonear.ensembl_requester import EnsemblRequest
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.load_gene import get_transcript_ids, construct_gene_object
from denovonear.site_specific_rates import SiteRates

KNOWN_PATH = "https://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz"

def get_options():
    '''
    '''
    
    parser = argparse.ArgumentParser('get mutation rates at all sites in ' \
        'dominant genes')
    parser.add_argument('--known', default=KNOWN_PATH,
        help='path or url to table of known genes')
    parser.add_argument('--output', default='dominant_rates.txt.gz',
        help='path write table of rates to')
    
    return parser.parse_args()

def load_dominant(url, include_possible=False):
    ''' get a set of known dominant DD genes
    
    Args:
        url: URL for online DDG2P
        include_possible: boolean for whether to include DDG2P genes where the
            confidence is only possible
    
    Returns:
        set of HGNC symbols for dominant DDG2P genes
    '''
    
    compress = None
    if url[-2:] == 'gz':
        compress = 'gzip'
    
    data = pandas.read_table(url, sep=",", compression=compress)
    
    data = data[~data['allelic requirement'].isnull()]
    data = data[~data['DDD category'].isnull()]
    
    # select the dominant DDG2P genes
    data = data[data['allelic requirement'].str.contains('monoallelic|x-linked dominant', regex=True)]
    
    string = 'confirmed|probable|both DD and IF'
    if include_possible:
        string += '|possible'
    
    data = data[data['DDD category'].str.contains(string, regex=True)]
    
    return set(data['gene symbol'])

def get_transcripts(symbol, ensembl):
    ''' get a list of Transcript objects for a gene
    
    Args:
        symbol: HGNC symbol for a gene
        ensembl: EnsemblRequest object, to retrieve gene data with
    
    Returns:
        list of Transcript objects (see denovonear), sorted by size (longest
        transcripts first)
    '''
    
    transcript_ids = get_transcript_ids(ensembl, symbol)
    
    transcripts = []
    for x in sorted(transcript_ids, key=transcript_ids.get, reverse=True):
        try:
            tx = construct_gene_object(ensembl, x)
            transcripts.append(tx)
        except ValueError:
            continue
    
    return transcripts

def rates_per_site(transcripts, mut_dict):
    ''' get table of mutation rates per site across all transcripts for a gene
    
    Args:
        transcripts: list of Transcript objects for a single gene (sorted by
            size, longest first)
        mut_dict:
    '''
    
    rates = []
    combined = None
    for tx in transcripts:
        sites = SiteRates(tx, mut_dict, masked_sites=combined)
        if combined is None:
            combined = tx
        combined += tx
        
        # for each consequence type, get all the sites for that consequence type,
        # along with the ref, alt and coordinates
        for cq in ['nonsense', 'missense', 'synonymous', 'splice_lof', 'splice_region']:
            for choice in sites[cq]:
                choice['pos'] = tx.get_position_on_chrom(choice['pos'], choice['offset'])
                choice['chrom'] = tx.get_chrom()
                choice['cq'] = cq
                rates.append(choice)
    
    return rates

def get_gene_rates(symbol, ensembl, mut_dict):
    ''' get per nucleotide mutation rates for all SNV alt alleles in a gene
    
    Args:
        symbol: HGNC symbol for gene
        ensembl: EnsemblRequest object, for extracting coordinates and sequence.
        mut_dict: list of lists of sequence context changes and associated
            mutation rates as [[initial, changed, rate], ...]
    
    Returns:
        pandas DataFrame of mutation rates at each possible SNV change within
        the coding sequence of a gene.
    '''
    
    transcripts = get_transcripts(symbol, ensembl)
    
    if len(transcripts) == 0:
        return pandas.DataFrame(columns=['symbol', 'chrom', 'pos', 'ref', 'alt',
            'cq', 'prob'])
    
    rates = rates_per_site(transcripts, mut_dict)
    
    # convert the list of dictionaries to a DataFrame, and track the gene
    rates = pandas.DataFrame(rates)
    rates['symbol'] = symbol
    
    rates.sort_values(by=['pos', 'alt'], inplace=True)
    
    return rates[['symbol', 'chrom', 'pos', 'ref', 'alt', 'cq', 'prob']]

def main():
    
    args = get_options()
    
    ensembl = EnsemblRequest('cache', 'grch37')
    mut_dict = load_mutation_rates()
    
    dominant = load_dominant(args.known)
    
    data = pandas.DataFrame(columns=['symbol', 'chrom', 'pos', 'ref', 'alt',
        'cq', 'prob'])
    data['pos'] = data['pos'].astype(int)
    for symbol in dominant:
        print(symbol)
        rates = get_gene_rates(symbol, ensembl, mut_dict)
        data = data.append(rates, ignore_index=True)
    
    with gzip.open(args.output, 'wt') as handle:
        data.to_csv(handle, sep='\t', index=False)

if __name__ == '__main__':
    main()
