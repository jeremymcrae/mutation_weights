
import argparse
import gzip
from itertools import groupby, count

import pysam
import pandas

rates_path = '/nfs/users/nfs_j/jm33/apps/mutation_weights/dominant_rates.txt.gz'
cadd_path = '/lustre/scratch115/projects/ddd/users/jm33/cadd/v1.0/whole_genome_SNVs.tsv.gz'
outpath = '/nfs/users/nfs_j/jm33/apps/mutation_weights/cadd_scores-1.0.txt.gz'

def get_options():
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--rates', default=rates_path,
        help='path to site specific rates per gene')
    parser.add_argument('--cadd', default=cadd_path,
        help='path to tabix indexed CADD scores for all possible SNVs (GRCh37)')
    parser.add_argument('--output', default=outpath)
    
    return parser.parse_args()

def load_cadd(cadd, chrom, start, end):
    ''' load cadd scores for a region
    '''
    def parse(line):
        chrom, pos, ref, alt, raw, scaled = line.split('\t')
        
        return {'chrom': chrom, 'pos': int(pos), 'ref': ref, 'alt': alt,
            'raw': float(raw), 'score': float(scaled)}
    
    return pandas.DataFrame([ parse(x) for x in cadd.fetch(chrom, start, end) ])

def load_rates(path):
    ''' get a DataFrame of mutation rates by site
    '''
    
    rates = pandas.read_table(path)
    rates['chrom'] = rates['chrom'].astype(str)
    
    return rates

def as_range(g):
    l = list(g)
    return l[0], l[-1]

def main():
    args = get_options()
    
    rates = load_rates(args.rates)
    cadd = pysam.TabixFile(args.cadd)
    
    cadd_data = pandas.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'raw', 'score'])
    for symbol, group in rates.groupby('symbol'):
        chrom = set(group['chrom'])
        assert len(chrom) == 1
        chrom = list(chrom)[0]
        
        print(symbol)
        # split the gene into contiguous blocks of positions relating to each,
        # exon, rather than loading cadd scores for all sites from gene start
        # to end which is slower and more memory intensive
        for _, g in groupby(sorted(set(group['pos'])), key=lambda n, c=count(): n-next(c)):
            start, end = as_range(g)
            temp = load_cadd(cadd, chrom, start-1, end)
            cadd_data = cadd_data.append(temp, ignore_index=True)
    
    cadd_data['pos'] = cadd_data['pos'].astype(int)
    
    with gzip.open(args.output, 'wt') as handle:
        cadd_data.to_csv(handle, sep='\t', index=False)

if __name__ == '__main__':
    main()
