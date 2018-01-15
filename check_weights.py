

import os
from itertools import accumulate

import pandas

from denovonear.ensembl_requester import EnsemblRequest
from denovonear.load_gene import construct_gene_object
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.site_specific_rates import SiteRates

from mupit.constants import LOF_CQ, MISSENSE_CQ

from weights.load_data import load_rates, load_cadd, count_trios, load_de_novos, \
    load_regional_constraint
from weights.chrX_correction import correct_for_x_chrom
from weights.plot_enrichment import plot_enrichment
from weights.constraint import get_constrained_positions

rates_path = '/nfs/users/nfs_j/jm33/apps/mutation_weights/dominant_rates.txt.gz'
cadd_path = '/nfs/users/nfs_j/jm33/apps/mutation_weights/cadd_scores-1.3.txt.gz'
cache_dir = '/nfs/users/nfs_j/jm33/apps/denovonear/scripts/cache'
de_novos_path = '/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.ddd_only.2015-11-24.txt'
validations_path = '/lustre/scratch113/projects/ddd/users/jm33/de_novos.validation_results.2015-11-24.txt'
trios_path = '/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/trios.txt'
families_path = '/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/family_relationships.txt'

constraint_path = '/nfs/users/nfs_j/jm33/apps/mutation_weights/genes_regional_data_cleaned_chisq_0.0_metrics_2016_10_05.txt.gz'

def get_expected_rates(rates, male, female):
    
    autosomal = 2 * (male + female)
    rates['prob'] = rates['prob'] * autosomal
    
    return correct_for_x_chrom(rates, male, female)

def merge_rates_and_cadd(rates, cadd):
    return rates.merge(cadd, how='left', on=['chrom', 'pos', 'ref', 'alt'])

def check_enrichment(missense, de_novos, min_threshold, max_threshold=None):
    ''' check enrichment of missense de novos within sites
    '''
    
    # select missense sites with CADD scores > threshold
    subset = missense[(missense['score'] >= min_threshold)]
    if max_threshold is not None:
        subset = subset[subset['score'] < max_threshold]
    
    expected = sum(subset['prob'])
    keys = set(zip(subset['chrom'], subset['pos'], subset['alt']))
    
    # count the observed missense SNVs within those sites
    de_novo_keys = set(zip(de_novos['chrom'], de_novos['start_pos'], de_novos['alt_allele']))
    observed = len(de_novo_keys & keys)
    
    # return the ratio of observed to expected
    return observed/expected, len(subset)

def annotate_constraint(data, constraint_path, threshold=1e-3, ratio=0.4):
    ''' annotate per-site rates by whether the site is under regional constraint
    '''
    # default to unconstrained
    data['constrained'] = False
    
    constraint = load_regional_constraint(constraint_path)
    mut_dict = load_mutation_rates()
    ensembl = EnsemblRequest(cache_dir, 'grch37')
    
    modified = []
    for symbol, group in data.groupby('symbol'):
        if symbol not in set(constraint['gene']):
            sites = set([])
        else:
            regional = constraint[constraint['gene'] == symbol]
            tx_id = list(regional['transcript'])[0]
            tx = construct_gene_object(ensembl, tx_id.split('.')[0])
            sites = get_constrained_positions(tx, regional, threshold, ratio)
        
        gene_constraint = group['constrained'].copy()
        gene_constraint.loc[group['pos'].isin(sites)] = True
        group['constrained'] = gene_constraint
        
        modified.append(group)
    
    return pandas.concat(modified)

rates = load_rates(rates_path)
cadd = load_cadd(cadd_path)

trios = count_trios(trios_path, families_path)
expected = get_expected_rates(rates, trios['male'], trios['female'])

data = merge_rates_and_cadd(expected, cadd)
data = annotate_constraint(data, constraint_path)

constrained = data[data['constrained']]
unconstrained = data[~data['constrained']]

dominant = set(rates['symbol'])
de_novos = load_de_novos(de_novos_path, validations_path)

# check the enrichment of PTV candidates within the PTV sites
ptv = de_novos[de_novos['hgnc'].isin(dominant) & de_novos['consequence'].isin(LOF_CQ)]
lof_enrich = len(ptv)/sum(rates['prob'][rates.cq.isin(['nonsense', 'splice_lof'])])

subsets = {'all': data, 'constrained': constrained, 'unconstrained': unconstrained}

for key in subsets:
    subset = subsets[key]
    missense = subset[subset['cq'] == 'missense'].copy()
    
    # get enrichment within CADD ranges (rather than in sites above a threshold)
    increment = 5
    thresh = list(range(0, 40, increment))
    enrich = [ check_enrichment(missense, de_novos, x, x + increment) for x in thresh ]
    cdf = [ x[1] for x in enrich ]
    cdf = [ x/sum(cdf) for x in accumulate(cdf) ]
    enrich = [ x[0] for x in enrich ]
    print(key)
    print(thresh)
    print(enrich)
    print(cdf)
    plot_enrichment(thresh, enrich, cdf, path='weights.v1.3.{}.pdf'.format(key))
