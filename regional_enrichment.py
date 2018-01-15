
from scipy.stats import chi2, poisson, fisher_exact
import pandas

from denovonear.ensembl_requester import EnsemblRequest
from denovonear.load_gene import construct_gene_object
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.site_specific_rates import SiteRates

from mupit.open_ddd_data import standardise_ddd_de_novos
from mupit.mutation_rates import get_expected_mutations
from mupit.count_de_novos import get_de_novo_counts

from weights.constraint import get_constrained_positions, classify_de_novos_by_constraint
from weights.load_data import load_de_novos, load_regional_constraint

# compare enrichment of de novo mutations in dominant genes in regions of high
# constraint vs regions without high constraint. Do PTV enrichment and PAV
# enrichment separately.

cache_dir = '/nfs/users/nfs_j/jm33/apps/denovonear/scripts/cache'
constraint_path = '/nfs/users/nfs_j/jm33/apps/mutation_weights/genes_regional_data_cleaned_chisq_0.0_metrics_2016_10_05.txt.gz'
de_novos_path = '/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.ddd_only.2015-11-24.txt'
validations_path = '/lustre/scratch113/projects/ddd/users/jm33/de_novos.validation_results.2015-11-24.txt'

def get_gene_rates(tx, sites, cqs, constrained_sites):
    gene_rates = {'constrained': dict(zip(cqs, [0.0] * len(cqs))),
        'unconstrained': dict(zip(cqs, [0.0] * len(cqs)))}
    for cq in cqs:
        for choice in sites[cq]:
            pos = tx.get_position_on_chrom(choice['pos'], choice['offset'])
            
            category = 'unconstrained'
            if pos in constrained_sites:
                category = 'constrained'
                
            gene_rates[category][cq] += choice['prob']
    
    return gene_rates

def prepare_rates(rates_dict):
    ''' convert a list of rates dicts to a DataFrame, and clean up the data
    '''
    
    rates = pandas.DataFrame(rates_dict)
    
    rename = {'synonymous': 'syn', 'missense': 'mis', 'nonsense': 'non',
        'splice_lof': 'splice_site', 'symbol': 'hgnc'}
    rates = rates.rename(columns=rename)
    return include_indel_rates(rates)

def include_indel_rates(rates):
    ''' add per-gene indel mutation rates to the output file
    '''
    
    # run through the file of rates to find the overall nonsense mutation rate,
    # and the total length of CDS regions in the file.
    nonsense_sum = sum(rates['non'])
    length_sum = sum(rates['length'])
    
    # add the frameshift rates to each gene
    rates['frameshift'] = (rates['length']/length_sum) * nonsense_sum * 1.25
    
    return rates

def get_rates_by_constraint(constraint, cache_dir, threshold=1e-4, ratio=1.0):
    ''' get mutation rates in and out of constrained regions
    '''
    
    rates = {'constrained': [], 'unconstrained': []}
    mut_dict = load_mutation_rates()
    ensembl = EnsemblRequest(cache_dir, 'grch37')
    for tx_id, group in constraint.groupby('transcript'):
        tx = construct_gene_object(ensembl, tx_id.split('.')[0])
        sites = SiteRates(tx, mut_dict)
        
        constrained_sites = get_constrained_positions(tx, group, threshold, ratio)
        
        cqs = ['nonsense', 'missense', 'synonymous', 'splice_lof', 'splice_region']
        gene_rates = get_gene_rates(tx, sites, cqs, constrained_sites)
        
        # now add the gene rates to the larger list of all genes
        for category in ['constrained', 'unconstrained']:
            gene_rates[category]['symbol'] = list(group['gene'])[0]
            gene_rates[category]['chrom'] = list(group['chr'])[0]
            gene_rates[category]['length'] = tx.chrom_pos_to_cds(tx.get_cds_end())['pos']
            
            rates[category].append(gene_rates[category])
    
    return rates

def enrichment(observed, expected):
    ''' assess enrichment of de novo mutations
    '''
    
    groups = {'PTV': ['lof_snv', 'lof_indel'],
        'PAV': ['missense_snv', 'missense_indel']}
    
    data = {}
    for x in groups:
        obs = sum([ observed[x].sum() for x in groups[x] ])
        exp = sum([ expected[x].sum() for x in groups[x] ])
        
        ratio = obs/exp
        p_value = poisson.sf(obs - 1, exp)
        data[x] = {'ratio': ratio, 'p_value': p_value, 'observed': obs,
            'expected': exp}
    
    return data

def compare_regions(constrained, unconstrained, category):
    ''' compare counts between constrained and unconstrained regions
    '''
    counts = [[constrained[category]['observed'], constrained[category]['expected']],
          [unconstrained[category]['observed'], unconstrained[category]['expected']]]
    
    return fisher_exact(counts)

def check_enrichment(constraint, de_novos, cache_dir, male, female, threshold, ratio):
    rates = get_rates_by_constraint(constraint, cache_dir, threshold, ratio)
    
    constrained_exp = get_expected_mutations(prepare_rates(rates['constrained']), male, female)
    unconstrained_exp = get_expected_mutations(prepare_rates(rates['unconstrained']), male, female)
    
    in_constraint = classify_de_novos_by_constraint(constraint, de_novos, cache_dir, threshold, ratio)
    
    constrained_obs = get_de_novo_counts(de_novos[in_constraint])
    unconstrained_obs = get_de_novo_counts(de_novos[[ not x for x in in_constraint ]])
    
    constrained_enrich = enrichment(constrained_obs, constrained_exp)
    unconstrained_enrich = enrichment(unconstrained_obs, unconstrained_exp)
    
    ptv_diff = compare_regions(constrained_enrich, unconstrained_enrich, 'PTV')
    pav_diff = compare_regions(constrained_enrich, unconstrained_enrich, 'PAV')
    
    return {'constrained_enrich': constrained_enrich,
        'unconstrained_enrich': unconstrained_enrich,
        'PTV diff': ptv_diff, 'PAV diff': pav_diff}

constraint = load_regional_constraint(constraint_path)
de_novos = load_de_novos(de_novos_path, validations_path, keep_indels=True)
de_novos = de_novos[de_novos['hgnc'].isin(constraint['gene'])]

trios = {'male': 2408, 'female': 1885}
for thresh in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
    for ratio in [0.2, 0.4, 0.6, 0.8, 1.0]:
        print(thresh, ratio)
        result = check_enrichment(constraint, de_novos, cache_dir, trios['male'],
            trios['female'], thresh, ratio)
        print(result)
