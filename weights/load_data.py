
import pandas

from mupit.open_ddd_data import standardise_ddd_de_novos

def load_rates(path):
    ''' get a DataFrame of mutation rates by site
    '''
    
    rates = pandas.read_table(path)
    rates['chrom'] = rates['chrom'].astype(str)
    
    # fix an issue with duplicated entries
    non_dups = rates[~rates[['chrom', 'pos', 'alt']].duplicated(keep=False)]
    dups = rates[rates[['chrom', 'pos', 'alt']].duplicated(keep=False)]
    dups = dups[~dups[['chrom', 'pos', 'alt']].duplicated(keep='first')]
    rates = non_dups.append(dups, ignore_index=True)
    
    return rates.sort_values(['symbol', 'chrom', 'pos'])

def load_cadd(path):
    ''' get cadd scores by sites
    '''
    cadd = pandas.read_table(path)
    cadd['chrom'] = cadd['chrom'].astype(str)
    cadd['pos'] = cadd['pos'].astype(int)
    
    return cadd

def count_trios(trios_path, families_path):
    ''' count the number of male and female trios in the cohort
    '''
    
    trios = pandas.read_table(trios_path)
    families = pandas.read_table(families_path)
    
    families = families[families['individual_id'].isin(trios['proband_stable_id'])]
    families['sex'] = families['sex'].map({'M': 'male', 'F': 'female'})
    
    return dict(families.sex.value_counts())

def load_de_novos(de_novos_path, validations_path, keep_indels=False):
    '''
    '''
    variants = standardise_ddd_de_novos(de_novos_path)
    
    validations = pandas.read_table(validations_path, sep="\t")
    variants = variants.merge(validations, how="left",
        on=["person_id", "chrom", "start_pos", "end_pos", "ref_allele",
            "alt_allele", "hgnc", "consequence"])
    
    variants = variants[~variants["status"].isin(["false_positive", "inherited"])]
    del variants["status"]
    
    if not keep_indels:
        # we only want SNVs for this
        variants = variants[variants['type'] == 'snv']
    
    return variants

def load_regional_constraint(path):
    ''' load in constraint regions, and determine chromosome coordinates
    '''
    data = pandas.read_table(path)
    data['gene'] = data['gene'].astype(str)
    return data
