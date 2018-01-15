
from scipy.stats import chi2

def aa_to_chrom(tx, region):
    ''' convert an amino acid region of a transcript to chromosomal coordinates
    
    Args:
        tx: Transcript object for a gene
        region: start and end amino acid positions (dash-separated) e.g. '1-260'
    
    Returns:
        tuple of start and end chromosomal coordinates
    '''
    start, end = region.split('-')
    start = (int(start) - 1) * 3
    end = ((int(end) - 1) * 3) + 2
    
    return tx.get_position_on_chrom(start), tx.get_position_on_chrom(end)

def get_constrained_positions(tx, group, threshold=1e-4, ratio_threshold=1.0):
    ''' get all the positions in the constrained regions
    '''
    constraint_sites = set([])
    for i, row in group.iterrows():
        row = row.squeeze()
        p_value = chi2.sf(row['chisq_diff_null'], df=1)
        if p_value > threshold:
            continue
        
        if row['obs_exp'] > ratio_threshold:
            continue
        
        region = row['amino_acids']
        start, end = aa_to_chrom(tx, region)
        if tx.get_strand() == '-':
            start, end = end, start
        
        for x in range(start, end + 1):
            if tx.in_coding_region(x):
                constraint_sites.add(x)
    
    return constraint_sites

def classify_de_novos_by_constraint(constraint, de_novos, cache_dir, threshold=1e-4, ratio=1.0):
    ''' determine whether de novos fall within constrained regions
    '''
    sites = set([])
    ensembl = EnsemblRequest(cache_dir, 'grch37')
    for tx_id, group in constraint.groupby('transcript'):
        tx = construct_gene_object(ensembl, tx_id.split('.')[0])
        constrained_sites = get_constrained_positions(tx, group, threshold, ratio)
        
        chrom = tx.get_chrom()
        sites |= set(( (chrom, x) for x in constrained_sites ))
    
    return [ x in sites for x in zip(de_novos['chrom'], de_novos['start_pos']) ]
