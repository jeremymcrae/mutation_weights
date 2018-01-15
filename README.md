I wanted to compare the enrichment of missense DNMs at various CADD thresholds.
My baseline for enrichment was the enrichment of truncating mutations within
dominant DDG2P genes. PTV DNMs in these genes are ~35-fold enriched compared to
numbers expected from null mutation rates.

To set this up I started by getting the mutation rates and CADD scores for all 
sites within dominant DDG2P genes. 

`get_rates.py` finds per nucleotide/alt null mutation rates for all sites within
dominant DDG2P genes.

`get_cadd_scores.py` finds the CADD scores for the same set of genes.

These just make it easier to work with subset of sites and CADD scores, rather
than every site in the exome or genome. Quicker to load files.

Given the expected mutation rates and CADD scores per site, I wanted to find the
enrichment of missense de novos within different CADD bins. This is the
`check_weights.py` script. Requires the rates, CADD scores, de novos (and 
validations), table of trios, and families (to count name of male and female
probands).
