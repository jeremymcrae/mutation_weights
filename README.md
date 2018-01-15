I wanted to find good scores to classify the severity of missense mutations.
I used CADD scores as a reasonable proxy of severity, but the point was to 
scale them to values that are more appropriate for our DD cohort, i.e to
find values more suitable for comparing the severity of different scores.

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
`check_weights.py` script. Requires rates, CADD scores, de novos (with failed
validations removed), and counts of male and female probands (determined from
tables of trios and families.

I used CADD score bins of 0-5, 5-10 etc, and found missense sites with CADD
scores in each window, then calculated the enrichment as observed over expected.

I also stratified by whether the sites were in regions with regional missense
constraint (Kaitlin's regional constraint, that is). Basically doing the same thing
as above, but once where I only included sites under regional constraint, and 
another time with only sites not under regional constraint.
