
import pandas

def correct_for_x_chrom(rates, male_n, female_n):
    
    autosomal = 2 * (male_n + female_n)
    female_transmit = male_n + female_n
    male_transmit = female_n
    
    # get scaling factors using the alpha from the most recent SFHS (Scottish
    # Family Health Study) phased de novo data.
    alpha = 3.4
    male_k = 2 / (1 + (1 / alpha))
    female_k = 2 / (1 + alpha)
    
    # correct the non-PAR chrX genes for fewer transmissions and lower rate
    # (dependent on alpha)
    chrX = rates['chrom'].isin(['X', 'chrX'])
    x_factor = ((male_transmit * male_k) + (female_transmit * female_k)) / autosomal
    x_factor = pandas.Series([x_factor] * len(chrX), index=rates.index)
    x_factor[~chrX] = 1
    
    rates['prob'] *= x_factor
    
    return rates
