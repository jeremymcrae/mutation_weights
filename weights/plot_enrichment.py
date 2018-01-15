
import os

import matplotlib
matplotlib.use('agg')

from matplotlib import pyplot

def plot_enrichment(threshold, enrichment, cdf, path='temp.pdf'):
    ''' plots enrichment of de novos at different CADD thresholds
    '''
    
    fig = pyplot.figure(figsize=(6, 6))
    ax = fig.gca()
    cdf_ax = ax.twinx()
    
    plot1 = ax.plot(threshold, enrichment, linestyle='None', marker='.',
        markersize=10, label='enrichment')
    plot2 = cdf_ax.plot(threshold, cdf, color='gray', linewidth=1.5, label='CDF')
    
    e = ax.set_ylim(0, max(enrichment) * 1.05)
    e = cdf_ax.set_ylim(0, max(cdf) * 1.05)
    
    e = ax.spines['top'].set_visible(False)
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['left'].set_linewidth(1.5)
    e = ax.spines['bottom'].set_linewidth(1.5)
    
    e = cdf_ax.spines['left'].set_visible(False)
    e = cdf_ax.spines['top'].set_visible(False)
    e = cdf_ax.spines['right'].set_linewidth(1.5)
    e = cdf_ax.spines['right'].set_color('gray')
    
    e = ax.xaxis.set_ticks_position('bottom')
    e = ax.yaxis.set_ticks_position('left')
    
    e = ax.tick_params(axis='both', length=7, width=1.5, which='major', labelsize=15)
    e = cdf_ax.tick_params(axis='both', length=7, width=1.5, which='major', labelsize=15, colors='gray')
    e = ax.set_xlabel('CADD threshold', fontsize=15)
    e = ax.set_ylabel('observed/expected', fontsize=15)
    e = cdf_ax.set_ylabel('Cumlative proportion of missense', fontsize=15, color='gray')
    e = cdf_ax.yaxis.label.set_color('gray')
    
    # add plot legend, merged across axes
    labs = [ l.get_label() for l in plot1 + plot2 ]
    e = ax.legend(plot1 + plot2, labs, frameon=False)
    
    filetype = os.path.splitext(path)[1][1:]
    
    fig.savefig(path, dpi=300, format=filetype, bbox_inches='tight', pad_inches=0,
        transparent=True)
