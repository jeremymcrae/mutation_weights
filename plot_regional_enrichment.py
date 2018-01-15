
import pandas

import matplotlib
matplotlib.use('agg')

from matplotlib import pyplot

def annotate_heatmap(ax, mesh, ages, increment):
    """Add textual labels with the value in each cell."""
    
    ages = [ x + increment/2 for x in ages[:-1] ]
    mesh.update_scalarmappable()
    xpos, ypos = numpy.meshgrid(ages, ages)
    for x, y, val, color in zip(xpos.flat, ypos.flat,
                                mesh.get_array(), mesh.get_facecolors()):
        if val is not numpy.ma.masked:
            l = relative_luminance(color)
            text_color = ".0" if l > .408 else "w"
            val = "{:.2f}".format(val)
            text_kwargs = dict(color=text_color, ha="center", va="center")
            ax.text(x, y, val, **text_kwargs)

# data collected from checking regional enrichment
data = pandas.DataFrame([[1.00E-02, 0.2, 4.01665026],
    [1.00E-02, 0.4, 3.509335743],
    [1.00E-02, 0.6, 2.033735364],
    [1.00E-02, 0.8, 1.355519331],
    [1.00E-02, 1, 2.098436074],
    [0.001, 0.2, 4.061536082],
    [0.001, 0.4, 3.532499489],
    [0.001, 0.6, 2.040980741],
    [0.001, 0.8, 1.369466278],
    [0.001, 1, 2.210569219],
    [0.0001, 0.2, 4.144171706],
    [0.0001, 0.4, 3.606730769],
    [0.0001, 0.6, 2.098480202],
    [0.0001, 0.8, 1.357861709],
    [0.0001, 1, 2.525013991],
    [1.00E-05, 0.2, 4.312007278],
    [1.00E-05, 0.4, 3.672265977],
    [1.00E-05, 0.6, 2.175069304],
    [1.00E-05, 0.8, 1.457008325],
    [1.00E-05, 1, 2.129166777],
    [1.00E-06, 0.2, 4.53175934],
    [1.00E-06, 0.4, 3.79919622],
    [1.00E-06, 0.6, 2.221172739],
    [1.00E-06, 0.8, 1.510796265],
    [1.00E-06, 1, 2.06468105]], columns=['p-value', 'ratio', 'enrichment'])

data = data.pivot(index='p-value', columns='ratio', values='enrichment')

fig = pyplot.figure(figsize=(6, 6))
ax = fig.gca()
mesh = ax.pcolor(data, cmap=pyplot.cm.Blues, alpha=0.8)
e = pyplot.colorbar(mesh, ticks=[1.0, 2.0, 3.0, 4.0])

e = ax.set_xticklabels(data.columns, minor=False, fontsize=15)
e = ax.set_yticklabels(data.index, minor=False, fontsize=15)

e = ax.set_xlabel('Observed/expected ratio', fontsize=15)
e = ax.set_ylabel('observed/expected p-value', fontsize=15)

fig.savefig('regional_constraint.heatmap.pdf', format='pdf', bbox_inches='tight', pad_inches=0,
    transparent=True)
