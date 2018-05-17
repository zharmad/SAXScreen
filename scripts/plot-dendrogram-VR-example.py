import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import is_valid_dm
import numpy as np
import general_scripts as gs

matplotlib.rcParams['font.family'] = ['sans-serif']
matplotlib.rcParams['font.size'] = 10
#matplotlib.rcParams['mathtext.default'] = ['regular']
matplotlib.rcParams['text.usetex'] = True
#params = {'mathtext.default': 'bf' }
#plt.rcParams.update(params)
matplotlib.rcParams['text.latex.preamble'] = [
    r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
    r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
    r'\usepackage{helvet}',    # set the normal font here
    r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
    r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]

label_font0 = {'fontname':'Nimbus Mono L', 'fontsize':6, 'fontweight':'normal' }
label_font1 = {'fontname':'Nimbus Mono L', 'fontsize':10, 'fontweight':'normal' }
label_font2 = {'fontname':'Nimbus Mono L', 'fontsize':14, 'fontweight':'bold' }

colorbar_min=0.0
colorbar_max=1.0
colorbar_max2=0.5
metric='V_R'
# Adopted from gnuplot schemes
#set palette defined (0 "black", 1 '#2233bb', 2 "yellow",  3 "red", 4 "pink", 5 "white")
# i.e. #00-00-00 , #22-33-bb, #FF-FF-00, #FF-00-00, #FF-C0-CB, #FF-FF-FF
cdict = {'red':   [(0.0, 0.000, 0.000),
                   (0.2, 0.133, 0.133),
                   (0.4, 1.000, 1.000),
                   (0.6, 0.800, 0.800),
                   (0.8, 1.000, 1.000),
                   (1.0, 1.000, 1.000)],
         'green': [(0.0, 0.000, 0.000),
                   (0.2, 0.199, 0.199),
                   (0.4, 1.000, 1.000),
                   (0.6, 0.000, 0.000),
                   (0.8, 0.750, 0.750),
                   (1.0, 1.000, 1.000)],
         'blue':  [(0.0, 0.000, 0.000),
                   (0.2, 0.730, 0.730),
                   (0.4, 0.000, 0.000),
                   (0.6, 0.000, 0.000),
                   (0.8, 0.793, 0.793),
                   (1.0, 1.000, 1.000)]
        }
#cdict = {'red': [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)], 'green': [(0.0, 1.0, 1.0), (1.0, 0.0, 0.0)], 'blue': [(0.0, 1.0, 1.0), (1.0, 0.0, 0.0)]}

cmap_segmented = matplotlib.colors.LinearSegmentedColormap('Gnuplot', cdict)

fig_dims=(5,5)
dpi=300

#box_style = dict(boxstyle='round', facecolor='#FFEEDD')
# place a text box in upper left in axes coords
#ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
#                verticalalignment='top', bbox=props)
rna_font = {'fontname':'Nimbus Mono L'}

#plt.close('all')
rna_labels=['A$_{10}$','U$_{10}$','U$_5$C$_5$','C$_5$U$_5$','C$_{10}$','U$_3$C$_4$U$_3$','C$_3$U$_4$C$_3$','U$_9$','U$_8$','U$_7$','U$_6$','U$_5$','UGU$_8$']
rna_labels_pos=range(len(rna_labels))

# Print matrix
#figC, (ax1, ax2) = plt.subplots( 1, 2, figsize=(6, 3), dpi=300 )
figC = plt.figure( figsize=(5, 3), dpi=300 )

Xaverage=[]
conc_list=['0.0','0.1','0.2','0.4','0.6','0.8','0.9','1.0']
#conc_list=['1.0']
numConc=len(conc_list)
bFirst=True

for e in range(numConc):
    in_file='./fitted_'+metric+'_'+conc_list[e]+'_matrix.dat'
    out_file='matrix_'+metric+'_'+conc_list[e]+'.pdf'

    X = gs.load_matrix(in_file)
    # Rescale and symmetrize matrix.
    Xp = np.maximum( np.zeros( X.shape ), X )
    Xp = 0.5*( Xp + Xp.T )
    print is_valid_dm(Xp)

    if bFirst:
        bFirst=False
        Xaverage=Xp
    else:
        Xaverage=Xaverage+Xp

    plt.clf()
    ax1  = plt.subplot2grid((1,7), (0,0), colspan=5)
    ax2  = plt.subplot2grid((1,7), (0,5), colspan=2)
    plt.subplots_adjust(left=0.04, right=0.99, bottom=0.03, top=0.76, wspace=3.0, hspace=None)
    plt.figtext(s='Prot:RNA-ratio', x=0.68, y=0.95,
        horizontalalignment='center', verticalalignment='center', **label_font1 )
    plt.figtext(s='1.0:%s' % conc_list[e], x=0.68, y=0.90,
        horizontalalignment='center', verticalalignment='center', **label_font2 )

    #fig1 = plt.figure(figsize=(5,5), dpi=300 )
    plt.sca(ax1)
    plt.cla()
    ax1.xaxis.tick_top()
    plt.xticks(rna_labels_pos, rna_labels, rotation='vertical', **rna_font)
    plt.yticks(rna_labels_pos, rna_labels, rotation='horizontal', **rna_font)

    graph = ax1.imshow(X,interpolation='none',cmap=cmap_segmented, vmin=colorbar_min, vmax=colorbar_max)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = figC.colorbar(graph, cax=cax)
    cbar.ax.set_ylabel('pairwise reduced-$\chi$')
    #plt.margins(1.0, tight=False)

    plt.sca(ax2)
    plt.cla()
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')
#    plt.xlabel('cluster distance')
    plt.xlabel('$d_{ab}=0.5(\chi_{ab}+\chi_{ba})$', **label_font0)
    ax2.spines['left'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    # generate the linkage matrix for the dendrogram.
    Z = linkage(Xp, 'average')
    dendrogram(Z,
        labels=rna_labels,
        color_threshold=0.5*max(Z[:,2]),
        orientation='right', #leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=10
    )
    figC.savefig(out_file)

    #print inconsistent(Z)

    print "= = Output for %s is complete." % conc_list[e]

# Do final plot with aggregate data
Xaverage /= numConc
out_file='matrix_'+metric+'_aggregate.pdf'

plt.clf()
ax1  = plt.subplot2grid((1,7), (0,0), colspan=5)
ax2  = plt.subplot2grid((1,7), (0,5), colspan=2)
plt.subplots_adjust(left=0.04, right=0.99, bottom=0.03, top=0.76, wspace=3.0, hspace=None)
plt.figtext(s='Aggregate', x=0.68, y=0.95,
    horizontalalignment='center', verticalalignment='center', **label_font1 )

#fig1 = plt.figure(figsize=(5,5), dpi=300 )
plt.sca(ax1)
plt.cla()
ax1.xaxis.tick_top()
plt.xticks(rna_labels_pos, rna_labels, rotation='vertical', **rna_font)
plt.yticks(rna_labels_pos, rna_labels, rotation='horizontal', **rna_font)

graph = ax1.imshow(Xaverage,interpolation='none',cmap=cmap_segmented, vmin=colorbar_min, vmax=colorbar_max2)
#plt.margins(1.0, tight=False)
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = figC.colorbar(graph, cax=cax)
cbar.ax.set_ylabel('pairwise reduced-$\chi$')

plt.sca(ax2)
plt.cla()
ax2.xaxis.tick_top()
ax2.xaxis.set_label_position('top')
plt.xlabel('cluster distance')
plt.xlabel('$d_{ab}=0.5(\chi_{ab}+\chi_{ba})$', **label_font0)
ax2.spines['left'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['right'].set_visible(False)
# generate the linkage matrix for the dendrogram.
Z = linkage(Xaverage, 'average')
dendrogram(Z,
    labels=rna_labels,
    color_threshold=0.5*max(Z[:,2]),
    orientation='right', #leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=10
)
figC.savefig(out_file)
