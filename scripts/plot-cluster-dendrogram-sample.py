import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.cluster.hierarchy import dendrogram, linkage, inconsistent
from scipy.spatial.distance import is_valid_dm
import numpy as np
import general_scripts as gs

matplotlib.rcParams['font.family'] = ['sans-serif']
matplotlib.rcParams['font.size'] = 12
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
label_font1 = {'fontname':'Nimbus Sans L', 'fontsize':16 , 'fontweight':'normal' }
label_font2 = {'fontname':'Nimbus Snas L', 'fontsize':20, 'fontweight':'bold' }
rna_font = {'fontname':'cmr', 'fontsize':20 , 'fontweight':'normal' }

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

#plt.close('all')
#rna_labels=['A$_{10}$','U$_{10}$','U$_5$C$_5$','C$_5$U$_5$','C$_{10}$','U$_3$C$_4$U$_3$','C$_3$U$_4$C$_3$','U$_9$','U$_8$','U$_7$','U$_6$','U$_5$','UGU$_8$']
rna_labels=['[G2]\,UGU$_{8}$' , '[G2]\,U$_{10}$', '[G2]\,UGU$_{3}$C$_{5}$' ,
            '[G2]\,U$_{2}$GU$_{7}$' , '[G2]\,U$_{3}$GU$_{6}$' , '[G2]\,U$_{4}$GU$_{5}$' ,
            '[G2]\,U$_{5}$GU$_{4}$' , '[G2]\,U$_{6}$GU$_{3}$' , '[G2]\,U$_{7}$GU$_{2}$' , '[G2]\,U$_{8}$GU' ,
            '[G2]\,UGU$_{7}$' , '[G2]\,UGU$_{6}$' , '[G2]\,UGU$_{5}$' ,'[G2]\,UGU$_{4}$' ,
            '[G2]\,UGU$_{3}$' , '[G2]\,UGU$_{2}$' , '[G2]\,(UGU$_{3}$)$_{2}$' ,
            '[G3]\,A$_{10}$' , '[G3]\,U$_{10}$' , '[G3]\,U$_{5}$C$_{5}$' , '[G3]\,C$_{5}$U$_{5}$' ,
            '[G3]\,C$_{10}$' , '[G3]\,U$_{3}$C$_{4}$U$_{3}$' , '[G3]\,C$_{3}$U$_{4}$C$_{3}$' ,
            '[G3]\,U$_{9}$' , '[G3]\,U$_{8}$' , '[G3]\,U$_{7}$' , '[G3]\,U$_{6}$' , '[G3]\,U$_{5}$' , '[G3]\,UGU$_{8}$' ]

rna_labels_pos=range(len(rna_labels))

# Print matrix
#figC, (ax1, ax2) = plt.subplots( 1, 2, figsize=(6, 3), dpi=300 )
figC = plt.figure( figsize=(6, 6), dpi=600 )

conc_list=['0.0','0.1','0.2','0.4','0.6','0.8','0.9','1.0']
#conc_list=['0.8']
bFirst=True
Xsum=[]
for e in range(len(conc_list)):
    #in_file='Grenoble3_SXL10GS/saxs-0.9-chiMatrix.dat'
    in_file='./saxs-'+conc_list[e]+'-chiMatrix-main.dat'
    X = gs.load_matrix(in_file)
    print type(X)
    # Rescale and symmetrize matrix.
    Xp = np.maximum( np.zeros( X.shape ), X-1.0 )
    Xp = 0.5*( Xp + Xp.T )
    print is_valid_dm(Xp)
    if bFirst:
        bFirst=False
        Xsum=Xp
    else:
        Xsum=Xsum+Xp

Xsum /= len(conc_list)

out_file='output-chimat-combined-avg-main-dendrogram.pdf'
plt.clf()
#ax1  = plt.subplot2grid((1,7), (0,0), colspan=5)
ax2  = plt.subplot2grid((1,1), (0,0) )
plt.subplots_adjust(left=0.18, right=0.99, bottom=0.02, top=0.88, hspace=None)
#plt.figtext(s='Prot:RNA-ratio', x=0.68, y=0.97,
#    horizontalalignment='center', verticalalignment='center', **label_font1 )
#plt.figtext(s='1.0:%s' % conc_list[e], x=0.68, y=0.92,
#    horizontalalignment='center', verticalalignment='center', **label_font2 )

#plt.sca(ax1)
#plt.cla()
#ax1.xaxis.tick_top()
#plt.xticks(rna_labels_pos, rna_labels, rotation='vertical', **rna_font)
#plt.yticks(rna_labels_pos, rna_labels, rotation='horizontal', **rna_font)

#graph = ax1.imshow(X,interpolation='none',cmap=cmap_segmented, vmin=0.0, vmax=5.0)
#if bFirst:
#    divider = make_axes_locatable(ax1)
#    cax = divider.append_axes("right", size="5%", pad=0.05)
#    cbar = figC.colorbar(graph, cax=cax)
#    cbar.ax.tick_params(labelsize=16.)
#    cbar.ax.set_ylabel('pairwise reduced-$\chi$', **label_font1)
#plt.margins(1.0, tight=False)

plt.sca(ax2)
plt.cla()
ax2.xaxis.tick_top()
ax2.xaxis.set_label_position('top')
ax2.tick_params(labelsize=16.)
plt.xlabel('mean cluster distance $\overline{d_{ab}}$', **label_font1)
#plt.xlabel('$d_{ab}=0.5(\chi_{ab}+\chi_{ba})-1$', **label_font1)
ax2.spines['left'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['right'].set_visible(False)
# generate the linkage matrix for the dendrogram.
Z = linkage(Xsum, 'average')
dendrogram(Z,
    labels=rna_labels,
    color_threshold=0.5*max(Z[:,2]),
    orientation='right', #leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=12.  # font size for the x axis labels
)
figC.savefig(out_file)

#print inconsistent(Z)

print "= = Output for %s is complete." % conc_list[e]
