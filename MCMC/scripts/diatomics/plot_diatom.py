import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import numpy as np

mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode'] = True

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times'

mpl.rcParams['figure.titlesize'] = 'xx-large'
mpl.rcParams['legend.fontsize'] = 'large'
mpl.rcParams['axes.labelsize'] = 'x-large'
mpl.rcParams['axes.titlesize'] = 'large'

mpl.rcParams['xtick.labelsize'] = 'large'
mpl.rcParams['ytick.labelsize'] = 'large'

def plot_distributions( dist1, dist2, bins_number = 300, alpha = 0.1, title = "" ):
    hist_1, _ = np.histogram( dist1, bins = bins_number, normed = True )
    hist_2, _ = np.histogram( dist2, bins = bins_number, normed = True ) 
    hist_diff = hist_1 - hist_2

    lim1 = np.arange(min(dist1), max(dist1), ( max(dist1) - min(dist1)) / bins_number )
    lim2 = np.arange(min(dist2), max(dist2), ( max(dist2) - min(dist2)) / bins_number )

    plt.bar( lim1, hist_1, color = '#e33054', alpha = alpha )
    plt.bar( lim2, hist_2, color = 'blue', alpha = alpha )

    plt.title( title )
    plt.grid( linestyle = ':', alpha = 0.7 )

#plt.subplot(1, 3, 1)
#plot_distributions( jx_mh, jx_d1, bins_number = 700, alpha = 0.2, title = r"J$_x$ distribution")

#plt.savefig('foo.png')

def read_file( filename ):
    with open(filename, mode = 'r') as inputfile:
        lines = inputfile.readlines()

    jx, jy, pR = [], [], []
    for line in lines:
        if '#' in line: 
            continue

        data = line.split()

        jx.append( float(data[0]) )
        jy.append( float(data[1]) )
        pR.append( float(data[2]) )

    return jx, jy, pR

jx_mh, jy_mh, pR_mh = read_file( 'diatomics_mh.txt' )
jx_d1, jy_d1, pR_d1 = read_file( 'diatomics_exact.txt' )

alpha = 0.5

fig, ax = plt.subplots( figsize=[8, 6] )

patch1 = mpatches.Patch( color='#e33054', label = 'MH', alpha = alpha )
patch2 = mpatches.Patch( color='#ff9a00', label = 'Exact', alpha = alpha ) 

plt.subplot(1, 3, 1)
plt.title(r'p$_R$ distribution')
N, bins, patches = plt.hist( pR_mh, bins = 500, normed = True, color = '#e33054', alpha = alpha )
N, bins, patches = plt.hist( pR_d1, bins = 500, normed = True, color = '#ff9a00', alpha = alpha )
plt.legend( handles = [patch1, patch2] ) 
plt.grid(linestyle = ':', alpha = 0.7)

plt.subplot(1, 3, 2)
plt.title(r'J$_x$ distribution')
N, bins, patches = plt.hist( jx_mh, bins = 500, normed = True, color = '#e33054', alpha = alpha )
N, bins, patches = plt.hist( jx_d1, bins = 500, normed = True, color = '#ff9a00', alpha = alpha )  
plt.legend( handles = [patch1, patch2] ) 
plt.grid(linestyle = ':', alpha = 0.7)

plt.subplot(1, 3, 3)
plt.title(r'J$_y$ distribution')
N, bins, patches = plt.hist( jy_mh, bins = 500, normed = True, color = '#e33054', alpha = alpha )
N, bins, patches = plt.hist( jy_d1, bins = 500, normed = True, color = '#ff9a00', alpha = alpha )
plt.legend( handles = [patch1, patch2] ) 
plt.grid(linestyle = ':', alpha = 0.7)

plt.tight_layout()
plt.show()


