import matplotlib as mpl
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

def read_file( filename, n ):
    with open(filename, mode = 'r') as inputfile:
        lines = inputfile.readlines()

    lists = [ [] for _ in range(n) ]

    for line in lines:
        if '#' in line:
            continue
        
        data = line.split()
        for i in range(len(data)):
            lists[i].append( float(data[i]) )

    return lists 

#pR_mh, jx_mh, jy_mh = read_file( '../../out/diatomics/diatomics_rconst_nonrep/out.txt', n = 3 )
pR_mh, jx_mh, jy_mh = read_file( '../../out/diatomics/diatomics_rconst_rep/out_12.txt', n = 3 )
pR_exact, jx_exact, jy_exact = read_file( '../../out/diatomics/diatomics_rconst_exact/out.txt', n = 3 )

alpha = 0.3
nbins = 150

fig, ax = plt.subplots( figsize=[8, 6] )

color1 = 'red' #'#e33054'
color2 = 'black' #'#ff9a00'

patch1 = mpatches.Patch( color = color1, label = 'MH', alpha = alpha )
patch2 = mpatches.Patch( color = color2, label = 'Exact', alpha = alpha ) 

x = np.linspace( 0, np.pi )

#-----------------------------------------------------------------------

plt.subplot(1, 3, 1)
plt.title(r'p$_R$ distribution')
N, bins, patches = plt.hist( pR_mh, bins = nbins, normed = True, color = color1, alpha = alpha )
N, bins, patches = plt.hist( pR_exact, bins = nbins, normed = True, color = color2, alpha = alpha )
plt.legend( handles = [patch1, patch2] ) 
plt.grid(linestyle = ':', alpha = 0.7)

plt.subplot(1, 3, 2)
plt.title(r'J$_x$ distribution')
N, bins, patches = plt.hist( jx_mh, bins = nbins, normed = True, color = color1, alpha = alpha)
N, bins, patches = plt.hist( jx_exact, bins = nbins, normed = True, color = color2, alpha = alpha )
plt.legend( handles = [patch1, patch2] )
plt.grid( linestyle = ':', alpha = 0.7 )

plt.subplot(1, 3, 3)
plt.title(r'J$_y$ distribution')
N, bins, patches = plt.hist( jy_mh, bins = nbins, normed = True, color = color1, alpha = alpha)
N, bins, patches = plt.hist( jy_exact, bins = nbins, normed = True, color = color2, alpha = alpha )
plt.legend( handles = [patch1, patch2] )
plt.grid( linestyle = ':', alpha = 0.7 )

#-----------------------------------------------------------------------

plt.show()
