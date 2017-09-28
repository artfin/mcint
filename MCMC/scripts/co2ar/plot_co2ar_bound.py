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

r_mh, theta_mh, pR_mh, pT_mh, jx_mh, jy_mh, jz_mh = read_file( '../../out/co2ar/co2ar_bound/out.txt', n = 7 )

alpha = 0.3

fig, ax = plt.subplots( figsize=[8, 6] )

patch1 = mpatches.Patch( color='#e33054', label = 'MH', alpha = alpha )

x = np.linspace( 0, np.pi )

# ---------------------------------------------------------------
#plt.title(r'R distribution')
#N, bins, patches = plt.hist( r_mh, bins = 500, normed = True, color = '#e33054', alpha = alpha, range = [5.5, 13.0])
#plt.legend( handles = [patch1] )
#plt.grid( linestyle = ':', alpha = 0.7 )
# ---------------------------------------------------------------

# ---------------------------------------------------------------
plt.subplot(1, 3, 1)
plt.title(r'$\Theta$ distribution')
N, bins, patches = plt.hist( theta_mh, bins = 500, normed = True, color = '#e33054', alpha = alpha )
plt.plot(x, 0.5 * np.sin(x), linestyle = 'dashed', color = 'r', alpha = 0.7 )
plt.legend( handles = [patch1] ) 
plt.grid(linestyle = ':', alpha = 0.7)

plt.subplot(1, 3, 2)
plt.title(r'p$_R$ distribution')
N, bins, patches = plt.hist( pR_mh, bins = 500, normed = True, color = '#e33054', alpha = alpha )
plt.legend( handles = [patch1] ) 
plt.grid(linestyle = ':', alpha = 0.7)

plt.subplot(1, 3, 3)
plt.title(r'p$_T$ distribution')
N, bins, patches = plt.hist( pT_mh, bins = 500, normed = True, color = '#e33054', alpha = alpha )
plt.legend( handles = [patch1] ) 
plt.grid(linestyle = ':', alpha = 0.7)
# ---------------------------------------------------------------

# ---------------------------------------------------------------
#plt.subplot(1, 3, 1)
#plt.title(r'J$_x$ distribution')
#N, bins, patches = plt.hist( jx_mh, bins = 500, normed = True, color = '#e33054', alpha = alpha)
#plt.legend( handles = [patch1] )
#plt.grid( linestyle = ':', alpha = 0.7 )

#plt.subplot(1, 3, 2)
#plt.title(r'J$_y$ distribution')
#N, bins, patches = plt.hist( jy_mh, bins = 500, normed = True, color = '#e33054', alpha = alpha)
#plt.legend( handles = [patch1] )
#plt.grid( linestyle = ':', alpha = 0.7 )

#plt.subplot(1, 3, 3)
#plt.title(r'J$_z$ distribution')
#N, bins, patches = plt.hist( jz_mh, bins = 500, normed = True, color = '#e33054', alpha = alpha)
#plt.legend( handles = [patch1] )
#plt.grid( linestyle = ':', alpha = 0.7 )
# ---------------------------------------------------------------

plt.show()
