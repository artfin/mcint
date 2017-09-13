import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

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

jx_mh, jy_mh, pR_mh = read_file( 'mh.txt' )
jx_d1, jy_d1, pR_d1 = read_file( 'diatomics_exact.txt' )

green_patch = mpatches.Patch(color = '#777777', label = 'Danila formulae')

fig, ax = plt.subplots( figsize=[8, 6] )

plt.subplot(1, 3, 1)
plt.title(r'p$_R$ distribution')
N_d1, bins_d1, patches_d1 = plt.hist( pR_d1, bins = 500, normed = True, color = '#e33054' )
plt.grid(linestyle = ':', alpha = 0.7)

plt.subplot(1, 3, 2)
plt.title(r'J$_x$ distribution')
N_d1, bins_d1, patches_d1 = plt.hist( jx_d1, bins = 500, normed = True, color = '#396161' )
plt.grid(linestyle = ':', alpha = 0.7)

plt.subplot(1, 3, 3)
plt.title(r'J$_y$ distribution')
N_d1, bins_d1, patches_d1 = plt.hist( jy_d1, bins = 500, normed = True, color = '#232458' )
plt.grid(linestyle = ':', alpha = 0.7)

plt.tight_layout()

plt.show()

