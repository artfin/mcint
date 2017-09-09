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
jx_d1, jy_d1, pR_d1 = read_file( 'd_with_s.txt' )
jx_d2, jy_d2, pR_d2 = read_file( 'd_without_s.txt' )

fig, ax = plt.subplots( figsize=[8, 6] )

ax.set_title('pR distribution')

# cumulative distribution 
n_mh, bins_mh, patches_mh = ax.hist( pR_mh, bins = 100, normed = True, color = 'red' )
N_d1, bins_d1, patches_d1 = ax.hist( pR_d1, bins = 100, normed = True, color = '#777777' )

red_patch = mpatches.Patch(color = 'red', label = 'Metropolis-Hastings')
green_patch = mpatches.Patch(color = '#777777', label = 'Danila formulae')
plt.legend(handles=[red_patch, green_patch])

plt.grid(linestyle = ':', alpha = 0.7)

# random walk picture
# plt.plot( jx, color = 'k', linewidth = 1.2 )

plt.show()

