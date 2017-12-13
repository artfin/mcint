import matplotlib as mpl
import matplotlib.pyplot as plt

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
    with open( filename, mode = 'r') as inputfile:
        lines = inputfile.readlines()

    x = []
    for line in lines:
        x.append( float(line) )

    return x

x = read_file( 'sine_samples.txt' )

bins = 500
alpha = 0.5

fig, ax = plt.subplots( figsize = [8, 6] )

N, bins, patches = ax.hist( x, bins = bins, normed = True, color = '#e33054', alpha = alpha )

plt.grid( linestyle = ':', alpha = 0.7 )

plt.show()
