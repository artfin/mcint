import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import numpy as np

mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode'] = True

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times'
#mpl.rcParams['text.latex.preamble'] = [
    #r"\usepackage{amsmath}"
#]

mpl.rcParams['figure.titlesize'] = "xx-large"
mpl.rcParams['legend.fontsize'] = "large"
mpl.rcParams['axes.labelsize'] = "x-large"
mpl.rcParams['axes.titlesize'] = "large"

mpl.rcParams['xtick.labelsize'] = "large"
mpl.rcParams['ytick.labelsize'] = "large"

def read_data( filename ):

    with open(filename, 'r') as inputfile:
        lines = inputfile.readlines()

    temps = []
    constants = []

    for line in lines:
        if len(line) > 0:
            data = line.split()

            temps.append(float(data[0]))
            constants.append(float(data[1]))

    return temps, constants

fig = plt.figure()

temps1, simple_constants = read_data('simple_constants.dat')
temps2, full_constants = read_data('parallel_full.dat')

simp = np.log( simple_constants )
full = np.log( full_constants )

temps1 = np.log( temps1 )
temps2 = np.log( temps2 ) 

lw = 1.75

#plt.title(r'\textbf{Temperature dependence of K$_p$(N$_2$-N$_2$)}')

#plt.xlabel(r'\textbf{T}, K')
#plt.ylabel(r'\textbf{K}, atm$^{-1}$')
plt.xlabel(r'Log \textbf{T}')
plt.ylabel(r'Log \textbf{K}')

l1, = plt.plot(temps1, simp, color = 'k', linewidth = lw)
l2, = plt.plot(temps2, full, color = 'r', linewidth = lw)

red_patch = mpatches.Patch( color = 'red', label = 'Full phase space' )
black_patch = mpatches.Patch( color = 'black', label = 'Simple')

plt.legend( handles = [red_patch, black_patch]) 

plt.grid(linestyle = ':', alpha = 0.7)

plt.show()

