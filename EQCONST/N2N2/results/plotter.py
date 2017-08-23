import matplotlib as mpl
import matplotlib.pyplot as plt

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

temperatures, simple_constants = read_data('simple_constants.dat')
temps_full, full_constants = read_data('full_constants.dat')

simp = np.log(simple_constants)
full = np.log(full_constants)



lw = 1.75

plt.title(r'\textbf{Temperature dependence of equilibrium constant of N$_2$-N$_2$ complex}')

plt.xlabel(r'\textbf{T}, K')
plt.ylabel(r'\textbf{Simple constant}, atm$^{-1}$')

#l1, = plt.plot(temperatures, simp, color = 'k', linewidth = lw)
#l2, = plt.plot(temps_full, full, color = 'r', linewidth = lw)

l, = plt.plot(temperatures, diff, color = 'k', linewidth = lw)

plt.grid(linestyle = ':', alpha = 0.7)

plt.show()

