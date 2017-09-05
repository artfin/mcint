import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np

def target( x ):
    return np.sin(x)**2 * np.sin(2*x)**2 * mlab.normpdf(x, 0, 1)

with open('vars.txt', 'r') as inputfile:
    lines = inputfile.readlines()

ys = []
for line in lines:
    y = float(line)
    if y > 0:
        ys.append( y )

fig, ax = plt.subplots(figsize=[8,6])

N, bins, patches = ax.hist(ys, bins = 100, color = '#777777', normed = True )

xs = np.linspace( min(ys), max(ys) )
targets = [ target(x) for x in xs ]
ax.plot(xs, targets, lw = 2, color = 'red')

plt.show()
