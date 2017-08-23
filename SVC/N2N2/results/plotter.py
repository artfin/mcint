import matplotlib as mpl
import matplotlib.pyplot as plt

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
    svcs = []

    for line in lines:
        if len(line) > 0:
            data = line.split()

            temps.append(float(data[0]))
            svcs.append(float(data[1]))

    return temps, svcs

avoird_temp = [75.0, 80.0, 90.0, 100.0, 110.0, 125.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0, 600.0, 700.0]
avoird_svc_class = [-280.2, -247.0, -196.8, -161.0, -134.1, -104.4, -71.77, -35.79, -16.58, -4.75, 8.86, 16.30, 20.87, 23.90]

fig = plt.figure()

temperatures, svcs = read_data('svc.dat')

lw = 1.75

plt.title(r'\textbf{Second Virial Coefficient of N$_2$}')

plt.xlabel(r'\textbf{T}, K')
plt.ylabel(r'\textbf{SVC}, cm$^{3}$ \ mol$^{-1}$')

l1, = plt.plot(temperatures, svcs, color = 'k', linewidth = lw)
plt.scatter(avoird_temp, avoird_svc_class, color = 'r', marker = 'x')

plt.grid(linestyle = ':', alpha = 0.7)

plt.show()


