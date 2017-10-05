import matplotlib.pyplot as plt

def read_file( filename ):
    with open( filename, mode = 'r' ) as inputfile:
        lines = inputfile.readlines()

    lists = [ [] for _ in range(3) ]
    for line in lines:
        data = line.split()

        for i in range(3):
            lists[i].append( float(data[i]) )

    return lists

lbounds, ubounds, contents = read_file( '../../out/co2ar/co2ar_rconst_lagrange_binning/theta.txt' )
#lbounds, ubounds, contents = read_file( '../../out/co2ar/co2ar_rconst_lagrange_binning/omegax.txt' )
center_of_bin = [ 0.5 * (l + u) for l, u in zip(lbounds, ubounds) ]

lw = 2.0
plt.plot( center_of_bin, contents, linewidth = lw, color = 'k' )
plt.grid( linestyle = ':', alpha = 0.7 )

plt.show()
