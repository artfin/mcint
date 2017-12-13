import matplotlib.pyplot as plt

def read_file( filename ):
    with open(filename, mode = 'r') as inputfile:
        lines = inputfile.readlines()

    x = []
    y = []
    for line in lines:
        data = line.split()
        x.append( float(data[0]) )
        y.append( float(data[1]) )

    return x, y

x, y = read_file( '2d-exp-sin.txt' )

fig, ax = plt.subplots( figsize = [8, 6] )

# n_x, bins_x, patches_x = ax.hist( x, bins = 100, normed = True, color = '#777777' )

n_y, bins_y, paches_y = ax.hist( y, bins = 100, normed = True, color = '#777777' )

plt.show()
