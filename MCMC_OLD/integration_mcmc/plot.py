import matplotlib.pyplot as plt

def read_file( filename ):
    with open(filename, mode = 'r') as inputfile:
        lines = inputfile.readlines()

    x = []
    for line in lines:
        if '#' in line:
            continue
        data = line.split()
        x.append( float(data[0]) )

    return x

x = read_file( 'out.dat' )
 
alpha = 0.5

fig, ax = plt.subplots( figsize = [8, 6] )

N, bins, patches = plt.hist( x, bins = 100, normed = True, color = '#e33054', alpha = alpha )

plt.show()
