import matplotlib.pyplot as plt

def read_file( filename ):
    with open(filename, mode = 'r') as inputfile:
        lines = inputfile.readlines()

    xs, ys = [], []
    for line in lines:
        data = line.split()
        xs.append( float(data[0]) )
        ys.append( float(data[1]) )

    return xs, ys

xs, ys = read_file("2d.txt")

plt.scatter( xs, ys, color = 'k', s = 1 )
plt.show()
