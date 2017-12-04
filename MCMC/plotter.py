import sys
import matplotlib.pyplot as plt

if ( len(sys.argv) != 2 ):
    print( "Not enough arguments." )
    sys.exit( 1 )

def read_file( filename ):
    with open( filename, mode = 'r' ) as inputfile:
        lines = inputfile.readlines()

    lb, ub, contents = [], [], []
    for line in lines:
        if ( len(line) > 0 ):
            data  = line.split()

            lb.append( float(data[0]) )
            ub.append( float(data[1]) )
            contents.append( float(data[2]) )

    return lb, ub, contents

lb, ub, contents = read_file( sys.argv[1] )

plt.scatter( lb, contents, color = 'k', s = 20 );
plt.show()
