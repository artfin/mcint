import sys
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode'] = True

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times'

mpl.rcParams['figure.titlesize'] = "xx-large"
mpl.rcParams['legend.fontsize'] = "large"
mpl.rcParams['axes.labelsize'] = "x-large"
mpl.rcParams['axes.titlesize'] = "large"

mpl.rcParams['xtick.labelsize'] = 15 
mpl.rcParams['ytick.labelsize'] = 15 

def read_file( filename ):
    with open( filename, mode = 'r' ) as inputfile:
        lines = inputfile.readlines()

    lbs, ubs, contents = [], [], []
    for line in lines:
        if len(line) > 0:
            data = line.split()
            lbs.append( float(data[0]) )
            ubs.append( float(data[1]) )
            contents.append( float(data[2]) )

    return lbs, ubs, contents

def mean( lbs, ubs ):
    return [ 0.5 * (lb + ub) for lb, ub in zip(lbs, ubs) ]

def modify( s ):
    return '\Large ' + s.replace("_", "-").split(".txt")[0]

def diff( arr1, arr2 ):
    return [ a - b for a, b in zip(arr1, arr2) ]

def plot_one( filename ):
    lbs, ubs, contents = read_file( filename )
    means = mean( lbs,ubs )

    plt.ylim(( 0.0, 1.1 * max(contents) ))

    plt.scatter( means, contents, color = 'k', s = 3, marker = 'o' )
    plt.grid( linestyle = ':', alpha = 0.7 )
    plt.show()

def plot_two( filename1, filename2 ):
    lbs1, ubs1, contents1 = read_file( filename1 )
    lbs2, ubs2, contents2 = read_file( filename2 )

    means = mean( lbs1, ubs1 )
    
    fig = plt.figure()
    
    plt.subplot(121)
    plt.title(r'\Large Exact and MCMC distribtions' )
    l1 = plt.scatter( means, contents1, color = 'k', s = 3, marker = 'o' )
    l2 = plt.scatter( means, contents2, color = 'r', s = 3, marker = '^' )
    plt.grid( linestyle = ':', alpha = 0.7 )
    
    plt.subplot(122)
    plt.title(r'\Large Difference')
    d = diff( contents1, contents2 ) 
    l3 = plt.scatter( means, d, color = 'g', s = 5, marker = 'x' )
    plt.grid( linestyle = ':', alpha = 0.7 )
    
    fig.legend( (l1, l2, l3), (modify(filename1), modify(filename2), r'\Large difference'), 'lower center', ncol = 3, fancybox = True, shadow = True )

    plt.tight_layout()
    plt.show()
    

if (len(sys.argv) != 2 and len(sys.argv) != 3):
    print(len(sys.argv))
    print "Usage: ./... (string) distirbution-filename"
    sys.exit( 1 )

if ( len(sys.argv) == 2 ):
    plot_one( sys.argv[1] )

if ( len(sys.argv) == 3 ):
    plot_two( sys.argv[1], sys.argv[2] )
