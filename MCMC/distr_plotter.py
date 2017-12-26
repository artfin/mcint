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

def read_file( filename, n ):
    with open( filename, mode = 'r' ) as inputfile:
        lines = inputfile.readlines()

    lists = [ [] for _ in range(n) ]

    for line in lines:
        if len(line) > 0:
            data = line.split()
            
            for i in range(n):
                lists[i].append( float(data[i]) )

    return lists 

def mean( lbs, ubs ):
    return [ 0.5 * (lb + ub) for lb, ub in zip(lbs, ubs) ]

def modify( s ):
    return '\Large ' + s.replace("_", "-").split(".txt")[0].split("/")[-1]

def diff( arr1, arr2 ):
    return [ a - b for a, b in zip(arr1, arr2) ]

def rel_diff( exact, arr ):
    res = []
    for e, a in zip(exact, arr):
        if a != 0:
            res.append( (a-e)/a )
        else:
            res.append( 0 )
    return res 

def rel_diff_unequal_size( means1, contents1, means2, contents2 ):
    res = []

    for i in range(len(contents2)):
        closest = min( means1, key = lambda x: abs(x - means2[i]) )
        
        if ( abs(closest - means2[i]) > 0.1 ):
            res.append( 0 )
            continue

        index_closest = means1.index( closest )
        
        if contents2[i] == 0:
            res.append( 0 )
        else:
            res.append( (contents2[i]-contents1[index_closest])/contents2[i] )
        
    return res

def plot_one( filename ):
    lbs, ubs, contents = read_file( filename, 3 )
    means = mean( lbs,ubs )

    plt.ylim(( 0.0, 1.1 * max(contents) ))

    plt.scatter( means, contents, color = 'k', s = 3, marker = 'o' )
    plt.grid( linestyle = ':', alpha = 0.7 )
    plt.show()

def plot_two( filename1, filename2 ):
    try:
        lbs1, ubs1, contents1 = read_file( filename1, 3 )
        means1 = mean( lbs1, ubs1 )
    except IndexError:
        means1, contents1 = read_file( filename1, 2 )

    try:
        lbs2, ubs2, contents2 = read_file( filename2, 3 )
        means2 = mean( lbs2, ubs2 )
    except IndexError:
        means2, contents2 = read_file( filename2, 2 ) 

    fig = plt.figure()
    
    plt.subplot(121)
    plt.title(r'\Large Exact and MCMC distribtions' )
    l1 = plt.scatter( means1, contents1, color = 'k', s = 3, marker = 'o' )
    l2 = plt.scatter( means2, contents2, color = 'r', s = 3, marker = '^' )
    plt.grid( linestyle = ':', alpha = 0.7 )
    
    plt.subplot(122)
    plt.title(r'\Large Relative difference')
    d = rel_diff_unequal_size( means1, contents1, means2, contents2 ) 
    l3 = plt.scatter( means2, d, color = 'g', s = 5, marker = 'x' )
    plt.grid( linestyle = ':', alpha = 0.7 )
    
    fig.legend( (l1, l2, l3), (modify(filename1), modify(filename2), r'\Large difference'), 'lower center', ncol = 3, fancybox = True, shadow = True )

    plt.tight_layout()
    plt.show()
   
def plot_three( filename1, filename2, filename3 ):
    lbs1, ubs1, contents1 = read_file(filename1, 3 ) # exact
    lbs2, ubs2, contents2 = read_file(filename2, 3 ) # mcmc
    means3, contents3 = read_file(filename3, 2 ) # integrated
    
    means1 = mean( lbs1, ubs1 )
    means2 = mean( lbs2, ubs2 )

    fig = plt.figure()

    plt.subplot(121)
    plt.title(r'\Large Exact, Integrated and MCMC distributions' )
    l1 = plt.scatter( means1, contents1, color = 'k', s = 3, marker = 'o' )
    l2 = plt.scatter( means2, contents2, color = 'r', s = 3, marker = 'o' )
    l3 = plt.scatter( means3, contents3, color = 'b', s = 3, marker = 'o' )
    plt.grid( linestyle = ':', alpha = 0.7 )

    plt.subplot(122)
    plt.title(r'\Large Relative difference (to exact distribution)')
    rel_mcmc = rel_diff( contents1, contents2 )
    rel_integrated = rel_diff_unequal_size( means1, contents1, means3, contents3 )
    l4 = plt.scatter( means1, rel_mcmc, color = 'r', s = 3, marker = 'x' )
    l5 = plt.scatter( means3, rel_integrated, color = 'b', s = 3, marker = 'x' )
    plt.ylim( (-1, 1) )
    plt.grid( linestyle = ':', alpha = 0.7 )

    fig.legend( (l1, l2, l3, l4), (modify(filename1), modify(filename2), modify(filename3), r'\Large relative diff'), 'lower center', ncol = 4, fancybox = True, shadow = True )

    plt.tight_layout()
    plt.show()

if (len(sys.argv) != 2 and len(sys.argv) != 3 and len(sys.argv) != 4):
    print(len(sys.argv))
    print "Usage: ./... (string) distribution-filename"
    print "       ./... (string) fname1 (string) fname2"
    print "       ./... (string) fname1 (string) fname2 (string) fname3"
    sys.exit( 1 )

if ( len(sys.argv) == 2 ):
    plot_one( sys.argv[1] )

if ( len(sys.argv) == 3 ):
    plot_two( sys.argv[1], sys.argv[2] )

if ( len(sys.argv) == 4 ):
    plot_three( sys.argv[1], sys.argv[2], sys.argv[3] )
