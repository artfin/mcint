import matplotlib.pyplot as plt
import numpy as np

def read_file( filename, n ):
    with open(filename, mode = 'r') as inputfile:
        lines = inputfile.readlines()

    lists = [ [] for _ in range(n) ]

    for line in lines:
        if '#' in line:
            continue
        
        data = line.split()
        for i in range(len(data)):
            lists[i].append( float(data[i]) )

    return lists 

def transform_theta( theta ):
    theta_i = []
    for t in theta:
        x = t
       
        if ( x < 0 ):
            while ( x < 2 * np.pi ):
                x += 2 * np.pi
       
        if ( x > 2 * np.pi ):
            while ( x > 0 ):
                x -= 2 * np.pi
 
        theta_i.append( x )
    
    return theta_i

theta, pR, pT, jx, jy, jz = read_file( 'co2ar.txt', n = 6 )
pr_d, pt_d, jx_d, jy_d, jz_d = read_file( 'distribution_arco2_danila.txt', n = 5 )
# theta = transform_theta( theta )

fig, ax = plt.subplots( figsize=[8, 6] )

ax.set_title("pR distribution")

n_mh, bins_mh, patches_mh = ax.hist( pR, bins = 100, normed = True, color = '#777777')
n_d, bins_d, patches_d = ax.hist( pr_d, bins = 100, normed = True, color = 'red' )

plt.grid( linestyle = ':', alpha = 0.7)
plt.show()
