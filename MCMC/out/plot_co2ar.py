import matplotlib.pyplot as plt
import numpy as np

def read_file( filename ):
    with open(filename, mode = 'r') as inputfile:
        lines = inputfile.readlines()

    theta = []
    pR = []
    pT = []
    jx = []
    jy = []
    jz = []

    for line in lines:
        if '#' in line:
            continue
        
        data = line.split()
        theta.append( float(data[0]) )
        pR.append( float(data[1]) )
        pT.append( float(data[2]) )
        jx.append( float(data[3]) )
        jy.append( float(data[4]) )
        jz.append( float(data[5]) )

    return theta, pR, pT, jx, jy, jz

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

theta, pR, pT, jx, jy, jz = read_file( 'co2ar.txt' )
# theta = transform_theta( theta )

fig, ax = plt.subplots( figsize=[8, 6] )

ax.set_title("pR distribution")

n_mh, bins_mh, patches_mh = ax.hist( pT, bins = 100, normed = True, color = '#777777')

plt.grid( linestyle = ':', alpha = 0.7)
plt.show()
