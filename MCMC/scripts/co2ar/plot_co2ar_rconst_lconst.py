import matplotlib.pyplot as plt
import numpy as np

def read_file ( filename ):
    with open( filename, mode = 'r' ) as inputfile:
        lines = inputfile.readlines()

    lists = [ [] for _ in range(3) ]
    for line in lines:

        if '#' in line:
            continue

        data = line.split()
        for i in range(3):
            lists[i].append( float(data[i]) )

    return lists

#-----------------------------------------------------------------------------------------
theta_lb, theta_rb, theta = read_file( '../../out/co2ar/co2ar_rconst_lconst/theta.txt' )
pr_lb, pr_rb, pr = read_file( '../../out/co2ar/co2ar_rconst_lconst/pr.txt' )
ptheta_lb, ptheta_rb, ptheta = read_file( '../../out/co2ar/co2ar_rconst_lconst/pt.txt' )
jx_lb, jx_rb, jx = read_file( '../../out/co2ar/co2ar_rconst_lconst/jx.txt' )
jy_lb, jy_rb, jy = read_file( '../../out/co2ar/co2ar_rconst_lconst/jy.txt' )
jz_lb, jz_rb, jz = read_file( '../../out/co2ar/co2ar_rconst_lconst/jz.txt' )
l2_lb, l2_rb, l2 = read_file( '../../out/co2ar/co2ar_rconst_lconst/l2.txt' )
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
theta_med = [ 0.5 * (lb + ub) for lb, ub in zip( theta_lb, theta_rb)]
pr_med = [ 0.5 * (lb + ub) for lb, ub in zip( pr_lb, pr_rb) ]
ptheta_med = [ 0.5 * (lb + ub) for lb, ub in zip( ptheta_lb, ptheta_rb) ]
jx_med = [ 0.5 * (lb + ub) for lb, ub in zip( jx_lb, jx_rb) ]
jy_med = [ 0.5 * (lb + ub) for lb, ub in zip( jy_lb, jy_rb) ] 
jz_med = [ 0.5 * (lb + ub) for lb, ub in zip( jz_lb, jz_rb) ]
l2_med = [ 0.5 * (lb + ub) for lb, ub in zip( l2_lb, l2_rb) ]
#-----------------------------------------------------------------------------------------

lw = 2.0
x = np.linspace( 0, np.pi, 100 )
#-----------------------------------------------------------------------------------------
#plt.subplot( 1, 3, 1 )
#plt.title(r'$\Theta$ distribution')
#plt.plot( theta_med, theta, linestyle = 'solid', color = 'k' )
#plt.plot(x, 0.5 * np.sin(x), linestyle = 'dashed', color = 'r', alpha = 0.7 )
#plt.grid(linestyle = ':', alpha = 0.7)
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#plt.subplot( 1, 3, 2 )
#plt.title(r'$p_R$ distribution')
#plt.plot( pr_med, pr, linestyle = 'solid', color = 'k' )
#plt.grid(linestyle = ':', alpha = 0.7)
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#plt.subplot(1, 3, 3 )
#plt.title( r'$L^2$ distribution' )
#plt.plot( l2_med, l2, linestyle = 'solid', color = 'k' )
#plt.grid( linestyle = ':', alpha = 0.7 )

#plt.subplot( 1, 3, 3 )
#plt.title(r'$p_\Theta$ distribution')
#plt.plot( ptheta_med, ptheta, linestyle = 'solid', color = 'k' )
#plt.xlim( (-10, 10) ) 
#plt.grid(linestyle = ':', alpha = 0.7)
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
plt.subplot( 1, 3, 1 )
plt.title(r'$J_x$ distribution')
plt.plot( jx_med, jx, linestyle = 'solid', color = 'k' )
plt.grid(linestyle = ':', alpha = 0.7)
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
plt.subplot( 1, 3, 2 )
plt.title(r'$J_y$ distribution')
plt.plot( jy_med, jy, linestyle = 'solid', color = 'k' )
plt.grid(linestyle = ':', alpha = 0.7)
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
plt.subplot( 1, 3, 3 )
plt.title(r'$J_z$ distribution')
plt.plot( jz_med, jz, linestyle = 'solid', color = 'k' )
plt.xlim( (-10, 10) ) 
plt.grid(linestyle = ':', alpha = 0.7)
#-----------------------------------------------------------------------------------------

plt.show()

