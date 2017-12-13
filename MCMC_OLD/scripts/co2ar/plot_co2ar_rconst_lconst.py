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
theta_lb_20, theta_rb_20, theta_20 = read_file( '../../out/co2ar/co2ar_rconst_lconst/20_binning/theta.txt' )
theta_lb_50, theta_rb_50, theta_50 = read_file( '../../out/co2ar/co2ar_rconst_lconst/50/theta.txt' )
theta_lb_100, theta_rb_100, theta_100 = read_file( '../../out/co2ar/co2ar_rconst_lconst/100/theta.txt' )

pr_lb_20, pr_rb_20, pr_20 = read_file( '../../out/co2ar/co2ar_rconst_lconst/20_binning/pr.txt' )
pr_lb_50, pr_rb_50, pr_50 = read_file( '../../out/co2ar/co2ar_rconst_lconst/50/pr.txt' )
pr_lb_100, pr_rb_100, pr_100 = read_file( '../../out/co2ar/co2ar_rconst_lconst/100/pr.txt' )

ptheta_lb_20, ptheta_rb_20, ptheta_20 = read_file( '../../out/co2ar/co2ar_rconst_lconst/20_binning/pt.txt' )
ptheta_lb_50, ptheta_rb_50, ptheta_50 = read_file( '../../out/co2ar/co2ar_rconst_lconst/50/pt.txt' )
ptheta_lb_100, ptheta_rb_100, ptheta_100 = read_file( '../../out/co2ar/co2ar_rconst_lconst/100/pt.txt' )

jx_lb_20, jx_rb_20, jx_20 = read_file( '../../out/co2ar/co2ar_rconst_lconst/20_binning/jx.txt' )
jx_lb_50, jx_rb_50, jx_50 = read_file( '../../out/co2ar/co2ar_rconst_lconst/50/jx.txt' )
jx_lb_100, jx_rb_100, jx_100 = read_file( '../../out/co2ar/co2ar_rconst_lconst/100/jx.txt' )

jy_lb_20, jy_rb_20, jy_20 = read_file( '../../out/co2ar/co2ar_rconst_lconst/20_binning/jy.txt' )
jy_lb_50, jy_rb_50, jy_50 = read_file( '../../out/co2ar/co2ar_rconst_lconst/50/jy.txt' )
jy_lb_100, jy_rb_100, jy_100 = read_file( '../../out/co2ar/co2ar_rconst_lconst/100/jy.txt' )

jz_lb_20, jz_rb_20, jz_20 = read_file( '../../out/co2ar/co2ar_rconst_lconst/20_binning/jz.txt' )
jz_lb_50, jz_rb_50, jz_50 = read_file( '../../out/co2ar/co2ar_rconst_lconst/50/jz.txt' )
jz_lb_100, jz_rb_100, jz_100 = read_file( '../../out/co2ar/co2ar_rconst_lconst/100/jz.txt' )
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
theta_med_20 = [ 0.5 * (lb + ub) for lb, ub in zip( theta_lb_20, theta_rb_20)]
theta_med_50 = [ 0.5 * (lb + ub) for lb, ub in zip( theta_lb_50, theta_rb_50)]
theta_med_100 = [ 0.5 * (lb + ub) for lb, ub in zip( theta_lb_100, theta_rb_100)]

pr_med_20 = [ 0.5 * (lb + ub) for lb, ub in zip( pr_lb_20, pr_rb_20) ]
pr_med_50 = [ 0.5 * (lb + ub) for lb, ub in zip( pr_lb_50, pr_rb_50) ]
pr_med_100 = [ 0.5 * (lb + ub) for lb, ub in zip( pr_lb_100, pr_rb_100) ]

ptheta_med_20 = [ 0.5 * (lb + ub) for lb, ub in zip( ptheta_lb_20, ptheta_rb_20) ]
ptheta_med_50 = [ 0.5 * (lb + ub) for lb, ub in zip( ptheta_lb_50, ptheta_rb_50) ]
ptheta_med_100 = [ 0.5 * (lb + ub) for lb, ub in zip( ptheta_lb_100, ptheta_rb_100) ]

jx_med_20 = [ 0.5 * (lb + ub) for lb, ub in zip( jx_lb_20, jx_rb_20) ]
jx_med_50 = [ 0.5 * (lb + ub) for lb, ub in zip( jx_lb_50, jx_rb_50) ]
jx_med_100 = [ 0.5 * (lb + ub) for lb, ub in zip( jx_lb_100, jx_rb_100) ]

jy_med_20 = [ 0.5 * (lb + ub) for lb, ub in zip( jy_lb_20, jy_rb_20) ] 
jy_med_50 = [ 0.5 * (lb + ub) for lb, ub in zip( jy_lb_50, jy_rb_50) ] 
jy_med_100 = [ 0.5 * (lb + ub) for lb, ub in zip( jy_lb_100, jy_rb_100) ] 

jz_med_20 = [ 0.5 * (lb + ub) for lb, ub in zip( jz_lb_20, jz_rb_20) ]
jz_med_50 = [ 0.5 * (lb + ub) for lb, ub in zip( jz_lb_50, jz_rb_50) ]
jz_med_100 = [ 0.5 * (lb + ub) for lb, ub in zip( jz_lb_100, jz_rb_100) ]
#-----------------------------------------------------------------------------------------

lw = 2.0
x = np.linspace( 0, np.pi, 100 )
#-----------------------------------------------------------------------------------------
#plt.subplot( 1, 3, 1 )
#plt.title(r'$\Theta$ distribution')
#plt.plot( theta_med_20, theta_20, linestyle = 'solid', color = 'r' )
#plt.plot( theta_med_50, theta_50, linestyle = 'solid', color = 'g' )
#plt.plot( theta_med_100, theta_100, linestyle = 'solid', color = 'b' )
#plt.plot(x, 0.5 * np.sin(x), linestyle = 'dashed', color = 'k', alpha = 0.7 )
#plt.grid(linestyle = ':', alpha = 0.7)
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#plt.subplot( 1, 3, 2 )
#plt.title(r'$p_R$ distribution')
#plt.plot( pr_med_20, pr_20, linestyle = 'solid', color = 'r' )
#plt.plot( pr_med_50, pr_50, linestyle = 'solid', color = 'g' )
#plt.plot( pr_med_100, pr_100, linestyle = 'solid', color = 'b' )
#plt.grid(linestyle = ':', alpha = 0.7)
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#plt.subplot( 1, 3, 3 )
#plt.title(r'$p_\Theta$ distribution')
#plt.plot( ptheta_med_20, ptheta_20, linestyle = 'solid', color = 'r' )
#plt.plot( ptheta_med_50, ptheta_50, linestyle = 'solid', color = 'g' )
#plt.plot( ptheta_med_100, ptheta_100, linestyle = 'solid', color = 'b' )
#plt.xlim( (-30, 30) ) 
#plt.grid(linestyle = ':', alpha = 0.7)
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
plt.subplot( 1, 3, 1 )
plt.title(r'$J_x$ distribution')
plt.plot( jx_med_20, jx_20, linestyle = 'solid', color = 'r' )
#plt.plot( jx_med_50, jx_50, linestyle = 'solid', color = 'g' )
#plt.plot( jx_med_100, jx_100, linestyle = 'solid', color = 'b' )
plt.grid(linestyle = ':', alpha = 0.7)
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
plt.subplot( 1, 3, 2 )
plt.title(r'$J_y$ distribution')
plt.plot( jy_med_20, jy_20, linestyle = 'solid', color = 'r' )
#plt.plot( jy_med_50, jy_50, linestyle = 'solid', color = 'g' )
#plt.plot( jy_med_100, jy_100, linestyle = 'solid', color = 'b' )
plt.grid(linestyle = ':', alpha = 0.7)
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
plt.subplot( 1, 3, 3 )
plt.title(r'$J_z$ distribution')
plt.plot( jz_med_20, jz_20, linestyle = 'solid', color = 'r' )
#plt.plot( jz_med_50, jz_50, linestyle = 'solid', color = 'g' )
#plt.plot( jz_med_100, jz_100, linestyle = 'solid', color = 'b' )
plt.xlim( (-50, 50) ) 
plt.grid(linestyle = ':', alpha = 0.7)
#-----------------------------------------------------------------------------------------

plt.show()

