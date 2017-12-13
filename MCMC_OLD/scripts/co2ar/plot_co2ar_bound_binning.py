import matplotlib.pyplot as plt

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
r_lb_50, r_ub_50, r_50 = read_file( '../../out/co2ar/co2ar_bound_binning/50/r.txt')
r_lb_75, r_ub_75, r_75 = read_file( '../../out/co2ar/co2ar_bound_binning/75/r.txt')
r_lb_100, r_ub_100, r_100 = read_file( '../../out/co2ar/co2ar_bound_binning/100/r.txt')
r_lb_150, r_ub_150, r_150 = read_file( '../../out/co2ar/co2ar_bound_binning/150/r.txt')
r_lb_300, r_ub_300, r_300 = read_file( '../../out/co2ar/co2ar_bound_binning/300/r.txt')
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
theta_lb_50, theta_rb_50, theta_50 = read_file( '../../out/co2ar/co2ar_bound_binning/50/theta.txt' )
theta_lb_75, theta_rb_75, theta_75 = read_file( '../../out/co2ar/co2ar_bound_binning/75/theta.txt' )
theta_lb_100, theta_rb_100, theta_100 = read_file( '../../out/co2ar/co2ar_bound_binning/100/theta.txt' )
theta_lb_150, theta_rb_150, theta_150 = read_file( '../../out/co2ar/co2ar_bound_binning/150/theta.txt' )
theta_lb_300, theta_rb_300, theta_300 = read_file( '../../out/co2ar/co2ar_bound_binning/300/theta.txt' )
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
pr_lb_50, pr_rb_50, pr_50 = read_file( '../../out/co2ar/co2ar_bound_binning/50/pr.txt' )
pr_lb_75, pr_rb_75, pr_75 = read_file( '../../out/co2ar/co2ar_bound_binning/75/pr.txt' )
pr_lb_100, pr_rb_100, pr_100 = read_file( '../../out/co2ar/co2ar_bound_binning/100/pr.txt' )
pr_lb_150, pr_rb_150, pr_150 = read_file( '../../out/co2ar/co2ar_bound_binning/150/pr.txt' )
pr_lb_300, pr_rb_300, pr_300 = read_file( '../../out/co2ar/co2ar_bound_binning/300/pr.txt' )
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
jx_lb_50, jx_rb_50, jx_50 = read_file( '../../out/co2ar/co2ar_bound_binning/50/jx.txt' )
jx_lb_75, jx_rb_75, jx_75 = read_file( '../../out/co2ar/co2ar_bound_binning/75/jx.txt' )
jx_lb_100, jx_rb_100, jx_100 = read_file( '../../out/co2ar/co2ar_bound_binning/100/jx.txt' )
jx_lb_150, jx_rb_150, jx_150 = read_file( '../../out/co2ar/co2ar_bound_binning/150/jx.txt' )
jx_lb_300, jx_rb_300, jx_300 = read_file( '../../out/co2ar/co2ar_bound_binning/300/jx.txt' )
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
jz_lb_50, jz_rb_50, jz_50 = read_file( '../../out/co2ar/co2ar_bound_binning/50/jz.txt' )
jz_lb_75, jz_rb_75, jz_75 = read_file( '../../out/co2ar/co2ar_bound_binning/75/jz.txt' )
jz_lb_100, jz_rb_100, jz_100 = read_file( '../../out/co2ar/co2ar_bound_binning/100/jz.txt' )
jz_lb_150, jz_rb_150, jz_150 = read_file( '../../out/co2ar/co2ar_bound_binning/150/jz.txt' )
jz_lb_300, jz_rb_300, jz_300 = read_file( '../../out/co2ar/co2ar_bound_binning/300/jz.txt' )
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
r_med_50 = [ 0.5 * (lb + ub) for lb, ub in zip( r_lb_50, r_ub_50 )]
r_med_75 = [ 0.5 * (lb + ub) for lb, ub in zip( r_lb_75, r_ub_75 )]
r_med_100 = [ 0.5 * (lb + ub) for lb, ub in zip( r_lb_100, r_ub_100 )]
r_med_150 = [ 0.5 * (lb + ub) for lb, ub in zip( r_lb_150, r_ub_150 )]
r_med_300 = [ 0.5 * (lb + ub) for lb, ub in zip( r_lb_300, r_ub_300 )]
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
theta_med_50 = [ 0.5 * (lb + ub) for lb, ub in zip( theta_lb_50, theta_rb_50)]
theta_med_75 = [ 0.5 * (lb + ub) for lb, ub in zip( theta_lb_75, theta_rb_75)]
theta_med_100 = [ 0.5 * (lb + ub) for lb, ub in zip( theta_lb_100, theta_rb_100)]
theta_med_150 = [ 0.5 * (lb + ub) for lb, ub in zip( theta_lb_150, theta_rb_150)]
theta_med_300 = [ 0.5 * (lb + ub) for lb, ub in zip( theta_lb_300, theta_rb_300)]
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
pr_med_50 = [ 0.5 * (lb + ub) for lb, ub in zip( pr_lb_50, pr_rb_50) ]
pr_med_75 = [ 0.5 * (lb + ub) for lb, ub in zip( pr_lb_75, pr_rb_75) ]
pr_med_100 = [ 0.5 * (lb + ub) for lb, ub in zip( pr_lb_100, pr_rb_100) ]
pr_med_150 = [ 0.5 * (lb + ub) for lb, ub in zip( pr_lb_150, pr_rb_150) ]
pr_med_300 = [ 0.5 * (lb + ub) for lb, ub in zip( pr_lb_300, pr_rb_300) ]
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
jx_med_50 = [ 0.5 * (lb + ub) for lb, ub in zip( jx_lb_50, jx_rb_50) ]
jx_med_75 = [ 0.5 * (lb + ub) for lb, ub in zip( jx_lb_75, jx_rb_75) ]
jx_med_100 = [ 0.5 * (lb + ub) for lb, ub in zip( jx_lb_100, jx_rb_100) ]
jx_med_150 = [ 0.5 * (lb + ub) for lb, ub in zip( jx_lb_150, jx_rb_150) ]
jx_med_300 = [ 0.5 * (lb + ub) for lb, ub in zip( jx_lb_300, jx_rb_300) ]
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
jz_med_50 = [ 0.5 * (lb + ub) for lb, ub in zip( jz_lb_50, jz_rb_50) ]
jz_med_75 = [ 0.5 * (lb + ub) for lb, ub in zip( jz_lb_75, jz_rb_75) ]
jz_med_100 = [ 0.5 * (lb + ub) for lb, ub in zip( jz_lb_100, jz_rb_100) ]
jz_med_150 = [ 0.5 * (lb + ub) for lb, ub in zip( jz_lb_150, jz_rb_150) ]
jz_med_300 = [ 0.5 * (lb + ub) for lb, ub in zip( jz_lb_300, jz_rb_300) ]
#-----------------------------------------------------------------------------------------

lw = 2.0
#-----------------------------------------------------------------------------------------
#plt.plot( r_med_50, r_50, linewidth = lw, color = 'k' )
#plt.plot( r_med_75, r_75, linewidth = lw, color = 'y' )
#plt.plot( r_med_100, r_100, linewidth = lw, color = 'b' )
#plt.plot( r_med_150, r_150, linewidth = lw, color = 'r' )
#plt.plot( r_med_300, r_300, linewidth = lw, color = 'g' )
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
#plt.plot( theta_med_50, theta_50, linewidth = lw, color = 'k' )
#plt.plot( theta_med_75, theta_75, linewidth = lw, color = 'y' )
#plt.plot( theta_med_100, theta_100, linewidth = lw, color = 'b' )
#plt.plot( theta_med_150, theta_150, linewidth = lw, color = 'r' )
#plt.plot( theta_med_300, theta_300, linewidth = lw, color = 'g' )
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
#plt.plot( pr_med_50, pr_50, linewidth = lw, color = 'k' )
#plt.plot( pr_med_75, pr_75, linewidth = lw, color = 'y' )
#plt.plot( pr_med_100, pr_100, linewidth = lw, color = 'b' )
#plt.plot( pr_med_150, pr_150, linewidth = lw, color = 'r' )
#plt.plot( pr_med_300, pr_300, linewidth = lw, color = 'g' )
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
#plt.plot( jx_med_50, jx_50, linewidth = lw, color = 'k' )
#plt.plot( jx_med_75, jx_75, linewidth = lw, color = 'y' )
#plt.plot( jx_med_100, jx_100, linewidth = lw, color = 'b' )
#plt.plot( jx_med_150, jx_150, linewidth = lw, color = 'r' )
#plt.plot( jx_med_300, jx_300, linewidth = lw, color = 'g' )
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
plt.plot( jz_med_50, jz_50, linewidth = lw, color = 'k' )
plt.plot( jz_med_75, jz_75, linewidth = lw, color = 'y' )
plt.plot( jz_med_100, jz_100, linewidth = lw, color = 'b' )
plt.plot( jz_med_150, jz_150, linewidth = lw, color = 'r' )
plt.plot( jz_med_300, jz_300, linewidth = lw, color = 'g' )
#-----------------------------------------------------------------------------------------

plt.grid( linestyle = ':', alpha = 0.7 )
plt.show()

