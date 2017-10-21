import matplotlib.pyplot as plt
import numpy as np

def read_file( filename ):
    
    with open( filename, mode = 'r' ) as inputfile:
        lines = inputfile.readlines()

    temps = []
    eqconsts = []

    for line in lines:
        if '#' in line:
            continue 

        data = line.split()

        temps.append( float(data[0]) )
        eqconsts.append( float(data[1]) )

    return temps, eqconsts

temp, eq025 = read_file( 'full_0_25.txt' )
temp, eq4 = read_file( 'full_4_0.txt' )

temp1, eq_full = read_file( '../parallel_full.dat' )
temp, eq_simple = read_file( '../simple_constants.dat' )

temp_log = np.log( temp )
temp1_log = np.log( temp1 )

eq025_log = np.log( eq025 )
eq4_log = np.log( eq4 )

eq_full_log = np.log( eq_full )
eq_simple_log = np.log( eq_simple )

lw = 2.0
plt.plot( temp_log, eq025_log, color = 'b', linewidth = lw )
plt.plot( temp_log, eq4_log, color = 'r', linewidth = lw )

plt.plot( temp1_log, eq_full_log, color = 'k', linewidth = lw, linestyle = 'solid' )
plt.plot( temp_log, eq_simple_log, color = 'k', linewidth = lw, linestyle = 'dashed' )

plt.show()

