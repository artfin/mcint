import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import numpy as np

def read_data( filename ):

    with open(filename, 'r') as inputfile:
        lines = inputfile.readlines()

    temps = []
    constants = []

    for line in lines:
        if len(line) > 0:
            data = line.split()

            temps.append(float(data[0]))
            constants.append(float(data[1]))

    return temps, constants

temp_simple, const_simple = read_data( '../simple_constants.dat' )
temp_simple_log = np.log( temp_simple )
const_simple_log = np.log( const_simple )

temp_025, const_025 = read_data( 'full_0_25.txt' )
temp, const_050 = read_data( 'full_0_5.txt' )
temp, const_2 = read_data( 'full_2_0.txt' )
temp, const_4 = read_data( 'full_4_0.txt' )
temp, const_100 = read_data( 'full_100.txt' )
temp, const_400 = read_data( 'full_400.txt' )

temp_log = np.log( temp )
temp_025_log = np.log( temp_025 )

const_025_log = np.log( const_025 )
const_050_log = np.log( const_050 )
const_2_log = np.log( const_2 )
const_4_log = np.log( const_4 )
const_100_log = np.log( const_100 )
const_400_log = np.log( const_400 )

plt.plot( temp_simple_log, const_simple_log, color = 'k' )

plt.plot( temp_025_log, const_025_log )
plt.plot( temp_log, const_050_log )
plt.plot( temp_log, const_2_log )
plt.plot( temp_log, const_4_log )
plt.plot( temp_log, const_100_log, color = 'r' )
plt.plot( temp_log, const_400_log, color = 'y' )

plt.grid( linestyle = ':', alpha = 0.7 )

plt.show()
