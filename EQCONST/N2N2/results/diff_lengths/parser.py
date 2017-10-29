with open( '400.out', mode = 'r' ) as inputfile:
    lines = inputfile.readlines()

temps = []
eqconsts = []

for line in lines:
    if 'Temperature' in line:
        data = line.split()

        temps.append(data[1].split(';')[0])
        eqconsts.append(float(data[-1]))

with open( 'result.out', mode = 'w') as out:
    for t, e in zip( temps, eqconsts ):
        out.write( str(t) + " " + str(e) + '\n')
