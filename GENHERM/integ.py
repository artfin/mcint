with open('output_w.txt', mode = 'r') as inputfile:
    lines_w = inputfile.readlines()

with open('output_x.txt', mode = 'r') as inputfile:
    lines_x = inputfile.readlines()

weights = []
for line in lines_w:
    try:
        weights.append( float(line) )
    except TypeError:
        print("Ouch")

coords = []
for line in lines_x:
    try:
        coords.append( float(line) )
    except TypeError:
        print("Ouch")

sum = 0
for w in weights:
    sum += w

print("Result value of sum is : {0}".format( sum ))
