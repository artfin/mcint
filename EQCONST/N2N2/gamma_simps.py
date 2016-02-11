import scipy.special as sp
from time import time

import numpy as np

print "gammainc( 3.5, 0.05 ): " + str(sp.gammainc( 3.5, 0.05 ))
print "gamma( 1.1 ): " + str(sp.gamma(1.1))

a = 2.5

def s_int(f, x0, x1, step):
    
    res = 0
    
    lb = x0
    rb = x0 + step

    while ( rb < x1 ):
        res += s_int_p( f, lb, rb )
    
        # print "lb: " + str(lb)
        # print "rb: " + str(rb)
        # print "res: " + str(res)

        lb += step
        rb += step
    
    return res
    
def s_int_p(f, x0, x1):
    h = x1 - x0

    f0 = f(x0)
    f1 = f(x1)
    f12 = (f0 + f1) / 2
    
    return h/6 * (f0 + 4 * f12 + f1)

def f(x):
    #return x**2
    return x**1.5 * np.exp(- x)

integ = s_int(f, 0.0, 50.0, 0.01)
print "integ: " + str(integ)

abs_error = sp.gamma(2.5) - integ
rel_error = abs_error / integ

print "abs_error: " + str(abs_error)
print "rel_error: " + str(rel_error)

