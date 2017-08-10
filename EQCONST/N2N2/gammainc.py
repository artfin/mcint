import scipy.special as sp
from time import time

start = time()
print sp.gammainc(1.0, 0.1)
print time() - start 

