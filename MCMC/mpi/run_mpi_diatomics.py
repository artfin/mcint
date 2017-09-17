#!/usr/bin/python
import subprocess

nsteps = 1000000
burnin = 1000000
alpha = 7.0

instances = 4

programs = {
    'main': ('mpi_diatomics', instances),
}

sys_call = '{0} -n {1} ./{2} {3} {4} {5}'.format('mpirun', programs['main'][1], programs['main'][0], nsteps, burnin, alpha)
print(sys_call)

subprocess.call([sys_call], shell = True)


