#!/usr/bin/python
# -*- encoding: utf8 -*-

""" Test the basic features of the library:
Setup of uniform grid
PV inversion of an analytical PV distribution
"""

#import qgsolver.qg as qg
from wd_solver.ml import ml_model
from wd_solver.io import write_nc


def uniform_grid_runs():
    ''' Tests with uniform grid, closed domains
    '''
    ml = ml_model(hgrid = {'Nx':150, 'Ny':100}, vgrid = {'Nz':3 },
            K = 0.e0)
    #
    ml.set_wd()
    ml.solve_uv()
    write_nc([ml.u, ml.v], ['u', 'v'], 'output.nc', ml)
    #
    
    return ml

if __name__ == "__main__":
    
    ml = uniform_grid_runs()
    
    if ml._verbose:
        print 'Test done \n'
