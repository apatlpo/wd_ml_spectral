#!/usr/bin/python
# -*- encoding: utf8 -*-

""" Test the basic features of the library:
Setup of uniform grid
PV inversion of an analytical PV distribution
"""

#import qgsolver.qg as qg
from solver.ml import ml_model
from solver.io import write_nc
import numpy as np

def uniform_grid_runs():
    ''' Tests with uniform grid, closed domains
    '''
    ml = ml_model(hgrid = {'Nx':100, 'Ny':200, 'Lx':1.e6, 'Ly':2.e6},
                  r = 1.e-5)
    #
    ml.set_wd()
    ml.set_ubar()
    ml.set_solver()
    #ml.set_U()
    #
    write_nc([ml.W], ['W'], 'output_wind.nc', ml)
    write_nc([ml.Ubar], ['Ubar'], 'output_ubar.nc', ml)
    #
    ml.solve_uv(domega = 1.e-6)
    #write_nc([ml.u, ml.v], ['u', 'v'], 'output.nc', ml)
    write_nc([ml.U], ['U'], 'output.nc', ml)

    omega = 10**np.linspace(-5,-3,200)
    domega = np.diff(omega)
    for dom in domega:
        print '\n'
        ml.solve_uv(domega=dom)
        write_nc([ml.U], ['U'], 'output.nc', ml, create=False)

    
    return ml

if __name__ == "__main__":
    
    ml = uniform_grid_runs()
    
    if ml._verbose:
        print 'Test done \n'
