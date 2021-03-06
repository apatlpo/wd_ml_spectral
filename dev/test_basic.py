#!/usr/bin/python
# -*- encoding: utf8 -*-

""" Test the basic features of the library:
Setup of uniform grid
PV inversion of an analytical PV distribution
"""

import sys
#import qgsolver.qg as qg
from solver.ml import ml_model
from solver.io import write_nc
import numpy as np


# midlatitude scenario
f0_midlatitude = 2. * 2. * np.pi / 86400. * np.sin(45. * np.pi / 180.)
beta0_midlatitude = 2. * 2. * np.pi / 86400. * np.cos(45. * np.pi / 180.) / (6371. * 1e3)


def uniform_grid_runs(prefix='data/output', analytical_ubar=0, r=5.e-6):
    ''' Tests with uniform grid, closed domains
    '''
    #
    sys.stdout = open(prefix+'.log', 'w')
    ml = ml_model(hgrid = {'Nx':200, 'Ny':400, 'Lx':1.e6, 'Ly':2.e6}, r = r, f0=f0_midlatitude)
    #
    ml.set_wd()
    ml.set_ubar(analytical_ubar=analytical_ubar)
    ml.set_solver()
    #ml.set_U()
    #
    write_nc([ml.W], ['W'], prefix+'_wind.nc', ml)
    write_nc([ml.Ubar], ['Ubar'], prefix+'_ubar.nc', ml)
    #
    ml.solve_uv(domega = 1.e-6)
    write_nc([ml.U], ['U'], prefix+'_uv.nc', ml)

    omega = 10**np.linspace(-5,-3,200)
    domega = np.diff(omega)
    for dom in domega:
        #print '\n'
        if ml._verbose: print '-- omega = %.1e' %(ml.omega+dom)
        ml.solve_uv(domega=dom)
        write_nc([ml.U], ['U'], prefix+'_uv.nc', ml, create=False)

    return ml


if __name__ == "__main__":


    # base case no eddy
    ml = uniform_grid_runs(prefix='data/base')

    # base case with r*10
    ml = uniform_grid_runs(prefix='data/base_10r', r=5.e-5)

    # base case with r/10
    ml = uniform_grid_runs(prefix='data/base_0.1r', r=5.e-7)

    # 1 eddy
    ml = uniform_grid_runs(prefix='data/eddy', analytical_ubar=1)

    # multiple eddies
    ml = uniform_grid_runs(prefix='data/meddies', analytical_ubar=2)


    if ml._verbose:
        print ' f0 = ', str(f0_midlatitude), ' [1/s]'
        print ' beta0 = ', str(beta0_midlatitude), ' [1/s/m]'
        print 'Test done \n'
