#!/usr/bin/python
# -*- encoding: utf8 -*-


import sys
#import petsc4py
#from petsc4py import PETSc

import numpy as np

#from .grid import *
from .set_L import *

#
#==================== Serial solver ============================================
#

class wdinversion():
    """ PV inversion, parallel
    """
    
    def __init__(self, ml):
        """ Setup the solver
        """
                
        self._verbose = ml._verbose
                        
        # create the operator
        self.L = ml.da.createMat()
        #self.I = ml.da.createMat()
        #
        if self._verbose>0:
            print 'Operator L declared'

        # Fill in operator values
        if ml.grid._flag_hgrid_uniform:
            set_L(self.L, ml)
        else:
            set_L_curv(self.L, ml)
            
        #
        if self._verbose>0:
            print 'Operator L filled'

        # global vector for solver
        self._W = ml.da.createGlobalVec()

        # local vectors
        #self._localQ  = qg.da.createLocalVec()
        #self._localPSI  = qg.da.createLocalVec()

        # create solver
        self.ksp = PETSc.KSP()
        self.ksp.create(PETSc.COMM_WORLD)
        self.ksp.setOperators(self.L)
        # use conjugate gradients
        #self.ksp.setType('cg')
        self.ksp.setType('gmres')
        self.ksp.setInitialGuessNonzero(True)
        # and incomplete Cholesky for preconditionning
        #self.ksp.getPC().setType('icc')
        # set tolerances
        #self.ksp.setTolerances(rtol=1e-10) # nope
        #
        #PETSc.Options().setValue('-ksp_view', None)
        #PETSc.Options().setValue('-ksp_monitor', None)
        #PETSc.Options().setValue('-ksp_converged_reason', None)        
        self.ksp.setFromOptions()
         

    def solve(self, omega, W, U, da):
        """ Compute the PV inversion
        """
        # copy Q into Qinv
        W.copy(self._W) 
        # fix boundaries
        #self.set_qinv_bdy(da)
        # add frequency component to operator:
        self.ksp.setOperators(self.L.shift(-np.sqrt(-1)*omega))
        # actually solves the pb
        self.ksp.solve(self._W, U)
        # tmp, test:
        #self.L.mult(PSI, Q)
        if self._verbose>1:
            print 'Inversion done'


    def set_qinv_bdy(self, da):
        """ Set bdy in order to implement boundary conditions
        Set q to 0 along boundaries for inversion, may be an issue
        for time stepping
        """ 
        # 
        q = da.getVecArray(self._Qinv)
        mx, my, mz = da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()        
        # bottom bdy
        if (zs==0):
            k=0
            for j in range(ys, ye):
                for i in range(xs, xe):
                    q[i, j, k] = 0.
        # upper bdy
        if (ze==mz):
            k=mz-1
            for j in range(ys, ye):
                for i in range(xs, xe):
                    q[i, j, k] = 0.


