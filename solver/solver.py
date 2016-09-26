#!/usr/bin/python
# -*- encoding: utf8 -*-


import sys
#import petsc4py
#from petsc4py import PETSc

#import petsc4py
#petsc4py.init(args=sys.argv)

import numpy as np

#from .grid import *
from .set_L import *

#
#==================== Parallel solver ============================================
#

class wdinversion():
    """ PV inversion, parallel
    """
    
    def __init__(self, ml):
        """ Setup the solver
        """
                
        self._verbose = ml._verbose

        #
        self.kx=ml.kx
        self.ky=ml.ky

        # create the operator
        self.L = ml.da.createMat()
        #
        if self._verbose>0:
            print 'Operator L declared'

        # Fill in operator values
        if ml.grid._flag_hgrid_uniform:
            set_L(self.L, ml)
        else:
            #set_L_curv(self.L, ml)
            pass
        #
        if self._verbose>0:
            print 'Operator L filled'
        if self._verbose > 1:
            print self.L.norm()

        # global vector for solver
        #self._W = ml.da.createGlobalVec()

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
        self.ksp.getPC().setType('none')
        #self.ksp.getPC().setReusePreconditioner(True)
        # set tolerances
        #self.ksp.setTolerances(rtol=1e-10) # nope
        #
        for opt in sys.argv[1:]:
            PETSc.Options().setValue(opt, None)
        #PETSc.Options().setValue('-ksp_view', None)
        #PETSc.Options().setValue('-ksp_monitor', None)
        #PETSc.Options().setValue('-ksp_converged_reason', None)        
        self.ksp.setFromOptions()
        #print self.ksp.getOptionsPrefix()
        if self._verbose > 0:
            print 'Solver is set up'
         

    def solve(self, domega, W, U, da):
        """ Compute the PV inversion
        """
        #
        ## copy W into _W
        #W.copy(self._W)
        #
        # add frequency component to operator:
        #
        self.L.shift(-domega*1j)
        self.ksp.setOperators(self.L)
        #
        # actually solves the pb
        #
        #self.ksp.solve(self._W, U)
        self.ksp.solve(W, U)
        #
        # log info or debug
        #
        #print self.ksp.getConvergenceHistory()
        #print self.ksp.getIterationNumber()
        #print str(W.norm())+' / '+str(U.norm())
        #self.ksp.solve(U, self._W)
        # tmp, test:
        #self.L.mult(U, W)
        #
        if self._verbose>0:
            print 'Inversion done'

