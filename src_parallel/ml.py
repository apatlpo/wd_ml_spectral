#!/usr/bin/python
# -*- encoding: utf8 -*-

from .grid import *
from .solver import *

import petsc4py
#from Cython.Compiler.Main import verbose
petsc4py.init(sys.argv)
from petsc4py import PETSc

import numpy as np
from .io import read_nc_petsc
#from netCDF4 import Dataset


class ml_model():
    """ QG object
    """
    
    def __init__(self,
                 hgrid = None,
                 f0 = 7e-5, K = 1.e2,
                 f0_file = None,
                 verbose = 1,
                 ):
        """ ML object creation
        Parameters:
        """

        #
        # Build grid object
        #        
        self.grid = grid(hgrid) 

        #
        # init petsc
        #
        
        #OptDB = PETSc.Options()

        # setup tiling
        #self.da = PETSc.DMDA().create([self.grid.Nx, self.grid.Ny, self.grid.Nz],
        #                              stencil_width=2)
        self.da = PETSc.DMDA().create(sizes = [self.grid.Nx, self.grid.Ny],
                                      proc_sizes = [2,4], dof=2,
                                      stencil_width = 2)
        # http://lists.mcs.anl.gov/pipermail/petsc-dev/2016-April/018889.html
        self.comm = self.da.getComm()
        self.rank = self.comm.getRank()
        # print tiling information
        if self.rank is 0 and verbose>0:
            print 'PETSc DMDA created'
            print 'The 2D grid is tiled according to (nproc_x, nproc_y) : '\
                    +str(self.da.proc_sizes) 
            #print 'rank='+str(self.rank)+' ranges='+str(self.da.ranges)
        
        # for lon/lat grids should load metric terms over tiles
        if not self.grid._flag_hgrid_uniform:
            self.grid.load_metric_terms(self.da, self.comm)

        # print out grid information
        if self.rank is 0 and verbose>0:
            self._verbose=verbose
        else:
            self._verbose=0
        #
        if self._verbose>0:
            # general information
            print 'A ML model object is being created'
            # print out grid parameters
            print self.grid

        #
        # Coriolis
        #
        if f0_file is not None:
            if self._verbose:
                print 'Reads f0 from '+f0_file
            #
            self.f0 = read_nc('f0', f0_file)
            #
            if self._verbose:
                print 'Reads Coriolis parameter f from '+f0_file
            self.grid.load_coriolis_parameter(f0_file, self.da, self.comm)
        else:
            self.f0 = f0
        #
        self.K = K

        #
        # declare petsc vectors
        #
        
        # wind vector
        self.W = self.da.createGlobalVec()
        # wind-driven current
        self.U = self.da.createGlobalVec()
        # background current
        self.Ubar = self.da.createGlobalVec()
        
        #
        # initiate pv inversion solver
        #
        self.wdinv = wdinversion(self)

        
    
        
    def set_psi(self, analytical_psi=True, file_psi=None):
        """ Set psi to a given value
        """
        if self._verbose:
            print 'Set psi to ... ?'
        #if analytical_psi:
            

    def set_q(self, analytical_q=True, file_q=None):
        """ Set q to a given value
        """
        #
        if file_q is not None:
            if self._verbose:
                print 'Set q from file '+file_q+' ...'
            read_nc_petsc(self.Q, 'q', file_q, self)
        elif analytical_q:
            if self._verbose:
                print 'Set q analytically '
            self.set_q_analytically()


    def set_q_analytically(self):
        """ Set q analytically
        """
        q = self.da.getVecArray(self.Q)
        mx, my, mz = self.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()
        #
        if self._verbose:
            pass
            #print 'Set q analytically \n'
        for k in range(zs, ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    q[i, j, k] = 1.e-5*np.exp(-((i/float(mx-1)-0.5)**2 
                                              + (j/float(my-1)-0.5)**2)/0.1**2)
                    q[i, j, k] *= np.sin(i/float(mx-1)*np.pi) 
                    q[i, j, k] *= np.sin(2*j/float(my-1)*np.pi)
            

    def set_q_bdy(self):
        """ Reset q at boundaries such that dq/dn=0 """
        #
        q = self.da.getVecArray(self.Q)
        #
        mx, my, mz = self.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()
        #
        # south bdy
        if (ys==0):
            j=0
            for k in range(zs, ze):
                for i in range(xs, xe):
                    q[i, j, k] = q[i,j+1,k]
        # north bdy
        if (ye==my):
            j=my-1
            for k in range(zs, ze):
                for i in range(xs, xe):
                    q[i, j, k] = q[i,j-1,k]
        # west bdy
        if (xs==0):
            i=0
            for k in range(zs, ze):
                for j in range(ys, ye):
                    q[i, j, k] = q[i+1,j,k]
        # east bdy
        if (xe==mx):
            i=mx-1
            for k in range(zs, ze):
                for j in range(ys, ye):
                    q[i, j, k] = q[i-1,j,k]                


    def solve_uv(self, omega=None):
        """ wrapper around solver solve method
        """
        self.wdinv.solve(omega, self.W, self.U, self.da)

    
    def get_uv(self):
        """ Compute horizontal velocities
        """



