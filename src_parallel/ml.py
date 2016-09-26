#!/usr/bin/python
# -*- encoding: utf8 -*-

from .grid import *
from .solver import *

#import petsc4py
#from Cython.Compiler.Main import verbose
#petsc4py.init(args=sys.argv)
#petsc4py.init(args=sys.argv[1:]) # does not seem to work
from petsc4py import PETSc

import numpy as np
from .io import read_nc_petsc
#from netCDF4 import Dataset


class ml_model():
    """ Wind driven mixed layer object
    """
    
    def __init__(self,
                 hgrid = None,
                 f0 = 7.e-5, K = 1.e2,
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

        #print sys.argv[1:]
        #petsc4py.init(sys.argv[1:])
        OptDB = PETSc.Options()
        #for o in OptDB:
        #    print o.prefix
        #print OptDB.prefix

        # test the complex build of PETSc
        #print(PETSc.ScalarType)

        # setup tiling
        #self.da = PETSc.DMDA().create([self.grid.Nx, self.grid.Ny, self.grid.Nz],
        #                              stencil_width=2)
        self.da = PETSc.DMDA().create(sizes = [self.grid.Nx, self.grid.Ny, 2],
                                      proc_sizes = [2,4,1], dof=1,
                                      stencil_width = 1)
        # http://lists.mcs.anl.gov/pipermail/petsc-dev/2016-April/018889.html
        self.comm = self.da.getComm()
        self.rank = self.comm.getRank()
        # print tiling information
        if self.rank is 0 and verbose>0:
            print 'PETSc DMDA created'
            print 'The 2D grid is tiled according to (nproc_x, nproc_y) : '\
                    +str(self.da.proc_sizes)
            # Note that the grid is 3D in fact (third dimension describes vectors)
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

        # indexes along third dimension corresponding to x and y directions
        self.kx=0; self.ky=1

        # wind vector
        self.W = self.da.createGlobalVec()
        self.W.zeroEntries()
        # wind-driven current
        self.U = self.da.createGlobalVec()
        self.U.zeroEntries()
        # background current
        self.Ubar = self.da.createGlobalVec()

        #
        # initiate pv inversion solver
        #
        self._wdinv = wdinversion(self)


    def set_wd(self, analytical_wd=True, file_wd=None):
        """ Set q to a given value
        """
        #
        if file_wd is not None:
            if self._verbose:
                print 'Set wind from file '+file_wd+' ...'
                print 'But not implemented yet !!!!'
            #read_nc_petsc(self.W, 'q', file_wd, self)
        elif analytical_wd:
            if self._verbose:
                print 'Set wind analytically '
            self.set_wd_analytically()

    def set_wd_analytically(self):
        """ Set wind analytically
        """
        w = self.da.getVecArray(self.W)
        mx, my, mz = self.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()
        #print zs, ze
        #
        for j in range(ys, ye):
            for i in range(xs, xe):
                # U component
                w[i, j, self.kx] = 1.e0*np.exp(-((i/float(mx-1)-0.5)**2
                                                 + (j/float(my-1)-0.5)**2)/0.1**2)
                #w[i, j, 0] *= np.sin(i/float(mx-1)*np.pi)
                #w[i, j, 0] *= np.sin(2*j/float(my-1)*np.pi)
                # V component
                w[i, j, self.ky] = 0.
                #w[i, j, 2] = 0.


    def solve_uv(self, domega=None):
        """ wrapper around solver solve method
        """
        self._wdinv.solve(domega, self.W, self.U, self.da)


    def set_U(self):
        """ Debug: set solution analytically
        """
        u = self.da.getVecArray(self.U)
        mx, my, mz = self.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()
        #
        for j in range(ys, ye):
            for i in range(xs, xe):
                # U component
                u[i, j, self.ky] = 1.e0
                # V component
                u[i, j, self.kx] = 0.

