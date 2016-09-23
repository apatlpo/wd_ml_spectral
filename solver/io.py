#!/usr/bin/python
# -*- encoding: utf8 -*-

from petsc4py import PETSc

import numpy as np
from netCDF4 import Dataset


def write_nc(V, vname, filename, ml, create=True):
    """ Write a variable to a netcdf file
    Parameters:
        V list of petsc vectors
        vname list of corresponding names
        filename
        ml object
    """

    # number of variables to be stored
    Nv=len(vname)
    # process rank
    rank = ml.rank

    if rank == 0 and create:

        ### create a netcdf file to store QG pv for inversion
        rootgrp = Dataset(filename, 'w',
                          format='NETCDF4_CLASSIC', clobber=True)

        # create dimensions
        rootgrp.createDimension('x', ml.grid.Nx)
        rootgrp.createDimension('y', ml.grid.Ny)
        rootgrp.createDimension('t', None)
        
        # create variables
        dtype='f8'
        nc_x = rootgrp.createVariable('x',dtype,('x'))
        nc_y = rootgrp.createVariable('y',dtype,('y'))
        nc_x[:], nc_y[:] = ml.grid.get_xy()
        # 2D variables but all are vectors
        nc_Vx=[]
        nc_Vy=[]
        for name in vname:
            nc_Vx.append(rootgrp.createVariable(name+'x',dtype,('t','y','x',)))
            nc_Vy.append(rootgrp.createVariable(name + 'y', dtype, ('t', 'y', 'x',)))

    elif rank == 0:
        ### open netcdf file
        rootgrp = Dataset(filename, 'a',
                          format='NETCDF4_CLASSIC')
        # 2D variables but all are vectors
        nc_Vx=[]
        nc_Vy=[]
        for name in vname:
            nc_Vx.append(rootgrp.variables[name+'x'])
            nc_Vy.append(rootgrp.variables[name+'y'])


    # loop around variables now and store them
    Vn = ml.da.createNaturalVec()
    for i in xrange(Nv):    
        ml.da.globalToNatural(V[i], Vn)
        scatter, Vn0 = PETSc.Scatter.toZero(Vn)
        scatter.scatter(Vn, Vn0, False, PETSc.Scatter.Mode.FORWARD)
        if rank == 0:
            Vf = Vn0[...].reshape(ml.da.sizes[::-1], order='c')
            if create:
                nc_Vx[i][:] = Vf[np.newaxis,0,...]
                nc_Vy[i][:] = Vf[np.newaxis,1,...]
            else:
                if i==0: it=nc_Vx[i].shape[0]
                nc_Vx[i][it,...] = Vf[0,:]
                nc_Vy[i][it, ...] = Vf[1, :]
        ml.comm.barrier()
      
    if rank == 0:
        # close the netcdf file
        rootgrp.close()
        


def read_nc_petsc(V, vname, filename, ml):
    """ Read a variable from a netcdf file and stores it in a petsc Vector
    Parameters:
        V one(!) petsc vector
        vname corresponding name in netcdf file
        filename
        ml object
    """
    v = qg.da.getVecArray(V)
    (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()
    rootgrp = Dataset(filename, 'r')
    ndim=len(rootgrp.variables['q'].shape)
    if ndim>3:
        for k in range(zs, ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    #v[i, j, k] = rootgrp.variables['q'][-1,k,j,i]
                    # line above does not work for early versions of netcdf4 python library
                    # print netCDF4.__version__  1.1.1 has a bug and one cannot call -1 for last index:
                    # https://github.com/Unidata/netcdf4-python/issues/306
                    v[i, j, k] = rootgrp.variables['q'][rootgrp.variables['q'].shape[0]-1,k,j,i]
    else:
        for k in range(zs, ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    v[i, j, k] = rootgrp.variables['q'][k,j,i]
    rootgrp.close()
    qg.comm.barrier()
    #if qg.rank ==0: print 'Variable '+vname+' read from '+filename
    if qg._verbose: print '... done'


        
def read_nc(vnames, filename):
    """ Read variables from a netcdf file
    Parameters:
        vnames list of variable names
        filename
    """

    # open netcdf file
    rootgrp = Dataset(filename, 'r')
    
    # loop around variables to load
    if isinstance(vnames, list):
        V=[]
        for name in vnames:
            V.append(rootgrp.variables[name][:])
    else:
        V = rootgrp.variables[vnames][:]
        
    # close the netcdf file
    rootgrp.close()
    
    return V


def read_hgrid_dimensions(hgrid_file):
    """ Reads grid dimension from netcdf file
    Could put dimension names as optional inputs ...
    """
    # open netcdf file
    rootgrp = Dataset(hgrid_file, 'r')
    Nx = len(rootgrp.dimensions['x'])
    Ny = len(rootgrp.dimensions['y'])    
    return Nx, Ny
    


        