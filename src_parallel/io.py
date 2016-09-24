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
                nc_Vx[i][:] = Vf[np.newaxis,ml.kx,...]
                nc_Vy[i][:] = Vf[np.newaxis,ml.ky,...]
            else:
                if i==0: it=nc_Vx[i].shape[0]
                nc_Vx[i][it,...] = Vf[ml.kx,:]
                nc_Vy[i][it, ...] = Vf[ml.ky, :]
        ml.comm.barrier()
      
    if rank == 0:
        # close the netcdf file
        rootgrp.close()
        



def read_nc():
    pass

def read_nc_petsc():
    pass

def read_hgrid_dimensions():
    pass