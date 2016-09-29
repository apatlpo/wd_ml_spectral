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
        # for complex support
        #rootgrp = Dataset(filename, 'w',
        #                  format='NETCDF4', clobber=True)

        # create dimensions
        rootgrp.createDimension('x', ml.grid.Nx)
        rootgrp.createDimension('y', ml.grid.Ny)
        rootgrp.createDimension('omega', None)
        
        # create variables
        dtype='f8'
        #complex128 = np.dtype([("real", np.float64), ("imag", np.float64)]) # create complex128 compound data type.
        #dtype = rootgrp.createCompoundType(complex128, "complex128")
        #
        nc_x = rootgrp.createVariable('x',dtype,('x'))
        nc_y = rootgrp.createVariable('y',dtype,('y'))
        nc_x[:], nc_y[:] = ml.grid.get_xy()
        # 2D variables but all are vectors
        nc_Vx_real=[]
        nc_Vx_imag=[]
        nc_Vy_real=[]
        nc_Vy_imag=[]
        for name in vname:
            nc_Vx_real.append(rootgrp.createVariable(name+'x_r',dtype,('omega','y','x',)))
            nc_Vx_imag.append(rootgrp.createVariable(name+'x_i',dtype,('omega','y','x',)))
            nc_Vy_real.append(rootgrp.createVariable(name + 'y_r', dtype, ('omega', 'y', 'x',)))
            nc_Vy_imag.append(rootgrp.createVariable(name + 'y_i', dtype, ('omega', 'y', 'x',)))

    elif rank == 0:
        ### open netcdf file
        #rootgrp = Dataset(filename, 'a',
        #                  format='NETCDF4_CLASSIC')
        rootgrp = Dataset(filename, 'a',
                          format='NETCDF4')
        # 2D variables but all are vectors
        nc_Vx_real=[]
        nc_Vx_imag=[]
        nc_Vy_real=[]
        nc_Vy_imag=[]
        for name in vname:
            nc_Vx_real.append(rootgrp.variables[name+'x_r'])
            nc_Vx_imag.append(rootgrp.variables[name+'x_i'])
            nc_Vy_real.append(rootgrp.variables[name+'y_r'])
            nc_Vy_imag.append(rootgrp.variables[name+'y_i'])


    # loop around variables now and store them
    Vn = ml.da.createNaturalVec()
    for i in xrange(Nv):    
        ml.da.globalToNatural(V[i], Vn)
        scatter, Vn0 = PETSc.Scatter.toZero(Vn)
        scatter.scatter(Vn, Vn0, False, PETSc.Scatter.Mode.FORWARD)
        if rank == 0:
            Vf = Vn0[...].reshape(ml.da.sizes[::-1], order='c')
            if create:
                nc_Vx_real[i][:] = Vf[np.newaxis,ml.kx,...].real
                nc_Vx_imag[i][:] = Vf[np.newaxis,ml.kx,...].imag
                nc_Vy_real[i][:] = Vf[np.newaxis,ml.ky,...].real
                nc_Vy_imag[i][:] = Vf[np.newaxis,ml.ky,...].imag
            else:
                if i==0: it=nc_Vx_real[i].shape[0]
                nc_Vx_real[i][it,...] = Vf[ml.kx,:].real
                nc_Vx_imag[i][it,...] = Vf[ml.kx,:].imag
                nc_Vy_real[i][it, ...] = Vf[ml.ky, :].real
                nc_Vy_imag[i][it, ...] = Vf[ml.ky, :].imag
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