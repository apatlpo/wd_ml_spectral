#!/usr/bin/python
# -*- encoding: utf8 -*-

import numpy as np
from netCDF4 import Dataset


""" Write a complex variable to a netcdf file
"""



### create a netcdf file to store QG pv for inversion
#rootgrp = Dataset(filename, 'w',
#                  format='NETCDF4_CLASSIC', clobber=True)
# for complex support
rootgrp = Dataset('test.nc', 'w', format='NETCDF4', clobber=True)

# create dimensions
rootgrp.createDimension('x', 100)
rootgrp.createDimension('y', 200)
rootgrp.createDimension('t', None)

# create variables
#dtype='f8'
complex128 = np.dtype([("real", np.float64), ("imag", np.float64)]) # create complex128 compound data type.
dtype = rootgrp.createCompoundType(complex128, "complex128")
#
#nc_x = rootgrp.createVariable('x','f8',('x'))
#nc_y = rootgrp.createVariable('y',dtype,('y'))
#nc_x[:], nc_y[:] = ml.grid.get_xy()
# 2D variables but all are vectors
nc_V = rootgrp.createVariable('V',dtype,('t','y','x',))

V = np.random.randn(100, 200)*1j

Vd = np.empty((100,200),complex128)
Vd['real'] = V.real
Vd['imag'] = V.imag

#nc_V[:] = V[np.newaxis,...] # does not work
nc_V[:] = Vd[np.newaxis,...]


# close the netcdf file
rootgrp.close()


