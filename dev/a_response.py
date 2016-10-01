#!/usr/bin/python
# -*- encoding: utf8 -*-

import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset


""" Load output and plot response function
"""

###
nc = Dataset('output.nc', 'r')

#
ur = nc.variables['Ux_r']
ui = nc.variables['Ux_i']
vr = nc.variables['Uy_r']
vi = nc.variables['Uy_i']

#
omega = nc.variables['omega'][:]/2./np.pi
x = nc.variables['x'][:]
y = nc.variables['y'][:]

# kinetic energy
ke = (ur[:]**2+ui[:]**2+vr[:]**2+vi[:]**2).mean(axis=2).mean(axis=1)

#
plt.figure(figsize=(6,4))
plt.ion()
plt.show()

plt.loglog(omega,ke,'k')

#plt.legend(loc=0, frameon=False, fontsize=6)
plt.xlabel('omega [cps]')
plt.ylabel('m^2/s^2 /(rad/s)')
#plt.axis([0, 4000, 0, 2.5*1e-3])
plt.grid()
plt.tight_layout()
plt.savefig("figs/ke.pdf")

sys.exit()


