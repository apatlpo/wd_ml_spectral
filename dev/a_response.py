#!/usr/bin/python
# -*- encoding: utf8 -*-

import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset

mpl.rcParams['font.size']= 9

# !!! f0 should be read from outputs ...
f0_midlatitude = 2. * 2. * np.pi / 86400. * np.sin(45. * np.pi / 180.)
beta0_midlatitude = 2. * 2. * np.pi / 86400. * np.cos(45. * np.pi / 180.) / (6371. * 1e3)

#f0 = 7.e-5
f0 = f0_midlatitude
beta0 = beta0_midlatitude

# normalized version for plots
f0n = f0/ 2. / np.pi


""" Load output and plot response function
"""


def plot_response(prefix):

    # for figure names
    pref=prefix.replace('data/','')+'_'

    ###
    nc = Dataset(prefix+'_uv.nc', 'r')
    nc_ubar = Dataset(prefix+'_ubar.nc', 'r')
    nc_wind = Dataset(prefix+'_wind.nc', 'r')

    # get response
    ur = nc.variables['Ux_r']
    ui = nc.variables['Ux_i']
    vr = nc.variables['Uy_r']
    vi = nc.variables['Uy_i']

    # get background flow
    ubar = nc_ubar.variables['Ubarx_r']
    vbar = nc_ubar.variables['Ubary_r']

    # get wind
    Wx = nc_wind.variables['Wx_r']
    Wy = nc_wind.variables['Wy_r']

    #
    omega = nc.variables['omega'][:]/2./np.pi
    x = nc.variables['x'][:]
    y = nc.variables['y'][:]

    ### kinetic energy
    ke = (ur[:]**2+ui[:]**2+vr[:]**2+vi[:]**2)
    ke_mean = ke.mean(axis=2).mean(axis=1)

    #
    plt.figure(figsize=(6,4))
    #plt.ion()
    #plt.show()

    #print omega.shape
    #print ke.shape
    #plt.loglog(omega,ke.reshape((ke.shape[0],ke.shape[1]*ke.shape[2])),color='0.7',lw=0.05)
    plt.loglog(omega/f0n,ke_mean,'k', label='mean')
    plt.loglog(omega/f0n,ke[:,ke.shape[1]/2,ke.shape[2]/2],'b', label='center')


    #plt.legend(loc=0, frameon=False, fontsize=6)
    plt.xlabel('omega/f0 [1]')
    plt.ylabel('m^2/s^2 /(rad/s)')
    plt.ylim((1e-6,1e0))
    #plt.axis([0, 4000, 0, 2.5*1e-3])
    plt.grid()
    plt.tight_layout()
    plt.title('KE, '+pref[:-1])
    plt.legend(frameon=False)

    plt.axvline(x=1.,color='r', lw=0.1)

    plt.savefig('figs/'+pref+'ke.pdf')


    ### plot wind
    plt.figure(figsize=(5,6))

    #plt.contourf(x,y,Wy[:].squeeze(),20)
    #plt.colorbar()
    plt.plot(Wy[0,:,0]*1e3,y/1.e3,'k')
    plt.xlabel('[Pa]')
    plt.ylabel('y [km]')
    plt.grid(True)

    #plt.axis('equal')
    plt.title('Wind stress y dir, '+pref[:-1])

    plt.savefig('figs/'+pref+'tauy.pdf')


    ### plot ubar, vbar
    plt.figure(figsize=(10, 6))

    plt.subplot(1, 2, 1)

    plt.contourf(x/1e3,y/1e3,ubar[:].squeeze(),20, cmap=plt.cm.RdBu_r)
    plt.colorbar()
    plt.xlabel('x [km]')
    plt.ylabel('y [km]')
    plt.grid(True)
    plt.xlim((x[0]/1e3,x[-1]/1e3))
    plt.ylim((y[0]/1e3,y[-1]/1e3))
    plt.axis('equal')
    plt.title('ubar, ' + pref[:-1])

    plt.subplot(1, 2, 2)
    plt.contourf(x/1e3, y/1e3, vbar[:].squeeze(), 20, cmap=plt.cm.RdBu_r)
    plt.colorbar()
    plt.xlabel('x [km]')
    #plt.ylabel('y [km]')
    plt.grid(True)
    plt.xlim((x[0]/1e3,x[-1]/1e3))
    plt.ylim((y[0]/1e3,y[-1]/1e3))
    plt.axis('equal')
    plt.title('vbar, ' + pref[:-1])

    plt.savefig('figs/' + pref + 'uvbar.pdf')




if __name__ == "__main__":

    P=[]

    # base case
    P.append('data/base')

    # base case with r*10
    P.append('data/base_10r')

    # base case with r/10
    P.append('data/base_0.1r')

    # 1 eddy
    P.append('data/eddy')

    # multiple eddies
    P.append('data/meddies')


    for prefix in P:
        plot_response(prefix)

