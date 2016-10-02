#!/usr/bin/python
# -*- encoding: utf8 -*-

import sys
import numpy as np
#import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset


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
    plt.loglog(omega,ke_mean,'k', label='mean')
    plt.loglog(omega,ke[:,ke.shape[1]/2,ke.shape[2]/2],'b', label='center')


    #plt.legend(loc=0, frameon=False, fontsize=6)
    plt.xlabel('omega [cps]')
    plt.ylabel('m^2/s^2 /(rad/s)')
    plt.ylim((1e-6,1e-2))
    #plt.axis([0, 4000, 0, 2.5*1e-3])
    plt.grid()
    plt.tight_layout()
    plt.title('KE, '+pref)

    plt.axvline(x=7e-5/2./np.pi,color='r')


    plt.savefig('figs/'+pref+'ke.pdf')


    ### plot wind
    plt.figure(figsize=(4,6))

    #plt.contourf(x,y,Wy[:].squeeze(),20)
    #plt.colorbar()
    plt.plot(Wy[0,:,0],y,'k')
    plt.xlabel('[Pa]')
    plt.ylabel('y [km]')

    #plt.axis('equal')
    plt.title('Wind stress y dir, '+pref[:-1])

    plt.savefig('figs/'+pref+'tauy.pdf')


    ### plot ubar, vbar
    plt.figure(figsize=(4, 6))

    plt.contourf(x,y,ubar[:].squeeze(),20)
    plt.colorbar()
    plt.xlabel('x [km]')
    plt.ylabel('y [km]')
    # plt.axis('equal')
    plt.title('ubar, ' + pref[:-1])

    plt.savefig('figs/' + pref + 'ubar.pdf')




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





