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

def load_response(prefix):

    # for figure names
    #pref=prefix.replace('data/','')+'_'

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


    ### spectral distribution
    for i, lomega in enumerate(omega):
        lEu, k = get_hspectrum(ur[i,...],ui[i,...],x,y)
        lEv, k = get_hspectrum(vr[i,...],vi[i,...],x,y)
        if i==0:
            E = (lEu+lEv)[None,:]
        else:
            E = np.concatenate((E, (lEu+lEv)[None,:]), axis=0)

    return ke, omega, E, k




def get_hspectrum(ur,ui,x,y):

    # wavenumbers
    kx = np.fft.fftfreq(x.size, x[1]-x[0])
    ky = np.fft.fftfreq(y.size, y[1]-y[0])
    kkx, kky = np.meshgrid(kx,ky)
    kk = np.sqrt(kkx**2+kky**2)
    # target  one
    k=kx[1:x.size/2]
    dk=k[1]-k[0]
    Nk=k.size
    subs = np.ceil((kk-k[0]+dk/2)/dk)-1
    #print subs.dtype  #float64 weird, int with python3.0 apparently
    ip = np.where( (subs>=0) & (subs<Nk) & (kkx>=0) )

    # takes a zonal average off
    #ur -= ur.mean(axis=1, keepdims=True)
    #ui -= ui.mean(axis=1, keepdims=True)

    # computes the tapering window
    win = np.hanning(y.size)[:, None]*np.hanning(x.size)[None,:]
    # apply the taper
    ur = ur * win
    ui = ui * win

    # computes the fft
    fur = np.fft.fft2(ur) / x.size / y.size
    fui = np.fft.fft2(ui) / x.size / y.size

    # average the spectrum
    Eur = np.bincount(subs[ip].astype(int), np.abs(fur[ip]) ** 2 / dk)
    Eui = np.bincount(subs[ip].astype(int), np.abs(fui[ip]) ** 2 / dk)

    return Eur+Eui, k



if __name__ == "__main__":

    pref='all_'

    P=[]; colors = []

    # base case
    P.append('data/base')
    colors.append('k')

    # base case with r*10
    P.append('data/base_10r')
    colors.append('0.5')

    # base case with r/10
    P.append('data/base_0.1r')
    colors.append('0.8')

    # 1 eddy
    P.append('data/eddy')
    colors.append('orange')

    # multiple eddies
    P.append('data/meddies')
    colors.append('cadetblue')

    ke = []; cname=[]
    E = [];
    for prefix in P:
        lke, omega, lE, k = load_response(prefix)
        ke.append(lke)
        E.append(lE)
        cname.append(prefix.replace('data/',''))


    if True:

        #
        # sensitivity to r
        #
        plt.figure(figsize=(6, 4))

        # plt.loglog(omega,ke.reshape((ke.shape[0],ke.shape[1]*ke.shape[2])),color='0.7',lw=0.05)
        for i,lke in enumerate(ke[:3]):
            #plt.loglog(omega/f0n, lke.mean(axis=2).mean(axis=1), '--', color=colors[i], label=cname[i]+' mean')
            plt.loglog(omega/f0n, lke[:, lke.shape[1] / 2, lke.shape[2] / 2], color=colors[i], label=cname[i]+' center')

        # plt.legend(loc=0, frameon=False, fontsize=6)
        plt.xlabel('omega/f0 [1]')
        plt.ylabel('m^2/s^2 /(rad/s)')
        plt.ylim((1e-6, 1e0))
        # plt.axis([0, 4000, 0, 2.5*1e-3])
        plt.grid()
        plt.tight_layout()
        plt.title('KE')
        plt.legend(loc=2, frameon=False)

        plt.axvline(x=1., color='r', lw=0.1)

        plt.savefig('figs/' + pref + 'ke_r.pdf')

        #
        # sensitivity to eddies
        #
        plt.figure(figsize=(6, 4))

        # plt.loglog(omega,ke.reshape((ke.shape[0],ke.shape[1]*ke.shape[2])),color='0.7',lw=0.05)
        i=0; lke=ke[i]
        plt.loglog(omega/f0n, lke.mean(axis=2).mean(axis=1), '--', color=colors[i], label=cname[i] + ' mean')
        plt.loglog(omega/f0n, lke[:, lke.shape[1] / 2, lke.shape[2] / 2], color=colors[i], label=cname[i] + ' center')

        for lke, col, lcname in zip(ke[3:],colors[3:], cname[3:]):
            plt.loglog(omega/f0n, lke.mean(axis=2).mean(axis=1), '--', color=col, label=lcname + ' mean')
            plt.loglog(omega/f0n, lke[:, lke.shape[1] / 2, lke.shape[2] / 2], color=col, label=lcname + ' center')

        # plt.legend(loc=0, frameon=False, fontsize=6)
        plt.xlabel('omega/f0 [1]')
        plt.ylabel('m^2/s^2 /(rad/s)')
        plt.ylim((1e-6, 1e0))
        # plt.axis([0, 4000, 0, 2.5*1e-3])
        plt.grid()
        plt.tight_layout()
        plt.title('KE')
        plt.legend(loc=2, frameon=False)

        plt.axvline(x = 1., color='r', lw=0.1)

        plt.savefig('figs/' + pref + 'ke.pdf')

        #
        # sensitivity to eddies, at omega=f0
        #
        plt.figure(figsize=(6, 4))

        iomega = np.where(omega/f0n > 1.)[0][0]
        # plt.loglog(omega,ke.reshape((ke.shape[0],ke.shape[1]*ke.shape[2])),color='0.7',lw=0.05)
        i = 0;
        lke = ke[i]
        plt.loglog(k, E[i][iomega,:], '-', color=colors[i], label=cname[i])

        for lE, col, lcname in zip(E[3:], colors[3:], cname[3:]):
            plt.loglog(k, lE[iomega,:], '--', color=col, label=lcname)

        # plt.legend(loc=0, frameon=False, fontsize=6)
        plt.xlabel('k [cpm]')
        plt.ylabel('m^2/s^2 /(rad/s)/(rad/m)')
        plt.ylim((1e-8, 1e4))
        plt.xlim((1e-6, 1e-4))
        # plt.axis([0, 4000, 0, 2.5*1e-3])
        plt.grid()
        plt.tight_layout()
        #plt.title('KE')
        plt.legend(frameon=False)

        #plt.axvline(x=1., color='r', lw=0.1)

        plt.savefig('figs/' + pref + 'ke_f0.pdf')

        #
        # sensitivity to eddies, at omega=f0+beta0*500km
        #
        plt.figure(figsize=(6, 4))

        iomega = np.where(omega / f0n > (f0 + beta0 * 500.e3) / f0)[0][0]
        # plt.loglog(omega,ke.reshape((ke.shape[0],ke.shape[1]*ke.shape[2])),color='0.7',lw=0.05)
        i = 0;
        lke = ke[i]
        plt.loglog(k, E[i][iomega, :], '-', color=colors[i], label=cname[i])

        for lE, col, lcname in zip(E[3:], colors[3:], cname[3:]):
            plt.loglog(k, lE[iomega, :], '-', color=col, label=lcname)

        # plt.legend(loc=0, frameon=False, fontsize=6)
        plt.xlabel('k [cpm]')
        plt.ylabel('m^2/s^2 /(rad/s)/(cpm)')
        plt.ylim((1e-8, 1e4))
        plt.xlim((1e-6, 1e-4))
        # plt.axis([0, 4000, 0, 2.5*1e-3])
        plt.grid()
        plt.tight_layout()
        # plt.title('KE')
        plt.legend(frameon=False)

        # add dispersion relation for a uniform stratification
        N = 1e-3;
        H = 4.e3;
        for i in xrange(4):
            c_mode = N * H / (float(i + 1) * np.pi)
            k_mode = np.sqrt( (f0 + beta0 * 500.e3)**2 - f0**2 )/c_mode /2/np.pi
            plt.axvline(x=k_mode, color='r', lw=0.1)

        # plt.axvline(x=1., color='r', lw=0.1)

        plt.savefig('figs/' + pref + 'ke_f0off.pdf')



    if False:
        #
        # plot the reponse in (omega,k) space
        #
        lvls = 20
        lvls = np.arange(-10.,5.,1.)
        for lE, lcname in zip(E,cname):

            plt.figure(figsize=(6, 4))

            toplt = lE
            toplt = np.log10(toplt)
            plt.contourf(k,omega/f0n,toplt, lvls, cmap=plt.cm.viridis)
            cbar = plt.colorbar()
            plt.ylabel('omega/f0 [1]')
            plt.xlabel('k [cpm]')
            plt.xscale('log')
            plt.yscale('log')
            plt.title('KE [(m/s)^2/cps/cpm], '+lcname)

            # Coriolis
            plt.axhline(y=1. , color='r', lw=0.1)
            plt.axhline(y=(f0+beta0*500.e3)/f0 , color='r', lw=0.1)

            # add dispersion relation for a uniform stratification
            N=1e-3; H=4.e3;
            plt.axhline(y=N/f0 , color='r', lw=0.1)
            for i in xrange(4):
                c_mode = N*H/(float(i+1)*np.pi)
                omega_mode = np.sqrt(f0**2 + (2.*np.pi*k)**2*c_mode**2)
                plt.plot(k,omega_mode/f0, '--', color='r',lw=0.1)

            plt.savefig('figs/' + lcname + '_E.pdf')



