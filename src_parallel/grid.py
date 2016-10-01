#!/usr/bin/python
# -*- encoding: utf8 -*-

import numpy as np
from .io import read_nc, read_hgrid_dimensions
# for curvilinear grids
from netCDF4 import Dataset

class grid(object):
    """ Grid object
    """



#
#==================== Builders ============================================
# 
    
    #
    # object init
    #
    def __init__(self, hgrid = None):
        
        #
        # horizontal global grids
        #
        hgrid_uniform_default = {'Lx':3.e2*1.e3, 'Ly':2e2*1.e3, 
                                 'Nx':150, 'Ny':100}
        self._flag_hgrid_uniform = False
        if hgrid is None or isinstance(hgrid,dict):
            # uniform grid
            self._flag_hgrid_uniform = True            
            #
            hgrid_input = hgrid_uniform_default
            for key, value in hgrid.items():
                hgrid_input[key]=value
            #
            self._build_hgrid_uniform(**hgrid_input)
        else:
            # curvilinear grid
            self._build_hgrid_curvilinear(hgrid)


    #
    # Uniform grids
    #
    def _build_hgrid_uniform(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        # compute metric terms
        self.dx=self.Lx/(self.Nx-1.)
        self.dy=self.Ly/(self.Ny-1.)
        
    
    #
    # Curvilinear horizontal grid
    #
    def _build_hgrid_curvilinear(self, hgrid_file):
        # store metric file but metric terms are loaded later
        self.hgrid_file = hgrid_file
        # loads dimensions for dmda creation
        self.Nx, self.Ny = read_hgrid_dimensions(self.hgrid_file)
        
    def load_metric_terms(self, da, comm):
        # create a 3D vector containing metric terms
        self.D = da.createGlobalVec()
        # load curvilinear metric terms
        v = da.getVecArray(self.D)
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()
        # indexes along the third dimension of 
        self._k_dx =zs
        self._k_dy =zs+1
        self._k_lon=zs+2
        self._k_lat=zs+3       
        # open and read netcdf file
        rootgrp = Dataset(self.hgrid_file, 'r')
        for j in range(ys, ye):
            for i in range(xs, xe):
                v[i, j, self._k_dx] = rootgrp.variables['e1'][j,i]
                v[i, j, self._k_dy] = rootgrp.variables['e2'][j,i]
                v[i, j, self._k_lon] = rootgrp.variables['lon'][j,i]
                v[i, j, self._k_lat] = rootgrp.variables['lat'][j,i]
        rootgrp.close()
        #
        comm.barrier()
        pass
    
    def load_coriolis_parameter(self, coriolis_file, da, comm):
        v = da.getVecArray(self.D)
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()
        # indexes along the third dimension 
        self._k_f=zs+4       
        # open and read netcdf file
        rootgrp = Dataset(coriolis_file, 'r')
        for j in range(ys, ye):
            for i in range(xs, xe):
                v[i, j, self._k_f] = rootgrp.variables['f'][j,i]                
        rootgrp.close()
        #
        comm.barrier()
        pass



#
#==================== Grid information ============================================
#             
              
    def __str__(self):
        
        if self._flag_hgrid_uniform:
            out = 'The horizontal grid is uniform with:\n' \
                + '  Nx = %i , Ny = %i \n' % (self.Nx, self.Ny) \
                + '  Lx = %e km , Ly = %e km \n' % (self.Lx/1e3, self.Ly/1e3) \
                + '  dx = %e , dy = %e \n' % (self.dx, self.dy)
        else:
            # get stats about metric terms
            # not trivial to implement as min/max needs to be taken across tiles ...
            out = 'The horizontal grid is curvlinear with:\n' \
                + '  Nx = %i , Ny = %i \n' % (self.Nx, self.Ny) 
                #+ '  min(dx) = %e , mean(dx) = %e, max(dx) = %e \n' % (np.min(self.dx), np.mean(self.dx), np.max(self.dx)) \
                #+ '  min(dy) = %e , mean(dy) = %e, max(dy) = %e \n' % (np.min(self.dy), np.mean(self.dy), np.max(self.dy))
                
        return out
      
                  
#
#==================== extract grid data ============================================
#

                  
    def get_xy(self, i=None, j=None):
        if i==None and j==None:
            if self._flag_hgrid_uniform:
                x=np.linspace(0,self.Lx,self.Nx)
                y=np.linspace(0,self.Ly,self.Ny)
            else:
                x=np.arange(0., float(self.Nx))
                y=np.arange(0., float(self.Ny))
        else:
            x = i * self.dx
            y = j * self.dy
        return x,y


    #def get_xy_from_ij(self,i,j):
    #    x = self.dx * i
    #    y = self.dy * j
    #    return x, y


    def get_dist_from_center(self, i, j):
        dx = self.dx * (i-self.Nx/2)
        dy = self.dy * (j-self.Ny/2)
        r = np.sqrt(dx**2+dy**2)
        return r



