#!/usr/bin/python
# -*- encoding: utf8 -*-

from petsc4py import PETSc

def set_L(L, ml):
    """ Builds the operator along with boundary conditions
        Horizontally uniform grid:
        ...
    """
    #
    mx, my, mz = ml.da.getSizes()
    dx, dy = ml.grid.dx, ml.grid.dy
    idx, idy = [1.0/dl for dl in [dx, dy]]
    #
    kx=ml.kx
    ky=ml.ky
    #
    (xs, xe), (ys, ye), (zs, ze) = ml.da.getRanges()
    #
    ### declare local vectors
    local_Ubar = ml.da.createLocalVec()
    ml.da.globalToLocal(ml.Ubar, local_Ubar)
    ubar = ml.da.getVecArray(local_Ubar)
    #
    L.zeroEntries()
    row = PETSc.Mat.Stencil()
    col = PETSc.Mat.Stencil()
    #
    for j in range(ys, ye):
        for i in range(xs, xe):
            #
            # U equation
            #
            row.index = (i, j, kx)
            row.field = 0
            #if (i==0    or j==0 or
            #    i==mx-1 or j==my-1):
            if False:
                # U points
                L.setValueStencil(row, row, 1.0)
            else: # interior points
                for index, value in [
                    ((i,j-1,kx), -ubar[i,j,ky] *idy *.5 ), # U points
                    ((i-1,j,kx), -ubar[i,j,kx] *idx *.5 ),
                    ((i, j, kx),  ml.r + (ubar[i+1,j, kx] - ubar[i-1,j,kx]) * idx * .5 ),
                    ((i+1,j,kx),  ubar[i,j,kx] *idx *.5),
                    ((i,j+1,kx),  ubar[i,j,ky] *idy *.5 ),
                    ((i,j-1,ky), 0.), # V points
                    ((i-1,j,ky), 0.),
                    ((i, j, ky), -ml.f0 \
                                + (ubar[i,j+1, kx] - ubar[i,j-1,kx]) * idy * .5),
                    ((i+1,j,ky), 0.),
                    ((i,j+1,ky), 0.)
                    ]:
                    col.index = index
                    col.field = 0
                    L.setValueStencil(row, col, value)
            #
            # V equation
            #
            row.index = (i, j, ky)
            row.field = 0
            #if (i == 0 or j == 0 or
            #    i == mx - 1 or j == my - 1):
            if False:
                # V points
                L.setValueStencil(row, row, 1.0)
            else: # interior points
                for index, value in [
                    ((i,j-1,kx), 0.), # U points
                    ((i-1,j,kx), 0.),
                    ((i, j, kx), ml.f0 + \
                                  (ubar[i+1, j, ky] - ubar[i-1, j, ky]) * idx * .5),
                    ((i+1,j,kx), 0.),
                    ((i,j+1,kx), 0.),
                    ((i,j-1,ky), -ubar[i,j,ky] *idy *.5), # V points
                    ((i-1,j,ky), -ubar[i,j,kx] *idx *.5),
                    ((i, j, ky), ml.r + (ubar[i,j+1, ky] - ubar[i+1,j,ky]) * idy * .5),
                    ((i+1,j,ky),  ubar[i,j,kx] *idx *.5),
                    ((i,j+1,ky),  ubar[i,j,ky] *idy *.5)
                    ]:
                    col.index = index
                    col.field = 0
                    L.setValueStencil(row, col, value)
    L.assemble()
    return



