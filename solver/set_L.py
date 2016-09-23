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
    (xs, xe), (ys, ye), (zs, ze) = ml.da.getRanges()
    #
    L.zeroEntries()
    row = PETSc.Mat.Stencil()
    col = PETSc.Mat.Stencil()
    #
    for j in range(ys, ye):
        for i in range(xs, xe):
            # U equation
            row.index = (i, j, 0)
            row.field = 0
            if (i==0    or j==0 or
                i==mx-1 or j==my-1):
                # U points
                L.setValueStencil(row, row, 1.0)
            else:
                # interior points
                for index, value in [
                    ((i,j-1,0), 0.), # U points
                    ((i-1,j,0), 0.),
                    ((i, j, 0), 0.),
                    ((i+1,j,0), 0.),
                    ((i,j+1,0), 0.),
                    ((i,j-1,1), 0.), # V points
                    ((i-1,j,1), 0.),
                    ((i, j, 1), -ml.f0),
                    ((i+1,j,1), 0.),
                    ((i,j+1,1), 0.)
                    ]:
                    col.index = index
                    col.field = 0
                    L.setValueStencil(row, col, value)
            #
            # V equation
            #
            row.index = (i, j, 1)
            row.field = 0
            if (i == 0 or j == 0 or
                i == mx - 1 or j == my - 1):
                # V points
                L.setValueStencil(row, row, 1.0)
            else:
                # interior points
                for index, value in [
                    ((i,j-1,0), 0.), # U points
                    ((i-1,j,0), 0.),
                    ((i, j, 0), ml.f0),
                    ((i+1,j,0), 0.),
                    ((i,j+1,0), 0.),
                    ((i,j-1,1), 0.), # V points
                    ((i-1,j,1), 0.),
                    ((i, j, 1), 0.),
                    ((i+1,j,1), 0.),
                    ((i,j+1,1), 0.)
                    ]:
                    col.index = index
                    col.field = 0
                    L.setValueStencil(row, col, value)
    L.assemble()
    return



