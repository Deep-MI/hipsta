"""
This module provides a function to create QC plots

"""

import os
import logging

import numpy as np
import plotly.graph_objects as go
    
from lapy import TriaMesh, io, plot as lpp
from plotly.subplots import make_subplots

from .getLevelsets import levelsetsTria

# ==============================================================================
# FUNCTIONS

def _sortLevelSets(LVL, dims, tol=1e-16):

    # create array of line segments
    tmpx = list()
    tmpy = list()

    for i in range(len(LVL[1][0])):
        tmpx.append((LVL[0][0][LVL[1][0][i][0] - 1][dims[0]], LVL[0][0][LVL[1][0][i][1] - 1][dims[0]]))
        tmpy.append((LVL[0][0][LVL[1][0][i][0] - 1][dims[1]], LVL[0][0][LVL[1][0][i][1] - 1][dims[1]]))

    tmpx = np.array(tmpx)
    tmpy = np.array(tmpy)

    # remove duplicate points
    tmpxy = np.unique(np.concatenate((tmpx, tmpy), axis=1), axis=0)
    tmpx = tmpxy[:, 0:2]
    tmpy = tmpxy[:, 2:4]

    # remove segments which are de-facto points
    tmpIdx = np.logical_or(np.abs(tmpx[:,0]-tmpx[:,1])>tol, np.abs(tmpy[:,0]-tmpy[:,1])>tol)
    tmpx = tmpx[tmpIdx, :]
    tmpy = tmpy[tmpIdx, :]

    # need to order array of line segments; whenever we encounter a
    # closed loop, we will already plot; otherwise, plot in the end
    sortIdx = np.array(range(0, len(tmpx)))

    tmpxSort = np.array(tmpx[sortIdx[0], ], ndmin=2)
    tmpySort = np.array(tmpy[sortIdx[0], ], ndmin=2)

    sortIdx = np.delete(sortIdx, sortIdx[0])

    while len(sortIdx) > 1:

        findIdx = np.array(np.where(np.logical_and(
            np.abs(tmpx[sortIdx, ] - tmpxSort[tmpxSort.shape[0]-1, 1]) < tol,
            np.abs(tmpy[sortIdx, ] - tmpySort[tmpySort.shape[0]-1, 1]) < tol)), ndmin=2).T

        # delete existing finds
        findIdxKeep = list()
        for k in range(findIdx.shape[0]):
            if not np.any(np.all(np.logical_or(tmpx[sortIdx[findIdx[k, 0]], 0] == tmpxSort, tmpx[sortIdx[findIdx[k, 0]], 1] == tmpxSort), axis=1)):
                findIdxKeep.append(k)
        findIdx = findIdx[findIdxKeep, ]

        if findIdx.shape[0] == 0:
            # reset (start new loop)
            tmpxSort = np.array(tmpx[sortIdx[0], ], ndmin=2)
            tmpySort = np.array(tmpy[sortIdx[0], ], ndmin=2)
            sortIdx = np.delete(sortIdx, 0)
        elif findIdx.shape[0] == 1:
            # add to current set
            if findIdx[0, 1] == 0:
                tmpxSort = np.append(tmpxSort, np.array(tmpx[sortIdx[findIdx[0, 0]], ::1], ndmin=2), axis=0)
                tmpySort = np.append(tmpySort, np.array(tmpy[sortIdx[findIdx[0, 0]], ::1], ndmin=2), axis=0)
            elif findIdx[0, 1] == 1:
                tmpxSort = np.append(tmpxSort, np.array(tmpx[sortIdx[findIdx[0, 0]], ::-1], ndmin=2), axis=0)
                tmpySort = np.append(tmpySort, np.array(tmpy[sortIdx[findIdx[0, 0]], ::-1], ndmin=2), axis=0)
            sortIdx = np.delete(sortIdx, findIdx[0, 0])
        elif findIdx.shape[0] > 1:
            # warning
            logging.warning("A problem occurred with the surface overlays")
            break

    return tmpxSort, tmpySort


def qcPlots(params, stage=None):

    # mesh
    if params.internal.noqc is False and stage=="mesh":

        triaMesh = TriaMesh.read_vtk(os.path.join(params.OUTDIR, params.HEMI + ".surf.vtk"))

        if params.HEMI =="lh":
            camera = dict(
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=-2, y=1, z=1))
        else:
            camera = dict(
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=2, y=-1, z=1))

        lpp.plot_tria_mesh(triaMesh, tcolor=[50,50,50], background_color="black", camera=camera, export_png=os.path.join(params.OUTDIR, 'qc', params.HEMI + '.mesh.png'), no_display=True, scale_png=0.5)

    # profile
    if params.internal.noqc is False and stage=="profile":

        triaMesh = TriaMesh.read_vtk(os.path.join(params.OUTDIR, 'tetra-cube', params.HEMI + ".rm.bnd.seam.rm.cut.vtk"))

        triaFunc = np.array(io.read_vfunc(os.path.join(params.OUTDIR, 'tetra-cube', params.HEMI + ".poisson1.rm.bnd.seam.rm.cut.psol")))

        #

        lvl2_0, lvl2i_0, lvl2j_0 = levelsetsTria(triaMesh.v, triaMesh.t, triaFunc, 0.05)
        lvl2_1, lvl2i_1, lvl2j_1 = levelsetsTria(triaMesh.v, triaMesh.t, triaFunc, 0.25)
        lvl2_2, lvl2i_2, lvl2j_2 = levelsetsTria(triaMesh.v, triaMesh.t, triaFunc, 0.5)
        lvl2_3, lvl2i_3, lvl2j_3 = levelsetsTria(triaMesh.v, triaMesh.t, triaFunc, 0.75)
        lvl2_4, lvl2i_4, lvl2j_4 = levelsetsTria(triaMesh.v, triaMesh.t, triaFunc, 0.95)

        #

        tmpxSort0, tmpySort0 = _sortLevelSets([lvl2_0, lvl2i_0, lvl2j_0], dims=[0, 2])
        tmpxSort1, tmpySort1 = _sortLevelSets([lvl2_1, lvl2i_1, lvl2j_1], dims=[0, 2])
        tmpxSort2, tmpySort2 = _sortLevelSets([lvl2_2, lvl2i_2, lvl2j_2], dims=[0, 2])
        tmpxSort3, tmpySort3 = _sortLevelSets([lvl2_3, lvl2i_3, lvl2j_3], dims=[0, 2])
        tmpxSort4, tmpySort4 = _sortLevelSets([lvl2_4, lvl2i_4, lvl2j_4], dims=[0, 2])

        #

        fig = make_subplots(rows=1, cols=5)

        fig.add_trace(go.Scatter(x=tmpxSort0[:,0], y=tmpySort0[:,0], mode="lines"), row=1, col=1)
        fig.add_trace(go.Scatter(x=tmpxSort1[:,0], y=tmpySort1[:,0], mode="lines"), row=1, col=2)
        fig.add_trace(go.Scatter(x=tmpxSort2[:,0], y=tmpySort2[:,0], mode="lines"), row=1, col=3)
        fig.add_trace(go.Scatter(x=tmpxSort3[:,0], y=tmpySort3[:,0], mode="lines"), row=1, col=4)
        fig.add_trace(go.Scatter(x=tmpxSort4[:,0], y=tmpySort4[:,0], mode="lines"), row=1, col=5)

        fig.update_layout(yaxis = dict(scaleanchor = 'x'),
                          yaxis2 = dict(scaleanchor = 'x'),
                          yaxis3 = dict(scaleanchor = 'x'),
                          yaxis4 = dict(scaleanchor = 'x'),
                          yaxis5 = dict(scaleanchor = 'x'))

        fig.write_image(os.path.join(params.OUTDIR, 'qc', params.HEMI + '.profile.png'))

    # hull
    if stage=="hull":

        triaMesh = TriaMesh.read_vtk(os.path.join(params.OUTDIR, "thickness", params.HEMI + ".hull.vtk"))

        if params.HEMI =="lh":
            camera = dict(
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=-2, y=1, z=1))
        else:
            camera = dict(
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=2, y=-1, z=1))

        lpp.plot_tria_mesh(triaMesh, tcolor=[50,50,50], background_color="black", camera=camera, export_png=os.path.join(params.OUTDIR, 'qc', params.HEMI + '.hull.png'), no_display=True, scale_png=0.5)

    # return
    return params
