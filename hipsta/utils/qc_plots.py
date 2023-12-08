"""
This module provides a function to create QC plots

"""

import logging
import os

import nibabel as nb
import numpy as np
import plotly.graph_objects as go
from lapy import TriaMesh, io
from lapy import plot as lpp
from plotly.subplots import make_subplots

from .get_levelsets import levelsetsTria

# ==============================================================================
# LOGGING

LOGGER = logging.getLogger(__name__)

# ==============================================================================
# FUNCTIONS


def _sortLevelSets(LVL, dims, tol=1e-16):
    # create array of line segments
    tmpx = list()
    tmpy = list()

    for i in range(len(LVL[1][0])):
        tmpx.append(
            (
                LVL[0][0][LVL[1][0][i][0] - 1][dims[0]],
                LVL[0][0][LVL[1][0][i][1] - 1][dims[0]],
            )
        )
        tmpy.append(
            (
                LVL[0][0][LVL[1][0][i][0] - 1][dims[1]],
                LVL[0][0][LVL[1][0][i][1] - 1][dims[1]],
            )
        )

    tmpx = np.array(tmpx)
    tmpy = np.array(tmpy)

    # remove duplicate points
    tmpxy = np.unique(np.concatenate((tmpx, tmpy), axis=1), axis=0)
    tmpx = tmpxy[:, 0:2]
    tmpy = tmpxy[:, 2:4]

    # remove segments which are de-facto points
    tmpIdx = np.logical_or(
        np.abs(tmpx[:, 0] - tmpx[:, 1]) > tol, np.abs(tmpy[:, 0] - tmpy[:, 1]) > tol
    )
    tmpx = tmpx[tmpIdx, :]
    tmpy = tmpy[tmpIdx, :]

    # need to order array of line segments; whenever we encounter a
    # closed loop, we will already plot; otherwise, plot in the end
    sortIdx = np.array(range(0, len(tmpx)))

    tmpxSort = np.array(tmpx[sortIdx[0],], ndmin=2)
    tmpySort = np.array(tmpy[sortIdx[0],], ndmin=2)

    sortIdx = np.delete(sortIdx, sortIdx[0])

    while len(sortIdx) > 1:
        findIdx = np.array(
            np.where(
                np.logical_and(
                    np.abs(tmpx[sortIdx,] - tmpxSort[tmpxSort.shape[0] - 1, 1]) < tol,
                    np.abs(tmpy[sortIdx,] - tmpySort[tmpySort.shape[0] - 1, 1]) < tol,
                )
            ),
            ndmin=2,
        ).T

        # delete existing finds
        findIdxKeep = list()
        for k in range(findIdx.shape[0]):
            if not np.any(
                np.all(
                    np.logical_or(
                        tmpx[sortIdx[findIdx[k, 0]], 0] == tmpxSort,
                        tmpx[sortIdx[findIdx[k, 0]], 1] == tmpxSort,
                    ),
                    axis=1,
                )
            ):
                findIdxKeep.append(k)
        findIdx = findIdx[findIdxKeep,]

        if findIdx.shape[0] == 0:
            # reset (start new loop)
            tmpxSort = np.array(tmpx[sortIdx[0],], ndmin=2)
            tmpySort = np.array(tmpy[sortIdx[0],], ndmin=2)
            sortIdx = np.delete(sortIdx, 0)
        elif findIdx.shape[0] == 1:
            # add to current set
            if findIdx[0, 1] == 0:
                tmpxSort = np.append(
                    tmpxSort,
                    np.array(tmpx[sortIdx[findIdx[0, 0]], ::1], ndmin=2),
                    axis=0,
                )
                tmpySort = np.append(
                    tmpySort,
                    np.array(tmpy[sortIdx[findIdx[0, 0]], ::1], ndmin=2),
                    axis=0,
                )
            elif findIdx[0, 1] == 1:
                tmpxSort = np.append(
                    tmpxSort,
                    np.array(tmpx[sortIdx[findIdx[0, 0]], ::-1], ndmin=2),
                    axis=0,
                )
                tmpySort = np.append(
                    tmpySort,
                    np.array(tmpy[sortIdx[findIdx[0, 0]], ::-1], ndmin=2),
                    axis=0,
                )
            sortIdx = np.delete(sortIdx, findIdx[0, 0])
        elif findIdx.shape[0] > 1:
            # warning
            LOGGER.warning("A problem occurred with the surface overlays")
            break

    return tmpxSort, tmpySort


def qcPlots(params, stage=None):
    # get axes
    img = nb.load(params.FILENAME)
    ras2ras_tkr = img.header.get_vox2ras_tkr() @ img.header.get_ras2vox()
    ornts = nb.orientations.io_orientation(np.linalg.inv(ras2ras_tkr))
    scale_factor = 2
    tilt_factor = 1
    if ornts[0, 0] == 2:
        up=dict(x=ornts[0, 1], y=0, z=0)
    elif ornts[1, 0] == 2:
        up=dict(x=0, y=ornts[0, 1], z=0)
    elif ornts[2, 0] == 2:
        up=dict(x=0, y=0, z=ornts[0, 1])
    if ornts[0, 0] == 0:
        if ornts[1, 0] == 1:
            # SCA
            if params.HEMI == "lh":
                eye = dict(x=scale_factor*ornts[0, 1], y=tilt_factor*ornts[0, 1], z=tilt_factor*ornts[0, 1])
            else:
                eye = dict(x=-scale_factor*ornts[0, 1], y=-tilt_factor*ornts[0, 1], z=tilt_factor*ornts[0, 1])
        else:
            # SAC
            if params.HEMI == "lh":
                eye = dict(x=scale_factor*ornts[0, 1], y=tilt_factor*ornts[0, 1], z=tilt_factor*ornts[0, 1])
            else:
                eye = dict(x=-scale_factor*ornts[0, 1], y=tilt_factor*ornts[0, 1], z=-tilt_factor*ornts[0, 1])
    if ornts[1, 0] == 0:
        if ornts[0, 0] == 1:
            # CSA
            if params.HEMI == "lh":
                eye = dict(x=tilt_factor*ornts[0, 1], y=scale_factor*ornts[0, 1], z=tilt_factor*ornts[0, 1])
            else:
                eye = dict(x=-tilt_factor*ornts[0, 1], y=-scale_factor*ornts[0, 1], z=tilt_factor*ornts[0, 1])
        else:
            # ASC
            if params.HEMI == "lh":
                eye = dict(x=tilt_factor*ornts[0, 1], y=scale_factor*ornts[0, 1], z=tilt_factor*ornts[0, 1])
            else:
                eye = dict(x=tilt_factor*ornts[0, 1], y=-scale_factor*ornts[0, 1], z=-tilt_factor*ornts[0, 1])
    if ornts[2, 0] == 0:
        if ornts[1, 0] == 1:
            # ACS
            if params.HEMI == "lh":
                eye = dict(x=tilt_factor*ornts[0, 1], y=tilt_factor*ornts[0, 1], z=scale_factor*ornts[0, 1])
            else:
                eye = dict(x=tilt_factor*ornts[0, 1], y=-tilt_factor*ornts[0, 1], z=-scale_factor*ornts[0, 1])
        else:
            # CAS
            if params.HEMI == "lh":
                eye = dict(x=tilt_factor*ornts[0, 1], y=tilt_factor*ornts[0, 1], z=scale_factor*ornts[0, 1])
            else:
                eye = dict(x=-tilt_factor*ornts[0, 1], y=tilt_factor*ornts[0, 1], z=-scale_factor*ornts[0, 1])

    # mesh
    if params.internal.no_qc is False and stage == "mesh":
        triaMesh = TriaMesh.read_vtk(
            os.path.join(params.OUTDIR, params.HEMI + ".surf.vtk")
        )

        camera = dict(
            up=up,
            center=dict(x=0, y=0, z=0),
            eye=eye,
            )
        lpp.plot_tria_mesh(
            triaMesh,
            tcolor=[25, 25, 25],
            background_color="white",
            camera=camera,
            export_png=os.path.join(params.OUTDIR, "qc", params.HEMI + ".mesh.png"),
            no_display=True,
            scale_png=0.5,
        )

    # profile
    if params.internal.no_qc is False and stage == "profile":
        triaMesh = TriaMesh.read_vtk(
            os.path.join(
                params.OUTDIR, "tetra-cube", params.HEMI + ".rm.bnd.seam.rm.cut.vtk"
            )
        )

        triaFunc = np.array(
            io.read_vfunc(
                os.path.join(
                    params.OUTDIR,
                    "tetra-cube",
                    params.HEMI + ".poisson1.rm.bnd.seam.rm.cut.psol",
                )
            )
        )

        #

        lvl2_0, lvl2i_0, lvl2j_0 = levelsetsTria(triaMesh.v, triaMesh.t, triaFunc, 0.10)
        lvl2_1, lvl2i_1, lvl2j_1 = levelsetsTria(triaMesh.v, triaMesh.t, triaFunc, 0.25)
        lvl2_2, lvl2i_2, lvl2j_2 = levelsetsTria(triaMesh.v, triaMesh.t, triaFunc, 0.50)
        lvl2_3, lvl2i_3, lvl2j_3 = levelsetsTria(triaMesh.v, triaMesh.t, triaFunc, 0.75)
        lvl2_4, lvl2i_4, lvl2j_4 = levelsetsTria(triaMesh.v, triaMesh.t, triaFunc, 0.90)

        #

        tmpxSort0, tmpySort0 = _sortLevelSets([lvl2_0, lvl2i_0, lvl2j_0], dims=[0, 2])
        tmpxSort1, tmpySort1 = _sortLevelSets([lvl2_1, lvl2i_1, lvl2j_1], dims=[0, 2])
        tmpxSort2, tmpySort2 = _sortLevelSets([lvl2_2, lvl2i_2, lvl2j_2], dims=[0, 2])
        tmpxSort3, tmpySort3 = _sortLevelSets([lvl2_3, lvl2i_3, lvl2j_3], dims=[0, 2])
        tmpxSort4, tmpySort4 = _sortLevelSets([lvl2_4, lvl2i_4, lvl2j_4], dims=[0, 2])

        #

        fig = make_subplots(rows=2, cols=3)

        fig.add_trace(
            go.Scatter(x=tmpxSort0[:, 0], y=tmpySort0[:, 0], mode="lines"), row=1, col=1
        )
        fig.add_trace(
            go.Scatter(x=tmpxSort1[:, 0], y=tmpySort1[:, 0], mode="lines"), row=1, col=2
        )
        fig.add_trace(
            go.Scatter(x=tmpxSort2[:, 0], y=tmpySort2[:, 0], mode="lines"), row=1, col=3
        )
        fig.add_trace(
            go.Scatter(x=tmpxSort3[:, 0], y=tmpySort3[:, 0], mode="lines"), row=2, col=1
        )
        fig.add_trace(
            go.Scatter(x=tmpxSort4[:, 0], y=tmpySort4[:, 0], mode="lines"), row=2, col=2
        )

        xmin = np.min((np.min(tmpxSort0[:, 0]), np.min(tmpxSort1[:, 0]), np.min(tmpxSort2[:, 0]), np.min(tmpxSort3[:, 0]), np.min(tmpxSort4[:, 0])))
        xmax = np.max((np.max(tmpxSort0[:, 0]), np.max(tmpxSort1[:, 0]), np.max(tmpxSort2[:, 0]), np.max(tmpxSort3[:, 0]), np.max(tmpxSort4[:, 0])))

        ymin = np.min((np.min(tmpySort0[:, 0]), np.min(tmpySort1[:, 0]), np.min(tmpySort2[:, 0]), np.min(tmpySort3[:, 0]), np.min(tmpySort4[:, 0])))
        ymax = np.max((np.max(tmpySort0[:, 0]), np.max(tmpySort1[:, 0]), np.max(tmpySort2[:, 0]), np.max(tmpySort3[:, 0]), np.max(tmpySort4[:, 0])))

        fig.update_layout(
            xaxis=dict(range=[xmin, xmax]),
            xaxis2=dict(range=[xmin, xmax]),
            xaxis3=dict(range=[xmin, xmax]),
            xaxis4=dict(range=[xmin, xmax]),
            xaxis5=dict(range=[xmin, xmax]),
            yaxis=dict(scaleanchor="x", range=[ymin, ymax]),
            yaxis2=dict(scaleanchor="x", range=[ymin, ymax]),
            yaxis3=dict(scaleanchor="x", range=[ymin, ymax]),
            yaxis4=dict(scaleanchor="x", range=[ymin, ymax]),
            yaxis5=dict(scaleanchor="x", range=[ymin, ymax]),
        )

        fig.update_layout(showlegend=False)

        fig.write_image(os.path.join(params.OUTDIR, "qc", params.HEMI + ".profile.png"))

    # hull
    if stage == "hull":
        triaMesh = TriaMesh.read_vtk(
            os.path.join(params.OUTDIR, "thickness", params.HEMI + ".hull.vtk")
        )

        camera = dict(
            up=up,
            center=dict(x=0, y=0, z=0),
            eye=eye,
        )
        lpp.plot_tria_mesh(
            triaMesh,
            tcolor=[25, 25, 25],
            background_color="white",
            camera=camera,
            export_png=os.path.join(params.OUTDIR, "qc", params.HEMI + ".hull.png"),
            no_display=True,
            scale_png=0.5,
        )

    # return
    return params
