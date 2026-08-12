"""
This module provides a function to create QC plots

"""

import logging
import os

import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np
import pyvista as pv
from lapy import TriaMesh, io

from .get_levelsets import levelsetsTria

# ==============================================================================
# LOGGING

LOGGER = logging.getLogger(__name__)

# ==============================================================================
# FUNCTIONS


def _sortLevelSets(LVL, dims, tol=1e-16):
    """
    Order unordered levelset line segments into continuous polylines for plotting.

    Parameters
    ----------
    LVL : list
        Levelset data as returned by :func:`levelsetsTria`, given as
        ``[vertices, segments, tria_indices]``, where ``vertices`` holds
        the interpolated levelset points and ``segments`` holds 1-based
        index pairs into ``vertices`` describing individual line
        segments.
    dims : list of int
        The two coordinate dimensions (0, 1, or 2) of the vertices to
        project onto for the 2D line plot.
    tol : float, default=1e-16
        Tolerance below which two coordinate values are considered
        identical, used to remove duplicate points and degenerate
        (zero-length) segments.

    Returns
    -------
    tmpxSort : numpy.ndarray
        Ordered array of segment start/end x-coordinates.
    tmpySort : numpy.ndarray
        Ordered array of segment start/end y-coordinates.

    """
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
    tmpIdx = np.logical_or(np.abs(tmpx[:, 0] - tmpx[:, 1]) > tol, np.abs(tmpy[:, 0] - tmpy[:, 1]) > tol)
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


def _rotateVector(v, axis, angle_deg):
    """
    Rotate a vector around an axis by a given angle, using Rodrigues' formula.

    Parameters
    ----------
    v : numpy.ndarray
        The 3D vector to rotate.
    axis : numpy.ndarray
        Unit-length 3D rotation axis.
    angle_deg : float
        Rotation angle in degrees.

    Returns
    -------
    numpy.ndarray
        The rotated 3D vector.

    """
    angle = np.radians(angle_deg)
    return v * np.cos(angle) + np.cross(axis, v) * np.sin(angle) + axis * np.dot(axis, v) * (1 - np.cos(angle))


def _renderTriaMesh(path, up, eye, out_path, scale=0.5, shape=(1, 4)):
    """
    Render a triangle mesh to a PNG file as a grid of rotated views.

    The `rows x cols` views are evenly spaced by rotating the camera
    around the superior/inferior (`up`) axis, giving QC screenshots of
    the mesh from all sides.

    Parameters
    ----------
    path : str
        Path to the mesh file (e.g. a ``.vtk`` file) to render.
    up : dict
        Superior/inferior axis, used both as the camera up vector and
        as the rotation axis for the views, given as
        ``{"x": .., "y": .., "z": ..}``.
    eye : dict
        Camera viewing direction (from the mesh center) for the first
        (0 degree) view, as ``{"x": .., "y": .., "z": ..}``. Only the
        direction is used; it is normalized and scaled by the mesh's
        bounding-box diagonal to place the camera at a suitable
        distance.
    out_path : str
        Destination path for the rendered PNG screenshot.
    scale : float, default=0.5
        Scale factor applied to a base 800x800 pixel window size per
        view to control output image resolution.
    shape : tuple of int, default=(1, 4)
        Number of ``(rows, cols)`` views to render, evenly spaced
        around a full rotation about the `up` axis.

    Returns
    -------
    None

    """
    m = TriaMesh.read_vtk(path)
    m.rm_free_vertices_()

    mesh = pv.PolyData.from_regular_faces(m.v, m.t)
    mesh = mesh.compute_normals(auto_orient_normals=True, consistent_normals=True)

    up_vec = np.array([up["x"], up["y"], up["z"]], dtype=float)
    up_vec /= np.linalg.norm(up_vec)
    eye_vec = np.array([eye["x"], eye["y"], eye["z"]], dtype=float)
    eye_vec /= np.linalg.norm(eye_vec)

    center = np.array(mesh.center)
    distance = mesh.length * 1.5

    rows, cols = shape
    plotter = pv.Plotter(
        shape=shape,
        off_screen=True,
        window_size=[int(800 * cols * scale), int(800 * rows * scale)],
    )

    angles = np.linspace(0, 360, rows * cols, endpoint=False)
    for i, angle in enumerate(angles):
        plotter.subplot(*divmod(i, cols))
        plotter.set_background("white")
        plotter.add_mesh(
            mesh,
            color=np.array([180, 180, 180]) / 255,
            show_edges=False,
            smooth_shading=True,
            ambient=0.3,
            diffuse=0.7,
            specular=0.05,
        )

        direction = _rotateVector(eye_vec, up_vec, angle)
        plotter.camera.focal_point = tuple(center)
        plotter.camera.position = tuple(center + direction * distance)
        plotter.camera.up = tuple(up_vec)

    plotter.screenshot(out_path)
    plotter.close()


def qcPlots(params, stage=None):
    """
    Create quality-control (QC) plots for a given processing stage.

    Depending on `stage`, this renders either a 3D screenshot of a
    triangle mesh (``"mesh"`` and ``"hull"`` stages) or a grid of 2D
    levelset profile plots (``"profile"`` stage), and writes the result
    as a PNG file to ``params.OUTDIR/qc``.

    Parameters
    ----------
    params : Namespace-like object
        Parameter object providing at least ``FILENAME``, ``OUTDIR``,
        ``HEMI``, and ``internal.no_qc`` attributes, as constructed
        elsewhere in the pipeline.
    stage : str, optional
        Which QC plot to create. One of ``"mesh"``, ``"profile"``, 
        `mid-surface`, or ``"hull"``. If ``None`` or unrecognized, no 
        plot is created.

    Returns
    -------
    params : Namespace-like object
        The unmodified `params` object, returned for chaining.

    """
    # get axes
    img = nb.load(params.FILENAME)
    ras2ras_tkr = img.header.get_vox2ras_tkr() @ img.header.get_ras2vox()
    ornts = nb.orientations.io_orientation(np.linalg.inv(ras2ras_tkr))
    scale_factor = 2
    tilt_factor = 1
    if ornts[0, 0] == 2:
        up = dict(x=ornts[0, 1], y=0, z=0)
    elif ornts[1, 0] == 2:
        up = dict(x=0, y=ornts[0, 1], z=0)
    elif ornts[2, 0] == 2:
        up = dict(x=0, y=0, z=ornts[0, 1])
    if ornts[0, 0] == 0:
        if ornts[1, 0] == 1:
            # SCA
            if params.HEMI == "lh":
                eye = dict(x=scale_factor * ornts[0, 1], y=tilt_factor * ornts[0, 1], z=tilt_factor * ornts[0, 1])
            else:
                eye = dict(x=-scale_factor * ornts[0, 1], y=-tilt_factor * ornts[0, 1], z=tilt_factor * ornts[0, 1])
        else:
            # SAC
            if params.HEMI == "lh":
                eye = dict(x=scale_factor * ornts[0, 1], y=tilt_factor * ornts[0, 1], z=tilt_factor * ornts[0, 1])
            else:
                eye = dict(x=-scale_factor * ornts[0, 1], y=tilt_factor * ornts[0, 1], z=-tilt_factor * ornts[0, 1])
    if ornts[1, 0] == 0:
        if ornts[0, 0] == 1:
            # CSA
            if params.HEMI == "lh":
                eye = dict(x=tilt_factor * ornts[0, 1], y=scale_factor * ornts[0, 1], z=tilt_factor * ornts[0, 1])
            else:
                eye = dict(x=-tilt_factor * ornts[0, 1], y=-scale_factor * ornts[0, 1], z=tilt_factor * ornts[0, 1])
        else:
            # ASC
            if params.HEMI == "lh":
                eye = dict(x=tilt_factor * ornts[0, 1], y=scale_factor * ornts[0, 1], z=tilt_factor * ornts[0, 1])
            else:
                eye = dict(x=tilt_factor * ornts[0, 1], y=-scale_factor * ornts[0, 1], z=-tilt_factor * ornts[0, 1])
    if ornts[2, 0] == 0:
        if ornts[1, 0] == 1:
            # ACS
            if params.HEMI == "lh":
                eye = dict(x=tilt_factor * ornts[0, 1], y=tilt_factor * ornts[0, 1], z=scale_factor * ornts[0, 1])
            else:
                eye = dict(x=tilt_factor * ornts[0, 1], y=-tilt_factor * ornts[0, 1], z=-scale_factor * ornts[0, 1])
        else:
            # CAS
            if params.HEMI == "lh":
                eye = dict(x=tilt_factor * ornts[0, 1], y=tilt_factor * ornts[0, 1], z=scale_factor * ornts[0, 1])
            else:
                eye = dict(x=-tilt_factor * ornts[0, 1], y=tilt_factor * ornts[0, 1], z=-scale_factor * ornts[0, 1])

    # mesh
    if params.internal.no_qc is False and stage == "mesh":
        try:
            _renderTriaMesh(
                os.path.join(params.OUTDIR, params.HEMI + ".surf.vtk"),
                up,
                eye,
                os.path.join(params.OUTDIR, "qc", params.HEMI + ".mesh.png"),
                scale=0.5,
            )
        except Exception as e:
            LOGGER.warning("Could not create QC plot for mesh stage: %s", e)

    # profile
    if params.internal.no_qc is False and stage == "profile":
        triaMesh = TriaMesh.read_vtk(os.path.join(params.OUTDIR, "tetra-cube", params.HEMI + ".rm.bnd.seam.rm.cut.vtk"))

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

        xmin = np.min(
            (
                np.min(tmpxSort0[:, 0]),
                np.min(tmpxSort1[:, 0]),
                np.min(tmpxSort2[:, 0]),
                np.min(tmpxSort3[:, 0]),
                np.min(tmpxSort4[:, 0]),
            )
        )
        xmax = np.max(
            (
                np.max(tmpxSort0[:, 0]),
                np.max(tmpxSort1[:, 0]),
                np.max(tmpxSort2[:, 0]),
                np.max(tmpxSort3[:, 0]),
                np.max(tmpxSort4[:, 0]),
            )
        )

        ymin = np.min(
            (
                np.min(tmpySort0[:, 0]),
                np.min(tmpySort1[:, 0]),
                np.min(tmpySort2[:, 0]),
                np.min(tmpySort3[:, 0]),
                np.min(tmpySort4[:, 0]),
            )
        )
        ymax = np.max(
            (
                np.max(tmpySort0[:, 0]),
                np.max(tmpySort1[:, 0]),
                np.max(tmpySort2[:, 0]),
                np.max(tmpySort3[:, 0]),
                np.max(tmpySort4[:, 0]),
            )
        )

        fig, axes = plt.subplots(1, 5)

        subplots = [
            (axes[0], tmpxSort0, tmpySort0),
            (axes[1], tmpxSort1, tmpySort1),
            (axes[2], tmpxSort2, tmpySort2),
            (axes[3], tmpxSort3, tmpySort3),
            (axes[4], tmpxSort4, tmpySort4),
        ]
        for ax, tmpxSort, tmpySort in subplots:
            ax.plot(tmpxSort[:, 0], tmpySort[:, 0])
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            ax.set_aspect("equal", adjustable="box")
            ax.axis("off")

        try:
            fig.savefig(os.path.join(params.OUTDIR, "qc", params.HEMI + ".profile.png"))
        except Exception as e:
            LOGGER.warning("Could not create QC plot for profile stage: %s", e)
        finally:
            plt.close(fig)

    # mid-surface
    if params.internal.no_qc is False and stage == "mid-surface":
        try:
            _renderTriaMesh(
                os.path.join(params.OUTDIR, "thickness", params.HEMI + ".mid-surface.vtk"),
                up,
                eye,
                os.path.join(params.OUTDIR, "qc", params.HEMI + ".mid-surface.png"),
                scale=0.5,
            )
        except Exception as e:
            LOGGER.warning("Could not create QC plot for mid-surface stage: %s", e)

    # hull
    if params.internal.no_qc is False and stage == "hull":
        try:
            _renderTriaMesh(
                os.path.join(params.OUTDIR, "thickness", params.HEMI + ".hull.vtk"),
                up,
                eye,
                os.path.join(params.OUTDIR, "qc", params.HEMI + ".hull.png"),
                scale=0.5,
            )
        except Exception as e:
            LOGGER.warning("Could not create QC plot for hull stage: %s", e)

    # return
    return params
