"""
This module provides a function to remove boundary tetras from tetrahedral meshes

"""

import os

import numpy as np
from lapy import TetMesh, io

# -----------------------------------------------------------------------------
# MAIN FUNCTION


def removeBoundaryMask(params):
    # -------------------------------------------------------------------------
    # message

    print()
    print("--------------------------------------------------------------------------------")
    print("Removing boundary tetras from mesh")
    print()

    # -------------------------------------------------------------------------
    # evaluate input

    VTKFile = os.path.join(params.OUTDIR, params.HEMI + ".tetra.vtk")
    PSOLFile = os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".tetra.psol")

    labelBndHead = params.LUTDICT["bndhead"]
    labelBndTail = params.LUTDICT["bndtail"]

    # -------------------------------------------------------------------------
    # load data

    tetMesh = TetMesh.read_vtk(VTKFile)

    lbl = io.read_vfunc(PSOLFile)
    lbl = np.array(lbl)

    vTail = tetMesh.v[lbl == labelBndTail,]
    vHead = tetMesh.v[lbl == labelBndHead,]

    # -------------------------------------------------------------------------
    # cutting surfaces based on point-cloud PCA

    vTaile, vTailc = np.linalg.eig(np.cov(vTail, rowvar=False))
    vHeade, vHeadc = np.linalg.eig(np.cov(vHead, rowvar=False))

    # need to order EVs

    vTailc = vTailc[:, np.flip(np.argsort(vTaile), axis=0)]
    vHeadc = vHeadc[:, np.flip(np.argsort(vHeade), axis=0)]
    vTaile = np.flip(np.sort(vTaile), axis=0)
    vHeade = np.flip(np.sort(vHeade), axis=0)

    # -------------------------------------------------------------------------
    # determine on which side of the plane a given point is

    # head

    # support vectors (move a little bit inwards, direction may differ); first,
    # determine if the richtungsvektor is pointing towards the center:

    if np.linalg.norm((np.mean(tetMesh.v, axis=0) - np.mean(vHead, axis=0))) > np.linalg.norm(
        np.mean(tetMesh.v, axis=0) - (np.mean(vHead, axis=0) + vHeadc[:, 2])
    ):  # yes
        sHead = np.mean(vHead, axis=0) + 1.0 * vHeadc[:, 2]
        dirHead = 1
    else:
        sHead = np.mean(vHead, axis=0) - 1.0 * vHeadc[:, 2]
        dirHead = 0

    # compute unit normal vector to plane

    uHead = vHeadc[:, 2] / np.linalg.norm(vHeadc[:, 2])

    # compute normals from each point to plane (could be v or vHead)

    dHead = np.zeros(tetMesh.v.shape[0])

    for i in range(0, tetMesh.v.shape[0]):
        # vector from point to support vector
        n = tetMesh.v[i, :] - sHead

        # distance from point to plane along normal
        dHead[i] = np.dot(n, uHead)

    # tail

    # support vectors (move a little bit inwards, direction may differ); first,
    # determine if the richtungsvektor is pointing towards the center:

    if np.linalg.norm((np.mean(tetMesh.v, axis=0) - np.mean(vTail, axis=0))) > np.linalg.norm(
        np.mean(tetMesh.v, axis=0) - (np.mean(vTail, axis=0) + vTailc[:, 2])
    ):
        sTail = np.mean(vTail, axis=0) + 1.0 * vTailc[:, 2]  # yes
        dirTail = 1
    else:
        sTail = np.mean(vTail, axis=0) - 1.0 * vTailc[:, 2]  # no
        dirTail = 0

    # compute unit normal vector to plane

    uTail = vTailc[:, 2] / np.linalg.norm(vTailc[:, 2])

    # compute normals from each point to plane (could be v or vTail)

    dTail = np.zeros(tetMesh.v.shape[0])

    for i in range(0, tetMesh.v.shape[0]):
        # vector from point to support vector
        n = tetMesh.v[i, :] - sTail

        # distance from point to plane along normal
        dTail[i] = np.dot(n, uTail)

    # -------------------------------------------------------------------------
    # remove triangles that contain 2 or 3 or 4 points

    if dirHead & dirTail:
        tcut = tetMesh.t[
            np.sum(np.reshape(np.in1d(tetMesh.t, np.where((dTail <= 0) | (dHead <= 0))), tetMesh.t.shape), axis=1) < 2,
            :,
        ]
        vcutIdxTail = np.where(dTail <= 0)[0]
        vcutIdxHead = np.where(dHead <= 0)[0]
    elif (not (dirHead)) & dirTail:
        tcut = tetMesh.t[
            np.sum(np.reshape(np.in1d(tetMesh.t, np.where((dTail <= 0) | (dHead > 0))), tetMesh.t.shape), axis=1) < 2, :
        ]
        vcutIdxTail = np.where(dTail <= 0)[0]
        vcutIdxHead = np.where(dHead > 0)[0]
    elif dirHead & (not (dirTail)):
        tcut = tetMesh.t[
            np.sum(np.reshape(np.in1d(tetMesh.t, np.where((dTail > 0) | (dHead <= 0))), tetMesh.t.shape), axis=1) < 2, :
        ]
        vcutIdxTail = np.where(dTail > 0)[0]
        vcutIdxHead = np.where(dHead <= 0)[0]
    elif (not (dirHead)) & (not (dirTail)):
        tcut = tetMesh.t[
            np.sum(np.reshape(np.in1d(tetMesh.t, np.where((dTail > 0) | (dHead > 0))), tetMesh.t.shape), axis=1) < 2, :
        ]
        vcutIdxTail = np.where(dTail > 0)[0]
        vcutIdxHead = np.where(dHead > 0)[0]

    # -------------------------------------------------------------------------
    # write PSOL

    vIdx = np.zeros(np.shape(tetMesh.v)[0])
    vIdx[vcutIdxHead] = labelBndTail
    vIdx[vcutIdxTail] = labelBndHead

    io.write_vfunc(os.path.join(params.OUTDIR, "tetra-cut", params.HEMI + ".tetra-remove_bnd.psol"), vIdx)

    # for visualization

    tetMeshBnd = tetMesh.boundary_tria()

    vRmBndKeep, vBndRmDel = tetMeshBnd.rm_free_vertices_()

    tetMeshBnd.orient_()

    tetMeshBnd.write_vtk(filename=os.path.join(params.OUTDIR, "tetra-cut", params.HEMI + ".rm.bnd.tetra.vtk"))

    io.write_vfunc(os.path.join(params.OUTDIR, "tetra-cut", params.HEMI + ".rm.bnd.tetra.psol"), vIdx)

    # -------------------------------------------------------------------------
    # return

    return params
