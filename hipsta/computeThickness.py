"""
This module provides a function to compute thickness

"""

import os
import logging

import numpy as np
import pandas as pd
import nibabel as nb

from lapy import TriaMesh, TetMesh, io
from scipy.sparse import csgraph as sc

from .utils.getLevelsets import levelsetsTetra

# ------------------------------------------------------------------------------
# MAIN FUNCTION
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# computeThickness

def computeThickness(params):

    # message

    print()
    print("--------------------------------------------------------------------------------")
    print("Computing thickness")
    print()

    # get data

    HEMI = params.HEMI
    VERBOSE = not(params.internal.CLEANUP)
    IN_MESH = os.path.join(params.OUTDIR, 'tetra-cube', HEMI + ".seam.rm.cut.vtk")
    IN_FUNC = os.path.join(params.OUTDIR, 'tetra-cube', HEMI + ".uvw.seam.rm.cut.vtk")
    OUT_DIR = os.path.join(params.OUTDIR, 'thickness')

    paramsTHXn = params.internal.THXn
    paramsTHXp = params.internal.THXp
    paramsTHXk = params.internal.THXk

    paramsTHYn = params.internal.THYn
    paramsTHYp = params.internal.THYp
    paramsTHYk = params.internal.THYk

    paramsTHZn = params.internal.THZn
    paramsTHZp = params.internal.THZp
    paramsTHZk = params.internal.THZk

    allowRagged = params.internal.allowRagged
    allowRaggedTriangles = params.internal.allowRaggedTriangles
    skipOrient = params.internal.skipOrient

    # load mesh

    tetMesh = TetMesh.read_vtk(IN_MESH)

    v4 = tetMesh.v
    t4 = tetMesh.t

    # load function

    tetMeshUVW = TetMesh.read_vtk(IN_FUNC)

    p4 = tetMeshUVW.v
    q4 = tetMeshUVW.t

    # determine levels

    lLVL4x = np.linspace(paramsTHXn, paramsTHXp, paramsTHXk)
    lLVL4y = np.linspace(paramsTHYn, paramsTHYp, paramsTHYk)
    lLVL4z = np.linspace(paramsTHZn, paramsTHZp, paramsTHZk)

    # compute levelsets in original space

    vLVL4x, tLVL4x, iLVL4x, jLVL4x, oLVL4x = levelsetsTetra(v4, t4, p4[:,0], lLVL4x)
    vLVL4y, tLVL4y, iLVL4y, jLVL4y, oLVL4y = levelsetsTetra(v4, t4, p4[:,1], lLVL4y)
    vLVL4z, tLVL4z, iLVL4z, jLVL4z, oLVL4z = levelsetsTetra(v4, t4, p4[:,2], lLVL4z)

    # compute levelsets in parameter space

    vLVL4px, tLVL4px, iLVL4px, jLVL4px, oLVL4px = levelsetsTetra(p4, t4, p4[:,0], lLVL4x)
    vLVL4py, tLVL4py, iLVL4py, jLVL4py, oLVL4py = levelsetsTetra(p4, t4, p4[:,1], lLVL4y)
    vLVL4pz, tLVL4pz, iLVL4pz, jLVL4pz, oLVL4pz = levelsetsTetra(p4, t4, p4[:,2], lLVL4z)

    # write out single slices in x, y, z directions

    if VERBOSE is True:

        for i in range(0,len(vLVL4px)):
            TriaMesh.write_vtk(TriaMesh(vLVL4x[i], tLVL4x[i]-1), os.path.join(OUT_DIR, HEMI + '.vlx.lvl'+str(i)+'.vtk'))
            TriaMesh.write_vtk(TriaMesh(vLVL4px[i], tLVL4px[i]-1), os.path.join(OUT_DIR, HEMI + '.uvw.vlx.lvl'+str(i)+'.vtk'))

        for i in range(0,len(vLVL4py)):
            TriaMesh.write_vtk(TriaMesh(vLVL4y[i], tLVL4y[i]-1), os.path.join(OUT_DIR, HEMI + '.vly.lvl'+str(i)+'.vtk'))
            TriaMesh.write_vtk(TriaMesh(vLVL4py[i], tLVL4py[i]-1), os.path.join(OUT_DIR, HEMI + '.uvw.vly.lvl'+str(i)+'.vtk'))

        for i in range(0,len(vLVL4pz)):
            TriaMesh.write_vtk(TriaMesh(vLVL4z[i], tLVL4z[i]-1), os.path.join(OUT_DIR, HEMI + '.vlz.lvl'+str(i)+'.vtk'))
            TriaMesh.write_vtk(TriaMesh(vLVL4pz[i], tLVL4pz[i]-1), os.path.join(OUT_DIR, HEMI + '.uvw.vlz.lvl'+str(i)+'.vtk'))

    # now loop across levelsets

    origV4 = np.empty((len(vLVL4px),len(vLVL4py),len(vLVL4pz),3))
    origV4[:] = np.nan

    origV4flat = list()

    ctr = np.empty((len(vLVL4px),len(vLVL4py),len(vLVL4pz)))
    ctr[:] = np.nan

    ctrflat = list()

    msgNoCmTet = 0
    msgCmTetNoInt = 0
    msgCmTetMnInt = 0

    for lx in range(0,len(vLVL4px)):

        for ly in range(0,len(vLVL4py)):

            for lz in range(0,len(vLVL4pz)):

                # which tetras are present in all three functions? note that
                # this can be more than one tetra, and we do not know the
                # number in advance.

                # variant based on parameter levelsets:
                tmp = np.intersect1d(iLVL4px[lx],np.intersect1d(iLVL4py[ly],iLVL4pz[lz]))

                if len(tmp) == 0:

                    msgNoCmTet += 1

                    origV4[lx,ly,lz,:] = [np.nan, np.nan, np.nan]
                    ctr[lx,ly,lz] = np.nan

                else:

                    # if using p4, which tetra holds point p4(lx,ly,lz)?

                    # create a barycentric coordinate system with the first
                    # vertex as the basis; point is inside the tetra if all
                    # barycentric oordinates are within the range [0,1]
                    # (sign matters!)

                    # then map from barycentric ctr back to xyz in v4 space

                    # Tp is the transformation matrix for the parameter space
                    # Tv is the transformation matrix for the vertex space

                    # b is the vector of barycentric coordinates

                    # ctr is the tetra within which the point is located

                    Tp = np.empty((3,3,len(tmp)))
                    Tp[:] = np.nan

                    b = np.empty((3,len(tmp)))
                    b[:] = np.nan

                    for i in range(0,len(tmp)):

                        Tp[:,:,i] = np.array([
                            p4[t4[tmp[i],1],:] - p4[t4[tmp[i], 0], :],
                            p4[t4[tmp[i],2],:] - p4[t4[tmp[i], 0], :],
                            p4[t4[tmp[i],3],:] - p4[t4[tmp[i], 0], :]
                            ]).transpose()

                        if np.linalg.matrix_rank(Tp[:,:,i]) == 3:

                            b[:,i] = np.matmul(
                                np.linalg.inv(Tp[:,:,i]),
                                np.array([[lLVL4x[lx],lLVL4y[ly],lLVL4z[lz]]-p4[t4[tmp[i],0],:]]).transpose()
                                ).squeeze()

                        else:

                            b[:,i] = np.nan

                    b = np.concatenate((b, (1 - np.sum(b, axis=0)).reshape((1,b.shape[1]))), axis=0)

                    tmpi = np.where(np.sum((b>=0)&(b<=1),axis=0)==4)[0]

                    if len(tmpi) == 0:
                        msgCmTetNoInt += 1
                    elif len(tmpi) > 1:
                        msgCmTetMnInt += 1
                    else:

                        Tv = np.array([
                            v4[t4[tmp[tmpi[0]],1],:] - v4[t4[tmp[tmpi[0]], 0], :],
                            v4[t4[tmp[tmpi[0]],2],:] - v4[t4[tmp[tmpi[0]], 0], :],
                            v4[t4[tmp[tmpi[0]],3],:] - v4[t4[tmp[tmpi[0]], 0], :]
                            ]).transpose()

                        origV4[lx,ly,lz,:] = np.matmul(Tv, b[0:3,tmpi[0]]) + v4[t4[tmp[tmpi[0]],0],:]

                        ctr[lx,ly,lz] = tmp[tmpi[0]]

                if any(np.isnan(origV4[lx,ly,lz,:])) == False:
                    origV4flat.append(np.hstack((np.array([lx,ly,lz]), origV4[lx,ly,lz,:])))

                if np.isnan(ctr[lx,ly,lz]) == False:
                    ctrflat.append(np.hstack((np.array([lx,ly,lz]), ctr[lx,ly,lz])))

    origV4flat = np.array(origV4flat)
    ctrflat = np.array(ctr.flat)

    logging.info('No common tetra found: ' + str(msgNoCmTet))
    logging.info('Common tetra, but no intersection: ' + str(msgCmTetNoInt))
    logging.info('Common tetra, but too many intersections: ' + str(msgCmTetMnInt))

    # -------------------------------------------------------------------------
    # lines

    # x->yz
    llx = [ [ list() for k in range(np.shape(origV4)[2]) ] for z in range(np.shape(origV4)[1]) ]
    for iy in range(0,np.shape(origV4)[1]):
        for iz in range(0,np.shape(origV4)[2]):
            llx[iy][iz] = origV4flat[np.logical_and(origV4flat[:, 1] == iy, origV4flat[:, 2] == iz), 0:6]

    # y->xz
    lly = [ [ list() for k in range(np.shape(origV4)[2]) ] for z in range(np.shape(origV4)[0]) ]
    for ix in range(0,np.shape(origV4)[0]):
        for iz in range(0,np.shape(origV4)[2]):
            lly[ix][iz] = origV4flat[np.logical_and(origV4flat[:, 0] == ix, origV4flat[:, 2] == iz), 0:6]

    # z->xy
    llz = [ [ list() for k in range(np.shape(origV4)[1]) ] for z in range(np.shape(origV4)[0]) ]
    for ix in range(0,np.shape(origV4)[0]):
        for iy in range(0,np.shape(origV4)[1]):
            llz[ix][iy] = origV4flat[np.logical_and(origV4flat[:, 0] == ix, origV4flat[:, 1] == iy), 0:6]

    # --------------------------------------------------------------------------
    # create fake streamlines (will also work with non-continuous lines)

    # x

    vlx = list()
    tlx = list()
    nlx = list()
    llxLgth = np.full((np.shape(origV4)[1], np.shape(origV4)[2]), np.nan)

    for ii in range(len(llx)):

        for ij in range(len(llx[ii])):

            if len(llx[ii][ij])>1:

                tmp = list()

                for ik in range(len(llx[ii][ij]) - 1):

                    if llx[ii][ij][ik, 0] + 1 == llx[ii][ij][ik + 1, 0]:

                        vlx.append(llx[ii][ij][ik, 3:6])

                        tlx.append([len(vlx), len(vlx), len(vlx) + 1])

                        nlx.append([ ii, ij, ik, np.linalg.norm(llx[ii][ij][ik + 1, 3:6] - llx[ii][ij][ik, 3:6]) ])

                        tmp.append(np.linalg.norm(llx[ii][ij][ik + 1, 3:6] - llx[ii][ij][ik, 3:6]))

                    elif llx[ii][ij][ik - 1, 0] + 1 == llx[ii][ij][ik, 0]:

                        # close a potentially open previous line
                        vlx.append(llx[ii][ij][ik, 3:6])

                if llx[ii][ij][ik + 1, 0] - 1 == llx[ii][ij][ik, 0]:

                    vlx.append(llx[ii][ij][ik + 1, 3:6])

                if np.sum(tmp): # empty lists are false
                    llxLgth[ii, ij] = np.sum(tmp)

    TriaMesh.write_vtk(TriaMesh(np.array(vlx), np.asarray(tlx)-1), os.path.join(OUT_DIR, HEMI + '.grid-lines-x.vtk'))

    dfx = pd.DataFrame(llxLgth)

    dz = zip(range(0, dfx.shape[1]), [ "z"+str(z) for z in range(0, dfx.shape[1]) ])
    dfx = dfx.rename(columns=dict(dz))

    dy = zip(range(0, dfx.shape[0]), [ "y"+str(z) for z in range(0, dfx.shape[0]) ])
    dfx = dfx.rename(index=dict(dy))

    dfx.to_csv(os.path.join(OUT_DIR, HEMI + '.grid-segments-x.csv'))

    # y

    vly = list()
    tly = list()
    nly = list()
    llyLgth = np.full((np.shape(origV4)[0], np.shape(origV4)[2]), np.nan)

    for ii in range(len(lly)):

        for ij in range(len(lly[ii])):

            if len(lly[ii][ij])>1:

                tmp = list()

                for ik in range(len(lly[ii][ij]) - 1):

                    if lly[ii][ij][ik, 1] + 1 == lly[ii][ij][ik + 1, 1]:

                        vly.append(lly[ii][ij][ik, 3:6])

                        tly.append([len(vly), len(vly), len(vly) + 1])

                        nly.append([ ii, ij, ik, np.linalg.norm(lly[ii][ij][ik + 1, 3:6] - lly[ii][ij][ik, 3:6]) ])

                        tmp.append(np.linalg.norm(lly[ii][ij][ik + 1, 3:6] - lly[ii][ij][ik, 3:6]))

                    elif lly[ii][ij][ik - 1, 1] + 1 == lly[ii][ij][ik, 1]:

                        # close a potentially open previous line
                        vly.append(lly[ii][ij][ik, 3:6])

                if lly[ii][ij][ik + 1, 1] - 1 == lly[ii][ij][ik, 1]:

                    vly.append(lly[ii][ij][ik + 1, 3:6])

                if np.sum(tmp): # empty lists are false
                    llyLgth[ii, ij] = np.sum(tmp)

    TriaMesh.write_vtk(TriaMesh(np.array(vly), np.asarray(tly)-1), os.path.join(OUT_DIR, HEMI + '.grid-lines-y.vtk'))

    dfy = pd.DataFrame(llyLgth)

    dz = zip(range(0, dfy.shape[1]), [ "z"+str(z) for z in range(0, dfy.shape[1]) ])
    dfy = dfy.rename(columns=dict(dz))

    dx = zip(range(0, dfy.shape[0]), [ "x"+str(z) for z in range(0, dfy.shape[0]) ])
    dfy = dfy.rename(index=dict(dx))

    dfy.to_csv(os.path.join(OUT_DIR, HEMI + '.grid-segments-y.csv'))

    # z

    vlz = list()
    tlz = list()
    nlz = list()
    llzLgth = np.full((np.shape(origV4)[0], np.shape(origV4)[1]), np.nan)

    for ii in range(len(llz)):

        for ij in range(len(llz[ii])):

            if len(llz[ii][ij])>1:

                tmp = list()

                for ik in range(len(llz[ii][ij]) - 1):

                    if llz[ii][ij][ik, 2] + 1 == llz[ii][ij][ik + 1, 2]:

                        vlz.append(llz[ii][ij][ik, 3:6])

                        tlz.append([len(vlz), len(vlz), len(vlz) + 1])

                        nlz.append([ ii, ij, ik, np.linalg.norm(llz[ii][ij][ik + 1, 3:6] - llz[ii][ij][ik, 3:6]) ])

                        tmp.append(np.linalg.norm(llz[ii][ij][ik + 1, 3:6] - llz[ii][ij][ik, 3:6]))

                    elif llz[ii][ij][ik - 1, 2] + 1 == llz[ii][ij][ik, 2]:

                        # close a potentially open previous line
                        vlz.append(llz[ii][ij][ik, 3:6])

                if llz[ii][ij][ik + 1, 2] - 1 == llz[ii][ij][ik, 2]:

                    vlz.append(llz[ii][ij][ik + 1, 3:6])

                if np.sum(tmp): # empty lists are false
                    llzLgth[ii, ij] = np.sum(tmp)

    TriaMesh.write_vtk(TriaMesh(np.array(vlz), np.asarray(tlz)-1), os.path.join(OUT_DIR, HEMI + '.grid-lines-z.vtk'))

    dfz = pd.DataFrame(llzLgth)

    dy = zip(range(0, dfz.shape[1]), [ "y"+str(z) for z in range(0, dfz.shape[1]) ])
    dfz = dfz.rename(columns=dict(dy))

    dx = zip(range(0, dfz.shape[0]), [ "x"+str(z) for z in range(0, dfz.shape[0]) ])
    dfz = dfz.rename(index=dict(dx))

    dfz.to_csv(os.path.join(OUT_DIR, HEMI + '.grid-segments-z.csv'))

    # --------------------------------------------------------------------------
    # export origV4flat

    pd.DataFrame(origV4flat).to_csv(os.path.join(OUT_DIR, HEMI + '.grid-lines.csv'), header=False, index=False)

    # --------------------------------------------------------------------------
    # create checkerboard (we only do it for one dimension: z)

    origV4flatTria = dict()
    origV4flatCol = dict()

    minIdx = 0
    midIdx = int(np.round((len(lLVL4z) - 1) / 2))
    maxIdx = len(lLVL4z)-1

    for k in [minIdx, midIdx, maxIdx]:

        origV4flatTria[k] = list()
        origV4flatCol[k] = list()

        for i in range(0, len(llz) - 1):

            if allowRagged or (
            (len(np.unique(origV4flat[np.logical_and(origV4flat[:, 0] == i, origV4flat[:, 2] == k), 1])) == paramsTHYk) and
            (len(np.unique(origV4flat[np.logical_and(origV4flat[:, 0] == i+1, origV4flat[:, 2] == k), 1])) == paramsTHYk)
            ):

                for j in range(0, len(llz[i])-1):

                    # result may be empty, but that's OK: nothing will be added then

                    # we are doing this per square, could be done per tria also (split
                    # into two if-clauses)

                    if allowRaggedTriangles is False:

                        if all((
                            (len(np.where(np.logical_and(origV4flat[:, 0] == i    , np.logical_and(origV4flat[:, 1] == j    , origV4flat[:, 2] == k)))[0])>0),
                            (len(np.where(np.logical_and(origV4flat[:, 0] == i + 1, np.logical_and(origV4flat[:, 1] == j    , origV4flat[:, 2] == k)))[0])>0),
                            (len(np.where(np.logical_and(origV4flat[:, 0] == i    , np.logical_and(origV4flat[:, 1] == j + 1, origV4flat[:, 2] == k)))[0])>0),
                            (len(np.where(np.logical_and(origV4flat[:, 0] == i + 1, np.logical_and(origV4flat[:, 1] == j + 1, origV4flat[:, 2] == k)))[0])>0)
                            )):

                            origV4flatTria[k].append(
                                [
                                    np.where(np.logical_and(origV4flat[:, 0] == i    , np.logical_and(origV4flat[:, 1] == j    , origV4flat[:, 2] == k)))[0][0],
                                    np.where(np.logical_and(origV4flat[:, 0] == i + 1, np.logical_and(origV4flat[:, 1] == j    , origV4flat[:, 2] == k)))[0][0],
                                    np.where(np.logical_and(origV4flat[:, 0] == i    , np.logical_and(origV4flat[:, 1] == j + 1, origV4flat[:, 2] == k)))[0][0]
                                ]
                            )

                            origV4flatTria[k].append(
                                [
                                    np.where(np.logical_and(origV4flat[:, 0] == i + 1, np.logical_and(origV4flat[:, 1] == j    , origV4flat[:, 2] == k)))[0][0],
                                    np.where(np.logical_and(origV4flat[:, 0] == i    , np.logical_and(origV4flat[:, 1] == j + 1, origV4flat[:, 2] == k)))[0][0],
                                    np.where(np.logical_and(origV4flat[:, 0] == i + 1, np.logical_and(origV4flat[:, 1] == j + 1, origV4flat[:, 2] == k)))[0][0]
                                ]
                            )

                            origV4flatCol[k].append(np.mod(i + j, 2))
                            origV4flatCol[k].append(np.mod(i + j, 2))

                    elif allowRaggedTriangles is True:

                        if all((
                            (len(np.where(np.logical_and(origV4flat[:, 0] == i    , np.logical_and(origV4flat[:, 1] == j    , origV4flat[:, 2] == k)))[0])>0),
                            (len(np.where(np.logical_and(origV4flat[:, 0] == i + 1, np.logical_and(origV4flat[:, 1] == j    , origV4flat[:, 2] == k)))[0])>0),
                            (len(np.where(np.logical_and(origV4flat[:, 0] == i    , np.logical_and(origV4flat[:, 1] == j + 1, origV4flat[:, 2] == k)))[0])>0)
                            )):

                            origV4flatTria[k].append(
                                [
                                    np.where(np.logical_and(origV4flat[:, 0] == i    , np.logical_and(origV4flat[:, 1] == j    , origV4flat[:, 2] == k)))[0][0],
                                    np.where(np.logical_and(origV4flat[:, 0] == i + 1, np.logical_and(origV4flat[:, 1] == j    , origV4flat[:, 2] == k)))[0][0],
                                    np.where(np.logical_and(origV4flat[:, 0] == i    , np.logical_and(origV4flat[:, 1] == j + 1, origV4flat[:, 2] == k)))[0][0]
                                ]
                            )

                            origV4flatCol[k].append(np.mod(i + j, 2))

                        if all((
                            (len(np.where(np.logical_and(origV4flat[:, 0] == i + 1, np.logical_and(origV4flat[:, 1] == j    , origV4flat[:, 2] == k)))[0])>0),
                            (len(np.where(np.logical_and(origV4flat[:, 0] == i    , np.logical_and(origV4flat[:, 1] == j + 1, origV4flat[:, 2] == k)))[0])>0),
                            (len(np.where(np.logical_and(origV4flat[:, 0] == i + 1, np.logical_and(origV4flat[:, 1] == j + 1, origV4flat[:, 2] == k)))[0])>0)
                            )):

                            origV4flatTria[k].append(
                                [
                                    np.where(np.logical_and(origV4flat[:, 0] == i + 1, np.logical_and(origV4flat[:, 1] == j    , origV4flat[:, 2] == k)))[0][0],
                                    np.where(np.logical_and(origV4flat[:, 0] == i    , np.logical_and(origV4flat[:, 1] == j + 1, origV4flat[:, 2] == k)))[0][0],
                                    np.where(np.logical_and(origV4flat[:, 0] == i + 1, np.logical_and(origV4flat[:, 1] == j + 1, origV4flat[:, 2] == k)))[0][0]
                                ]
                            )

                            origV4flatCol[k].append(np.mod(i + j, 2))

    #

    vMin = origV4flat[:, 3:6]
    vMid = origV4flat[:, 3:6]
    vMax = origV4flat[:, 3:6]

    tMin = np.array(origV4flatTria[minIdx])
    tMid = np.array(origV4flatTria[midIdx])
    tMax = np.array(origV4flatTria[maxIdx])

    #

    cMin = origV4flat[:, 0:3]
    cMid = origV4flat[:, 0:3]
    cMax = origV4flat[:, 0:3]

    triaMin = TriaMesh(vMin, tMin)
    triaMid = TriaMesh(vMid, tMid)
    triaMax = TriaMesh(vMax, tMax)

    vMinKeep, vMinDel = triaMin.rm_free_vertices_()
    vMidKeep, vMidDel = triaMid.rm_free_vertices_()
    vMaxKeep, vMaxDel = triaMax.rm_free_vertices_()

    cMinRm = cMid[vMinKeep,:]
    cMidRm = cMid[vMidKeep,:]
    cMaxRm = cMid[vMaxKeep,:]

    triaMinRm = TriaMesh(triaMin.v, triaMin.t)
    triaMidRm = TriaMesh(triaMid.v, triaMid.t)
    triaMaxRm = TriaMesh(triaMax.v, triaMax.t)

    #

    if skipOrient is False:
        if sc.connected_components(triaMinRm._construct_adj_sym())[0]==1:
            triaMinRm.orient_()
            skipOrientMin = False
        else:
            logging.warning("Warning: exterior surface contains more than one component, not orienting. Check your results.")
            skipOrientMin = True
        if sc.connected_components(triaMidRm._construct_adj_sym())[0]==1:
            triaMidRm.orient_()
            skipOrientMid = False
        else:
            logging.warning("Warning: mid-surface contains more than one component, not orienting. Check your results.")
            skipOrientMid = True
        if sc.connected_components(triaMaxRm._construct_adj_sym())[0]==1:
            triaMaxRm.orient_()
            skipOrientMax = False
        else:
            logging.warning("Warning: interior surface contains more than one component, not orienting. Check your results.")
            skipOrientMax = True
    else:
        skipOrientMin = True
        skipOrientMid = True
        skipOrientMax = True

    #

    TriaMesh.write_vtk(triaMinRm, os.path.join(OUT_DIR, HEMI + '.ext-surface.vtk'))
    TriaMesh.write_vtk(triaMidRm, os.path.join(OUT_DIR, HEMI + '.mid-surface.vtk'))
    TriaMesh.write_vtk(triaMaxRm, os.path.join(OUT_DIR, HEMI + '.int-surface.vtk'))

    pd.DataFrame(np.concatenate((cMinRm, triaMinRm.v), axis=1)).to_csv(os.path.join(OUT_DIR, HEMI + '.ext-surface.csv'), header=False, index=False)
    pd.DataFrame(np.concatenate((cMidRm, triaMidRm.v), axis=1)).to_csv(os.path.join(OUT_DIR, HEMI + '.mid-surface.csv'), header=False, index=False)
    pd.DataFrame(np.concatenate((cMaxRm, triaMaxRm.v), axis=1)).to_csv(os.path.join(OUT_DIR, HEMI + '.int-surface.csv'), header=False, index=False)

    # --------------------------------------------------------------------------
    # export thickness

    thickness = list()

    for i in range(0, len(cMidRm)):
        thickness.append(np.array(dfz)[int(cMidRm[i,0]), int(cMidRm[i,1])])

    io.write_vfunc(os.path.join(OUT_DIR, HEMI + '.mid-surface.thickness.psol'), np.array(thickness))

    nb.freesurfer.save(nb.freesurfer.MGHImage(dataobj=np.array(thickness).astype("float32"), affine=None), filename=os.path.join(OUT_DIR, HEMI + '.mid-surface.thickness.mgh'))

    # --------------------------------------------------------------------------
    # export curvature

    if skipOrientMin is False:

        u1Min, u2Min, c1Min, c2Min = triaMinRm.curvature_tria()

        mean_curvMin= (c1Min + c2Min) / 2

        gauss_curvMin = c1Min * c2Min

        mean_curvMin = triaMinRm.map_tfunc_to_vfunc(mean_curvMin)

        gauss_curvMin = triaMinRm.map_tfunc_to_vfunc(gauss_curvMin)

        io.write_vfunc(os.path.join(OUT_DIR, HEMI + '.ext-surface.mean-curv.psol'), mean_curvMin)

        io.write_vfunc(os.path.join(OUT_DIR, HEMI + '.ext-surface.gauss-curv.psol'), gauss_curvMin)

        nb.freesurfer.save(nb.freesurfer.MGHImage(dataobj=mean_curvMin.astype("float32"), affine=None), filename=os.path.join(OUT_DIR, HEMI + '.ext-surface.mean-curv.mgh'))

        nb.freesurfer.save(nb.freesurfer.MGHImage(dataobj=gauss_curvMin.astype("float32"), affine=None), filename=os.path.join(OUT_DIR, HEMI + '.ext-surface.gauss-curv.mgh'))

        pd.DataFrame(np.concatenate((cMinRm, np.array(mean_curvMin, ndmin=2).transpose()), axis=1)).to_csv(os.path.join(OUT_DIR, HEMI + '.ext-surface.mean-curv.csv'), header=False, index=False)

        pd.DataFrame(np.concatenate((cMinRm, np.array(gauss_curvMin, ndmin=2).transpose()), axis=1)).to_csv(os.path.join(OUT_DIR, HEMI + '.ext-surface.gauss-curv.csv'), header=False, index=False)

    else:

        logging.warning("Warning: Not computing curvature for exterior surface, check surface.")

    if skipOrientMid is False:

        u1Mid, u2Mid, c1Mid, c2Mid = triaMidRm.curvature_tria()

        mean_curvMid= (c1Mid + c2Mid) / 2

        gauss_curvMid = c1Mid * c2Mid

        mean_curvMid = triaMidRm.map_tfunc_to_vfunc(mean_curvMid)

        gauss_curvMid = triaMidRm.map_tfunc_to_vfunc(gauss_curvMid)

        io.write_vfunc(os.path.join(OUT_DIR, HEMI + '.mid-surface.mean-curv.psol'), mean_curvMid)

        io.write_vfunc(os.path.join(OUT_DIR, HEMI + '.mid-surface.gauss-curv.psol'), gauss_curvMid)

        nb.freesurfer.save(nb.freesurfer.MGHImage(dataobj=mean_curvMid.astype("float32"), affine=None), filename=os.path.join(OUT_DIR, HEMI + '.mid-surface.mean-curv.mgh'))

        nb.freesurfer.save(nb.freesurfer.MGHImage(dataobj=gauss_curvMid.astype("float32"), affine=None), filename=os.path.join(OUT_DIR, HEMI + '.mid-surface.gauss-curv.mgh'))

        pd.DataFrame(np.concatenate((cMidRm, np.array(mean_curvMid, ndmin=2).transpose()), axis=1)).to_csv(os.path.join(OUT_DIR, HEMI + '.mid-surface.mean-curv.csv'), header=False, index=False)

        pd.DataFrame(np.concatenate((cMidRm, np.array(gauss_curvMid, ndmin=2).transpose()), axis=1)).to_csv(os.path.join(OUT_DIR, HEMI + '.mid-surface.gauss-curv.csv'), header=False, index=False)

    else:

        logging.warning("Warning: Not computing curvature for mid-surface, check surface.")

    if skipOrientMax is False:

        u1Max, u2Max, c1Max, c2Max = triaMaxRm.curvature_tria()

        mean_curvMax= (c1Max + c2Max) / 2

        gauss_curvMax = c1Max * c2Max

        mean_curvMax = triaMaxRm.map_tfunc_to_vfunc(mean_curvMax)

        gauss_curvMax = triaMaxRm.map_tfunc_to_vfunc(gauss_curvMax)

        io.write_vfunc(os.path.join(OUT_DIR, HEMI + '.int-surface.mean-curv.psol'), mean_curvMax)

        io.write_vfunc(os.path.join(OUT_DIR, HEMI + '.int-surface.gauss-curv.psol'), gauss_curvMax)

        nb.freesurfer.save(nb.freesurfer.MGHImage(dataobj=mean_curvMax.astype("float32"), affine=None), filename=os.path.join(OUT_DIR, HEMI + '.int-surface.mean-curv.mgh'))

        nb.freesurfer.save(nb.freesurfer.MGHImage(dataobj=gauss_curvMax.astype("float32"), affine=None), filename=os.path.join(OUT_DIR, HEMI + '.int-surface.gauss-curv.mgh'))

        pd.DataFrame(np.concatenate((cMaxRm, np.array(mean_curvMax, ndmin=2).transpose()), axis=1)).to_csv(os.path.join(OUT_DIR, HEMI + '.int-surface.mean-curv.csv'), header=False, index=False)

        pd.DataFrame(np.concatenate((cMaxRm, np.array(gauss_curvMax, ndmin=2).transpose()), axis=1)).to_csv(os.path.join(OUT_DIR, HEMI + '.int-surface.gauss-curv.csv'), header=False, index=False)

    else:

        logging.warning("Warning: Not computing curvature for interior surface, check surface.")

    # -----------------------------------------------------------------------------
    # create hull

    # here: fast --> slow: d1 --> d2 --> d3

    f = origV4flat

    d1 = paramsTHXk
    d2 = paramsTHYk
    d3 = paramsTHZk

    #

    vert = np.full((d1, d2, d3, 3), np.nan)
    tria = list()

    for l in range(0, len(f)):

        # get indices
        i, j, k = f[l, 0:3].astype(int)

        # assign coordinates
        vert[i,j,k,:] = f[l, 3:6]

        # get trias

        # x==0|x==d1-1
        if (i==0) or i==(d1-1):

            # x:yz
            if (np.mod(j+1, d2)>0) and (np.mod(k+1, d3)>0):
                tria.append([i*d2*d3+j*d3+k, i*d2*d3+j*d3+k+1, i*d2*d3+j*d3+k+d3+1])
                tria.append([i*d2*d3+j*d3+k, i*d2*d3+j*d3+k+d3+1, i*d2*d3+j*d3+k+d3])

        # y==0|y==d2-1
        if (j==0) or (j==(d2-1)):

            # y:xz
            if (np.mod(i+1, d1)>0) and (np.mod(k+1, d3)>0):
                tria.append([i*d2*d3+j*d3+k, i*d2*d3+j*d3+k+d3*d2, i*d2*d3+j*d3+k+1])
                tria.append([i*d2*d3+j*d3+k+1, i*d2*d3+j*d3+k+d3*d2, i*d2*d3+j*d3+k+d3*d2+1])

        # z==0|z==d3-1
        if (k==0) or (k==(d3-1)):

            # z:xy
            if (np.mod(i+1, d1)>0) and (np.mod(j+1, d2)>0):
                tria.append([i*d2*d3+j*d3+k, i*d2*d3+j*d3+k+d2*d3, i*d2*d3+j*d3+k+d3])
                tria.append([i*d2*d3+j*d3+k+d3, i*d2*d3+j*d3+k+d2*d3, i*d2*d3+j*d3+k+d2*d3+d3])


    vertList = np.reshape(vert, (d1*d2*d3, 3))
    triaNoNan = np.array(tria)[np.sum(np.sum(np.isnan(vertList[tria,:]), axis=2), axis=1)==0,:]

    TriaMesh.write_vtk(TriaMesh(vertList, triaNoNan), os.path.join(OUT_DIR, HEMI + '.hull.vtk'))

    # --------------------------------------------------------------------------
    # return

    return(params)


# ------------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# MAIN PART

if __name__ == "__main__":

    # Command line options and error checking

    options = options_parse()

    # Run analysis

    computeThickness(params=None,
        IN_MESH=options.IN_MESH,
        IN_FUNC=options.IN_FUNC,
        OUT_DIR=options.OUT_DIR,
        HEMI=options.HEMI)
