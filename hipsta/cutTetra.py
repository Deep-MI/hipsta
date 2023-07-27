"""
This module provides a function to cut tetrahedral meshes

"""

import os

import numpy as np
import scipy.sparse as sp
import scipy.linalg as la
import itertools as it

from lapy import TriaMesh, TetMesh, io, TetMesh, Solver

from .utils.get_levelsets import levelsetsTetra

# ------------------------------------------------------------------------------
# AUXILIARY FUNCTIONS
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# _importData()

def _importData(tetraFile, tetraIdxFile):

    # read tetra mesh

    tetMesh = TetMesh.read_vtk(tetraFile)
    v4 = tetMesh.v
    t4 = tetMesh.t

    # read tetra index file

    l4 = np.array(io.read_vfunc(tetraIdxFile))

    # return

    return v4, t4, l4


# ------------------------------------------------------------------------------
# _exportData()

def _exportData(tetraCutOutFile, tetraCutOutFileFunc, tetraCutOutDir, hemi, v4c, t4c, i4c):

    # write out cut mesh and overlay

    tetMesh = TetMesh(v=v4c, t=t4c)

    TetMesh.write_vtk(tetMesh, tetraCutOutFile)

    io.write_vfunc(tetraCutOutFileFunc, i4c)

    # create boundary mesh and output

    triaMesh4cBnd = tetMesh.boundary_tria()

    triaMesh4cBnd.orient_()

    TriaMesh.write_vtk(triaMesh4cBnd, os.path.join(tetraCutOutDir, hemi + '.bnd.cut.vtk'))

    # remove all trias that share a vertex at the cut planes, because we want
    # to have an open boundary mesh; also remove free vertices

    t4cBndOpen = triaMesh4cBnd.t[np.sum(i4c[triaMesh4cBnd.t]==1, axis=1)>=1, :]

    triaMesh4cBndOpen = TriaMesh(v=v4c, t=t4cBndOpen)

    triaMesh4cBndOpen.orient_()

    TriaMesh.write_vtk(triaMesh4cBndOpen, os.path.join(tetraCutOutDir, hemi + '.open.bnd.cut.vtk'))

    # remove free vertices from boundary mesh and output

    triaMesh4cBnd.rm_free_vertices_()

    TriaMesh.write_vtk(triaMesh4cBnd, os.path.join(tetraCutOutDir, hemi + '.rm.bnd.cut.vtk'))

    # remove free vertices from open mesh and output ...

    triaMesh4cBndOpen.rm_free_vertices_()

    TriaMesh.write_vtk(triaMesh4cBndOpen, os.path.join(tetraCutOutDir, hemi + '.rm.open.bnd.cut.vtk'))

    # write out mapping between open bnd and v4c

    with open(os.path.join(tetraCutOutDir, hemi + '.rm.open.bnd.cut.lst'), "w") as f:

        [ print(str(x), file=f) for x in np.unique(t4cBndOpen) ]


# ------------------------------------------------------------------------------
# _preprocessData()

def _preprocessData(v4, t4, l4):

    tetMesh = TetMesh(v=v4, t=t4)

    tetMeshSolver = Solver(tetMesh)

    A = tetMeshSolver.stiffness

    M = tetMeshSolver.mass

    P = tetMeshSolver.poisson(dtup=(np.where(l4!=0)[0],l4[np.where(l4!=0)[0]]))

    return A, M, P


# ------------------------------------------------------------------------------
# _getTetra()

def _getTetra(aa, ii):

    # we compute which subset of three of neighbors of x has an adj mtx of
    # ones exclusively

    tt = list()

    for iii in range(len(ii)):

        n = np.where(aa[:, ii[iii]])[0]

        n = n[np.where(n!=ii[iii])[0]]

        if len(n)>=3:

            nt = list(it.combinations(n, 3))

            for ij in range(len(nt)):

                if aa[np.ix_(nt[ij], nt[ij])].all():

                    tt.append(np.sort(np.hstack((ii[iii], nt[ij]))))

    if len(tt)>0:

        tt = np.unique(tt, axis=0)

    return tt


# ------------------------------------------------------------------------------
# _tetra26()

def _tetra26(tmp1,tmp2,tmp12,ind1,ind2,ind3,ind4,casevar):

    # create a local adj matrix: we know that tmp1, tmp2(1,:) and
    # tmp2(2,:) are connected within themselves
    ind = np.unique(
        np.vstack((
            list(it.combinations(np.where(np.isin(tmp12,tmp1))[0], 2)),
            list(it.combinations(np.where(np.isin(tmp12,tmp2[0,:]))[0], 2)),
            list(it.combinations(np.where(np.isin(tmp12,tmp2[1,:]))[0], 2))
            )),
        axis=0
        )
    a = np.zeros((tmp12.shape[0], tmp12.shape[0]))
    a[ind[:,0], ind[:,1]] = 1
    a = a + a.transpose()
    a = a + np.diag(np.ones(a.shape[0]))

    # add what we know
    if casevar==1:

        # 1-3, 1-4, 2-5, 2-6
        a[0, 2] = 1
        a[2, 0] = 1
        a[0, 3] = 1
        a[3, 0] = 1
        a[1, 4] = 1
        a[4, 1] = 1
        a[1, 5] = 1
        a[5, 1] = 1

    elif casevar==2:

        # 1-3, 1-5, 2-4, 2-6
        a[0, 2] = 1
        a[2, 0] = 1
        a[0, 4] = 1
        a[4, 0] = 1
        a[1, 3] = 1
        a[3, 1] = 1
        a[1, 5] = 1
        a[5, 1] = 1

    # add indices supplied by arguments
    a[ind1, ind2] = 1
    a[ind2, ind1] = 1
    a[ind3, ind4] = 1
    a[ind4, ind3] = 1

    # rewiring of the cut plane necessary if 1-6, 2-3
    if ind1==1 and ind2==2 and ind3==0 and ind4==5:
        a[3, 4] = 0
        a[4, 3] = 0
        a[2, 5] = 1
        a[5, 2] = 1

    # get tetra
    tetra = tmp12[_getTetra(a, range(0, 6))]

    # return
    return tetra, a


# ------------------------------------------------------------------------------
# _tetra33()

def _tetra33(tmp1,tmp2,tmp12,ind1,ind2,ind3,v,A):

    # create a local adj matrix: we know that tmp1 and tmp2 are
    # connected within themselves, and we know from the levelsetsTetra
    # code that tmp1(1) is connected with tmp2(1), tmp1(2) is connected
    # with tmp2(2), and tmp1(3) is connected with tmp2(3)
    a = la.block_diag(np.ones((tmp1.shape[0],tmp1.shape[0])), np.ones((tmp1.shape[0],tmp1.shape[0]))) + np.diag(np.ones(tmp1.shape[0]), k=3) + np.diag(np.ones(tmp1.shape[0]), k=-3)
    a = a + np.diag(np.ones(a.shape[0]))

    # now add those edges as specified by the arguments
    a[ind1[0], ind1[1]] = 1
    a[ind1[1], ind1[0]] = 1
    a[ind2[0], ind2[1]] = 1
    a[ind2[1], ind2[0]] = 1
    a[ind3[0], ind3[1]] = 1
    a[ind3[1], ind3[0]] = 1

    # get all tetras
    tetra = tmp12[_getTetra(a, range(0, 6))]

    # need to add a point in two special cases:
    if (ind1[0]==0 and ind1[1]==5 and ind2[0]==2 and ind2[1]==4 and ind3[0]==1 and ind3[1]==3) or (ind1[0]==2 and ind1[1]==3 and ind2[0]==1 and ind2[1]==5 and ind3[0]==0 and ind3[1]==4):

        # compute midpoints on each edge
        p1 = v[tmp1[0], :] + np.divide([ v[tmp2[0], :] - v[tmp1[0], :] ], 2)
        p2 = v[tmp1[1], :] + np.divide([ v[tmp2[1], :] - v[tmp1[1], :] ], 2)
        p3 = v[tmp1[2], :] + np.divide([ v[tmp2[2], :] - v[tmp1[2], :] ], 2)

        # compute center of triangle (p1, p2, p3)
        pc = np.divide((p1 + p2 + p3 ), 3)

        # add point to the list of vertices
        v = np.concatenate((v, pc), axis=0)

        # create new list of points
        tmp123 = np.hstack((tmp1, tmp2, np.shape(v)[0] - 1))

        # append a
        a = np.vstack((a,np.ones((1,np.shape(a)[1]))))
        a = np.hstack((a,np.ones((np.shape(a)[0],1))))

        # get tetra
        tetra = tmp123[_getTetra(a, range(0, 7))]

        # add one row/col to A and set diag elem to 1
        A = sp.coo_matrix(
            (np.ones(np.shape(A.nonzero())[1] + 1), (np.hstack((A.nonzero()[0], np.shape(v)[0] - 1)), np.hstack((A.nonzero()[1], np.shape(v)[0] - 1)))),
            shape=(np.shape(v)[0],np.shape(v)[0]))


    # return
    return tetra, a, v, A


# ------------------------------------------------------------------------------
# MAIN FUNCTION
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# cutTetra()

def cutTetra(params):
    """
    This is a function to cut tetrahedral meshes

    """

    # -------------------------------------------------------------------------
    # message
    # -------------------------------------------------------------------------

    print()
    print("--------------------------------------------------------------------------------")
    print("Cutting tetrahedral mesh")
    print()

    # -------------------------------------------------------------------------
    # evaluate input
    # -------------------------------------------------------------------------

    tetraFile = os.path.join(params.OUTDIR, params.HEMI + ".tetra.vtk")

    tetraIdxFile = os.path.join(params.OUTDIR, 'tetra-cut', params.HEMI + ".tetra-remove_bnd.psol")

    tetraCutOutFile = os.path.join(params.OUTDIR, params.HEMI + ".cut.vtk")

    tetraCutOutFileFunc = os.path.join(params.OUTDIR, "tetra-cut", params.HEMI + ".cut.psol")

    tetraCutOutDir = os.path.join(params.OUTDIR, "tetra-cut")

    labelHead = params.LUTDICT['jointhead']
    labelTail = params.LUTDICT['jointtail']
    labelBndHead = params.LUTDICT['bndhead']
    labelBndTail = params.LUTDICT['bndtail']
    labelBndCA4 = params.LUTDICT['bndca4']

    cutRange = params.internal.cut_range

    hemi = params.HEMI

    # -------------------------------------------------------------------------
    # get data
    # -------------------------------------------------------------------------

    v4, t4, l4idx = _importData(tetraFile, tetraIdxFile)

    l4idx[np.where(l4idx<labelBndTail)[0]] = 0
    l4idx[np.where(l4idx==labelBndCA4)[0]] = 0
    l4idx[np.where(l4idx==labelBndHead)[0]] = -1
    l4idx[np.where(l4idx==labelBndTail)[0]] = 1

    # -------------------------------------------------------------------------
    # preprocess data
    # -------------------------------------------------------------------------

    Ap, Mp, l4 = _preprocessData(v4, t4, l4idx)

    # -------------------------------------------------------------------------
    # create a smooth boundary
    # -------------------------------------------------------------------------

    # compute level sets
    vlTail, llTail, ilTail, jlTail, olTail = levelsetsTetra(v4, t4, l4, cutRange[0])
    vlHead, llHead, ilHead, jlHead, olHead = levelsetsTetra(v4, t4, l4, cutRange[1])

    # find all tetras that do not exceed the cutting criteria
    t4c1 = t4[np.where(np.sum(np.isin(t4, np.where((l4>cutRange[0]) & (l4<cutRange[1]))), axis=1)==4), :][0]

    # add new points to v and generate new triangles
    v4c = np.concatenate((v4, vlTail[0], vlHead[0]), axis=0)
    t4c = t4c1

    # create adjacency matrix
    idx = np.concatenate(
        (t4c[:,[0, 1]], t4c[:,[0, 2]], t4c[:,[0, 3]], t4c[:,[1, 2]], t4c[:,[1, 3]], t4c[:,[2, 3]]),
        axis=0)
    A = sp.eye(np.shape(v4c)[0])
    A = A + sp.coo_matrix(
        (np.ones((np.shape(idx)[0])), (idx[:,0],idx[:,1])),
        shape=(np.shape(v4c)[0],np.shape(v4c)[0])
        )
    A = A + A.transpose()
    A[A>0] = 1
    A.tolil()

    # -------------------------------------------------------------------------
    # iterate through all tetras on the level set labelTail
    # -------------------------------------------------------------------------

    for i in range(len(ilTail[0])):

        # check which points of the original tetra are not exceeding the
        # plane, can be 1 or 2 or 3 points
        tmp1 = t4[ilTail[0][i], olTail[0][i][:]]

        # find points on the level set, can be 1x3 or 2x3 points, where 1x3
        # only occurs with tmp1==1 and tmp1==3, and 2x3 only with tmp1==2. If
        # there are 2x3 points, only 4 of them are unique.
        # Note that we need to change from 1-based to 0-based indices here,
        # therefore we add '-1' at the end.
        tmp2 = (llTail[0][np.where(jlTail[0]==ilTail[0][i])[0][:], :] + np.shape(v4)[0] - 1).squeeze()

        # the general strategy is to test out all degrees of freedom that we
        # have to see if there is something fitting. That are 1 deg for 1|3,
        # 2 deg for 2|6, and 3*2 deg for 3|3. this does not include any
        # rewiring or adding points.

        if (np.size(tmp1)==1) & (np.size(tmp2)==3):

            # no check against A should be necessary here, because we have no
            # degrees of freedom in constructing the tetra. in other words,
            # this tetra should always fit.

            # concatenate tmp1 and tmp2. this is already the tetra.
            tetra = np.hstack((tmp1, tmp2))[np.newaxis, :]

        elif (np.size(tmp1)==2) & (np.size(tmp2)==6):

            # concatenate tmp1 and tmp2. drop duplicate points.
            idx = np.unique(tmp2.flatten(), return_index=True)[1]
            tmp12 = np.hstack((tmp1, [ tmp2.flatten()[idxi] for idxi in sorted(idx) ]))

            # compute tetras and local adjacency matrix for various edge
            # insertion options
            tetra21, A21 = _tetra26(tmp1, tmp2, tmp12, 0, 3, 0, 5, 2)
            tetra22, A22 = _tetra26(tmp1, tmp2, tmp12, 0, 3, 1, 4, 2)
            tetra23, A23 = _tetra26(tmp1, tmp2, tmp12, 1, 2, 0, 5, 2)
            tetra24, A24 = _tetra26(tmp1, tmp2, tmp12, 1, 2, 1, 4, 2)

            # check which variant is compatible with A: all edges from the
            # relevant subset of A must be present in the local adj mtx.
            # However, the updated local adj mtx may contain new edges that
            # will be added to A.
            if not ((A21-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra21
            elif not ((A22-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra22
            elif not ((A23-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra23
            elif not ((A24-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra24

        elif (np.size(tmp1)==3) & (np.size(tmp2)==3):

            # concatenate tmp1 and tmp2
            tmp12 = np.hstack((tmp1, tmp2))

            # compute tetras and local adjacency matrix for various edge
            # insertion options
            tetra111, a111, v4c111, A111 = _tetra33(tmp1,tmp2,tmp12,[0,5],[2,4],[0,4],v4c,A)
            tetra112, a112, v4c112, A112 = _tetra33(tmp1,tmp2,tmp12,[0,5],[2,4],[1,3],v4c,A)
            tetra121, a121, v4c121, A121 = _tetra33(tmp1,tmp2,tmp12,[0,5],[1,5],[0,4],v4c,A)
            tetra122, a122, v4c122, A122 = _tetra33(tmp1,tmp2,tmp12,[0,5],[1,5],[1,3],v4c,A)
            tetra211, a211, v4c211, A211 = _tetra33(tmp1,tmp2,tmp12,[2,3],[2,4],[0,4],v4c,A)
            tetra212, a212, v4c212, A212 = _tetra33(tmp1,tmp2,tmp12,[2,3],[2,4],[1,3],v4c,A)
            tetra221, a221, v4c221, A221 = _tetra33(tmp1,tmp2,tmp12,[2,3],[1,5],[0,4],v4c,A)
            tetra222, a222, v4c222, A222 = _tetra33(tmp1,tmp2,tmp12,[2,3],[1,5],[1,3],v4c,A)

            # check which variant is compatible with A: all edges from the
            # relevant subset of A must be present in the local adj mtx.
            # However, the updated local adj mtx may contain new edges that
            # will be added to A.
            if not ((a111[0:6,0:6]-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra111
                v4c = v4c111
                A = A111.tolil()
            elif not ((a112[0:6,0:6]-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra112
                v4c = v4c112
                A = A112.tolil()
            elif not ((a121[0:6,0:6]-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra121
                v4c = v4c121
                A = A121.tolil()
            elif not ((a122[0:6,0:6]-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra122
                v4c = v4c122
                A = A122.tolil()
            elif not ((a211[0:6,0:6]-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra211
                v4c = v4c211
                A = A211.tolil()
            elif not ((a212[0:6,0:6]-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra212
                v4c = v4c212
                A = A212.tolil()
            elif not ((a221[0:6,0:6]-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra221
                v4c = v4c221
                A = A221.tolil()
            elif not ((a222[0:6,0:6]-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra222
                v4c = v4c222
                A = A222.tolil()

        # update adjacency matrix

        for k in range(np.shape(tetra)[0]):

            A[tetra[k,0], tetra[k,1]] = 1
            A[tetra[k,1], tetra[k,0]] = 1

            A[tetra[k,0], tetra[k,2]] = 1
            A[tetra[k,2], tetra[k,0]] = 1

            A[tetra[k,0], tetra[k,3]] = 1
            A[tetra[k,3], tetra[k,0]] = 1

            A[tetra[k,1], tetra[k,2]] = 1
            A[tetra[k,2], tetra[k,1]] = 1

            A[tetra[k,1], tetra[k,3]] = 1
            A[tetra[k,3], tetra[k,1]] = 1

            A[tetra[k,2], tetra[k,3]] = 1
            A[tetra[k,3], tetra[k,2]] = 1

        # update list of tetras

        t4c = np.concatenate((t4c, tetra), axis=0)

    # -------------------------------------------------------------------------
    # iterate through all tetras on the level set labelHead
    # -------------------------------------------------------------------------

    for i in range(len(ilHead[0])):

        # check which points of the original tetra are not exceeding the
        # plane, can be 1 or 2 or 3 points; note the setdiff --> outlying
        # points are defined from the other direction as compared to vTail.
        tmp1 = t4[ ilHead[0][i], np.setdiff1d((0, 1, 2, 3), olHead[0][i][:]) ]

        # find points on the level set, can be 1x3 or 2x3 points, where 1x3
        # only occurs with tmp1==1 and tmp1==3, and 2x3 only with tmp1==2. If
        # there are 2x3 points, only 4 of them are unique.
        # Note that we need to change from 1-based to 0-based indices here,
        # therefore we add '-1' at the end.
        tmp2 = (llHead[0][np.where(jlHead[0]==ilHead[0][i])[0][:], :] + np.shape(v4)[0] + np.shape(vlTail[0])[0] - 1).squeeze()

        # the general strategy is to test out all degrees of freedom that we
        # have to see if there is something fitting. That are 1 deg for 1|3,
        # 2 deg for 2|6, and 3*2 deg for 3|3. this does not include any
        # rewiring or adding points.

        if (np.size(tmp1)==1) & (np.size(tmp2)==3):

            # no check against A should be necessary here, because we have no
            # degrees of freedom in constructing the tetra. in other words,
            # this tetra should always fit.

            # concatenate tmp1 and tmp2. this is already the tetra.
            tetra = np.hstack((tmp1, tmp2))[np.newaxis, :]

        elif (np.size(tmp1)==2) & (np.size(tmp2)==6):

            # concatenate tmp1 and tmp2. drop duplicate points.
            idx = np.unique(tmp2.flatten(), return_index=True)[1]
            tmp12 = np.hstack((tmp1, [ tmp2.flatten()[idxi] for idxi in sorted(idx) ]))

            # compute tetras and local adjacency matrix for various edge
            # insertion options
            tetra11, A11 = _tetra26(tmp1, tmp2, tmp12, 0, 4, 0, 5, 1)
            tetra12, A12 = _tetra26(tmp1, tmp2, tmp12, 0, 4, 1, 3, 1)
            tetra13, A13 = _tetra26(tmp1, tmp2, tmp12, 1, 2, 0, 5, 1)
            tetra14, A14 = _tetra26(tmp1, tmp2, tmp12, 1, 2, 1, 3, 1)

            # check which variant is compatible with A: all edges from the
            # relevant subset of A must be present in the local adj mtx.
            # However, the updated local adj mtx may contain new edges that
            # will be added to A.

            if not ((A11-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra11
            elif not ((A12-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra12
            elif not ((A13-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra13
            elif not ((A14-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra14

        elif (np.size(tmp1)==3) & (np.size(tmp2)==3):

            # concatenate tmp1 and tmp2
            tmp12 = np.hstack((tmp1, tmp2))

            # compute tetras and local adjacency matrix for various edge
            # insertion options
            tetra111, a111, v4c111, A111 = _tetra33(tmp1,tmp2,tmp12,[0,5],[2,4],[0,4],v4c,A)
            tetra112, a112, v4c112, A112 = _tetra33(tmp1,tmp2,tmp12,[0,5],[2,4],[1,3],v4c,A)
            tetra121, a121, v4c121, A121 = _tetra33(tmp1,tmp2,tmp12,[0,5],[1,5],[0,4],v4c,A)
            tetra122, a122, v4c122, A122 = _tetra33(tmp1,tmp2,tmp12,[0,5],[1,5],[1,3],v4c,A)
            tetra211, a211, v4c211, A211 = _tetra33(tmp1,tmp2,tmp12,[2,3],[2,4],[0,4],v4c,A)
            tetra212, a212, v4c212, A212 = _tetra33(tmp1,tmp2,tmp12,[2,3],[2,4],[1,3],v4c,A)
            tetra221, a221, v4c221, A221 = _tetra33(tmp1,tmp2,tmp12,[2,3],[1,5],[0,4],v4c,A)
            tetra222, a222, v4c222, A222 = _tetra33(tmp1,tmp2,tmp12,[2,3],[1,5],[1,3],v4c,A)

            # check which variant is compatible with A: all edges from the
            # relevant subset of A must be present in the local adj mtx.
            # However, the updated local adj mtx may contain new edges that
            # will be added to A.

            if not ((a111[0:6,0:6]-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra111
                v4c = v4c111
                A = A111.tolil()
            elif not ((a112[0:6,0:6]-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra112
                v4c = v4c112
                A = A112.tolil()
            elif not ((a121[0:6,0:6]-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra121
                v4c = v4c121
                A = A121.tolil()
            elif not ((a122[0:6,0:6]-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra122
                v4c = v4c122
                A = A122.tolil()
            elif not ((a211[0:6,0:6]-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra211
                v4c = v4c211
                A = A211.tolil()
            elif not ((a212[0:6,0:6]-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra212
                v4c = v4c212
                A = A212.tolil()
            elif not ((a221[0:6,0:6]-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra221
                v4c = v4c221
                A = A221.tolil()
            elif not ((a222[0:6,0:6]-A[tmp12,:][:,tmp12].todense())<0).any():
                tetra = tetra222
                v4c = v4c222
                A = A222.tolil()

        # update adjacency matrix

        for k in range(np.shape(tetra)[0]):

            A[tetra[k,0], tetra[k,1]] = 1
            A[tetra[k,1], tetra[k,0]] = 1

            A[tetra[k,0], tetra[k,2]] = 1
            A[tetra[k,2], tetra[k,0]] = 1

            A[tetra[k,0], tetra[k,3]] = 1
            A[tetra[k,3], tetra[k,0]] = 1

            A[tetra[k,1], tetra[k,2]] = 1
            A[tetra[k,2], tetra[k,1]] = 1

            A[tetra[k,1], tetra[k,3]] = 1
            A[tetra[k,3], tetra[k,1]] = 1

            A[tetra[k,2], tetra[k,3]] = 1
            A[tetra[k,3], tetra[k,2]] = 1

        # update list of tetras

        t4c = np.concatenate((t4c, tetra), axis=0)

    #
    t4c = np.unique(t4c, axis=0)

    #
    i4c = np.concatenate((
        1 * np.ones(np.shape(v4)[0]), # original points (includes exterior and interior points)
        labelTail * np.ones(np.shape(vlTail[0])[0]), # labelTail points
        labelHead * np.ones(np.shape(vlHead[0])[0]), # labelHead points
        2 * np.ones(np.shape(v4c)[0] - np.shape(v4)[0]- np.shape(vlTail[0])[0] - np.shape(vlHead[0])[0]) # new points
    ))

    # -------------------------------------------------------------------------
    # export data
    # -------------------------------------------------------------------------

    _exportData(tetraCutOutFile, tetraCutOutFileFunc, tetraCutOutDir, hemi, v4c, t4c, i4c)

    # -------------------------------------------------------------------------
    # return
    # -------------------------------------------------------------------------

    return params
