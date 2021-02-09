"""
This module provides a function to compute a cube parametrization

"""

# ------------------------------------------------------------------------------
# AUXILIARY FUNCTIONS
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# getSeam()

def getSeam(v4c, t4c, i4c, k4c, v4cBndOpenKeep, t4cBndOpen, anisoLaplEvec):

    # --------------------------------------------------------------------------
    # getSeam subfunctions
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    # getSeam subfunction #1

    # case 1: one nan in the tetra
    def getSeamCase1(v4c, t4c, i4c, k4c, vfuncXEv1, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn, i):

        # do we have one or two neg values?
        # case 1a: one neg value
        if np.sum(vfuncXEv1[t4c[i, :]]<0) == 1:

            # get indices
            idxPos1 = np.where(vfuncXEv1[t4c[i, :]]>0)[0][0]
            idxPos2 = np.where(vfuncXEv1[t4c[i, :]]>0)[0][1]
            idxNeg1 = np.where(vfuncXEv1[t4c[i, :]]<0)[0][0]
            idxNan1 = np.where(np.isnan(vfuncXEv1[t4c[i, :]]))[0][0]

            # --------------------------------------------------------------
            # get new point 1 (we go neg-->pos)

            # scaling factor
            sf1 = np.abs(vfuncXEv1[t4c[i, idxNeg1]]) / (vfuncXEv1[t4c[i, idxPos1]] - vfuncXEv1[t4c[i, idxNeg1]])

            # difference vector
            dv1 = v4c[t4c[i, idxPos1], :] - v4c[t4c[i, idxNeg1], :]

            # compute new point
            pt1 = sf1 * dv1 + v4c[t4c[i, idxNeg1], :]

            # --------------------------------------------------------------
            # get new point 2

            # scaling factor
            sf2 = np.abs(vfuncXEv1[t4c[i, idxNeg1]]) / (vfuncXEv1[t4c[i, idxPos2]] - vfuncXEv1[t4c[i, idxNeg1]])

            # difference vector
            dv2 = v4c[t4c[i, idxPos2], :] - v4c[t4c[i, idxNeg1], :]

            # compute new point
            pt2 = sf2 * dv2 + v4c[t4c[i, idxNeg1], :]

            # --------------------------------------------------------------
            # store and connect

            # check for duplicate points (pt1)
            if newVtcsAdj[t4c[i, idxPos1], t4c[i, idxNeg1]] == 0:
                # add to list
                v4c = np.concatenate((v4c, np.array(pt1, ndmin=2)))
                # keep index
                idxPt1 = len(v4c) - 1
                # update adj mtx
                newVtcsAdj[t4c[i, idxPos1], t4c[i, idxNeg1]] = idxPt1
                newVtcsAdj[t4c[i, idxNeg1], t4c[i, idxPos1]] = idxPt1
                # also append i4c and k4c
                i4c = np.concatenate((i4c, st.mode(i4c[t4c[i,:]])[0]))
                k4c = np.concatenate((k4c, st.mode(k4c[t4c[i,:]])[0]))
                # let's remember the new indices, and whether or not it's a
                # 'positive' or 'negative' tetra (in terms of the 2nd
                # eigenfunction)
                newVtcs = np.concatenate((newVtcs, [idxPt1]))
                tmpSgn = np.unique(vfuncXEv2[t4c[i,:]][np.logical_not(np.isnan(vfuncXEv2[t4c[i,:]]))])
                if len(tmpSgn)==1:
                    newVtcsSgn = np.concatenate((newVtcsSgn, tmpSgn))
                else:
                    print("Inconsistency while creating seam, exiting.")
                    sys.exit(1)
            else:
                idxPt1 = newVtcsAdj[t4c[i, idxPos1], t4c[i, idxNeg1]]

            # check for duplicate points (pt2)
            if newVtcsAdj[t4c[i, idxPos2], t4c[i, idxNeg1]] == 0:
                # add to list
                v4c = np.concatenate((v4c, np.array(pt2, ndmin=2)))
                # keep index
                idxPt2 = len(v4c) - 1
                # update adj mtx
                newVtcsAdj[t4c[i, idxPos2], t4c[i, idxNeg1]] = idxPt2
                newVtcsAdj[t4c[i, idxNeg1], t4c[i, idxPos2]] = idxPt2
                # also append i4c and k4c
                i4c = np.concatenate((i4c, st.mode(i4c[t4c[i,:]])[0]))
                k4c = np.concatenate((k4c, st.mode(k4c[t4c[i,:]])[0]))
                # let's remember the new indices, and whether or not it's a
                # 'positive' or 'negative' tetra (in terms of the 2nd
                # eigenfunction)
                newVtcs = np.concatenate((newVtcs, [idxPt2]))
                tmpSgn = np.unique(vfuncXEv2[t4c[i,:]][np.logical_not(np.isnan(vfuncXEv2[t4c[i,:]]))])
                if len(tmpSgn)==1:
                    newVtcsSgn = np.concatenate((newVtcsSgn, tmpSgn))
                else:
                    print("Inconsistency while creating seam, exiting.")
                    sys.exit(1)
            else:
                idxPt2 = newVtcsAdj[t4c[i, idxPos2], t4c[i, idxNeg1]]

            # create new tetras
            # Neg1, pt1, pt2, nan
            # pt1, pt2, Pos1, nan
            # t2, Pos1, Pos2, nan

            newTetra = np.concatenate((newTetra, np.array((t4c[i, idxNeg1], idxPt1, idxPt2, t4c[i, idxNan1]), ndmin=2)))
            newTetra = np.concatenate((newTetra, np.array((idxPt1, idxPt2, t4c[i, idxPos1], t4c[i, idxNan1]), ndmin=2)))
            newTetra = np.concatenate((newTetra, np.array((idxPt2, t4c[i, idxPos1], t4c[i, idxPos2], t4c[i, idxNan1]), ndmin=2)))

        # case 1b: two neg value
        elif np.sum(vfuncXEv1[t4c[i, :]]<0) == 2:

            # get indices
            idxPos1 = np.where(vfuncXEv1[t4c[i, :]]>0)[0][0]
            idxNeg1 = np.where(vfuncXEv1[t4c[i, :]]<0)[0][0]
            idxNeg2 = np.where(vfuncXEv1[t4c[i, :]]<0)[0][1]
            idxNan1 = np.where(np.isnan(vfuncXEv1[t4c[i, :]]))[0][0]

            # --------------------------------------------------------------
            # get new point 1

            # scaling factor
            sf1 = np.abs(vfuncXEv1[t4c[i, idxNeg1]]) / (vfuncXEv1[t4c[i, idxPos1]] - vfuncXEv1[t4c[i, idxNeg1]])

            # difference vector
            dv1 = v4c[t4c[i, idxPos1], :] - v4c[t4c[i, idxNeg1], :]

            # compute new point
            pt1 = sf1 * dv1 + v4c[t4c[i, idxNeg1], :]

            # --------------------------------------------------------------
            # get new point 2

            # scaling factor
            sf2 = np.abs(vfuncXEv1[t4c[i, idxNeg2]]) / (vfuncXEv1[t4c[i, idxPos1]] - vfuncXEv1[t4c[i, idxNeg2]])

            # difference vector
            dv2 = v4c[t4c[i, idxPos1], :] - v4c[t4c[i, idxNeg2], :]

            # compute new point
            pt2 = sf2 * dv2 + v4c[t4c[i, idxNeg2], :]

            # --------------------------------------------------------------
            # store and connect

            # check for duplicate points (pt1)
            if newVtcsAdj[t4c[i, idxPos1], t4c[i, idxNeg1]] == 0:
                # add to list
                v4c = np.concatenate((v4c, np.array(pt1, ndmin=2)))
                # keep index
                idxPt1 = len(v4c) - 1
                # update adj mtx
                newVtcsAdj[t4c[i, idxPos1], t4c[i, idxNeg1]] = idxPt1
                newVtcsAdj[t4c[i, idxNeg1], t4c[i, idxPos1]] = idxPt1
                # also append i4c and k4c
                i4c = np.concatenate((i4c, st.mode(i4c[t4c[i,:]])[0]))
                k4c = np.concatenate((k4c, st.mode(k4c[t4c[i,:]])[0]))
                # let's remember the new indices, and whether or not it's a
                # 'positive' or 'negative' tetra (in terms of the 2nd
                # eigenfunction)
                newVtcs = np.concatenate((newVtcs, [idxPt1]))
                tmpSgn = np.unique(vfuncXEv2[t4c[i,:]][np.logical_not(np.isnan(vfuncXEv2[t4c[i,:]]))])
                if len(tmpSgn)==1:
                    newVtcsSgn = np.concatenate((newVtcsSgn, tmpSgn))
                else:
                    print("Inconsistency while creating seam, exiting.")
                    sys.exit(1)
            else:
                idxPt1 = newVtcsAdj[t4c[i, idxPos1], t4c[i, idxNeg1]]

            # check for duplicate points (pt2)
            if newVtcsAdj[t4c[i, idxPos1], t4c[i, idxNeg2]] == 0:
                # add to list
                v4c = np.concatenate((v4c, np.array(pt2, ndmin=2)))
                # keep index
                idxPt2 = len(v4c) - 1
                # update adj mtx
                newVtcsAdj[t4c[i, idxPos1], t4c[i, idxNeg2]] = idxPt2
                newVtcsAdj[t4c[i, idxNeg2], t4c[i, idxPos1]] = idxPt2
                # also append i4c and k4c
                i4c = np.concatenate((i4c, st.mode(i4c[t4c[i,:]])[0]))
                k4c = np.concatenate((k4c, st.mode(k4c[t4c[i,:]])[0]))
                # let's remember the new indices, and whether or not it's a
                # 'positive' or 'negative' tetra (in terms of the 2nd
                # eigenfunction)
                newVtcs = np.concatenate((newVtcs, [idxPt2]))
                tmpSgn = np.unique(vfuncXEv2[t4c[i,:]][np.logical_not(np.isnan(vfuncXEv2[t4c[i,:]]))])
                if len(tmpSgn)==1:
                    newVtcsSgn = np.concatenate((newVtcsSgn, tmpSgn))
                else:
                    print("Inconsistency while creating seam, exiting.")
                    sys.exit(1)
            else:
                idxPt2 = newVtcsAdj[t4c[i, idxPos1], t4c[i, idxNeg2]]

            # create new tetras
            # pos1, pt1, pt2, nan
            # pt1, pt2, neg1, nan
            # pt2, neg1, neg2, nan

            newTetra = np.concatenate((newTetra, np.array((t4c[i, idxPos1], idxPt1, idxPt2, t4c[i, idxNan1]), ndmin=2)))
            newTetra = np.concatenate((newTetra, np.array((idxPt1, idxPt2, t4c[i, idxNeg1], t4c[i, idxNan1]), ndmin=2)))
            newTetra = np.concatenate((newTetra, np.array((idxPt2, t4c[i, idxNeg1], t4c[i, idxNeg2], t4c[i, idxNan1]), ndmin=2)))

        else:
            print("Inconsistency while creating seam (incorrect number of negative values), exiting.")
            sys.exit(1)

        # return
        return v4c, t4c, i4c, k4c, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn

    # --------------------------------------------------------------------------
    # getSeam subfunction #2

    # case 2: two nans in tetra
    def getSeamCase2(v4c, t4c, i4c, k4c, vfuncXEv1, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn, i):

        # get indices
        idxPos1 = np.where(vfuncXEv1[t4c[i, :]]>0)[0][0]
        idxNeg1 = np.where(vfuncXEv1[t4c[i, :]]<0)[0][0]
        idxNan1 = np.where(np.isnan(vfuncXEv1[t4c[i, :]]))[0][0]
        idxNan2 = np.where(np.isnan(vfuncXEv1[t4c[i, :]]))[0][1]

        # --------------------------------------------------------------
        # get new point
        np.where(vfuncXEv1[t4c[i, :]]<0)[0][0]
        np.where(vfuncXEv1[t4c[i, :]]>0)[0][0]

        # scaling factor
        sf1 = np.abs(vfuncXEv1[t4c[i, idxNeg1]]) / (vfuncXEv1[t4c[i, idxPos1]] - vfuncXEv1[t4c[i, idxNeg1]])

        # difference vector
        dv1 = v4c[t4c[i, idxPos1], :] - v4c[t4c[i, idxNeg1], :]

        # compute new point
        pt1 = sf1 * dv1 + v4c[t4c[i, idxNeg1], :]

        # --------------------------------------------------------------
        # store and connect

        # check for duplicate points (pt1)
        if newVtcsAdj[t4c[i, idxPos1], t4c[i, idxNeg1]] == 0:
            # add to list
            v4c = np.concatenate((v4c, np.array(pt1, ndmin=2)))
            # keep index
            idxPt1 = len(v4c) - 1
            # update adj mtx
            newVtcsAdj[t4c[i, idxPos1], t4c[i, idxNeg1]] = idxPt1
            newVtcsAdj[t4c[i, idxNeg1], t4c[i, idxPos1]] = idxPt1
            # also append i4c and k4c
            i4c = np.concatenate((i4c, st.mode(i4c[t4c[i,:]])[0]))
            k4c = np.concatenate((k4c, st.mode(k4c[t4c[i,:]])[0]))
            # let's remember the new indices, and whether or not it's a
            # 'positive' or 'negative' tetra (in terms of the 2nd
            # eigenfunction)
            newVtcs = np.concatenate((newVtcs, [idxPt1]))
            tmpSgn = np.unique(vfuncXEv2[t4c[i,:]][np.logical_not(np.isnan(vfuncXEv2[t4c[i,:]]))])
            if len(tmpSgn)==1:
                newVtcsSgn = np.concatenate((newVtcsSgn, tmpSgn))
            else:
                print("Inconsistency while creating seam, exiting.")
                sys.exit(1)
        else:
            idxPt1 = newVtcsAdj[t4c[i, idxPos1], t4c[i, idxNeg1]]

        # create new tetras
        # pos1, pt1, nan1, nan2
        # neg1, pt1, nan1, nan2
        newTetra = np.concatenate((newTetra, np.array((t4c[i, idxPos1], idxPt1, t4c[i, idxNan1], t4c[i, idxNan2]), ndmin=2)))
        newTetra = np.concatenate((newTetra, np.array((t4c[i, idxNeg1], idxPt1, t4c[i, idxNan1], t4c[i, idxNan2]), ndmin=2)))

        # return
        return v4c, t4c, i4c, k4c, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn

    # --------------------------------------------------------------------------
    # getSeam subfunction #3

    # case 3: zero nans in the tetra (could happen for boundary tetras
    # towards two sides)
    def getSeamCase3(v4c, t4c, i4c, k4c, vfuncXEv1, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn, i):

        # determine which edges of the current tetra belong to the
        # boundary surface
        numEdgeBnd = np.where((
            np.logical_or(t4c[i, 0] == np.unique(e4cBndOpen, axis=0), t4c[i, 1] == np.unique(e4cBndOpen, axis=0)).all(axis=1).any(),
            np.logical_or(t4c[i, 0] == np.unique(e4cBndOpen, axis=0), t4c[i, 2] == np.unique(e4cBndOpen, axis=0)).all(axis=1).any(),
            np.logical_or(t4c[i, 0] == np.unique(e4cBndOpen, axis=0), t4c[i, 3] == np.unique(e4cBndOpen, axis=0)).all(axis=1).any(),
            np.logical_or(t4c[i, 1] == np.unique(e4cBndOpen, axis=0), t4c[i, 2] == np.unique(e4cBndOpen, axis=0)).all(axis=1).any(),
            np.logical_or(t4c[i, 1] == np.unique(e4cBndOpen, axis=0), t4c[i, 3] == np.unique(e4cBndOpen, axis=0)).all(axis=1).any(),
            np.logical_or(t4c[i, 2] == np.unique(e4cBndOpen, axis=0), t4c[i, 3] == np.unique(e4cBndOpen, axis=0)).all(axis=1).any()))[0]

        # number of trias that belong to the surface (and which ones)
        numTriaBnd = np.where((
            np.logical_or(np.logical_or(t4c[i, 0] == t4cBndOpen, t4c[i, 1] == t4cBndOpen), t4c[i, 2] == t4cBndOpen).all(axis=1).any(),
            np.logical_or(np.logical_or(t4c[i, 0] == t4cBndOpen, t4c[i, 1] == t4cBndOpen), t4c[i, 3] == t4cBndOpen).all(axis=1).any(),
            np.logical_or(np.logical_or(t4c[i, 0] == t4cBndOpen, t4c[i, 2] == t4cBndOpen), t4c[i, 3] == t4cBndOpen).all(axis=1).any(),
            np.logical_or(np.logical_or(t4c[i, 1] == t4cBndOpen, t4c[i, 2] == t4cBndOpen), t4c[i, 3] == t4cBndOpen).all(axis=1).any()))[0]

        # examine values: we can only cut trias or edges that have at
        # least one pos and one neg value

        # number of edges that have pos / neg values (and which ones)
        numEdgeCross = np.where(np.sum(np.sign(vfuncXEv1[t4c[i, ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))]]), axis=1) == 0)[0]

        # number of trias that have pos / neg values (and which ones)
        numTriaCross = np.where(np.abs(np.sum(np.sign(vfuncXEv1[t4c[i, ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))]]), axis=1)) == 1)[0]

        #
        #print(
        #    "surface edges: " + str(len(numEdgeBnd)) + ", " +
        #    "surface trias: " + str(len(numTriaBnd)) + ", " +
        #    "pos/neg edges: " + str(len(numEdgeCross)) + ", " +
        #    "pos/neg trias: " + str(len(numTriaCross)) + ", " +
        #    "pos/neg surface edges: " + str(len(np.intersect1d(numEdgeBnd, numEdgeCross))) + ", " +
        #    "pos/neg surface trias: " + str(len(np.intersect1d(numTriaBnd, numTriaCross))))

        # we now analyze, based on the above info, which case the current one
        # actually is, i.e. ignoring the vertex or edge that is at the opposite
        # site etc.

        # distinguish cases
        if len(numEdgeBnd) == 2 and len(numTriaBnd) == 0 and len(numEdgeCross) == 3 and len(numTriaCross) == 3 and len(np.intersect1d(numEdgeBnd, numEdgeCross)) == 1 and len(np.intersect1d(numTriaBnd, numTriaCross)) == 0:
            # surface edges: 2, surface trias: 0, pos/neg edges: 3, pos/neg trias: 3, pos/neg surface edges: 1, pos/neg surface trias: 0
            # should be case 2: two nans in tetra, due to only one single pos/neg surface edge
            # now find which edge is on the surface, and set the other points
            # temporarily to nan
            tmpvfuncXEv1 = vfuncXEv1.copy()
            tmpvfuncXEv1[t4c[i, np.setdiff1d((0, 1, 2, 3), ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))[np.intersect1d(numEdgeBnd, numEdgeCross).item()])]] = np.nan
            v4c, t4c, i4c, k4c, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn = getSeamCase2(v4c, t4c, i4c, k4c, tmpvfuncXEv1, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn, i)
        elif len(numEdgeBnd) == 3 and len(numTriaBnd) == 1 and len(numEdgeCross) == 3 and len(numTriaCross) == 3 and len(np.intersect1d(numEdgeBnd, numEdgeCross)) == 2 and len(np.intersect1d(numTriaBnd, numTriaCross)) == 1:
            # surface edges: 3, surface trias: 1, pos/neg edges: 3, pos/neg trias: 3, pos/neg surface edges: 2, pos/neg surface trias: 1
            # should be case 1: one nan in tetra, due to one single surface tria
            # now find which tria is on the surface, and set the remaining vtx
            # temporarily to nan
            tmpvfuncXEv1 = vfuncXEv1.copy()
            tmpvfuncXEv1[t4c[i, np.setdiff1d((0, 1, 2, 3), ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))[np.intersect1d(numTriaBnd, numTriaCross).item()])]] = np.nan
            v4c, t4c, i4c, k4c, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn = getSeamCase1(v4c, t4c, i4c, k4c, tmpvfuncXEv1, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn, i)
        elif len(numEdgeBnd) == 3 and len(numTriaBnd) == 1 and len(numEdgeCross) == 4 and len(numTriaCross) == 4 and len(np.intersect1d(numEdgeBnd, numEdgeCross)) == 2 and len(np.intersect1d(numTriaBnd, numTriaCross)) == 1:
            # surface edges: 3, surface trias: 1, pos/neg edges: 4, pos/neg trias: 4, pos/neg surface edges: 2, pos/neg surface trias: 1
            # should be case 1: one nan in tetra, due to one single surface tria
            # now find which tria is on the surface, and set the remaining vtx
            # temporarily to nan - essentially the same case as above, but kept
            # separate for clarity
            tmpvfuncXEv1 = vfuncXEv1.copy()
            tmpvfuncXEv1[t4c[i, np.setdiff1d((0, 1, 2, 3), ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))[np.intersect1d(numTriaBnd, numTriaCross).item()])]] = np.nan
            v4c, t4c, i4c, k4c, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn = getSeamCase1(v4c, t4c, i4c, k4c, tmpvfuncXEv1, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn, i)
        else:
            print("Unknown getSeam() case, exiting.")
            sys.exit(1)

        # return
        return v4c, t4c, i4c, k4c, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn

    # --------------------------------------------------------------------------
    # main part of function
    # --------------------------------------------------------------------------

    # imports
    import os
    import sys

    import lapy as lp
    import numpy as np

    from shapetools.triaUtils import writeVTK, writePSOL
    from shapetools.triaUtils import tetra_get_boundary_tria, tria_rm_free_vertices

    from scipy import sparse as sp
    from scipy import stats as st

    # initialize
    newVtcs = np.empty(shape=(0), dtype=np.int)
    newVtcsSgn = np.empty(shape=(0), dtype=np.int)
    newVtcsAdj = sp.lil_matrix((len(v4c), len(v4c)), dtype=np.int)
    newTetra = np.empty(shape=(0, 4), dtype=np.int)

    # assign first eigenfunction
    vfuncXEv1 = np.zeros(np.shape(v4c)[0]) * np.nan
    vfuncXEv1[v4cBndOpenKeep] = anisoLaplEvec[:, 1]

    # assign second eigenfunction
    vfuncXEv2 = np.zeros(np.shape(v4c)[0]) * np.nan
    vfuncXEv2[v4cBndOpenKeep] = np.sign(anisoLaplEvec[:,2])

    # find tetra with zeros: first, transform boundary trias to edges; here,
    # the important trick is to start working with the boundary trias from
    # t4cBndOpen, which is to avoid trias that have e.g. 2 points at the upper
    # surface and 1 point on the lower surface (which can happen). Then
    # transform to edges and identify those with pos and neg values. We use
    # edges because these need to be cut also, not only trias, and because a
    # tetra can also just share an edge with the surface, not necessarily a
    # tria. Finally select all tetras that include one such edge.
    e4cBndOpen = np.concatenate((t4cBndOpen[:,(0,1)], t4cBndOpen[:,(0,2)], t4cBndOpen[:,(1,2)]), axis=0)
    e4cBndOpen = np.unique(np.sort(e4cBndOpen, axis=1), axis=1)
    zeroEdgeIdx = e4cBndOpen[np.where(np.abs(np.sum(np.sign(vfuncXEv1[e4cBndOpen]), axis=1))==0)[0]]
    zeroTetraIdx = list()
    for i in zeroEdgeIdx:
        zeroTetraIdx.extend(np.where(np.sum(np.isin(t4c, i), axis=1)==2)[0].tolist())
    zeroTetraIdx = np.unique(np.array(zeroTetraIdx))

    # from lapy import Plot as lpp
    #
    # tetMesh = lp.TetMesh(v4c, t4c)
    #
    # vfunc = np.zeros(np.shape(v4c)[0])
    # vfunc[np.unique(t4c[zeroTetraIdx])] = 1
    #
    # lpp.plot_tet_mesh(tetMesh, html_output=True, vfunc=vfunc, plot_edges=True, width=1200, height=1200)

    # loop across tetras with zeros (locally disable warning that np.nan < 0 is
    # False and np.nan > 0 is also False, which is intended)
    with np.errstate(invalid='ignore'):

        for i in zeroTetraIdx:

            # we can typically have one or two nans in the tetra. we may
            # sometimes have zero nans in the tetra.
            if np.sum(np.isnan(vfuncXEv1[t4c[i]])) == 1:
                # if there is one nan, we need to make sure that the three vtcs form a surface triangle; alternatively we check if it is at least a surface edge
                if np.isin(t4cBndOpen, t4c[i, np.logical_not(np.isnan(vfuncXEv1[t4c[i]]))]).all(axis=1).any():
                    v4c, t4c, i4c, k4c, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn = getSeamCase1(
                    v4c, t4c, i4c, k4c, vfuncXEv1, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn, i)
                elif np.isin(e4cBndOpen, t4c[i, np.logical_not(np.isnan(vfuncXEv1[t4c[i]]))]).all(axis=1).any():
                    # need to temporarily set one vtx to nan
                    tmpEdge = (
                        np.isin(e4cBndOpen, t4c[i, (0, 1)]).all(axis=1).any(),
                        np.isin(e4cBndOpen, t4c[i, (0, 2)]).all(axis=1).any(),
                        np.isin(e4cBndOpen, t4c[i, (0, 3)]).all(axis=1).any(),
                        np.isin(e4cBndOpen, t4c[i, (1, 2)]).all(axis=1).any(),
                        np.isin(e4cBndOpen, t4c[i, (1, 3)]).all(axis=1).any(),
                        np.isin(e4cBndOpen, t4c[i, (2, 3)]).all(axis=1).any())
                    if len(np.where(tmpEdge)) > 1:
                        print("Cave: found " + str(len(np.where(tmpEdge))) + " edges, exiting.")
                        sys.exit(1)
                    tmpEdgeNot = np.setdiff1d((0, 1, 2, 3),  ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))[np.where(tmpEdge)[0].item()])
                    tmpvfuncXEv1 = vfuncXEv1.copy()
                    tmpvfuncXEv1[t4c[i, tmpEdgeNot]] = np.nan
                    v4c, t4c, i4c, k4c, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn = getSeamCase2(
                    v4c, t4c, i4c, k4c, tmpvfuncXEv1, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn, i)
                else:
                    print("Vtx " + str(i) + " not in triangle or edge, skipping")
                    sys.exit(1)

            elif np.sum(np.isnan(vfuncXEv1[t4c[i]])) == 2:
                # if there are two nans, we need to make sure that they are linked by a surface edge; this should have happenend already.
                tmpEdge = (
                    np.isin(e4cBndOpen, t4c[i, (0, 1)]).all(axis=1).any(),
                    np.isin(e4cBndOpen, t4c[i, (0, 2)]).all(axis=1).any(),
                    np.isin(e4cBndOpen, t4c[i, (0, 3)]).all(axis=1).any(),
                    np.isin(e4cBndOpen, t4c[i, (1, 2)]).all(axis=1).any(),
                    np.isin(e4cBndOpen, t4c[i, (1, 3)]).all(axis=1).any(),
                    np.isin(e4cBndOpen, t4c[i, (2, 3)]).all(axis=1).any())
                if len(np.where(tmpEdge)) > 1:
                    print("Cave: found " + str(len(np.where(tmpEdge))) + " edges, exiting.")
                    sys.exit(1)
                v4c, t4c, i4c, k4c, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn = getSeamCase2(
                v4c, t4c, i4c, k4c, vfuncXEv1, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn, i)

            elif np.sum(np.isnan(vfuncXEv1[t4c[i]])) == 0:
                # if there are no nans, the checking is done within the getSeamCase3 function
                v4c, t4c, i4c, k4c, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn = getSeamCase3(
                v4c, t4c, i4c, k4c, vfuncXEv1, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn, i)

            else:
                # we should never have four nans
                print("Inconsistency while creating seam (incorrect number of nans), exiting.")
                sys.exit(1)

    # remove tetras
    t4c = np.delete(t4c, zeroTetraIdx, axis=0)

    # add new tetras
    t4c = np.concatenate((t4c, newTetra))

    # return
    return v4c, t4c, i4c, k4c, newVtcs, newVtcsSgn, newTetra

# ------------------------------------------------------------------------------
# my_fem_tria_aniso()

def my_fem_tria_aniso(tria, u1, u2, c1, c2, mean_curv=None, wtMeanCurv=0.0, alpha=0.0, lump=False, setMinCurvToZero=False, setMaxCurvToZero=False):
    """
    computeABtria(v,t) computes the two sparse symmetric matrices representing
           the Laplace Beltrami Operator for a given triangle mesh using
           the linear finite element method (assuming a closed mesh or
           the Neumann boundary condition).

    Inputs:   v  - vertices : list of lists of 3 floats
              t  - triangles: N list of lists of 3 int of indices (>=0) into v array
              u1 - min curv:  min curvature direction per triangle (Nx3 floats)
              u2 - max curv:  max curvature direction per triangle (Nx3 floats)
              aniso_mat  - anisotropy matrix: diagonal elements in u1,u2 basis per
                              triangle (Nx2 floats)

    Outputs:  A  - sparse sym. (n x n) positive semi definite numpy matrix
              B  - sparse sym. (n x n) positive definite numpy matrix (inner product)

    Can be used to solve sparse generalized Eigenvalue problem: A x = lambda B x
    or to solve Poisson equation: A x = B f (where f is function on mesh vertices)
    or to solve Laplace equation: A x = 0
    or to model the operator's action on a vector x:   y = B\(Ax)
    """
    import sys
    import numpy as np
    from scipy import sparse
    # Compute vertex coordinates and a difference vector for each triangle:
    t1 = tria.t[:, 0]
    t2 = tria.t[:, 1]
    t3 = tria.t[:, 2]
    v1 = tria.v[t1, :]
    v2 = tria.v[t2, :]
    v3 = tria.v[t3, :]
    v2mv1 = v2 - v1
    v3mv2 = v3 - v2
    v1mv3 = v1 - v3
    # transform edge e to basis U = (U1,U2) via U^T * e
    # Ui is n x 3, e is n x 1, result is n x 2
    uv2mv1 = np.column_stack((np.sum(u1 * v2mv1, axis=1), np.sum(u2 * v2mv1, axis=1)))
    uv3mv2 = np.column_stack((np.sum(u1 * v3mv2, axis=1), np.sum(u2 * v3mv2, axis=1)))
    uv1mv3 = np.column_stack((np.sum(u1 * v1mv3, axis=1), np.sum(u2 * v1mv3, axis=1)))
    # combined aniso mat
    if wtMeanCurv > 0.0:
        # normalize curv to -1 ... +1 interval
        norm_mean_curv  = mean_curv / np.max(np.abs((mean_curv)))
        # compute new weighting parameter
        beta = 1 + wtMeanCurv * norm_mean_curv
        # compute aniso
        aniso_mat = np.empty((tria.t.shape[0], 2))
        if setMinCurvToZero:
            aniso_mat[:, 1] = 1
        else:
            aniso_mat[:, 1] = np.exp(-alpha * np.abs(c1) - wtMeanCurv * mean_curv)
        if setMaxCurvToZero:
            aniso_mat[:, 0] = 1
        else:
            aniso_mat[:, 0] = np.exp(-alpha * np.abs(c2) - wtMeanCurv * mean_curv)
    else:
        # compute aniso
        aniso_mat = np.empty((tria.t.shape[0], 2))
        if setMinCurvToZero:
            aniso_mat[:, 1] = 1
        else:
            aniso_mat[:, 1] = np.exp(-alpha * np.abs(c1))
        if setMaxCurvToZero:
            aniso_mat[:, 0] = 1
        else:
            aniso_mat[:, 0] = np.exp(-alpha * np.abs(c2))
    # Compute cross product and 4*vol for each triangle:
    cr = np.cross(v3mv2, v1mv3)
    vol = 2 * np.sqrt(np.sum(cr * cr, axis=1))
    # zero vol will cause division by zero below, so set to small value:
    vol_mean = 0.0001 * np.mean(vol)
    vol[vol < sys.float_info.epsilon] = vol_mean
    # compute cotangents for A
    # using that v2mv1 = - (v3mv2 + v1mv3) this can also be seen by
    # summing the local matrix entries in the old algorithm
    # Also: here aniso_mat is the two diagonal entries, not full matrices
    a12 = np.sum(uv3mv2 * aniso_mat * uv1mv3, axis=1) / vol
    a23 = np.sum(uv1mv3 * aniso_mat * uv2mv1, axis=1) / vol
    a31 = np.sum(uv2mv1 * aniso_mat * uv3mv2, axis=1) / vol
    # compute diagonals (from row sum = 0)
    a11 = -a12 - a31
    a22 = -a12 - a23
    a33 = -a31 - a23
    # stack columns to assemble data
    local_a = np.column_stack((a12, a12, a23, a23, a31, a31, a11, a22, a33)).reshape(-1)
    i = np.column_stack((t1, t2, t2, t3, t3, t1, t1, t2, t3)).reshape(-1)
    j = np.column_stack((t2, t1, t3, t2, t1, t3, t1, t2, t3)).reshape(-1)
    # Construct sparse matrix:
    # a = sparse.csr_matrix((local_a, (i, j)))
    a = sparse.csc_matrix((local_a, (i, j)), dtype=np.float32)
    if not lump:
        # create b matrix data (account for that vol is 4 times area)
        b_ii = vol / 24
        b_ij = vol / 48
        local_b = np.column_stack((b_ij, b_ij, b_ij, b_ij, b_ij, b_ij, b_ii, b_ii, b_ii)).reshape(-1)
        b = sparse.csc_matrix((local_b, (i, j)), dtype=np.float32)
    else:
        # when lumping put all onto diagonal  (area/3 for each vertex)
        b_ii = vol / 12
        local_b = np.column_stack((b_ii, b_ii, b_ii)).reshape(-1)
        i = np.column_stack((t1, t2, t3)).reshape(-1)
        b = sparse.csc_matrix((local_b, (i, i)), dtype=np.float32)
    return a, b

# ------------------------------------------------------------------------------
# MAIN FUNCTION
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# computeCubeParam

def computeCubeParam(params, cutTetraMeshDir=None, cutTetraMeshFile=None,
    cutTetraIndexFile=None, openBndCutTetraMeshFile=None, tetraIndexFile=None,
    hemi=None, outputDir=None, paramAddNoise = False, paramLegacy = True,
    paramAnisoSmooth = 3, paramAnisoAlpha = None, paramWeightMeanCurv = 0.0,
    labelPrsbc=234, labelSbc=236, labelCA1=238, labelCA3=240, labelBndCA4=2420,
    labelTail=226, labelHead=232):

    # --------------------------------------------------------------------------
    # computeCubeParam subfunctions
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    # options_parse()

    def options_parse():
        """
        Command Line Options Parser:
        initiate the option parser and return the parsed object
        """

        # imports

        import sys
        import optparse

        # define helptext
        HELPTEXT = """
        SUMMARY

        This is an auxiliary script for the shapetools.py script and is usually
        called from within that script.

        The script requires four arguments:

        --cutTetraMeshDir <cutTetraMeshDir>
        --cutTetraMeshFile <cutTetraMeshFile>
        --cutTetraIndexFile <cutTetraIndexFile>
        --openBndCutTetraMeshFile <openBndCutTetraMeshFile>
        --tetraIndexFile <tetraIndexFile>
        --hemi <lh|rh>

        """

        # message
        print("\nReading input options ...")

        # initialize
        parser = optparse.OptionParser(usage=HELPTEXT)

        # help text
        h_cutTetraMeshDir = 'cutTetraMeshDir'
        h_cutTetraMeshFile = 'cutTetraMeshFile'
        h_cutTetraIndexFile = 'cutTetraIndexFile'
        h_openBndCutTetraMeshFile = 'openBndCutTetraMeshFile'
        h_tetraIndexFile = 'tetraIndexFile'
        h_hemi = 'hemisphere'

        # specify inputs
        group = optparse.OptionGroup(parser, "Required Options:", "...")
        group.add_option('--cutTetraMeshDir', dest='cutTetraMeshDir', help=h_cutTetraMeshDir)
        group.add_option('--cutTetraMeshFile', dest='cutTetraMeshFile', help=h_cutTetraMeshFile)
        group.add_option('--cutTetraIndexFile', dest='cutTetraIndexFile', help=h_cutTetraIndexFile)
        group.add_option('--openBndCutTetraMeshFile', dest='openBndCutTetraMeshFile', help=h_openBndCutTetraMeshFile)
        group.add_option('--tetraIndexFile', dest='tetraIndexFile', help=h_tetraIndexFile)
        group.add_option('--hemi', dest='hemi', help=h_hemi)
        parser.add_option_group(group)

        # parse arguments
        (options, args) = parser.parse_args()

        # check if there are any inputs
        if len(sys.argv) != 13:
            print(HELPTEXT)
            sys.exit(0)

        # check if cutTetraMeshDir file is given
        if options.cutTetraMeshDir is None:
            print('\nERROR: Specify --cutTetraMeshDir\n')
            sys.exit(1)
        else:
            print('... Found cutTetraMeshDir file ' + options.cutTetraMeshDir)

        # check if cutTetraMeshFile file is given
        if options.cutTetraMeshFile is None:
            print('\nERROR: Specify --cutTetraMeshFile\n')
            sys.exit(1)
        else:
            print('... Found cutTetraMeshFile file ' + options.cutTetraMeshFile)

        # check if cutTetraIndexFile file is given
        if options.cutTetraIndexFile is None:
            print('\nERROR: Specify --cutTetraIndexFile\n')
            sys.exit(1)
        else:
            print('... Found cutTetraIndexFile file ' + options.cutTetraIndexFile)

        # check if openBndCutTetraMeshFile file is given
        if options.openBndCutTetraMeshFile is None:
            print('\nERROR: Specify --openBndCutTetraMeshFile\n')
            sys.exit(1)
        else:
            print('... Found openBndCutTetraMeshFile file ' + options.openBndCutTetraMeshFile)

        # check if tetraIndexFile file is given
        if options.tetraIndexFile is None:
            print('\nERROR: Specify --tetraIndexFile\n')
            sys.exit(1)
        else:
            print('... Found tetraIndexFile file ' + options.tetraIndexFile)

        # check if hemi is given
        if options.hemi is None:
            print('\nERROR: Specify --hemi\n')
            sys.exit(1)
        else:
            print('... Found hemisphere ' + options.hemi)

        # return
        return options

    # --------------------------------------------------------------------------
    # main function part
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    # import

    import os, sys

    import numpy as np

    import lapy as lp
    import lapy.DiffGeo as DiffGeo

    from scipy.sparse.linalg import LinearOperator, eigsh, splu

    from shapetools.triaUtils import readVTK, writeVTK, readPSOL, writePSOL
    from shapetools.triaUtils import tria_rm_free_vertices, tetra_get_boundary_tria

    # --------------------------------------------------------------------------
    # message

    print()
    print("-------------------------------------------------------------------------")
    print()
    print("Computing cube parametrization")
    print()
    print("-------------------------------------------------------------------------")
    print()

    # --------------------------------------------------------------------------
    # get data

    # evaluate params

    if params is not None:

        cutTetraMeshDir = os.path.join(params.OUTDIR, 'tetra-cube')
        cutTetraMeshFile = os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_09 + ".vtk")
        cutTetraIndexFile = os.path.join(params.OUTDIR, 'tetra-cut', params.HEMI + "." + params.internal.HSFLABEL_09 + ".psol")
        openBndCutTetraMeshFile = os.path.join(params.OUTDIR, 'tetra-cut', params.HEMI + ".open.bnd.cut.tetra.vtk")
        tetraIndexFile = os.path.join(params.OUTDIR, 'tetra-labels', params.HEMI + "." + params.internal.HSFLABEL_07 + ".psol")

        hemi = params.HEMI

        paramAnisoAlpha = params.internal.anisoAlpha
        paramAnisoSmooth = params.internal.anisoSmooth
        paramWeightMeanCurv = params.internal.weightMeanCurv
        paramAddNoise = params.internal.cubeAddNoise
        paramLegacy = params.internal.cubeWriteLegacyVTK

        if params.LUT == "fs711":
            labelPrsbc = 234
            labelSbc = 236
            labelCA1 = 238
            labelCA3 = 240
            labelBndCA4 = 2420
            labelTail = 226
            labelHead = 232

        elif params.LUT == "ashs":
            labelPrsbc = 8
            labelSbc = 8
            labelCA1 = 1
            labelCA3 = 4
            labelBndCA4 = 4
            labelTail = 254
            labelHead = 255

    # load cut tetra mesh

    v4c, t4c = readVTK(cutTetraMeshFile)

    # load indices of cutting planes

    i4c = readPSOL(cutTetraIndexFile)

    # load subfields indices

    j4c = readPSOL(tetraIndexFile)

    # append subfield indices

    k4c = np.concatenate((j4c, np.zeros(len(i4c)-len(j4c))))

    # --------------------------------------------------------------------------
    # get a tria mesh that is cut open at its ends, but still consistent with
    # the tetra mesh

    v4cBndOpen, t4cBndOpen = readVTK(openBndCutTetraMeshFile)

    # remove free vertices and orient

    v4cBndOpenRm, t4cBndOpenRm, v4cBndOpenKeep, v4cBndOpenDel = tria_rm_free_vertices(v4cBndOpen, t4cBndOpen)

    triaMesh4cBndOpenRm = lp.TriaMesh(v4cBndOpenRm, t4cBndOpenRm)

    triaMesh4cBndOpenRm.orient_()

    t4cBndOpenRm = triaMesh4cBndOpenRm.t

    # --------------------------------------------------------------------------
    # add noise (to avoid some numerical problems, especially for artifical
    # data; option is disabled by default)

    if paramAddNoise is True:

        v4cBndOpenRm = v4cBndOpenRm + (np.random.random_sample(v4cBndOpenRm.shape)-0.5)/1000

    # --------------------------------------------------------------------------
    # compute anisotropic laplace

    if paramWeightMeanCurv == 0.0:

        # compute first eigenfunction

        fem = lp.Solver(triaMesh4cBndOpenRm, lump=True, aniso=paramAnisoAlpha, aniso_smooth=paramAnisoSmooth)

        anisoLaplEval, anisoLaplEvec = fem.eigs(k=3)

    elif paramWeightMeanCurv > 0.0:

        # compute first eigenfunction with modified anisotropy matrix

        u1, u2, c1, c2 = triaMesh4cBndOpenRm.curvature_tria(smoothit=paramAnisoSmooth)

        mean_curv = (c1 + c2) / 2
        gauss_curv = c1 * c2

        writePSOL(os.path.join(cutTetraMeshDir, hemi + '.lapy.mean-curv.psol'), triaMesh4cBndOpenRm.map_tfunc_to_vfunc(mean_curv))
        writePSOL(os.path.join(cutTetraMeshDir, hemi + '.lapy.gauss-curv.psol'), triaMesh4cBndOpenRm.map_tfunc_to_vfunc(gauss_curv))

        #

        fem = lp.Solver(triaMesh4cBndOpenRm)

        fem.stiffness, fem.mass = my_fem_tria_aniso(triaMesh4cBndOpenRm, u1, u2, c1, c2, mean_curv=mean_curv, wtMeanCurv=paramWeightMeanCurv, alpha=paramAnisoAlpha, lump=True)

        fem.geotype = type(triaMesh4cBndOpenRm)

        anisoLaplEval, anisoLaplEvec = fem.eigs(k=3)

    # --------------------------------------------------------------------------
    # get mapping of open boundary vertices to subfields

    with open(os.path.join(os.path.dirname(openBndCutTetraMeshFile), hemi + '.rm.open.bnd.cut.tetra.lst'), "r") as f:

        hsfList = f.read().splitlines()

    hsfList = np.array(hsfList).astype(np.int)

    writePSOL(os.path.join(os.path.dirname(openBndCutTetraMeshFile), hemi + '.hsf.rm.open.bnd.cut.tetra.psol'), k4c[hsfList])

    # --------------------------------------------------------------------------
    # post-process eigenfunction (order, flipping)

    # decide whether or not to flip anisoLaplEvec[:, 1] and anisoLaplEvec[:, 2]:
    # anisoLaplEvec[:, 1] should be zero at the extrema of 234 and 240.
    # anisoLaplEvec[:, 2] should have extremal values at 234 and 240.
    # We therefore check if the zeros of Evec1 are located in 234 or 240/2420
    # rather than 236 or 238. If not, we change order.

    # assign first eigenfunction
    vfuncXEv1 = np.zeros(np.shape(v4c)[0]) * np.nan
    vfuncXEv1[v4cBndOpenKeep] = anisoLaplEvec[:, 1]
    # get edges that contain zeros, and turn into list of vertices
    e4cBndOpen = np.concatenate((t4cBndOpen[:,(0,1)], t4cBndOpen[:,(0,2)], t4cBndOpen[:,(1,2)]), axis=0)
    e4cBndOpen = np.unique(np.sort(e4cBndOpen, axis=1), axis=1)
    zeroEdgeIdx = e4cBndOpen[np.where(np.abs(np.sum(np.sign(vfuncXEv1[e4cBndOpen]), axis=1))==0)[0]]
    zeroVtxIdx = np.unique(zeroEdgeIdx)

    if np.sum(k4c[zeroVtxIdx]==labelPrsbc)>0 and np.sum(np.logical_or(k4c[zeroVtxIdx]==labelCA3, k4c[zeroVtxIdx]==labelBndCA4))>0 and (np.sum(np.logical_or(k4c[zeroVtxIdx]==labelCA3, k4c[zeroVtxIdx]==labelBndCA4)) > np.sum(np.logical_or(k4c[zeroVtxIdx]==labelSbc, k4c[zeroVtxIdx]==labelCA1))):
         print("Not necessary to change order of EV1 and EV2")
    else:
         print("Changing order of EV1 and EV2")
         #import pdb; pdb.set_trace()
         anisoLaplEvec = anisoLaplEvec[:, (0, 2, 1)]

    # decide whether or not to flip anisoLaplEvec[:, 1]  (should be inf -> sup)

    # check EV1 for flips; this is done indirectly: we get the locations
    # of highest and lowest values in 234 (PrSbc); highest values should have a
    # higher value on the z-axis than lower values; note: could do the same
    #  (i.e., reverse) for 240 as well

    if np.median(v4cBndOpenRm[np.where(np.logical_and(anisoLaplEvec[:,1]> 0, k4c[hsfList]==labelPrsbc))[0], 2]) > np.median(v4cBndOpenRm[np.where(np.logical_and(anisoLaplEvec[:,1]<0, k4c[hsfList]==labelPrsbc))[0], 2]):
        ("No flip necessary for EV1")
    elif np.median(v4cBndOpenRm[np.where(np.logical_and(anisoLaplEvec[:,1]> 0, k4c[hsfList]==labelPrsbc))[0], 2]) < np.median(v4cBndOpenRm[np.where(np.logical_and(anisoLaplEvec[:,1]<0, k4c[hsfList]==labelPrsbc))[0], 2]):
        print("Flipping EV1")
        anisoLaplEvec[:, 1] = -anisoLaplEvec[:, 1]
    else:
        print("Inconsistency detected for EV1, exiting.")
        #import pdb; pdb.set_trace()
        sys.exit(1)

    # decide whether or not to flip anisoLaplEvec[:, 2] (should be 234 -> 240)

    if np.median(anisoLaplEvec[np.where(k4c[hsfList]==labelPrsbc)[0], 2])<0 and np.median(anisoLaplEvec[np.where(np.logical_or(k4c[hsfList]==labelCA3, k4c[hsfList]==labelBndCA4))[0], 2])>0:
        print("No flip necessary for EV2")
    elif np.median(anisoLaplEvec[np.where(k4c[hsfList]==labelPrsbc)[0], 2])>0 and np.median(anisoLaplEvec[np.where(np.logical_or(k4c[hsfList]==labelCA3, k4c[hsfList]==labelBndCA4))[0], 2])<0:
        print("Flipping EV2")
        anisoLaplEvec[:, 2] = -anisoLaplEvec[:, 2]
    else:
        print("Inconsistency detected for EV2, exiting.")
        #import pdb; pdb.set_trace()
        sys.exit(1)

    writePSOL(os.path.join(cutTetraMeshDir, hemi + '.lapy.aLBO.EV1.psol'), anisoLaplEvec[:,1])
    writePSOL(os.path.join(cutTetraMeshDir, hemi + '.lapy.aLBO.EV2.psol'), anisoLaplEvec[:,2])

    # --------------------------------------------------------------------------
    # define boundary conditions (note that we changed dimensions and directions!)

    # get seam

    v4c, t4c, i4c, k4c, newVtcs, newVtcsSgn, newTetra = getSeam(v4c, t4c, i4c, k4c, v4cBndOpenKeep, t4cBndOpen, anisoLaplEvec)

    # x: 234 -> 240

    vfuncX = np.zeros(np.shape(v4c)[0])
    srcX = np.zeros(np.shape(v4c)[0])
    srcX[newVtcs[newVtcsSgn>0]] = 1
    snkX = np.zeros(np.shape(v4c)[0])
    snkX[newVtcs[newVtcsSgn<0]] = -1
    vfuncX = srcX + snkX

    # y: Tail (226) -> Head (232) (post --> ant)

    vfuncY = np.zeros(np.shape(v4c)[0])
    vfuncY[np.where(i4c==labelTail)[0]] = -1
    vfuncY[np.where(i4c==labelHead)[0]] = 1

    # z: inf->sup

    vfuncZ = np.zeros(np.shape(v4c)[0])
    srcZ = np.zeros(np.shape(anisoLaplEvec[:,1])[0])
    srcZ[anisoLaplEvec[:,1]<0] = -1
    snkZ = np.zeros(np.shape(anisoLaplEvec[:,1])[0])
    snkZ[anisoLaplEvec[:,1]>0] = 1
    vfuncZ[v4cBndOpenKeep] = srcZ + snkZ

    # --------------------------------------------------------------------------
    # compute ABtetra and Poisson functions

    # first, need to remove free vtcs

    v4cRm, t4cRm, v4cRmKeep, v4cRmDel = tria_rm_free_vertices(v4c,t4c)

    writeVTK(outfile=os.path.join(cutTetraMeshDir, hemi + '.seam.rm.cut.tetra.vtk'), v=v4cRm, t=t4cRm)

    #

    vfuncXRm = vfuncX[v4cRmKeep]
    vfuncYRm = vfuncY[v4cRmKeep]
    vfuncZRm = vfuncZ[v4cRmKeep]

    #

    bndXRm = (np.where(vfuncXRm)[0], vfuncXRm[np.where(vfuncXRm)])
    bndYRm = (np.where(vfuncYRm)[0], vfuncYRm[np.where(vfuncYRm)])
    bndZRm = (np.where(vfuncZRm)[0], vfuncZRm[np.where(vfuncZRm)])

    #

    tetMesh4cRm = lp.TetMesh(v4cRm, t4cRm)

    fem = lp.Solver(tetMesh4cRm)

    P0 = fem.poisson(dtup=bndXRm)
    P1 = fem.poisson(dtup=bndYRm)
    P2 = fem.poisson(dtup=bndZRm)

    # --------------------------------------------------------------------------
    # output

    # write out functions

    writePSOL(os.path.join(cutTetraMeshDir, hemi + '.vfuncX.seam.rm.cut.tetra.psol'), vfuncXRm)
    writePSOL(os.path.join(cutTetraMeshDir, hemi + '.vfuncY.seam.rm.cut.tetra.psol'), vfuncYRm)
    writePSOL(os.path.join(cutTetraMeshDir, hemi + '.vfuncZ.seam.rm.cut.tetra.psol'), vfuncZRm)

    writePSOL(os.path.join(cutTetraMeshDir, hemi + '.poisson0.seam.rm.cut.tetra.psol'), P0)
    writePSOL(os.path.join(cutTetraMeshDir, hemi + '.poisson1.seam.rm.cut.tetra.psol'), P1)
    writePSOL(os.path.join(cutTetraMeshDir, hemi + '.poisson2.seam.rm.cut.tetra.psol'), P2)

    # also write out boundary files for visualization

    t4cRmBnd = tetra_get_boundary_tria(v4cRm, t4cRm)
    v4cRmBndRm, t4cRmBndRm, v4cRmBndRmKeep, v4cRmBndRmDel = tria_rm_free_vertices(v4cRm,t4cRmBnd)

    writeVTK(outfile=os.path.join(cutTetraMeshDir, hemi + '.rm.bnd.seam.rm.cut.tetra.vtk'), v=v4cRmBndRm, t=t4cRmBndRm)

    writePSOL(os.path.join(cutTetraMeshDir, hemi + '.vfuncX.rm.bnd.seam.rm.cut.tetra.psol'), vfuncXRm[v4cRmBndRmKeep])
    writePSOL(os.path.join(cutTetraMeshDir, hemi + '.vfuncY.rm.bnd.seam.rm.cut.tetra.psol'), vfuncYRm[v4cRmBndRmKeep])
    writePSOL(os.path.join(cutTetraMeshDir, hemi + '.vfuncZ.rm.bnd.seam.rm.cut.tetra.psol'), vfuncZRm[v4cRmBndRmKeep])

    writePSOL(os.path.join(cutTetraMeshDir, hemi + '.poisson0.rm.bnd.seam.rm.cut.tetra.psol'), P0[v4cRmBndRmKeep])
    writePSOL(os.path.join(cutTetraMeshDir, hemi + '.poisson1.rm.bnd.seam.rm.cut.tetra.psol'), P1[v4cRmBndRmKeep])
    writePSOL(os.path.join(cutTetraMeshDir, hemi + '.poisson2.rm.bnd.seam.rm.cut.tetra.psol'), P2[v4cRmBndRmKeep])

    # write out cube

    w4cRm = np.vstack((P0,P1,P2)).transpose()

    writeVTK(outfile=os.path.join(cutTetraMeshDir, hemi + '.uvw.seam.rm.cut.tetra.vtk'), v=w4cRm, t=t4cRm)

    #

    return(params)


# ------------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# MAIN PART

if __name__ == "__main__":

    # Command line options and error checking

    options = options_parse()

    # Run analysis

    computeCubeParam(params=None,
        cutTetraMeshDir=options.cutTetraMeshDir,
        cutTetraMeshFile=options.cutTetraMeshFile,
        cutTetraIndexFile=options.cutTetraIndexFile,
        openBndCutTetraMeshFile=options.openBndCutTetraMeshFile,
        tetraIndexFile=options.tetraIndexFile, hemi=options.hemi,
        outputDir=options.outputdir, paramAnisoSmooth=options.anisoSmooth,
        paramAnisoAlpha=options.anisoAlpha,
        paramWeightMeanCurv=options.weightMeanCurv, paramAddNoise=False,
        paramLegacy=False)
