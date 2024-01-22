"""
This module provides a function to compute a cube parametrization

"""

import logging
import os

import nibabel as nb
import numpy as np
from lapy import Solver, TetMesh, TriaMesh, io
from scipy import sparse as sp
from scipy import stats as st
from sklearn.decomposition import PCA

# ==============================================================================
# LOGGING

LOGGER = logging.getLogger(__name__)

# ==============================================================================
# FUNCTIONS

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
        if np.sum(vfuncXEv1[t4c[i, :]] < 0) == 1:
            # get indices
            idxPos1 = np.where(vfuncXEv1[t4c[i, :]] > 0)[0][0]
            idxPos2 = np.where(vfuncXEv1[t4c[i, :]] > 0)[0][1]
            idxNeg1 = np.where(vfuncXEv1[t4c[i, :]] < 0)[0][0]
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
                i4c = np.concatenate((i4c, st.mode(i4c[t4c[i, :]])[0]))
                k4c = np.concatenate((k4c, st.mode(k4c[t4c[i, :]])[0]))
                # let's remember the new indices, and whether or not it's a
                # 'positive' or 'negative' tetra (in terms of the 2nd
                # eigenfunction)
                newVtcs = np.concatenate((newVtcs, [idxPt1]))
                tmpSgn = np.unique(vfuncXEv2[t4c[i, :]][np.logical_not(np.isnan(vfuncXEv2[t4c[i, :]]))])
                if len(tmpSgn) == 1:
                    newVtcsSgn = np.concatenate((newVtcsSgn, tmpSgn))
                else:
                    newVtcsSgn = np.concatenate((newVtcsSgn, np.array([0])))
                    raise RuntimeError("Inconsistency while creating seam 1, exiting.")
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
                i4c = np.concatenate((i4c, st.mode(i4c[t4c[i, :]])[0]))
                k4c = np.concatenate((k4c, st.mode(k4c[t4c[i, :]])[0]))
                # let's remember the new indices, and whether or not it's a
                # 'positive' or 'negative' tetra (in terms of the 2nd
                # eigenfunction)
                newVtcs = np.concatenate((newVtcs, [idxPt2]))
                tmpSgn = np.unique(vfuncXEv2[t4c[i, :]][np.logical_not(np.isnan(vfuncXEv2[t4c[i, :]]))])
                if len(tmpSgn) == 1:
                    newVtcsSgn = np.concatenate((newVtcsSgn, tmpSgn))
                else:
                    newVtcsSgn = np.concatenate((newVtcsSgn, np.array([0])))
                    raise RuntimeError("Inconsistency while creating seam 2, exiting.")
            else:
                idxPt2 = newVtcsAdj[t4c[i, idxPos2], t4c[i, idxNeg1]]

            # create new tetras
            # Neg1, pt1, pt2, nan
            # pt1, pt2, Pos1, nan
            # t2, Pos1, Pos2, nan

            newTetra = np.concatenate((newTetra, np.array((t4c[i, idxNeg1], idxPt1, idxPt2, t4c[i, idxNan1]), ndmin=2)))
            newTetra = np.concatenate((newTetra, np.array((idxPt1, idxPt2, t4c[i, idxPos1], t4c[i, idxNan1]), ndmin=2)))
            newTetra = np.concatenate(
                (newTetra, np.array((idxPt2, t4c[i, idxPos1], t4c[i, idxPos2], t4c[i, idxNan1]), ndmin=2))
            )

        # case 1b: two neg value
        elif np.sum(vfuncXEv1[t4c[i, :]] < 0) == 2:
            # get indices
            idxPos1 = np.where(vfuncXEv1[t4c[i, :]] > 0)[0][0]
            idxNeg1 = np.where(vfuncXEv1[t4c[i, :]] < 0)[0][0]
            idxNeg2 = np.where(vfuncXEv1[t4c[i, :]] < 0)[0][1]
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
                i4c = np.concatenate((i4c, st.mode(i4c[t4c[i, :]])[0]))
                k4c = np.concatenate((k4c, st.mode(k4c[t4c[i, :]])[0]))
                # let's remember the new indices, and whether or not it's a
                # 'positive' or 'negative' tetra (in terms of the 2nd
                # eigenfunction)
                newVtcs = np.concatenate((newVtcs, [idxPt1]))
                tmpSgn = np.unique(vfuncXEv2[t4c[i, :]][np.logical_not(np.isnan(vfuncXEv2[t4c[i, :]]))])
                if len(tmpSgn) == 1:
                    newVtcsSgn = np.concatenate((newVtcsSgn, tmpSgn))
                else:
                    newVtcsSgn = np.concatenate((newVtcsSgn, np.array([0])))
                    raise RuntimeError("Inconsistency while creating seam 3, exiting.")
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
                i4c = np.concatenate((i4c, st.mode(i4c[t4c[i, :]])[0]))
                k4c = np.concatenate((k4c, st.mode(k4c[t4c[i, :]])[0]))
                # let's remember the new indices, and whether or not it's a
                # 'positive' or 'negative' tetra (in terms of the 2nd
                # eigenfunction)
                newVtcs = np.concatenate((newVtcs, [idxPt2]))
                tmpSgn = np.unique(vfuncXEv2[t4c[i, :]][np.logical_not(np.isnan(vfuncXEv2[t4c[i, :]]))])
                if len(tmpSgn) == 1:
                    newVtcsSgn = np.concatenate((newVtcsSgn, tmpSgn))
                else:
                    newVtcsSgn = np.concatenate((newVtcsSgn, np.array([0])))
                    raise RuntimeError("Inconsistency while creating seam 4, exiting.")
            else:
                idxPt2 = newVtcsAdj[t4c[i, idxPos1], t4c[i, idxNeg2]]

            # create new tetras
            # pos1, pt1, pt2, nan
            # pt1, pt2, neg1, nan
            # pt2, neg1, neg2, nan

            newTetra = np.concatenate((newTetra, np.array((t4c[i, idxPos1], idxPt1, idxPt2, t4c[i, idxNan1]), ndmin=2)))
            newTetra = np.concatenate((newTetra, np.array((idxPt1, idxPt2, t4c[i, idxNeg1], t4c[i, idxNan1]), ndmin=2)))
            newTetra = np.concatenate(
                (newTetra, np.array((idxPt2, t4c[i, idxNeg1], t4c[i, idxNeg2], t4c[i, idxNan1]), ndmin=2))
            )

        else:
            raise RuntimeError("Inconsistency while creating seam (incorrect number of negative values), exiting.")

        # return
        return v4c, t4c, i4c, k4c, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn

    # --------------------------------------------------------------------------
    # getSeam subfunction #2

    # case 2: two nans in tetra
    def getSeamCase2(v4c, t4c, i4c, k4c, vfuncXEv1, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn, i):
        # get indices
        idxPos1 = np.where(vfuncXEv1[t4c[i, :]] > 0)[0][0]
        idxNeg1 = np.where(vfuncXEv1[t4c[i, :]] < 0)[0][0]
        idxNan1 = np.where(np.isnan(vfuncXEv1[t4c[i, :]]))[0][0]
        idxNan2 = np.where(np.isnan(vfuncXEv1[t4c[i, :]]))[0][1]

        # --------------------------------------------------------------
        # get new point
        np.where(vfuncXEv1[t4c[i, :]] < 0)[0][0]
        np.where(vfuncXEv1[t4c[i, :]] > 0)[0][0]

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
            i4c = np.concatenate((i4c, st.mode(i4c[t4c[i, :]])[0]))
            k4c = np.concatenate((k4c, st.mode(k4c[t4c[i, :]])[0]))
            # let's remember the new indices, and whether or not it's a
            # 'positive' or 'negative' tetra (in terms of the 2nd
            # eigenfunction)
            newVtcs = np.concatenate((newVtcs, [idxPt1]))
            tmpSgn = np.unique(vfuncXEv2[t4c[i, :]][np.logical_not(np.isnan(vfuncXEv2[t4c[i, :]]))])
            if len(tmpSgn) == 1:
                newVtcsSgn = np.concatenate((newVtcsSgn, tmpSgn))
            else:
                newVtcsSgn = np.concatenate((newVtcsSgn, np.array([0])))
                raise RuntimeError("Inconsistency while creating seam 5, exiting.")
        else:
            idxPt1 = newVtcsAdj[t4c[i, idxPos1], t4c[i, idxNeg1]]

        # create new tetras
        # pos1, pt1, nan1, nan2
        # neg1, pt1, nan1, nan2
        newTetra = np.concatenate(
            (newTetra, np.array((t4c[i, idxPos1], idxPt1, t4c[i, idxNan1], t4c[i, idxNan2]), ndmin=2))
        )
        newTetra = np.concatenate(
            (newTetra, np.array((t4c[i, idxNeg1], idxPt1, t4c[i, idxNan1], t4c[i, idxNan2]), ndmin=2))
        )

        # return
        return v4c, t4c, i4c, k4c, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn

    # --------------------------------------------------------------------------
    # getSeam subfunction #3

    # case 3: zero nans in the tetra (could happen for boundary tetras
    # towards two sides)
    def getSeamCase3(v4c, t4c, i4c, k4c, vfuncXEv1, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn, i):
        # determine which edges of the current tetra belong to the
        # boundary surface
        numEdgeBnd = np.where(
            (
                np.logical_or(t4c[i, 0] == np.unique(e4cBndOpen, axis=0), t4c[i, 1] == np.unique(e4cBndOpen, axis=0))
                .all(axis=1)
                .any(),
                np.logical_or(t4c[i, 0] == np.unique(e4cBndOpen, axis=0), t4c[i, 2] == np.unique(e4cBndOpen, axis=0))
                .all(axis=1)
                .any(),
                np.logical_or(t4c[i, 0] == np.unique(e4cBndOpen, axis=0), t4c[i, 3] == np.unique(e4cBndOpen, axis=0))
                .all(axis=1)
                .any(),
                np.logical_or(t4c[i, 1] == np.unique(e4cBndOpen, axis=0), t4c[i, 2] == np.unique(e4cBndOpen, axis=0))
                .all(axis=1)
                .any(),
                np.logical_or(t4c[i, 1] == np.unique(e4cBndOpen, axis=0), t4c[i, 3] == np.unique(e4cBndOpen, axis=0))
                .all(axis=1)
                .any(),
                np.logical_or(t4c[i, 2] == np.unique(e4cBndOpen, axis=0), t4c[i, 3] == np.unique(e4cBndOpen, axis=0))
                .all(axis=1)
                .any(),
            )
        )[0]

        # number of trias that belong to the surface (and which ones)
        numTriaBnd = np.where(
            (
                np.logical_or(np.logical_or(t4c[i, 0] == t4cBndOpen, t4c[i, 1] == t4cBndOpen), t4c[i, 2] == t4cBndOpen)
                .all(axis=1)
                .any(),
                np.logical_or(np.logical_or(t4c[i, 0] == t4cBndOpen, t4c[i, 1] == t4cBndOpen), t4c[i, 3] == t4cBndOpen)
                .all(axis=1)
                .any(),
                np.logical_or(np.logical_or(t4c[i, 0] == t4cBndOpen, t4c[i, 2] == t4cBndOpen), t4c[i, 3] == t4cBndOpen)
                .all(axis=1)
                .any(),
                np.logical_or(np.logical_or(t4c[i, 1] == t4cBndOpen, t4c[i, 2] == t4cBndOpen), t4c[i, 3] == t4cBndOpen)
                .all(axis=1)
                .any(),
            )
        )[0]

        # examine values: we can only cut trias or edges that have at
        # least one pos and one neg value

        # number of edges that have pos / neg values (and which ones)
        numEdgeCross = np.where(
            np.sum(np.sign(vfuncXEv1[t4c[i, ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))]]), axis=1) == 0
        )[0]

        # number of trias that have pos / neg values (and which ones)
        numTriaCross = np.where(
            np.abs(np.sum(np.sign(vfuncXEv1[t4c[i, ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))]]), axis=1)) == 1
        )[0]

        #
        # print(
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
        if (
            len(numEdgeBnd) == 2
            and len(numTriaBnd) == 0
            and len(numEdgeCross) == 3
            and len(numTriaCross) == 3
            and len(np.intersect1d(numEdgeBnd, numEdgeCross)) == 1
            and len(np.intersect1d(numTriaBnd, numTriaCross)) == 0
        ):
            # surface edges: 2, surface trias: 0, pos/neg edges: 3, pos/neg trias: 3, pos/neg surface edges: 1, pos/neg surface trias: 0
            # should be case 2: two nans in tetra, due to only one single pos/neg surface edge
            # now find which edge is on the surface, and set the other points
            # temporarily to nan
            tmpvfuncXEv1 = vfuncXEv1.copy()
            tmpvfuncXEv1[
                t4c[
                    i,
                    np.setdiff1d(
                        (0, 1, 2, 3),
                        ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))[
                            np.intersect1d(numEdgeBnd, numEdgeCross).item()
                        ],
                    ),
                ]
            ] = np.nan
            v4c, t4c, i4c, k4c, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn = getSeamCase2(
                v4c, t4c, i4c, k4c, tmpvfuncXEv1, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn, i
            )
        elif (
            len(numEdgeBnd) == 3
            and len(numTriaBnd) == 1
            and len(numEdgeCross) == 3
            and len(numTriaCross) == 3
            and len(np.intersect1d(numEdgeBnd, numEdgeCross)) == 2
            and len(np.intersect1d(numTriaBnd, numTriaCross)) == 1
        ):
            # surface edges: 3, surface trias: 1, pos/neg edges: 3, pos/neg trias: 3, pos/neg surface edges: 2, pos/neg surface trias: 1
            # should be case 1: one nan in tetra, due to one single surface tria
            # now find which tria is on the surface, and set the remaining vtx
            # temporarily to nan
            tmpvfuncXEv1 = vfuncXEv1.copy()
            tmpvfuncXEv1[
                t4c[
                    i,
                    np.setdiff1d(
                        (0, 1, 2, 3),
                        ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))[np.intersect1d(numTriaBnd, numTriaCross).item()],
                    ),
                ]
            ] = np.nan
            v4c, t4c, i4c, k4c, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn = getSeamCase1(
                v4c, t4c, i4c, k4c, tmpvfuncXEv1, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn, i
            )
        elif (
            len(numEdgeBnd) == 3
            and len(numTriaBnd) == 1
            and len(numEdgeCross) == 4
            and len(numTriaCross) == 4
            and len(np.intersect1d(numEdgeBnd, numEdgeCross)) == 2
            and len(np.intersect1d(numTriaBnd, numTriaCross)) == 1
        ):
            # surface edges: 3, surface trias: 1, pos/neg edges: 4, pos/neg trias: 4, pos/neg surface edges: 2, pos/neg surface trias: 1
            # should be case 1: one nan in tetra, due to one single surface tria
            # now find which tria is on the surface, and set the remaining vtx
            # temporarily to nan - essentially the same case as above, but kept
            # separate for clarity
            tmpvfuncXEv1 = vfuncXEv1.copy()
            tmpvfuncXEv1[
                t4c[
                    i,
                    np.setdiff1d(
                        (0, 1, 2, 3),
                        ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))[np.intersect1d(numTriaBnd, numTriaCross).item()],
                    ),
                ]
            ] = np.nan
            v4c, t4c, i4c, k4c, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn = getSeamCase1(
                v4c, t4c, i4c, k4c, tmpvfuncXEv1, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn, i
            )
        else:
            raise RuntimeError("Unknown getSeam() case, exiting.")

        # return
        return v4c, t4c, i4c, k4c, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn

    # --------------------------------------------------------------------------
    # main part of function
    # --------------------------------------------------------------------------

    # initialize
    newVtcs = np.empty(shape=(0), dtype=int)
    newVtcsSgn = np.empty(shape=(0), dtype=int)
    newVtcsAdj = sp.lil_matrix((len(v4c), len(v4c)), dtype=int)
    newTetra = np.empty(shape=(0, 4), dtype=int)

    # assign first eigenfunction
    vfuncXEv1 = np.zeros(np.shape(v4c)[0]) * np.nan
    vfuncXEv1[v4cBndOpenKeep] = anisoLaplEvec[:, 1]

    # assign second eigenfunction
    vfuncXEv2 = np.zeros(np.shape(v4c)[0]) * np.nan
    vfuncXEv2[v4cBndOpenKeep] = np.sign(anisoLaplEvec[:, 2])

    # find tetra with zeros: first, transform boundary trias to edges; here,
    # the important trick is to start working with the boundary trias from
    # t4cBndOpen, which is to avoid trias that have e.g. 2 points at the upper
    # surface and 1 point on the lower surface (which can happen). Then
    # transform to edges and identify those with pos and neg values. We use
    # edges because these need to be cut also, not only trias, and because a
    # tetra can also just share an edge with the surface, not necessarily a
    # tria. Finally select all tetras that include one such edge.
    e4cBndOpen = np.concatenate((t4cBndOpen[:, (0, 1)], t4cBndOpen[:, (0, 2)], t4cBndOpen[:, (1, 2)]), axis=0)
    e4cBndOpen = np.unique(np.sort(e4cBndOpen, axis=1), axis=1)
    zeroEdgeIdx = e4cBndOpen[np.where(np.abs(np.sum(np.sign(vfuncXEv1[e4cBndOpen]), axis=1)) == 0)[0]]
    zeroTetraIdx = list()
    for i in zeroEdgeIdx:
        zeroTetraIdx.extend(np.where(np.sum(np.isin(t4c, i), axis=1) == 2)[0].tolist())
    zeroTetraIdx = np.unique(np.array(zeroTetraIdx))

    # from lapy import plot as lpp
    #
    # tetMesh = TetMesh(v4c, t4c)
    #
    # vfunc = np.zeros(np.shape(v4c)[0])
    # vfunc[np.unique(t4c[zeroTetraIdx])] = 1
    #
    # lpp.plot_tet_mesh(tetMesh, html_output=True, vfunc=vfunc, plot_edges=True, width=1200, height=1200)

    # loop across tetras with zeros (locally disable warning that np.nan < 0 is
    # False and np.nan > 0 is also False, which is intended)
    with np.errstate(invalid="ignore"):
        for i in zeroTetraIdx:
            # we can typically have one or two nans in the tetra. we may
            # sometimes have zero nans in the tetra.
            if np.sum(np.isnan(vfuncXEv1[t4c[i]])) == 1:
                # if there is one nan, we need to make sure that the three vtcs form a surface triangle; alternatively we check if it is at least a surface edge
                if np.isin(t4cBndOpen, t4c[i, np.logical_not(np.isnan(vfuncXEv1[t4c[i]]))]).all(axis=1).any():
                    v4c, t4c, i4c, k4c, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn = getSeamCase1(
                        v4c, t4c, i4c, k4c, vfuncXEv1, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn, i
                    )
                elif np.isin(e4cBndOpen, t4c[i, np.logical_not(np.isnan(vfuncXEv1[t4c[i]]))]).all(axis=1).any():
                    # need to temporarily set one vtx to nan
                    tmpEdge = (
                        np.isin(e4cBndOpen, t4c[i, (0, 1)]).all(axis=1).any(),
                        np.isin(e4cBndOpen, t4c[i, (0, 2)]).all(axis=1).any(),
                        np.isin(e4cBndOpen, t4c[i, (0, 3)]).all(axis=1).any(),
                        np.isin(e4cBndOpen, t4c[i, (1, 2)]).all(axis=1).any(),
                        np.isin(e4cBndOpen, t4c[i, (1, 3)]).all(axis=1).any(),
                        np.isin(e4cBndOpen, t4c[i, (2, 3)]).all(axis=1).any(),
                    )
                    if len(np.where(tmpEdge)) > 1:
                        raise RuntimeError("Found " + str(len(np.where(tmpEdge))) + " edges, exiting.")
                    tmpEdgeNot = np.setdiff1d(
                        (0, 1, 2, 3), ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))[np.where(tmpEdge)[0].item()]
                    )
                    tmpvfuncXEv1 = vfuncXEv1.copy()
                    tmpvfuncXEv1[t4c[i, tmpEdgeNot]] = np.nan
                    v4c, t4c, i4c, k4c, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn = getSeamCase2(
                        v4c, t4c, i4c, k4c, tmpvfuncXEv1, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn, i
                    )
                else:
                    raise RuntimeError("Vtx " + str(i) + " not in triangle or edge, skipping")

            elif np.sum(np.isnan(vfuncXEv1[t4c[i]])) == 2:
                # if there are two nans, we need to make sure that they are linked by a surface edge; this should have happenend already.
                tmpEdge = (
                    np.isin(e4cBndOpen, t4c[i, (0, 1)]).all(axis=1).any(),
                    np.isin(e4cBndOpen, t4c[i, (0, 2)]).all(axis=1).any(),
                    np.isin(e4cBndOpen, t4c[i, (0, 3)]).all(axis=1).any(),
                    np.isin(e4cBndOpen, t4c[i, (1, 2)]).all(axis=1).any(),
                    np.isin(e4cBndOpen, t4c[i, (1, 3)]).all(axis=1).any(),
                    np.isin(e4cBndOpen, t4c[i, (2, 3)]).all(axis=1).any(),
                )
                if len(np.where(tmpEdge)) > 1:
                    raise RuntimeError("Found " + str(len(np.where(tmpEdge))) + " edges, exiting.")
                v4c, t4c, i4c, k4c, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn = getSeamCase2(
                    v4c, t4c, i4c, k4c, vfuncXEv1, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn, i
                )

            elif np.sum(np.isnan(vfuncXEv1[t4c[i]])) == 0:
                # if there are no nans, the checking is done within the getSeamCase3 function
                v4c, t4c, i4c, k4c, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn = getSeamCase3(
                    v4c, t4c, i4c, k4c, vfuncXEv1, e4cBndOpen, newTetra, newVtcs, newVtcsAdj, newVtcsSgn, i
                )

            else:
                # we should never have four nans
                raise RuntimeError("Inconsistency while creating seam (incorrect number of nans), exiting.")

    # remove tetras
    t4c = np.delete(t4c, zeroTetraIdx, axis=0)

    # add new tetras
    t4c = np.concatenate((t4c, newTetra))

    # return
    return v4c, t4c, i4c, k4c, newVtcs, newVtcsSgn, newTetra


# ------------------------------------------------------------------------------
# MAIN FUNCTION
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# computeCubeParam


def computeCubeParam(params):
    # --------------------------------------------------------------------------
    # message

    print()
    print("--------------------------------------------------------------------------------")
    print("Computing cube parametrization")
    print()

    # --------------------------------------------------------------------------
    # get data

    # evaluate params

    cutTetraMeshDir = os.path.join(params.OUTDIR, "tetra-cube")
    cutTetraMeshFile = os.path.join(params.OUTDIR, params.HEMI + ".cut.vtk")
    cutTetraIndexFile = os.path.join(params.OUTDIR, "tetra-cut", params.HEMI + ".cut.psol")
    openBndCutTetraMeshFile = os.path.join(params.OUTDIR, "tetra-cut", params.HEMI + ".open.bnd.cut.vtk")
    tetraIndexFile = os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".tetra.psol")

    if len(params.internal.aniso_alpha) == 1:
        paramAnisoAlpha = params.internal.aniso_alpha[0]
    else:
        paramAnisoAlpha = params.internal.aniso_alpha

    # load cut tetra mesh

    tetMeshCut = TetMesh.read_vtk(cutTetraMeshFile)
    v4c = tetMeshCut.v
    t4c = tetMeshCut.t

    # load indices of cutting planes

    i4c = io.read_vfunc(cutTetraIndexFile)
    i4c = np.array(i4c)

    # load subfields indices

    j4c = io.read_vfunc(tetraIndexFile)
    j4c = np.array(j4c)

    # append subfield indices

    k4c = np.concatenate((j4c, np.zeros(len(i4c) - len(j4c))))

    # --------------------------------------------------------------------------
    # get a tria mesh that is cut open at its ends, but still consistent with
    # the tetra mesh

    triaMesh4cBndOpen = TriaMesh.read_vtk(openBndCutTetraMeshFile)

    v4cBndOpen = triaMesh4cBndOpen.v
    t4cBndOpen = triaMesh4cBndOpen.t

    # remove free vertices and orient

    v4cBndOpenKeep, v4cBndOpenDel = triaMesh4cBndOpen.rm_free_vertices_()

    triaMesh4cBndOpenRm = TriaMesh(v=triaMesh4cBndOpen.v, t=triaMesh4cBndOpen.t)

    triaMesh4cBndOpenRm.orient_()

    v4cBndOpenRm = triaMesh4cBndOpenRm.v
    t4cBndOpenRm = triaMesh4cBndOpenRm.t

    # --------------------------------------------------------------------------
    # compute anisotropic laplace

    # compute first eigenfunction

    fem = Solver(triaMesh4cBndOpenRm, lump=True, aniso=paramAnisoAlpha, aniso_smooth=params.internal.aniso_smooth)

    anisoLaplEval, anisoLaplEvec = fem.eigs(k=3)

    # --------------------------------------------------------------------------
    # get mapping of open boundary vertices to subfields

    with open(os.path.join(os.path.dirname(openBndCutTetraMeshFile), params.HEMI + ".rm.open.bnd.cut.lst"), "r") as f:
        hsfList = f.read().splitlines()

    hsfList = np.array(hsfList).astype(int)

    io.write_vfunc(
        os.path.join(os.path.dirname(openBndCutTetraMeshFile), params.HEMI + ".hsf.rm.open.bnd.cut.psol"), k4c[hsfList]
    )

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
    e4cBndOpen = np.concatenate((t4cBndOpen[:, (0, 1)], t4cBndOpen[:, (0, 2)], t4cBndOpen[:, (1, 2)]), axis=0)
    e4cBndOpen = np.unique(np.sort(e4cBndOpen, axis=1), axis=1)
    zeroEdgeIdx = e4cBndOpen[np.where(np.abs(np.sum(np.sign(vfuncXEv1[e4cBndOpen]), axis=1)) == 0)[0]]
    zeroVtxIdx = np.unique(zeroEdgeIdx)

    # decide whether or not to flip anisoLaplEvec[:, 1]  (should be inf -> sup)

    # check EV1 for flips; this is done indirectly: we get the locations
    # of highest and lowest values in 234 (PrSbc); highest values should have a
    # higher value on the z-axis than lower values; note: could do the same
    #  (i.e., reverse) for 240 as well

    if np.median(
        v4cBndOpenRm[
            np.where(np.logical_and(anisoLaplEvec[:, 1] > 0, k4c[hsfList] == params.LUTDICT["presubiculum"]))[0], 2
        ]
    ) > np.median(
        v4cBndOpenRm[
            np.where(np.logical_and(anisoLaplEvec[:, 1] < 0, k4c[hsfList] == params.LUTDICT["presubiculum"]))[0], 2
        ]
    ):
        pass
    elif np.median(
        v4cBndOpenRm[
            np.where(np.logical_and(anisoLaplEvec[:, 1] > 0, k4c[hsfList] == params.LUTDICT["presubiculum"]))[0], 2
        ]
    ) < np.median(
        v4cBndOpenRm[
            np.where(np.logical_and(anisoLaplEvec[:, 1] < 0, k4c[hsfList] == params.LUTDICT["presubiculum"]))[0], 2
        ]
    ):
        LOGGER.info("Flipping EV1")
        anisoLaplEvec[:, 1] = -anisoLaplEvec[:, 1]
    else:
        io.write_vfunc(os.path.join(cutTetraMeshDir, params.HEMI + ".lapy.aLBO.EV1.psol"), anisoLaplEvec[:, 1])
        io.write_vfunc(os.path.join(cutTetraMeshDir, params.HEMI + ".lapy.aLBO.EV2.psol"), anisoLaplEvec[:, 2])
        raise RuntimeError("Inconsistency detected for EV1, exiting.")

    # decide whether or not to flip anisoLaplEvec[:, 2] (should be 234 -> 240)

    if (
        np.median(anisoLaplEvec[np.where(k4c[hsfList] == params.LUTDICT["presubiculum"])[0], 2]) < 0
        and np.median(
            anisoLaplEvec[
                np.where(
                    np.logical_or(k4c[hsfList] == params.LUTDICT["ca3"], k4c[hsfList] == params.LUTDICT["bndca4"])
                )[0],
                2,
            ]
        )
        > 0
    ):
        pass
    elif (
        np.median(anisoLaplEvec[np.where(k4c[hsfList] == params.LUTDICT["presubiculum"])[0], 2]) > 0
        and np.median(
            anisoLaplEvec[
                np.where(
                    np.logical_or(k4c[hsfList] == params.LUTDICT["ca3"], k4c[hsfList] == params.LUTDICT["bndca4"])
                )[0],
                2,
            ]
        )
        < 0
    ):
        LOGGER.info("Flipping EV2")
        anisoLaplEvec[:, 2] = -anisoLaplEvec[:, 2]
    else:
        io.write_vfunc(os.path.join(cutTetraMeshDir, params.HEMI + ".lapy.aLBO.EV1.psol"), anisoLaplEvec[:, 1])
        io.write_vfunc(os.path.join(cutTetraMeshDir, params.HEMI + ".lapy.aLBO.EV2.psol"), anisoLaplEvec[:, 2])
        raise RuntimeError("Inconsistency detected for EV2, exiting.")

    io.write_vfunc(os.path.join(cutTetraMeshDir, params.HEMI + ".lapy.aLBO.EV1.psol"), anisoLaplEvec[:, 1])
    io.write_vfunc(os.path.join(cutTetraMeshDir, params.HEMI + ".lapy.aLBO.EV2.psol"), anisoLaplEvec[:, 2])

    # --------------------------------------------------------------------------
    # define boundary conditions (note that we changed dimensions and directions!)

    # get seam

    v4c, t4c, i4c, k4c, newVtcs, newVtcsSgn, newTetra = getSeam(
        v4c, t4c, i4c, k4c, v4cBndOpenKeep, t4cBndOpen, anisoLaplEvec
    )

    # get 2nd principal axis, note that direction may be arbitrary, will
    # possibly be switched below; note that this replaces EV2-based newVtcsSgn,
    # which we may remove in the future

    pca = PCA(n_components=3)
    pca_scores = pca.fit_transform(v4c)

    # x: medial -> lateral

    vfuncX = np.zeros(np.shape(v4c)[0])
    srcX = np.zeros(np.shape(v4c)[0])
    srcX[newVtcs[v4c[newVtcs, 0] < np.mean(v4c[newVtcs, 0])]] = -1
    ###srcX[newVtcs[pca_scores[newVtcs,1]<0]] = -1
    ###srcX[newVtcs[v4c[newVtcs,0]>0]] = 1
    ###srcX[newVtcs[newVtcsSgn>0]] = 1
    snkX = np.zeros(np.shape(v4c)[0])
    snkX[newVtcs[v4c[newVtcs, 0] >= np.mean(v4c[newVtcs, 0])]] = 1
    ###snkX[newVtcs[pca_scores[newVtcs,1]>0]] = 1
    ###snkX[newVtcs[v4c[newVtcs,0]<0]] = -1
    ###snkX[newVtcs[newVtcsSgn<0]] = -1
    vfuncX = srcX + snkX

    # y: Tail (226) -> Head (232) (post --> ant)

    vfuncY = np.zeros(np.shape(v4c)[0])
    vfuncY[np.where(i4c == params.LUTDICT["jointtail"])[0]] = -1
    vfuncY[np.where(i4c == params.LUTDICT["jointhead"])[0]] = 1

    # z: inf->sup

    vfuncZ = np.zeros(np.shape(v4c)[0])
    srcZ = np.zeros(np.shape(anisoLaplEvec[:, 1])[0])
    srcZ[anisoLaplEvec[:, 1] < 0] = -1
    snkZ = np.zeros(np.shape(anisoLaplEvec[:, 1])[0])
    snkZ[anisoLaplEvec[:, 1] > 0] = 1
    vfuncZ[v4cBndOpenKeep] = srcZ + snkZ

    # determine left / right

    img = nb.load(os.path.join(params.OUTDIR, params.HEMI + ".image.mgz"))
    mat = np.linalg.inv(img.header.get_vox2ras_tkr())
    pts = np.concatenate((v4c, np.ones((v4c.shape[0], 1))), axis=1)
    ras = np.matmul(np.matmul(pts, mat.transpose()), img.header.get_vox2ras().transpose())

    if params.HEMI == "lh":
        if (ras[:, 0] < 0.0).all():
            # lateral (snk, +1) should be more negative than medial (src, -1), otherwise flip
            if np.mean(ras[np.where(vfuncX > 0)[0], 0]) > np.mean(ras[np.where(vfuncX < 0)[0], 0]):
                vfuncX = -vfuncX
        else:
            raise RuntimeError("Ambiguous hemisphere info, exiting.")
    elif params.HEMI == "rh":
        if (ras[:, 0] > 0.0).all():
            # lateral (snk, +1) should be more positive than medial (src, -1), otherwise flip
            if np.mean(ras[np.where(vfuncX > 0)[0], 0]) < np.mean(ras[np.where(vfuncX < 0)[0], 0]):
                vfuncX = -vfuncX
        else:
            raise RuntimeError("Ambiguous hemisphere info, exiting.")

    # --------------------------------------------------------------------------
    # compute ABtetra and Poisson functions

    # first, need to remove free vtcs

    tetMesh4c = TetMesh(v4c, t4c)

    v4cRmKeep, v4cRmDel = tetMesh4c.rm_free_vertices_()

    v4cRm = tetMesh4c.v
    t4cRm = tetMesh4c.t

    TetMesh.write_vtk(tetMesh4c, os.path.join(cutTetraMeshDir, params.HEMI + ".seam.rm.cut.vtk"))

    #

    vfuncXRm = vfuncX[v4cRmKeep]
    vfuncYRm = vfuncY[v4cRmKeep]
    vfuncZRm = vfuncZ[v4cRmKeep]

    #

    bndXRm = (np.where(vfuncXRm)[0], vfuncXRm[np.where(vfuncXRm)])
    bndYRm = (np.where(vfuncYRm)[0], vfuncYRm[np.where(vfuncYRm)])
    bndZRm = (np.where(vfuncZRm)[0], vfuncZRm[np.where(vfuncZRm)])

    #

    tetMesh4cRm = TetMesh(v4cRm, t4cRm)

    fem = Solver(tetMesh4cRm)

    P0 = fem.poisson(dtup=bndXRm)
    P1 = fem.poisson(dtup=bndYRm)
    P2 = fem.poisson(dtup=bndZRm)

    # --------------------------------------------------------------------------
    # output

    # write out functions

    io.write_vfunc(os.path.join(cutTetraMeshDir, params.HEMI + ".vfuncX.seam.rm.cut.psol"), vfuncXRm)
    io.write_vfunc(os.path.join(cutTetraMeshDir, params.HEMI + ".vfuncY.seam.rm.cut.psol"), vfuncYRm)
    io.write_vfunc(os.path.join(cutTetraMeshDir, params.HEMI + ".vfuncZ.seam.rm.cut.psol"), vfuncZRm)

    io.write_vfunc(os.path.join(cutTetraMeshDir, params.HEMI + ".poisson0.seam.rm.cut.psol"), P0)
    io.write_vfunc(os.path.join(cutTetraMeshDir, params.HEMI + ".poisson1.seam.rm.cut.psol"), P1)
    io.write_vfunc(os.path.join(cutTetraMeshDir, params.HEMI + ".poisson2.seam.rm.cut.psol"), P2)

    # also write out boundary files for visualization

    tetMesh4cRmBnd = tetMesh4cRm.boundary_tria()

    v4cRmBndRmKeep, v4cRmBndRmDel = tetMesh4cRmBnd.rm_free_vertices_()

    v4cRmBndRm = tetMesh4cRmBnd.v
    t4cRmBndRm = tetMesh4cRmBnd.t

    TriaMesh.write_vtk(tetMesh4cRmBnd, filename=os.path.join(cutTetraMeshDir, params.HEMI + ".rm.bnd.seam.rm.cut.vtk"))

    io.write_vfunc(
        os.path.join(cutTetraMeshDir, params.HEMI + ".vfuncX.rm.bnd.seam.rm.cut.psol"), vfuncXRm[v4cRmBndRmKeep]
    )
    io.write_vfunc(
        os.path.join(cutTetraMeshDir, params.HEMI + ".vfuncY.rm.bnd.seam.rm.cut.psol"), vfuncYRm[v4cRmBndRmKeep]
    )
    io.write_vfunc(
        os.path.join(cutTetraMeshDir, params.HEMI + ".vfuncZ.rm.bnd.seam.rm.cut.psol"), vfuncZRm[v4cRmBndRmKeep]
    )

    io.write_vfunc(os.path.join(cutTetraMeshDir, params.HEMI + ".poisson0.rm.bnd.seam.rm.cut.psol"), P0[v4cRmBndRmKeep])
    io.write_vfunc(os.path.join(cutTetraMeshDir, params.HEMI + ".poisson1.rm.bnd.seam.rm.cut.psol"), P1[v4cRmBndRmKeep])
    io.write_vfunc(os.path.join(cutTetraMeshDir, params.HEMI + ".poisson2.rm.bnd.seam.rm.cut.psol"), P2[v4cRmBndRmKeep])

    # write out cube

    w4cRm = np.vstack((P0, P1, P2)).transpose()

    TetMesh.write_vtk(TetMesh(v=w4cRm, t=t4cRm), os.path.join(cutTetraMeshDir, params.HEMI + ".uvw.seam.rm.cut.vtk"))

    #

    return params
