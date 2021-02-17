"""
This module contains a set of auxiliary functions for mesh processing,
specifically:

- createSTL(filename, v, t)
- levelsetsTria(v, t, p, levelsets)
- levelsetsTetra(v, t, p, levelsets)

"""


# ------------------------------------------------------------------------------
# createSTL

def createSTL(filename, v, t):
    """
    createSTL(filename, v, t)

    A function to write STL files

    """

    import numpy as np

    with open(filename, 'w') as f:

        print("solid lap-data", file=f)

        for i in range(0, len(t)):

            v1mv0 = v[t[i, 1],:] - v[t[i, 0],:]
            v2mv0 = v[t[i, 2],:] - v[t[i, 0],:]

            normal = np.cross(v1mv0, v2mv0) / np.linalg.norm(np.cross(v1mv0, v2mv0))

            print("  facet normal", file=f)

            print("    outer loop", file=f)

            print("      vertex %f %f %f " % (v[t[i, 0], 0], v[t[i, 0], 1],v[t[i, 0], 2]), file=f)
            print("      vertex %f %f %f " % (v[t[i, 1], 0], v[t[i, 1], 1],v[t[i, 1], 2]), file=f)
            print("      vertex %f %f %f " % (v[t[i, 2], 0], v[t[i, 2], 1],v[t[i, 2], 2]), file=f)

            print("    endloop", file=f)

            print("  endfacet", file=f)

        print("end solid lap-data", file=f)


# ------------------------------------------------------------------------------
#

def levelsetsTria(v, t, p, levelsets):

    import numpy as np
    from scipy.sparse import csr_matrix, lil_matrix

    vLVL = list()
    lLVL = list()
    iLVL = list()

    levelsets = (np.array(levelsets, ndmin=2))

    for l in range(len(levelsets)):

        A = lil_matrix((np.shape(v)[0],np.shape(v)[0]))

        lvl = levelsets[l]

        nlvl = p[t] > lvl

        n = np.where(np.logical_or(np.sum(nlvl, axis=1) == 1 , np.sum(nlvl, axis=1) == 2))[0]

        # interpolate points

        ti = list()
        vi = list()

        for i in range(len(n)):

            # which are the outlying points in the current tria?
            oi = np.where(nlvl[n[i],:])[0]

            #  convert 2 --> 1
            if len(oi) == 2:
                oi = np.setdiff1d((0,1,2), oi)

            # find the two non - outyling points
            oix = np.setdiff1d((0, 1, 2), oi)

            # check if we have interpolated for one or both of these points before

            if np.count_nonzero(A[t[n[i], oi], t[n[i], oix[0]]]) ==0 :

                # compute difference vectors between outlying point and other points

                d10 = v[ t[n[i], oix[0]], :] - v[ t[n[i], oi], :]

                # compute differences of all points to lvl to get interpolation factors

                s10 = (lvl - p[t[n[i], oi]]) / (p[t[n[i], oix[0]]] - p[t[n[i], oi]])

                # compute new points

                v10 = s10 * d10 + v[ t[n[i], oi], :]

                # update vi and index(order matters)

                vi.append(v10.tolist()[0])

                ti10 = len(vi)

                # store between which two points we are interpolating (to avoid having duplicate points)

                A[ t[n[i], oi], t[n[i], oix[0]] ] = ti10
                A[ t[n[i], oix[0]], t[n[i], oi] ] = ti10

            else:

                ti10 = int(A[ t[n[i], oi], t[n[i], oix[0]] ].toarray().item())

            # essentially the same as above, just for oix[1]

            if np.count_nonzero(A[t[n[i], oi], t[n[i], oix[1]]]) == 0:

                d20 = v[ t[n[i], oix[1]], :] - v[ t[n[i], oi], :]

                s20 = (lvl - p[t[n[i], oi]]) / (p[t[n[i], oix[1]]] - p[t[n[i], oi]])

                v20 = s20 * d20 + v[ t[n[i], oi], :]

                # update vi and index(order matters)

                vi.append(v20.tolist()[0])

                ti20 = len(vi)

                A[ t[n[i], oi], t[n[i], oix[1]] ] = ti20
                A[ t[n[i], oix[1]], t[n[i], oi] ] = ti20

            else:

                ti20 = int(A[ t[n[i], oi], t[n[i], oix[1]] ].toarray().item())

            # store new indices

            ti.append((ti10,ti20))

            # clean up

            # clear oi oix d10 d20 s10 s20 v10 v20 t10 t20

        # store

        vLVL.append(vi)
        lLVL.append(ti)
        iLVL.append(n)

    return vLVL, lLVL, iLVL


# ------------------------------------------------------------------------------
# levelsetsTetra

def levelsetsTetra(v, t, p, levelsets):

    import numpy as np
    from scipy.sparse import csr_matrix, lil_matrix

    vLVL = list()
    tLVL = list()
    iLVL = list()
    jLVL = list()
    oLVL = list()

    levelsets = (np.array(levelsets, ndmin=2))

    for l in range(len(levelsets[0])):

        # matrix to store previous interpolation

        A = lil_matrix((np.shape(v)[0], np.shape(v)[0]))

        # current levelset

        lvl = levelsets[0][l]

        # for each tetra, which of its four points are outlying? nlvl is yes/no, n contains indices

        nlvl = p[t] > lvl

        n = np.where(np.logical_or(np.sum(nlvl, axis=1) == 1 , np.logical_or(np.sum(nlvl, axis=1) == 2, np.sum(nlvl, axis=1) == 3)))[0]

        # interpolate points

        ti = list()
        vi = list()
        m = list()
        o = list()

        for i in range(len(n)):

            # which are the outlying points in the current tetra?

            oi = np.where(nlvl[n[i],:])[0]

            # store oi prior to converting 3 --> 1

            o.append(oi)

            #  convert 3 --> 1

            if len(oi) == 3:
                oi = np.setdiff1d((0,1,2,3), oi)

            # interpolate

            if len(oi) == 1:

                # find the three non - outyling points

                oix = np.setdiff1d((0, 1, 2, 3), oi)

                # check if we have interpolated for one or more of these points before

                if np.sum(np.count_nonzero(A[t[n[i], oi], t[n[i], oix[0]]])) == 0:

                    # compute difference vectors between outlying point and other points

                    d10 = v[ t[n[i], oix[0]], :] - v[ t[n[i], oi], :]

                    # compute differences of all points to lvl to get interpolation factors

                    s10 = (lvl - p[t[n[i], oi]]) / (p[t[n[i], oix[0]]] - p[t[n[i], oi]])

                    # compute new points

                    v10 = s10 * d10 + v[ t[n[i], oi], :]

                    # update vi and index(order matters)

                    vi.append(v10.tolist()[0])

                    ti10 = len(vi)

                    # store between which two points we are interpolating (to avoid having duplicate points)

                    A[ t[n[i], oi], t[n[i], oix[0]] ] = ti10
                    A[ t[n[i], oix[0]], t[n[i], oi] ] = ti10

                else:

                    ti10 = int(A[ t[n[i], oi], t[n[i], oix[0]] ].toarray().item())

                # essentially the same as above, just for oix[1]

                if np.sum(np.count_nonzero(A[t[n[i], oi], t[n[i], oix[1]]])) == 0:

                    d20 = v[ t[n[i], oix[1]], :] - v[ t[n[i], oi], :]

                    s20 = (lvl - p[t[n[i], oi]]) / (p[t[n[i], oix[1]]] - p[t[n[i], oi]])

                    v20 = s20 * d20 + v[ t[n[i], oi], :]

                    # update vi and index(order matters)

                    vi.append(v20.tolist()[0])

                    ti20 = len(vi)

                    A[ t[n[i], oi], t[n[i], oix[1]] ] = ti20
                    A[ t[n[i], oix[1]], t[n[i], oi] ] = ti20

                else:

                    ti20 = int(A[ t[n[i], oi], t[n[i], oix[1]] ].toarray().item())

                # essentially the same as above, just for oix[2]

                if np.sum(np.count_nonzero(A[t[n[i], oi], t[n[i], oix[2]]])) == 0:

                    d30 = v[ t[n[i], oix[2]], :] - v[ t[n[i], oi], :]

                    s30 = (lvl - p[t[n[i], oi]]) / (p[t[n[i], oix[2]]] - p[t[n[i], oi]])

                    v30 = s30 * d30 + v[ t[n[i], oi], :]

                    # update vi and index(order matters)

                    vi.append(v30.tolist()[0])

                    ti30 = len(vi)

                    A[ t[n[i], oi], t[n[i], oix[2]] ] = ti30
                    A[ t[n[i], oix[2]], t[n[i], oi] ] = ti30

                else:

                    ti30 = int(A[ t[n[i], oi], t[n[i], oix[2]] ].toarray().item())

                # store new indices

                ti.append((ti10,ti20,ti30))

                m.append(n[i])

            elif len(oi) == 2:

                # find the two non-outyling points

                oix = np.setdiff1d((0, 1, 2, 3), oi)

                # check if we have interpolated for one or more of these points
                # before

                if np.sum(np.count_nonzero(A[t[n[i], oi[0]], t[n[i], oix[0]]])) == 0:

                    # compute difference vectors between outlying point and other points

                    d20 = v[ t[n[i], oix[0]], :] - v[ t[n[i], oi[0]], :]

                    # compute differences of all points to lvl to get interpolation
                    # factors

                    s20 = (lvl - p[t[n[i], oi[0]]]) / (p[t[n[i], oix[0]]] - p[t[n[i], oi[0]]])

                    # compute new points

                    v20 = s20 * d20 + v[ t[n[i], oi[0]], :]

                    # update vi and index (order matters)

                    vi.append(v20.tolist())

                    ti20 = len(vi)

                    # store between which two points we are interpolating (to
                    # avoid having duplicate points)

                    A[ t[n[i], oi[0]], t[n[i], oix[0]] ] = ti20
                    A[ t[n[i], oix[0]], t[n[i], oi[0]] ] = ti20

                else:

                    ti20 = int(A[ t[n[i], oi[0]], t[n[i], oix[0]] ])

                # same as above for oi[0] and oix[1]

                if np.sum(np.count_nonzero(A[t[n[i], oi[0]], t[n[i], oix[1]]])) == 0:

                    # compute difference vectors between outlying point and other points

                    d30 = v[ t[n[i], oix[1]], :] - v[ t[n[i], oi[0]], :]

                    # compute differences of all points to lvl to get interpolation
                    # factors

                    s30 = (lvl - p[t[n[i], oi[0]]]) / (p[t[n[i], oix[1]]] - p[t[n[i], oi[0]]])

                    # compute new points

                    v30 = s30 * d30 + v[ t[n[i], oi[0]], :]

                    # update vi and index (order matters)

                    vi.append(v30.tolist())

                    ti30 = len(vi)

                    # store between which two points we are interpolating (to
                    # avoid having duplicate points)

                    A[ t[n[i], oi[0]], t[n[i], oix[1]] ] = ti30
                    A[ t[n[i], oix[1]], t[n[i], oi[0]] ] = ti30

                else:

                    ti30 = int(A[ t[n[i], oi[0]], t[n[i], oix[1]] ])

                # same as above for oi[1] and oix[0]

                if np.sum(np.count_nonzero(A[t[n[i], oi[1]], t[n[i], oix[0]]])) == 0:

                    # compute difference vectors between outlying point and other points

                    d21 = v[ t[n[i], oix[0]], :] - v[ t[n[i], oi[1]], :]

                    # compute differences of all points to lvl to get interpolation
                    # factors

                    s21 = (lvl - p[t[n[i], oi[1]]]) / (p[t[n[i], oix[0]]] - p[t[n[i], oi[1]]])

                    # compute new points

                    v21 = s21 * d21 + v[ t[n[i], oi[1]], :]

                    # update vi and index (order matters)

                    vi.append(v21.tolist())

                    ti21 = len(vi)

                    # store between which two points we are interpolating (to
                    # avoid having duplicate points)

                    A[ t[n[i], oi[1]], t[n[i], oix[0]] ] = ti21
                    A[ t[n[i], oix[0]], t[n[i], oi[1]] ] = ti21

                else:

                    ti21 = int(A[ t[n[i], oi[1]], t[n[i], oix[0]] ])

                # same as above for oi[1] and oix[1]

                if np.sum(np.count_nonzero(A[t[n[i], oi[1]], t[n[i], oix[1]]])) == 0:

                    # compute difference vectors between outlying point and other points

                    d31 = v[ t[n[i], oix[1]], :] - v[ t[n[i], oi[1]], :]

                    # compute differences of all points to lvl to get interpolation
                    # factors

                    s31 = (lvl - p[t[n[i], oi[1]]]) / (p[t[n[i], oix[1]]] - p[t[n[i], oi[1]]])

                    # compute new points

                    v31 = s31 * d31 + v[ t[n[i], oi[1]], :]

                    # update vi and index (order matters)

                    vi.append(v31.tolist())

                    ti31 = len(vi)

                    # store between which two points we are interpolating (to
                    # avoid having duplicate points)

                    A[ t[n[i], oi[1]], t[n[i], oix[1]] ] = ti31
                    A[ t[n[i], oix[1]], t[n[i], oi[1]] ] = ti31

                else:

                    ti31 = int(A[ t[n[i], oi[1]], t[n[i], oix[1]] ])

                # store new indices

                ti.append((ti20,ti21,ti30))
                ti.append((ti21,ti30,ti31))

                m.append(n[i])
                m.append(n[i])

        # store

        vLVL.append(np.array(vi)) # interpolated points
        tLVL.append(np.array(ti)) # triangles of interpolated points (1-based indices)
        iLVL.append(n) # indices of tetras with 1-3 outlying points
        jLVL.append(m) # for each interpolated point, store its tetra
        oLVL.append(o) # outlying points in the current tetra

    return vLVL, tLVL, iLVL, jLVL, oLVL
