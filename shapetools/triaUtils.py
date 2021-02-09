"""
This module contains a set of auxiliary functions for mesh processing,
specifically:

- smooth(v,t,f=None,n=1)
- get_vertex_areas(v,t)
- get_adj(t,directed=False,store_tidx=False)

- readVTK(infile)
- writeVTK(outfile, v, t)
- writePSOL(filename,p)
- readPSOL(filename)

- createSTL(filename, v, t)

- levelsetsTria(v, t, p, levelsets)
- levelsetsTetra(v, t, p, levelsets)

- poissonSolver
- computeABtetra
- computeABtria
- computeABtriaAniso

- tria_rm_free_vertices
- tria_is_oriented
- tria_get_adj
- tetra_get_boundary_tria
- orientTria

The following deprecated functions will be removed in the future:

- importVTKDeprecated(infile)
- exportVTKDeprecated(v,t,outfile)
- readPSOLDeprecated(filename)

"""

# ------------------------------------------------------------------------------
# smooth

def smooth(v,t,f=None,n=1):
    """
    Smoothes vector float function on the mesh iteratively

    Inputs:   v - vertices   List of lists of 3 float coordinates
              t - trias      List of lists of 3 int of indices (>=0) into v array
              f - function   Float vector of values at vertices,
                             if empty, use vertex coordiantes
              n              Number of iterations for smoothing

    Output:   f              Smoothed surface vertex function
    """

    import numpy as np
    from scipy import sparse

    v = np.array(v)
    t = np.array(t)
    if f is None:
        f = v
    f = np.array(f)
    if (f.shape[0] != v.shape[0]):
        raise ValueError('Error: length of vfunc needs to match number of vertices')
    areas = get_vertex_areas(v,t)[:,np.newaxis]
    adj = get_adj(t)
    # binarize:
    adj.data = np.ones(adj.data.shape)
    # adjust columns to contain areas of vertex i
    adj2 = adj.multiply(areas)
    # rowsum is the area sum of 1-ring neighbors
    rowsum = np.sum(adj2,axis=1)
    # normalize rows to sum = 1
    adj2 = adj2.multiply(1.0/rowsum)
    # apply multiple times
    adj2 = adj2.__pow__(n)
    # and avearge values
    return adj2 * f


# ------------------------------------------------------------------------------
# get_vertex_areas

def get_vertex_areas(v,t):
    """
    Computes the area associated to each vertex (1/3 of one-ring trias)

    Inputs:   v - vertices   List of lists of 3 float coordinates
              t - trias      List of lists of 3 int of indices (>=0) into v array

    Output:   varea          Array of vertex areas
    """

    import numpy as np

    v1 = v[t[:,0],:]
    v2 = v[t[:,1],:]
    v3 = v[t[:,2],:]
    v2mv1 = v2-v1
    v3mv1 = v3-v1
    cr = np.cross(v2mv1,v3mv1)
    area = 0.5* np.sqrt(np.sum(cr*cr,axis=1))
    area3 = np.repeat(area[:,np.newaxis], 3, 1)
    #varea = accumarray(t(:),area3(:))./3;
    varea = np.bincount(t.flatten(),area3.flatten())/3.0
    return varea


# ------------------------------------------------------------------------------
# get_adj

def get_adj(t,directed=False,store_tidx=False):
    """
    Compute the adjacency matrix (edge graph) of the triangle mesh t

    Inputs:   t - trias      List of lists of 3 int of indices (>=0) into v array
                             Ordering is important: All triangles should be
                             oriented the same way (counter-clockwise, when
                             looking from above)

              directed       False (default): The non-directed adjacency matrix
                             will be symmetric. Each inner edge (i,j) will have
                             the number of triangles that contain this edge.
                             Inner edges usually 2, boundary edges 1. Higher
                             numbers can occur when there are non-manifold triangles.

                             True: The directed adjacency matrix is not symmetric.
                             For manifold meshes, there are only entries with
                             value 1. Symmetric entries are inner edges. Non-symmetric
                             are boundary edges. The direction prescribes a direction
                             on the boundary loops.

              store_tidx     If true, stores the tria idx +1 instead of one
                             in the matrix (allows lookup of vertex to tria).
                             Works only on directed adj matrix, and manifold.

    Output:   adj            The sparse adjacency matrix in CSC format
                             Note: for non-directed, adj is symmetric and
                             contains a 2 at inner edges and 1 at boundaries.
                             It can be binarized via:
                                adj.data = np.ones(adj.data.shape)
                             For directed adj it contains 1 for each directed
                             edge and is nonsymmetric if boundaries exist.
                             Adding this to its transpose creates the non-
                             directed version.

    """

    import numpy as np
    from scipy import sparse

    t1 = t[:,0]
    t2 = t[:,1]
    t3 = t[:,2]
    if directed:
        I = np.column_stack((t1, t2, t3))
        J = np.column_stack((t2, t3, t1))
    else:
        if store_tidx:
            raise ValueError('store_tidx requires directed adjacency matrix')
        I = np.column_stack((t1, t2, t2, t3, t3, t1))
        J = np.column_stack((t2, t1, t3, t2, t1, t3))
    I = I.flatten()
    J = J.flatten()
    dat = np.ones(I.shape)
    if store_tidx:
        # store tria idx +1  (zero means no edge here)
        dat = np.repeat(np.arange(1,t.shape[0]+1),3)
    adj = sparse.csc_matrix((dat, (I, J)))
    return adj


# ------------------------------------------------------------------------------
# importVTKDeprecated

def importVTKDeprecated(infile):
    """
    Load VTK triangle mesh
    """
    import numpy as np

    verbose = 1
    if (verbose > 0):
        print("--> VTK format         ... ")

    try:
        f = open(infile,'r')
    except IOError:
        print("[file not found or not readable]\n")
        return

    # skip comments
    line = f.readline()
    while line[0] == '#':
        line = f.readline()

    # search for ASCII keyword in first 5 lines:
    count = 0
    while count < 5 and not line.startswith("ASCII"):
        line = f.readline()
        #print line
        count = count+1
    if not line.startswith("ASCII"):
        print("[ASCII keyword not found] --> FAILED\n")
        return

    # expect Dataset Polydata line after ASCII:
    line = f.readline()
    if not line.startswith("DATASET POLYDATA"):
        print("[read: "+ line+" expected DATASET POLYDATA] --> FAILED\n")
        return

    # read number of points
    line = f.readline()
    larr = line.split()
    if larr[0]!="POINTS" or larr[2] != "float":
        print("[read: " + line + " expected POINTS # float] --> FAILED\n")
        return
    pnum = int(larr[1])
    # read points as chunch
    v=np.fromfile(f,'float32',3*pnum,' ')
    v.shape = (pnum, 3)

    # expect polygon or tria_strip line
    line = f.readline()
    larr = line.split()

    if larr[0]=="POLYGONS":
        tnum = int(larr[1])
        ttnum = int(larr[2])
        npt = float(ttnum) / tnum;
        if (npt != 4.0) :
            print("[having: " + str(npt)+ " data per tria, expected trias 3+1] --> FAILED\n")
            return
        t = np.fromfile(f,'int',ttnum,' ')
        t.shape = (tnum, 4)
        if t[tnum-1][0] != 3:
            print("[can only read triangles] --> FAILED\n")
            return
        t = np.delete(t,0,1)

    elif larr[0]=="TRIANGLE_STRIPS":
        tnum = int(larr[1])
        ttnum = int(larr[2])
        tt = []
        for i in xrange(tnum):
            larr = f.readline().split()
            if len(larr)==0:
                print("[error reading triangle strip (i)] --> FAILED\n")
                return
            n = larr[0]
            if len(larr)!=n+1:
                print("[error reading triangle strip (ii)] --> FAILED\n")
                return
            # create triangles from strip
            # note that larr tria info starts at index 1
            for ii in range(2,n):
                if (ii%2 == 0):
                    tria = [larr[ii-1], larr[ii], larr[ii+1]]
                else:
                    tria = [larr[ii], larr[ii-1], larr[ii+1]]
                tt.append(tria)
        t = np.array(tt)

    else:
        print("[read: "+line+ " expected POLYGONS or TRIANGLE_STRIPS] --> FAILED\n")
        return

    f.close()

    print(" --> DONE ( V: " + str(v.shape[0]) + " , T: " + str(t.shape[0]) + " )\n")

    return v, t


# ------------------------------------------------------------------------------
# exportVTKDeprecated

def exportVTKDeprecated(v,t,outfile):
    """
    Save VTK file

    usage: exportVTK(vertices,triangles,outfile)

    """

    # imports
    import numpy as np

    # open file
    try:
        f = open(outfile,'w')
    except IOError:
        print("[File "+outfile+" not writable]")
        return

    # check data structure

    # ...

    #
    f.write('# vtk DataFile Version 1.0\n')
    f.write('vtk output\n')
    f.write('ASCII\n')
    f.write('DATASET POLYDATA\n')
    f.write('POINTS '+str(np.shape(v)[0])+' float\n')

    for i in range(np.shape(v)[0]):
        f.write(' '.join(map(str,v[i,:])))
        f.write('\n')

    f.write('POLYGONS '+str(np.shape(t)[0])+' '+str(4*np.shape(t)[0])+'\n')

    for i in range(np.shape(t)[0]):
        f.write(' '.join(map(str,np.append(3,t[i,:]))))
        f.write('\n')

    # close file
    f.close()


# ------------------------------------------------------------------------------
# readVTK

def readVTK(infile):
    """
    v, t = readVTK(filename)

    Read VTK mesh
    : can be tria or tetra
    : can be vtk1 or vtk2

    """

    # imports
    import numpy as np

    # open file
    try:
        f = open(infile, 'r')
    except IOError:
        print("[file not found or not readable]\n")
        return

    # skip comments
    line = f.readline()
    while line[0] == '#':
        line = f.readline()

    # search for ASCII keyword in first 5 lines:
    count = 0
    while count < 5 and not line.startswith("ASCII"):
        line = f.readline()
        count = count + 1
    if not line.startswith("ASCII"):
        print("[ASCII keyword not found] --> FAILED\n")
        return

    # expect Dataset Polydata or Unstructured Grid line after ASCII:
    line = f.readline()
    if (not line.startswith("DATASET POLYDATA")) and (not line.startswith("DATASET UNSTRUCTURED_GRID")):
        print("[read: " + line + " expected DATASET POLYDATA OR UNSTRUCTURED GRID] --> FAILED\n")
        return

    # read number of points
    line = f.readline()
    larr = line.split()
    if larr[0] != "POINTS" or (larr[2] != "float" and larr[2] != "double"):
        print("[read: " + line + " expected POINTS # float|double] --> FAILED\n")
        return
    pnum = int(larr[1])

    # read points as chunch
    v = np.fromfile(f, 'float32', 3 * pnum, ' ')
    v.shape = (pnum, 3)

    # expect polygon or cells line
    line = f.readline()
    larr = line.split()

    # read trias / tetras
    if larr[0] == "POLYGONS" or larr[0] == "CELLS":
        tnum = int(larr[1])
        ttnum = int(larr[2])
        npt = float(ttnum) / tnum
        if (npt != 4.0) and (npt != 5.0):
            print("[having: " + str(npt) + " data per tria/tetra, expected 3+1 or 4+1] --> FAILED\n")
            return
        t = np.fromfile(f, 'int', ttnum, ' ')
        if (npt == 4.0):
            t.shape = (tnum, 4)
        elif (npt == 5.0):
            t.shape = (tnum, 5)
        t = np.delete(t, 0, 1)
    else:
        print("[read: " + line + " expected POLYGONS or CELLS] --> FAILED\n")
        return

    # close file
    f.close()

    # return
    return v, t


# ------------------------------------------------------------------------------
# writeVTK

def writeVTK(outfile, v, t, legacy=False):
    """
    writeVTK(filename, v, t)

    Write VTK mesh
    : can be tria or tetra
    : can only be vtk1

    """

    # imports
    import numpy as np

    # open file
    try:
        f = open(outfile, 'w')
    except IOError:
        print("[File " + outfile + " not writable]")
        return

    #
    f.write('# vtk DataFile Version 1.0\n')
    f.write('vtk output\n')
    f.write('ASCII\n')
    f.write('DATASET POLYDATA\n')
    f.write('POINTS ' + str(np.shape(v)[0]) + ' float\n')

    for i in range(np.shape(v)[0]):
        f.write(' '.join(map(str, v[i, :])))
        f.write('\n')

    if legacy is False:
        f.write('POLYGONS ' + str(np.shape(t)[0]) + ' ' + str((np.shape(t)[1] + 1) * np.shape(t)[0]) + '\n')
    else:
        # this is to allow coping with a bug in old LapGUI versions
        f.write('POLYGONS ' + str(np.shape(t)[0]) + ' ' + str(4 * np.shape(t)[0]) + '\n')

    for i in range(np.shape(t)[0]):
        f.write(' '.join(map(str, np.append(np.shape(t)[1], t[i, :]))) + '\n')

    f.write('\n')

    # close file
    f.close()


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
# writePSOL

def writePSOL(filename,p):
    """
    writePSOL(filename, p)

    A function to export PSOL files

    """

    with open(filename, 'w') as f:
        f.write('Solution:\n')
        f.write("("+",".join(p.astype(str))+")")


# ------------------------------------------------------------------------------
# readPSOLDeprecated

def readPSOLDeprecated(filename):
    """
    p = readPSOL(filename)

    A function to import PSOL files

    """

    import re
    import numpy as np

    with open(filename) as f:
        txt = f.readlines()

    i=0
    while i<len(txt):
        if "Solution:" in txt[i]:
            i=i+1
            tmp1=list()
            while i<len(txt):
                if txt[i].isspace():
                    break
                tmp1.append(txt[i].strip())
                i=i+1
            vals=list()
            for tmp2 in re.split('[;,]',re.sub('[{()}]','',''.join(tmp1))):
                vals.append(float(tmp2))
            del(tmp1,tmp2)
        i=i+1

    vals = np.ndarray(vals)

    return vals


# ------------------------------------------------------------------------------
# readPSOL

def readPSOL(filename):
    """
    p = readPSOL(filename)

    A function to import PSOL files

    """

    import re
    import numpy as np

    with open(filename) as f:
        txt = f.readlines()

    txt = [ x.strip() for x in txt ]

    txt.remove("Solution:")

    txt = [ re.sub("[{()}]",'', x) for x in txt ]

    if len(txt) == 1:

        txt = [ re.split("[,;]", x) for x in txt ][0]

    txt = [ np.float(x) for x in txt ]

    txt = np.array(txt)

    return txt


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


# ------------------------------------------------------------------------------
# computeABtetra

def computeABtetra(v, t, lump = False):
    """
    computeABtetra(v,t) computes the two sparse symmetric matrices representing
           the Laplace Beltrami Operator for a given tetrahedral mesh using
           the linear finite element method (Neumann boundary condition).

    Inputs:   v - vertices : list of lists of 3 floats
              t - tetras   : list of lists of 4 int of indices (>=0) into v array
                             Ordering is important: first three vertices for
                             triangle counter-clockwise when looking at it from
                             the inside of the tetrahedron

    Outputs:  A - sparse sym. (n x n) positive semi definite numpy matrix
              B - sparse sym. (n x n) positive definite numpy matrix (inner product)

    Can be used to solve sparse generalized Eigenvalue problem: A x = lambda B x
    or to solve Poisson equation: A x = B f (where f is function on mesh vertices)
    or to solve Laplace equation: A x = 0
    or to model the operator's action on a vector x:   y = B\(Ax)
    """

    # imports
    import numpy as np
    from scipy import sparse

    # initialize
    v = np.array(v)
    t = np.array(t)

    # Compute vertex coordinates and a difference vector for each triangle:
    t1 = t[:, 0]
    t2 = t[:, 1]
    t3 = t[:, 2]
    t4 = t[:, 3]
    v1 = v[t1, :]
    v2 = v[t2, :]
    v3 = v[t3, :]
    v4 = v[t4, :]
    e1 = v2 - v1
    e2 = v3 - v2
    e3 = v1 - v3
    e4 = v4 - v1
    e5 = v4 - v2
    e6 = v4 - v3

    # Compute cross product and 6 * vol for each triangle:
    cr  = np.cross(e1,e3)
    vol = np.abs(np.sum(e4*cr, axis=1))
    # zero vol will cause division by zero below, so set to small value:
    vol_mean = 0.0001*np.mean(vol)
    vol[vol==0] = vol_mean

    # compute dot products of edge vectors
    e11 = np.sum(e1*e1,axis=1)
    e22 = np.sum(e2*e2,axis=1)
    e33 = np.sum(e3*e3,axis=1)
    e44 = np.sum(e4*e4,axis=1)
    e55 = np.sum(e5*e5,axis=1)
    e66 = np.sum(e6*e6,axis=1)
    e12 = np.sum(e1*e2,axis=1)
    e13 = np.sum(e1*e3,axis=1)
    e14 = np.sum(e1*e4,axis=1)
    e15 = np.sum(e1*e5,axis=1)
    e23 = np.sum(e2*e3,axis=1)
    e25 = np.sum(e2*e5,axis=1)
    e26 = np.sum(e2*e6,axis=1)
    e34 = np.sum(e3*e4,axis=1)
    e36 = np.sum(e3*e6,axis=1)

    # compute entries for A (negations occur when one edge direction is flipped)
    # these can be computed multiple ways
    # basically for ij, take opposing edge (call it Ek) and two edges from the starting
    # point of Ek to point i (=El) and to point j (=Em), then these are of the
    # scheme:   (El * Ek)  (Em * Ek) - (El * Em) (Ek * Ek)
    # where * is vector dot product
    A12 = (- e36 * e26 + e23 * e66)/vol
    A13 = (- e15 * e25 + e12 * e55)/vol
    A14 = (e23 * e26 - e36 * e22)/vol
    A23 = (- e14 * e34 + e13 * e44)/vol
    A24 = (e13 * e34 - e14 * e33)/vol
    A34 = (- e14 * e13 + e11 * e34)/vol

    # compute diagonals (from row sum = 0)
    A11 = -A12-A13-A14
    A22 = -A12-A23-A24
    A33 = -A13-A23-A34
    A44 = -A14-A24-A34

    # stack columns to assemble data
    localA = np.column_stack((A12, A12, A23, A23, A13, A13, A14, A14, A24, A24, A34, A34, A11, A22, A33, A44))
    I = np.column_stack((t1, t2, t2, t3, t3, t1, t1, t4, t2, t4, t3, t4, t1, t2, t3, t4))
    J = np.column_stack((t2, t1, t3, t2, t1, t3, t4, t1, t4, t2, t4, t3, t1, t2, t3, t4))
    # Flatten arrays:
    I = I.flatten()
    J = J.flatten()
    localA = localA.flatten()/6.0
    # Construct sparse matrix:
    #A = sparse.csr_matrix((localA, (I, J)))
    A = sparse.csc_matrix((localA, (I, J)))

    if not lump:
        # create b matrix data (account for that vol is 6 times tet volume)
        Bii = vol / 60.0
        Bij = vol / 120.0
        localB = np.column_stack((Bij, Bij, Bij, Bij, Bij, Bij, Bij, Bij, Bij, Bij, Bij, Bij, Bii, Bii, Bii, Bii))
        localB = localB.flatten()
        B = sparse.csc_matrix((localB, (I, J)))
    else:
        # when lumping put all onto diagonal (volume/4 for each vertex)
        Bii = vol / 24.0
        localB = np.column_stack((Bii, Bii, Bii, Bii))
        I = np.column_stack((t1, t2, t3, t4))
        I = I.flatten()
        localB = localB.flatten()
        B = sparse.csc_matrix((localB, (I, I)))

    return A, B


# ------------------------------------------------------------------------------
# poissonSolver

def poissonSolver(A, B, h=0.0, dtup=(), ntup=(), useCholmod=True):
    """
    poissonSolver solves the poisson equation with boundary conditions
           based on the A and B Laplace matrices:  A x = B h

    Inputs:   A - sparse sym. (n x n) positive semi definite numpy matrix
              B - sparse sym. (n x n) positive definite numpy matrix (inner product)
                  A and B are obtained via computeAB for either triangle or tet mesh
    Optional:
              h - right hand side, can be constant or array with vertex values
                  Default: 0.0 (Laplace equation A x = 0)
              dtup - Dirichlet boundary condition as a tuple.
                     Tuple contains index and data arrays of same length.
                     Default: no Dirichlet condition
              ntup - Neumann boundary condition as a tuple.
                     Tuple contains index and data arrays of same length.
                     Default: Neumann on all boundaries
              useCholmod - default True: try to find Cholesky decomp from sksparse
                           falls back to scipy sparse splu if not found.

    Outputs:  x - array with vertex values of solution
    """

    # imports
    import numpy as np
    from scipy import sparse
    from scipy.sparse.linalg import splu
    if useCholmod:
        try:
            from sksparse.cholmod import cholesky
        except ImportError:
            useCholmod = False
    # check matrices
    dim = A.shape[0]
    if A.shape != B.shape or A.shape[1] != dim:
        raise ValueError('Error: Square input matrices should have same number of rows and columns')
    # create vector h
    if np.isscalar(h):
        h = np.full((dim,1),h,dtype='float64')
    elif (not np.isscalar(h)) and h.size != dim:
        raise ValueError('h should be either scalar or column vector with row num of A')
    # create vector d
    didx = []
    if dtup:
        if len(dtup) != 2:
            raise ValueError('dtup should contain index and data arrays')
        didx = dtup[0]
        ddat = dtup[1]
        if np.unique(didx).size != len(didx):
            raise ValueError('dtup indices need to be unique')
        if not (len(didx) > 0 and len(didx) == len(ddat)):
            raise ValueError('dtup should contain index and data arrays (same lengths > 0)')
        dvec = sparse.csc_matrix((ddat, (didx, np.zeros(len(didx),dtype=np.uint32))),(dim,1))
    # create vector n
    nvec=0
    if ntup:
        if len(ntup) != 2:
            raise ValueError('ntup should contain index and data arrays')
        nidx = ntup[0]
        ndat = ntup[1]
        if not (len(nidx) > 0 and len(nidx) == len(ndat)):
            raise ValueError('dtup should contain index and data arrays (same lengths > 0)')
        nvec = sparse.csc_matrix((ndat, (nidx, np.zeros(len(nidx),dtype=np.uint32))),(dim,1))
    # compute right hand side
    b = B * (h-nvec)
    if len(didx) >0:
        b= b - A * dvec
    # remove Dirichlet Nodes
    if len(didx)>0:
        mask=np.full(dim,True,dtype=bool)
        mask[didx] = False
        b = b[mask]
        # we need to keep A sparse and do col and row slicing
        # only on the right format:
        if A.getformat() == 'csc':
            A = A[:,mask].tocsr()
            A = A[mask,:]
            A = A.tocsc()
        elif A.getformat() == 'csr':
            A = A[mask,:].tocrc()
            A = A[:,mask]
        else:
            raise ValueError('A matrix needs to be sparse CSC or CSR')
    # solve A x = b
    print("Matrix Format now: "+A.getformat())
    if useCholmod:
        print("Solver: cholesky decomp - performance optimal ...")
        chol = cholesky(A)
        x = chol(b)
    else:
        print("Package scikit-sparse not found (Cholesky decomp)")
        print("Solver: spsolve (LU decomp) - performance not optimal ...")
        lu = splu(A)
        x = np.array(lu.solve(b.astype("float32"))).squeeze()
    # pad Dirichlet nodes
    if len(didx)>0:
        xfull = np.zeros(dim)
        xfull[mask] = x
        xfull[didx] = ddat
        return xfull
    return x

# ------------------------------------------------------------------------------
# computeABtria()

def computeABtria(v, t, lump = False):
    """
    computeABtria(v,t) computes the two sparse symmetric matrices representing
           the Laplace Beltrami Operator for a given triangle mesh using
           the linear finite element method (assuming a closed mesh or
           the Neumann boundary condition).

    Inputs:   v - vertices : list of lists of 3 floats
              t - triangles: list of lists of 3 int of indices (>=0) into v array

    Outputs:  A - sparse sym. (n x n) positive semi definite numpy matrix
              B - sparse sym. (n x n) positive definite numpy matrix (inner product)

    Can be used to solve sparse generalized Eigenvalue problem: A x = lambda B x
    or to solve Poisson equation: A x = B f (where f is function on mesh vertices)
    or to solve Laplace equation: A x = 0
    or to model the operator's action on a vector x:   y = B\(Ax)
    """

    # Imports
    import sys
    import numpy as np
    from scipy import sparse

    #
    v = np.array(v)
    t = np.array(t)

    # Compute vertex coordinates and a difference vector for each triangle:
    t1 = t[:, 0]
    t2 = t[:, 1]
    t3 = t[:, 2]
    v1 = v[t1, :]
    v2 = v[t2, :]
    v3 = v[t3, :]
    v2mv1 = v2 - v1
    v3mv2 = v3 - v2
    v1mv3 = v1 - v3

    # Compute cross product and 4*vol for each triangle:
    cr  = np.cross(v3mv2,v1mv3)
    vol = 2 * np.sqrt(np.sum(cr*cr, axis=1))
    # zero vol will cause division by zero below, so set to small value:
    vol_mean = 0.0001*np.mean(vol)
    vol[vol<sys.float_info.epsilon] = vol_mean

    # compute cotangents for A
    # using that v2mv1 = - (v3mv2 + v1mv3) this can also be seen by
    # summing the local matrix entries in the old algorithm
    A12 = np.sum(v3mv2*v1mv3,axis=1)/vol
    A23 = np.sum(v1mv3*v2mv1,axis=1)/vol
    A31 = np.sum(v2mv1*v3mv2,axis=1)/vol

    # compute diagonals (from row sum = 0)
    A11 = -A12-A31
    A22 = -A12-A23
    A33 = -A31-A23

    # stack columns to assemble data
    localA = np.column_stack((A12, A12, A23, A23, A31, A31, A11, A22, A33))
    I = np.column_stack((t1, t2, t2, t3, t3, t1, t1, t2, t3))
    J = np.column_stack((t2, t1, t3, t2, t1, t3, t1, t2, t3))

    # Flatten arrays
    I = I.flatten()
    J = J.flatten()
    localA = localA.flatten()

    # Construct sparse matrix:
    A = sparse.csc_matrix((localA, (I, J)))

    #
    if not lump:
        # create b matrix data (account for that vol is 4 times area)
        Bii = vol / 24
        Bij = vol / 48
        localB = np.column_stack((Bij, Bij, Bij, Bij, Bij, Bij, Bii, Bii, Bii))
        localB = localB.flatten()
        B = sparse.csc_matrix((localB, (I, J)))
    else:
        # when lumping put all onto diagonal  (area/3 for each vertex)
        Bii = vol / 12
        localB = np.column_stack((Bii, Bii, Bii))
        I = np.column_stack((t1, t2, t3))
        I = I.flatten()
        localB = localB.flatten()
        B = sparse.csc_matrix((localB, (I, I)))

    #
    return A, B


# ------------------------------------------------------------------------------
# computeABtriaAniso

def computeABtriaAniso(v, t, u1, u2, D, lump = False):
    """
    computeABtria(v,t) computes the two sparse symmetric matrices representing
           the Laplace Beltrami Operator for a given triangle mesh using
           the linear finite element method (assuming a closed mesh or
           the Neumann boundary condition).

    Inputs:   v  - vertices : list of lists of 3 floats
              t  - triangles: N list of lists of 3 int of indices (>=0) into v array
              u1 - min curv:  min curvature direction per triangle (Nx3 floats)
              u2 - max curv:  max curvature direction per triangle (Nx3 floats)
              D  - anisotropy matrix: diagonal elelments in u1,u2 basis per
                              triangle (Nx2 floats)

    Outputs:  A  - sparse sym. (n x n) positive semi definite numpy matrix
              B  - sparse sym. (n x n) positive definite numpy matrix (inner product)

    Can be used to solve sparse generalized Eigenvalue problem: A x = lambda B x
    or to solve Poisson equation: A x = B f (where f is function on mesh vertices)
    or to solve Laplace equation: A x = 0
    or to model the operator's action on a vector x:   y = B\(Ax)
    """

    # Imports
    import numpy as np
    from scipy import sparse
    import sys

    #
    v = np.array(v)
    t = np.array(t)

    # Compute vertex coordinates and a difference vector for each triangle:
    t1 = t[:, 0]
    t2 = t[:, 1]
    t3 = t[:, 2]
    v1 = v[t1, :]
    v2 = v[t2, :]
    v3 = v[t3, :]
    v2mv1 = v2 - v1
    v3mv2 = v3 - v2
    v1mv3 = v1 - v3

    # transform edge e to basis U = (U1,U2) via U^T * e
    # Ui is n x 3, e is n x 1, result is n x 2
    uv2mv1 = np.column_stack((np.sum(u1*v2mv1,axis=1),np.sum(u2*v2mv1,axis=1)))
    uv3mv2 = np.column_stack((np.sum(u1*v3mv2,axis=1),np.sum(u2*v3mv2,axis=1)))
    uv1mv3 = np.column_stack((np.sum(u1*v1mv3,axis=1),np.sum(u2*v1mv3,axis=1)))

    # Compute cross product and 4*vol for each triangle:
    cr  = np.cross(v3mv2,v1mv3)
    vol = 2 * np.sqrt(np.sum(cr*cr, axis=1))
    # zero vol will cause division by zero below, so set to small value:
    vol_mean = 0.0001*np.mean(vol)
    vol[vol<sys.float_info.epsilon] = vol_mean

    # compute cotangents for A
    # using that v2mv1 = - (v3mv2 + v1mv3) this can also be seen by
    # summing the local matrix entries in the old algorithm
    # Also: here D is the two diagonal entries, not full matrices
    A12 = np.sum(uv3mv2*D*uv1mv3,axis=1)/vol
    A23 = np.sum(uv1mv3*D*uv2mv1,axis=1)/vol
    A31 = np.sum(uv2mv1*D*uv3mv2,axis=1)/vol

    # compute diagonals (from row sum = 0)
    A11 = -A12-A31
    A22 = -A12-A23
    A33 = -A31-A23

    # stack columns to assemble data
    localA = np.column_stack((A12, A12, A23, A23, A31, A31, A11, A22, A33))
    I = np.column_stack((t1, t2, t2, t3, t3, t1, t1, t2, t3))
    J = np.column_stack((t2, t1, t3, t2, t1, t3, t1, t2, t3))

    # Flatten arrays:
    I = I.flatten()
    J = J.flatten()
    localA = localA.flatten()

    # Construct sparse matrix:
    A = sparse.csc_matrix((localA, (I, J)))

    #
    if not lump:
        # create b matrix data (account for that vol is 4 times area)
        Bii = vol / 24
        Bij = vol / 48
        localB = np.column_stack((Bij, Bij, Bij, Bij, Bij, Bij, Bii, Bii, Bii))
        localB = localB.flatten()
        B = sparse.csc_matrix((localB, (I, J)))
    else:
        # when lumping put all onto diagonal  (area/3 for each vertex)
        Bii = vol / 12
        localB = np.column_stack((Bii, Bii, Bii))
        I = np.column_stack((t1, t2, t3))
        I = I.flatten()
        localB = localB.flatten()
        B = sparse.csc_matrix((localB, (I, I)))

    #
    return A, B

# ------------------------------------------------------------------------------
# tria_rm_free_vertices

def tria_rm_free_vertices(v,t):
    """
    Remove unused (free) vertices from v and t. These are vertices that are not
    used in any triangle. They can produce problems when constructing, e.g.,
    Laplace matrices.

    Inputs:   v - vertices   List of lists of 3 float coordinates
              t              List of lists of 3 or 4 int of indices (>=0)
                             into v array (works for trias and tetras)

    Output:   vnew           New list of vertices
              tnew           New list of trias/tetras
              vkeep          Indices (from original list) of kept vertices
              vdel           Indices of deleted (unused) vertices
    """
    # import
    import numpy as np
    #
    tflat = t.flatten()
    vnum = np.max(v.shape)
    if (np.max(tflat) >= vnum):
        raise ValueError('Max index exceeds number of vertices')
    # determine which vertices to keep
    vkeep = np.full(vnum,False,dtype=bool)
    vkeep[tflat] = True
    # list of deleted vertices (old indices)
    vdel = np.nonzero(~vkeep)[0]
    # if nothing to delete return
    if len(vdel) == 0:
        return v,t , np.arange(vnum), []
    # delete unused vertices
    vnew = v[vkeep,:]
    # create lookup table
    tlookup = np.cumsum(vkeep)-1
    # reindex tria
    tnew = tlookup[t]
    # convert vkeep to index list
    vkeep = np.nonzero(vkeep)[0]
    return vnew, tnew, vkeep, vdel


# ------------------------------------------------------------------------------
# tria_is_oriented

def tria_is_oriented(t):
    """
    Check if triangle mesh is oriented. True if all triangles are oriented
    counter-clockwise, when looking from above.

    Inputs:   t - trias      List of lists of 3 int of indices (>=0) into v array

    Output:   oriented       bool True if max(adj_directed)=1
    """
    # import
    import numpy as np
    #
    adj = tria_get_adj(t,directed=True)
    oriented = (np.max(adj.data) == 1)
    return oriented

# ------------------------------------------------------------------------------
# tria_get_adj

def tria_get_adj(t,directed=False,store_tidx=False):
    """
    Compute the adjacency matrix (edge graph) of the triangle mesh t

    Inputs:   t - trias      List of lists of 3 int of indices (>=0) into v array
                             Ordering is important: All triangles should be
                             oriented the same way (counter-clockwise, when
                             looking from above)

              directed       False (default): The non-directed adjacency matrix
                             will be symmetric. Each inner edge (i,j) will have
                             the number of triangles that contain this edge.
                             Inner edges usually 2, boundary edges 1. Higher
                             numbers can occur when there are non-manifold triangles.

                             True: The directed adjacency matrix is not symmetric.
                             For manifold meshes, there are only entries with
                             value 1. Symmetric entries are inner edges. Non-symmetric
                             are boundary edges. The direction prescribes a direction
                             on the boundary loops.

              store_tidx     If true, stores the tria idx +1 instead of one
                             in the matrix (allows lookup of vertex to tria).
                             Works only on directed adj matrix, and manifold.

    Output:   adj            The sparse adjacency matrix in CSC format
                             Note: for non-directed, adj is symmetric and
                             contains a 2 at inner edges and 1 at boundaries.
                             It can be binarized via:
                                adj.data = np.ones(adj.data.shape)
                             For directed adj it contains 1 for each directed
                             edge and is nonsymmetric if boundaries exist.
                             Adding this to its transpose creates the non-
                             directed version.

    """

    # import
    import numpy as np
    from scipy import sparse
    #
    t1 = t[:,0]
    t2 = t[:,1]
    t3 = t[:,2]
    if directed:
        I = np.column_stack((t1, t2, t3))
        J = np.column_stack((t2, t3, t1))
    else:
        if store_tidx:
            raise ValueError('store_tidx requires directed adjacency matrix')
        I = np.column_stack((t1, t2, t2, t3, t3, t1))
        J = np.column_stack((t2, t1, t3, t2, t1, t3))
    I = I.flatten()
    J = J.flatten()
    dat = np.ones(I.shape)
    if store_tidx:
        # store tria idx +1  (zero means no edge here)
        dat = np.repeat(np.arange(1,t.shape[0]+1),3)
    adj = sparse.csc_matrix((dat, (I, J)))
    return adj

# ------------------------------------------------------------------------------
# tetra_get_boundary_tria

def tetra_get_boundary_tria(v,t,tetfunc = None):
    """
    Get boundary triangle mesh of tetrahedra (can have multiple connected components).

    Inputs:   v - vertices   List of lists of 3 float coordinates
              t              List of lists of 4 int of indices (>=0) into v array
              tetfunc        List of tetra function values

    Output:   tria           List of lists of 3 int: triangle mesh with idx into v
              triafunc       List of tria function values (mapped from tetfunc)
    """
    # import
    import numpy as np
    #
    v = np.array(v)
    t = np.array(t)
    # get all triangles
    allT = np.vstack( (t[:,np.array([3, 1, 2])] , t[:,np.array([2, 0, 3])], t[:,np.array([1, 3, 0])], t[:,np.array([0, 2, 1])]) )
    # sort rows so that faces are reorder in ascending order of indices
    allT.sort(axis=1)
    # find unique trias without a neighbor
    tria, indices, count = np.unique(allT, axis=0, return_index=True, return_counts=True)
    tria = tria[count == 1, :]
    print('Found '+str(np.size(tria,0))+' triangles on boundary.')
    # if we have tetra function, map these to the boundary triangles
    if tetfunc is not None:
        allTidx = np.tile(np.arange(t.shape[0]),4)
        Tidx = allTidx[indices[count == 1]]
        triafunc = tetfunc[Tidx]
        return tria, triafunc
    #
    return tria

# ------------------------------------------------------------------------------
# orientTria()

def orientTria(v, t, start=0, up=(1,1,1)):

    # imports

    import sys
    import numpy as np

    # subfunctions

    def getTriaNeighbors(t, idx):

        nb = np.where(np.sum(np.isin(t, t[idx,]),axis=1)==2)[0]
        return(nb)

    def alignTria(tt, tri):

        # restrict list to trias that have no nans (i.e. have been visited already)
        tt = tt[np.sum(np.isnan(tt), axis=1)==0,:]

        # restrict list to trias that share an edge with tri
        tt = tt[np.where(np.sum(np.isin(tt, tri), axis=1)==2)[0],:]

        # initialise list of inversions
        doInvert = np.zeros(len(tt))

        # loop through tt
        for i in range(0, len(tt)):
            # get shared edge
            e = tt[i, np.isin(tt[i,:], tri)]
            # compute positions of vertices of shared edge
            tt_i = (np.where(tt[i,:] == e[0])[0][0], np.where(tt[i,:] == e[1])[0][0])
            tri_i = (np.where(tri == e[0])[0][0], np.where(tri == e[1])[0][0])
            # invert tri if positions follow same direction
            if (tt_i == (0,1)) | (tt_i == (1,2)) | (tt_i == (2,0)):
                if (tri_i == (0,1)) | (tri_i == (1,2)) | (tri_i == (2,0)):
                    doInvert[i] = 1
            elif (tt_i == (1,0)) | (tt_i == (2,1)) | (tt_i == (0,2)):
                if (tri_i == (1,0)) | (tri_i == (2,1)) | (tri_i == (0,2)):
                    doInvert[i] = 1

        # invert
        if all(doInvert):
            tri = np.array((tri[0],tri[2],tri[1]))
        if any(doInvert) and not all(doInvert):
            print("Inversion problem")
            sys.exit(1)

        # return
        return(tri)

    # create array of nans for new trias

    newTria = np.empty(shape=np.shape(t), dtype="int64") * np.nan

    # create list of visited trias

    triaVisited = list()

    # create tria stack (list)

    triaStack = list()

    # add initial tria to visited tria array

    triaVisited.append(start)

    # check if orientation of initial tria conforms to 'up', if not, invert

    v1mv0 = v[t[start][1],] - v[t[start][0],]
    v2mv0 = v[t[start][2],] - v[t[start][0],]

    n = np.cross(v1mv0, v2mv0)

    if np.dot(n, up) < 0:
        newTria[start] = t[start,(0,2,1)]
    else:
        newTria[start] = t[start]

    # get neighbors of initial triangle

    nb = getTriaNeighbors(t, start)

    # push neighbors of initial triangle to tria stack

    triaStack.extend(nb)

    # unless tria stack is empty, continue

    while len(triaStack) > 0:

        # take first tria from tria stack and align with already visited trias

        tri = alignTria(newTria, t[triaStack[0],:])

        # add (possibly inverted) tria to newTria

        newTria[triaStack[0],:] = tri

        # add neighbors of current tria to stack (exclude those that are
        # already on the stack and those that have previously been removed from
        # the stack; this means that list of neighbors can be empty)

        nb = getTriaNeighbors(t, triaStack[0])

        triaStack.extend(np.setdiff1d(np.setdiff1d(nb, triaVisited), triaStack))

        # add current tria to visited tria

        triaVisited.append(triaStack[0])

        # remove current tria from stack

        triaStack.pop(0)

    # check

    if np.isnan(newTria).any():
        print("NaNs in newTria")
        sys.exit(1)

    # convert to int (was float due to np.nan)

    newTria = newTria.astype(np.int)

    # return

    return(newTria)
