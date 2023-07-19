"""
This module provides a function to create a label files for a tetrahedral mesh

"""

# -----------------------------------------------------------------------------
# AUXILIARY FUNCTIONS

# createBoundaryMask()

def createBoundaryMask(hbtFile, mskFile, outFile, label, bndlabel):
    """
    createBoundaryMask(hbtFile, mskFile)

    """

    # imports
    import numpy as np
    import nibabel as nb

    # create a new label for boundary voxels
    print('\nCreating a new label for boundary voxels ...')
    print('... Reading ' + hbtFile)

    # load segmentations
    hbt = nb.load(hbtFile)
    msk = nb.load(mskFile)

    # define new data structures
    bnd = msk
    bnd_data = bnd.get_fdata(caching='fill')
    hbt_data_flat = hbt.get_fdata().flatten()

    # find nonzero elements of msk
    mskx, msky, mskz = np.where(msk.get_fdata())

    # cycle through each nonzero element of msk
    print('... Processing mask')
    for i in range(0, len(mskx)):

        # create indices xyz of local neighborhood (one more +1 due to python indexing)
        mx, my, mz = np.meshgrid(range(mskx[i] - 1, mskx[i] + 1 + 1), range(msky[i] - 1, msky[i] + 1 + 1), range(mskz[i] - 1, mskz[i] + 1 + 1))

        # convert xyz to i
        # sel = range(0,27) # point conn
        # sel = range(0,len(mx.flatten())) # use this if increasing the neighborhood size above, e.g. to 3 or 6
        sel = np.array((5, 11, 13, 14, 15, 17, 23)) - 1  # face conn
        m = np.ravel_multi_index([mx.flatten()[sel], my.flatten()[sel], mz.flatten()[sel]], msk.shape)

        # get values from hbt and set bnd to 2260 if it contains any tail voxel
        # (226) and to 2320 if it contains any head voxel (232)

        for j in range(0, len(label)):
            if any(np.isin(hbt_data_flat[m], label[j])):
                bnd_data[mskx[i], msky[i], mskz[i]] = bndlabel[j]

    # write boundary mask; since we have used 'caching='fill'' while loading
    # the data, the contents of bnd will be automatically updated via bnd_data
    out = nb.MGHImage(bnd_data, affine=bnd.affine, header=bnd.header)
    nb.freesurfer.save(out, outFile)


# ------------------------------------------------------------------------------
# MAIN FUNCTION

# createTetraLabels

def createTetraLabels(params):
    """

    """

    # imports

    import os
    import subprocess

    import numpy as np
    import nibabel as nb

    from lapy import TriaMesh, TetMesh, io

    # message

    print()
    print("-------------------------------------------------------------------")
    print()
    print("Creating label files for tetrahedral meshes")
    print()
    print("-------------------------------------------------------------------")
    print()

    # determine files

    MSKFILE = os.path.join(params.OUTDIR, params.HEMI + ".labels.mgz")

    if params.internal.AUTOMASK_HEAD is True or params.internal.AUTOMASK_TAIL is True:
        HBTFILE = os.path.join(params.OUTDIR, "labels", params.HEMI + ".automask_image.mgz")
    else:
        HBTFILE = os.path.join(params.OUTDIR, params.HEMI + ".image.mgz")

    # create new masks that includes boundary labels to head and tail, and to CA4

    createBoundaryMask(hbtFile=HBTFILE, mskFile=MSKFILE, outFile=os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + "labels-bnd.mgz"), label=[params.LUTDICT['tail'], params.LUTDICT['head']], bndlabel=[params.LUTDICT["bndtail"], params.LUTDICT["bndhead"]])
    createBoundaryMask(hbtFile=HBTFILE, mskFile=MSKFILE, outFile=os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + "labels-ca4.mgz"), label=[params.LUTDICT['ca4']], bndlabel=[params.LUTDICT["bndca4"]])

    # merge bnd and ca4 masks (add ca4 to bnd)

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mris_calc") + " " \
        + "--output " \
        + os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + "labels-bnd.mgz") + " " \
        + os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + "labels-bnd.mgz") + " " \
        + "lowerlimit " \
        + os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + "labels-ca4.mgz")

    print(cmd)

    subprocess.run(cmd.split())

    # load bnd image
    
    bndImage = nb.freesurfer.load(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".labels-bnd.mgz"))

    # load tet mesh

    tetMesh = TetMesh.read_vtk(os.path.join(params.OUTDIR, params.HEMI + ".tetra.vtk"))    

    # convert vertex indices to voxel indices

    vtxVxIdx = np.round(np.matmul(np.linalg.inv(bndImage.header.get_vox2ras_tkr()), np.concatenate((tetMesh.v, np.ones((tetMesh.v.shape[0], 1))), axis=1).T).T).astype(int)[:, 0:3]

    # lookup labels of converted vertex indices in bnd image

    vtxLabels = bndImage.get_fdata()[vtxVxIdx[:, 0], vtxVxIdx[:, 1], vtxVxIdx[:, 2]]
                       
    # write out psol file
 
    with open(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".tetra.psol"), "w") as f:

        print("Solution:", file=f)
        [ print(str(int(x)), file=f) for x in vtxLabels ]

    # write out visualization files

    tetMeshBnd = tetMesh.boundary_tria()
    vRmBndKeep, vBndRmDel = tetMeshBnd.rm_free_vertices_()
    tetMeshBnd.orient_()
    TriaMesh.write_vtk(tetMeshBnd, os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".rm.bnd.tetra.vtk"))
    io.write_vfunc(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".rm.bnd.tetra.psol"), np.array([ float(x) for x in vtxLabels ])[vRmBndKeep])
    nb.freesurfer.save(nb.freesurfer.MGHImage(dataobj=np.array([ float(x) for x in vtxLabels ])[vRmBndKeep].astype("float32"), affine=None), filename=os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".rm.bnd.label.mgh"))

    # return

    return(params)


def createTetraLabelsRemove(params):
    """

    """

    # imports

    import os
    import subprocess

    import numpy as np
    import nibabel as nb

    from lapy import TetMesh, io

    # message

    print()
    print("-------------------------------------------------------------------")
    print()
    print("Creating label files for tetrahedral meshes")
    print()
    print("-------------------------------------------------------------------")
    print()

    # settings

    jointtail = 226
    jointhead = 232
    bndtail = 2260
    bndhead = 2320
    bndca4 = 2420

    # determine files

    MSKFILE = os.path.join(params.OUTDIR, params.HEMI + ".labels.mgz")
    HBTFILE = os.path.join(params.OUTDIR, params.HEMI + ".image.mgz")

    # create vertex labels 0a: create new masks that includes boundary labels to head and tail, and to CA4

    createBoundaryMask(hbtFile=HBTFILE, mskFile=MSKFILE, outFile=os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + "labels-bnd.mgz"), label=[params.LUTDICT['tail'], params.LUTDICT['head']], bndlabel=[bndtail, bndhead])
    createBoundaryMask(hbtFile=HBTFILE, mskFile=MSKFILE, outFile=os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + "labels-ca4.mgz"), label=[params.LUTDICT['ca4']], bndlabel=[bndca4])

    # create vertex labels 0b: merge bnd and ca4 masks (add ca4 to bnd)

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mris_calc") + " " \
        + "--output " \
        + os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + "labels-bnd.mgz") + " " \
        + os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + "labels-bnd.mgz") + " " \
        + "lowerlimit " \
        + os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + "labels-ca4.mgz")

    print(cmd)

    subprocess.run(cmd.split())

    # create vertex labels 1a: create label file from mask file

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_vol2label") + " " \
        + "--i " + os.path.join(params.OUTDIR, params.HEMI + ".mask.mgz")  + " " \
        + "--id 1 " \
        + "--l " + os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".mask.label")

    print(cmd)

    subprocess.run(cmd.split())

    # create vertex labels 1b: rewrite label file (update header line)

    with open(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".mask.label")) as g:

        lab = g.readlines()

    with open(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".mask.label"), "w") as f:

        print("#!ascii label " + params.HEMI + ".mask.label", file=f)

        [ print(x.strip(), file=f) for x in lab[1:] ]

    # create vertex labels 1c: create intermediate label file from bnd image file (uses voxel coordinates)

    assigned_bnd = nb.freesurfer.load(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".labels-bnd.mgz"))
    assigned_bnd_v2r_tkr = assigned_bnd.header.get_vox2ras_tkr()
    assigned_bnd_data = assigned_bnd.get_fdata()

    labNum = [ x.rstrip().split() for x in lab[2:] ]
    labNum = np.array(labNum).astype(float)[:,1:]

    labMat = np.concatenate((labNum[:,0:3], np.ones((len(labNum), 1))), axis=1)
    labMatXYZ = np.round(np.matmul(np.linalg.inv(assigned_bnd_v2r_tkr), labMat.transpose()).transpose()[:,0:3]).astype(int)

    assigned_bnd_labels = np.array(assigned_bnd_data[labMatXYZ[:,0],labMatXYZ[:,1],labMatXYZ[:,2]], ndmin=2).transpose()
    labMatList = np.concatenate((labMatXYZ, assigned_bnd_labels), axis=1).astype(int)

    with open(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".labels-bnd.list"), "w") as f:

        [ print(str(x[0])+" "+str(x[1])+" "+str(x[2])+" "+str(x[3]), file=f) for x in labMatList ]

    # create vertex labels 2: create label file from tetra file (does not include labels yet)

    tetMesh = TetMesh.read_vtk(os.path.join(params.OUTDIR, params.HEMI + ".tetra.vtk"))

    with open(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".tetra.label"), "w") as f:

        print("#!ascii label " + params.HEMI + ".tetra.label", file=f)
        print(str(np.shape(tetMesh.v)[0]), file=f)
        [ print("-1\t" + str(x).replace("[","").replace("]","") + "\t0", file=f) for x in tetMesh.v ]

    # create vertex labels 3a: 

    labFile = os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".tetra.label")
    lstFile = os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + "labels-bnd.list")

    lst = np.loadtxt(lstFile, dtype="float")
    lab = np.loadtxt(labFile, dtype="float", skiprows=2)

    crs = np.array(lst)  # to avoid mere referencing
    crs[:, 3] = 1

    ras = np.matmul(assigned_bnd_v2r_tkr, crs.transpose()).transpose()
    ras = np.append(ras[:, 0:3], lst[:, 3:4], axis=1)

    vtx = np.zeros(shape=np.shape(lab)[0])
    for i in range(np.shape(vtx)[0]):
        tmp = np.linalg.norm(ras[:, 0:3] - np.repeat(lab[i:(i + 1), 1:4], np.shape(ras)[0], axis=0), ord=2, axis=1)
        vtx[i] = lst[np.where(tmp == tmp.min())[0][0], 3]  # note: we added [0][0] here

    # the following lines will produce zero-mean, step-one labels
    # key = np.array(range(0, len(np.unique(vtx)))) - np.mean(range(0, len(np.unique(vtx))))
    # index = np.digitize(vtx, np.unique(vtx), right=True)
    # np.savetxt(fname=labFile.replace(".label", ".asc"), X=key[index])

    # the following line will keep original labels
    np.savetxt(fname=labFile.replace(".label", ".asc"), X=vtx)

    # create vertex labels 3b: create psol overlay file

    with open(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".tetra.asc"), "r") as g:

        asc = g.read().splitlines()

    with open(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".tetra.psol"), "w") as f:

        print("Solution:", file=f)
        [ print(str(float(x.strip())).replace("[","(").replace("]",")"), file=f) for x in asc ]

    # write out visualization files

    tetMeshBnd = tetMesh.boundary_tria()
    vRmBndKeep, vBndRmDel = tetMeshBnd.rm_free_vertices_()
    tetMeshBnd.orient_()
    TetMesh.write_vtk(tetMeshBnd, os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".rm.bnd.tetra.vtk"))
    io.write_vfunc(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".rm.bnd.tetra.psol"), np.array([ float(x) for x in asc ])[vRmBndKeep])
    nb.freesurfer.save(nb.freesurfer.MGHImage(dataobj=np.array([ float(x) for x in asc ])[vRmBndKeep].astype("float32"), affine=None),
        filename=os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".rm.bnd.label.mgh"))

    # clean up

    if params.skipCLEANUP is False:

        os.remove(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + "labels-bnd.list"))
        os.remove(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".tetra.asc"))

    # return

    return(params)
