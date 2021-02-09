"""
This module provides a function to create a label files for a tetrahedral mesh

"""

# -----------------------------------------------------------------------------
# AUXILIARY FUNCTIONS

# createBoundaryMask()

def createBoundaryMask(hbtFile, mskFile, suffix="bnd", label=[226, 232]):
    """
    createCylMsk(hbtFile, mskFile)

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
        mx, my, mz = np.meshgrid(range(mskx[i] - 1, mskx[i] + 1 + 1), range(msky[i] - 1, msky[i] + 1 + 1),
                                 range(mskz[i] - 1, mskz[i] + 1 + 1))

        # convert xyz to i
        # sel = range(0,27) # point conn
        # sel = range(0,len(mx.flatten())) # use this if increasing the neighborhood size above, e.g. to 3 or 6
        sel = np.array((5, 11, 13, 14, 15, 17, 23)) - 1  # face conn
        m = np.ravel_multi_index([mx.flatten()[sel], my.flatten()[sel], mz.flatten()[sel]], msk.shape)

        # get values from hbt and set bnd to 2260 if it contains any tail voxel
        # (226) and to 2320 if it contains any head voxel (232) (or to k*10 for
        # custom label k)

        for j in label:
            if any(hbt_data_flat[m] == j):
                bnd_data[mskx[i], msky[i], mskz[i]] = j * 10

        # TODO: evtl noch auf largest connected component einschraenken ...

    # write boundary mask; since we have used 'caching='fill'' while loading
    # the data, the contents of bnd will be automatically updated via bnd_data
    print('... Writing ' + mskFile.replace('.mgz', '-' + suffix + '.mgz'))
    out = nb.MGHImage(bnd_data, affine=bnd.affine, header=bnd.header)
    nb.freesurfer.save(out, mskFile.replace('.mgz', '-' + suffix + '.mgz'))


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

    from shapetools.createVertexLabels import createVertexLabels
    from shapetools.triaUtils import readVTK, writeVTK, writePSOL
    from shapetools.triaUtils import tetra_get_boundary_tria, tria_rm_free_vertices

    # message

    print()
    print("-------------------------------------------------------------------")
    print()
    print("Creating label files for tetrahedral meshes")
    print()
    print("-------------------------------------------------------------------")
    print()

    # determine files

    if params.LUT == "fs711":

        HBTFILE = os.path.join(params.SUBJDIR, params.SUBJID, "mri", params.HEMI + "." + params.SFX + ".HBT.mgz")
        MSKFILE = os.path.join(params.OUTDIR, "mask", params.HEMI + "." + params.internal.HSFLABEL_02 + "_assigned.mgz")
        CA4FILE = os.path.join(params.SUBJDIR, params.SUBJID, "mri", params.HEMI + "." + params.SFX + ".mgz")

        print("Using " + HBTFILE + " as head-body-tail file.")
        print("Using " + MSKFILE + " as mask file.")
        print()

        # create vertex labels 0a: create new mask that includes boundary labels

        createBoundaryMask(hbtFile=HBTFILE, mskFile=MSKFILE, suffix="bnd", label=[226, 232])

        createBoundaryMask(hbtFile=CA4FILE, mskFile=MSKFILE, suffix="ca4", label=[242])

    elif params.LUT == "ashs":

        HBTFILE = os.path.join(params.SUBJDIR, params.SUBJID, "mri", params.HEMI + "." + params.SFX + ".mgz")
        MSKFILE = os.path.join(params.OUTDIR, "mask", params.HEMI + "." + params.internal.HSFLABEL_02 + "_assigned.mgz")
        CA4FILE = os.path.join(params.SUBJDIR, params.SUBJID, "mri", params.HEMI + "." + params.SFX + ".mgz")

        createBoundaryMask(hbtFile=HBTFILE, mskFile=MSKFILE, suffix="bnd", label=[254, 255])

        createBoundaryMask(hbtFile=CA4FILE, mskFile=MSKFILE, suffix="ca4", label=[3])

    # create vertex labels 0b:

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mris_calc") + " " \
        + "--output " \
        + MSKFILE.replace('.mgz', '-bnd.mgz') + " " \
        + MSKFILE.replace('.mgz', '-bnd.mgz') + " " \
        + "lowerlimit " \
        + MSKFILE.replace('.mgz', '-ca4.mgz')

    print(cmd)

    subprocess.run(cmd.split())

    # create vertex labels 1a:

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_vol2label") + " " \
        + "--i " + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_02 + ".mgz")  + " " \
        + "--id 1 " \
        + "--l " + os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + params.internal.HSFLABEL_02 + "_assigned-bnd.tmp.label")

    print(cmd)

    subprocess.run(cmd.split())

    # create vertex labels 1b:

    with open(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + params.internal.HSFLABEL_02 + "_assigned-bnd.tmp.label")) as g:

        lab = g.readlines()

    with open(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + params.internal.HSFLABEL_02 + "_assigned-bnd.label"), "w") as f:

        print("#!ascii label " + params.HEMI + "." + params.internal.HSFLABEL_02 + "_assigned-bnd.label", file=f)

        [ print(x.strip(), file=f) for x in lab[1:] ]

    # create vertex labels 1c (only works with FS6.0):

    #cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_vol2roi") + " " \
    #    + "--srcvol " + os.path.join(params.OUTDIR, "mask", params.HEMI + "." + params.internal.HSFLABEL_02 + "_assigned-bnd.mgz") + " " \
    #    + "--label " + os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + params.internal.HSFLABEL_02 + "_assigned-bnd.label") + " " \
    #    + "--list " + os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + params.internal.HSFLABEL_02 + "_assigned-bnd.list")

    #print(cmd)

    #subprocess.run(cmd.split())

    # create vertex labels 1c (works with FS6.0 and FS7.11):

    assigned_bnd = nb.freesurfer.load(os.path.join(params.OUTDIR, "mask", params.HEMI + "." + params.internal.HSFLABEL_02 + "_assigned-bnd.mgz"))
    assigned_bnd_v2r_tkr = assigned_bnd.header.get_vox2ras_tkr()
    assigned_bnd_data = assigned_bnd.get_fdata()

    labNum = [ x.rstrip().split() for x in lab[2:] ]
    labNum = np.array(labNum).astype(float)[:,1:]

    labMat = np.concatenate((labNum[:,0:3], np.ones((len(labNum), 1))), axis=1)
    labMatXYZ = np.round(np.matmul(np.linalg.inv(assigned_bnd_v2r_tkr), labMat.transpose()).transpose()[:,0:3]).astype(int)

    assigned_bnd_labels = np.array(assigned_bnd_data[labMatXYZ[:,0],labMatXYZ[:,1],labMatXYZ[:,2]], ndmin=2).transpose()
    labMatList = np.concatenate((labMatXYZ, assigned_bnd_labels), axis=1).astype(int)

    with open(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + params.internal.HSFLABEL_02 + "_assigned-bnd.list"), "w") as f:

        [ print(str(x[0])+" "+str(x[1])+" "+str(x[2])+" "+str(x[3]), file=f) for x in labMatList ]

    # create vertex labels 1d:

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_info") + " " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + "." + params.internal.HSFLABEL_02 + "_assigned-bnd.mgz") + " " \
        + "--vox2ras-tkr "

    print(cmd)

    vox2ras = subprocess.check_output(cmd.split(), universal_newlines=True)

    with open(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + params.internal.HSFLABEL_02 + "_assigned-bnd.vox2ras"), "w") as f:

        print(vox2ras.strip(), file=f)

#    # create vertex labels 1e:
#
#    with open(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + params.internal.HSFLABEL_02 + "_assigned-bnd.list"), "r") as g:
#
#        lst = g.readlines()
#
#    with open(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + params.internal.HSFLABEL_02 + "_assigned-bnd.label"), "w") as f:
#
#        print(str(len(lst)), file=f)
#        [ print(x.strip(), file=f) for x in lab[1:] ]

    # create vertex labels 2

    v, t = readVTK(os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_07 + ".vtk"))

    with open(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + params.internal.HSFLABEL_07 + ".label"), "w") as f:

        print("#!ascii label " + params.HEMI + "." + params.internal.HSFLABEL_07 + ".label", file=f)
        print(str(np.shape(v)[0]), file=f)
        [ print("-1\t" + str(x).replace("[","").replace("]","") + "\t0", file=f) for x in v ]

    # create vertex labels 3a

    createVertexLabels(
        labFile=os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + params.internal.HSFLABEL_07 + ".label"),
        matFile=os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + params.internal.HSFLABEL_02 + "_assigned-bnd.vox2ras"),
        lstFile=os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + params.internal.HSFLABEL_02 + "_assigned-bnd.list"))

    # create vertex labels 3b: create psol overlay file

    with open(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + params.internal.HSFLABEL_07 + ".asc"), "r") as g:

        asc = g.read().splitlines()

    with open(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + params.internal.HSFLABEL_07 + ".psol"), "w") as f:

        print("Solution:", file=f)
        [ print(str(float(x.strip())).replace("[","(").replace("]",")"), file=f) for x in asc ]

    # write out visualization files

    tBnd = tetra_get_boundary_tria(v, t)

    vBndRm, tBndRm, vRmBndKeep, vBndRmDel = tria_rm_free_vertices(v,tBnd)

    writeVTK(outfile=os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".rm.bnd." + params.internal.HSFLABEL_07 + ".vtk"), v=vBndRm, t=tBndRm)

    writePSOL(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".rm.bnd." + params.internal.HSFLABEL_07 + ".psol"), np.array([ float(x) for x in asc ])[vRmBndKeep])

    # clean up

    if params.skipCLEANUP is False:

        os.remove(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + params.internal.HSFLABEL_02 + "_assigned-bnd.tmp.label"))
        os.remove(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + params.internal.HSFLABEL_02 + "_assigned-bnd.list"))
        os.remove(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + params.internal.HSFLABEL_02 + "_assigned-bnd.vox2ras"))
        os.remove(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + params.internal.HSFLABEL_07 + ".asc"))

    # return

    return(params)
