"""
This module provides a function to create a label files for a tetrahedral mesh

"""

import os
import subprocess

import nibabel as nb
import numpy as np
from lapy import TetMesh, TriaMesh, io

# -----------------------------------------------------------------------------
# AUXILIARY FUNCTIONS

# createBoundaryMask()


def _createBoundaryMask(hbtFile, mskFile, outFile, label, bndlabel):
    """
    createBoundaryMask(hbtFile, mskFile)

    """

    # load segmentations
    hbt = nb.load(hbtFile)
    msk = nb.load(mskFile)

    # define new data structures
    bnd = msk
    bnd_data = bnd.get_fdata(caching="fill")
    hbt_data_flat = hbt.get_fdata().flatten()

    # find nonzero elements of msk
    mskx, msky, mskz = np.where(msk.get_fdata())

    # cycle through each nonzero element of msk
    for i in range(0, len(mskx)):
        # create indices xyz of local neighborhood (one more +1 due to python indexing)
        mx, my, mz = np.meshgrid(
            range(mskx[i] - 1, mskx[i] + 1 + 1),
            range(msky[i] - 1, msky[i] + 1 + 1),
            range(mskz[i] - 1, mskz[i] + 1 + 1),
        )

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
    """ """

    # message

    print()
    print("--------------------------------------------------------------------------------")
    print("Creating label files for tetrahedral meshes")
    print()

    # determine files

    MSKFILE = os.path.join(params.OUTDIR, params.HEMI + ".labels.mgz")

    if params.internal.AUTOMASK_HEAD is True or params.internal.AUTOMASK_TAIL is True:
        HBTFILE = os.path.join(params.OUTDIR, "labels", params.HEMI + ".automask_image.mgz")
    else:
        HBTFILE = os.path.join(params.OUTDIR, params.HEMI + ".image.mgz")

    # create new masks that includes boundary labels to head and tail, and to CA4

    _createBoundaryMask(
        hbtFile=HBTFILE,
        mskFile=MSKFILE,
        outFile=os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + "labels-bnd.mgz"),
        label=[params.LUTDICT["tail"], params.LUTDICT["head"]],
        bndlabel=[params.LUTDICT["bndtail"], params.LUTDICT["bndhead"]],
    )
    _createBoundaryMask(
        hbtFile=HBTFILE,
        mskFile=MSKFILE,
        outFile=os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + "labels-ca4.mgz"),
        label=[params.LUTDICT["ca4"]],
        bndlabel=[params.LUTDICT["bndca4"]],
    )

    # merge bnd and ca4 masks (add ca4 to bnd)

    cmd = (
        os.path.join(os.environ.get("FREESURFER_HOME"), "bin", "mris_calc")
        + " "
        + "--output "
        + os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + "labels-bnd.mgz")
        + " "
        + os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + "labels-bnd.mgz")
        + " "
        + "lowerlimit "
        + os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + "." + "labels-ca4.mgz")
    )

    print(cmd)

    subprocess.run(cmd.split(), capture_output=True)

    # load bnd image

    bndImage = nb.freesurfer.load(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".labels-bnd.mgz"))

    # load tet mesh

    tetMesh = TetMesh.read_vtk(os.path.join(params.OUTDIR, params.HEMI + ".tetra.vtk"))

    # convert vertex indices to voxel indices

    vtxVxIdx = np.round(
        np.matmul(
            np.linalg.inv(bndImage.header.get_vox2ras_tkr()),
            np.concatenate((tetMesh.v, np.ones((tetMesh.v.shape[0], 1))), axis=1).T,
        ).T
    ).astype(int)[:, 0:3]

    # lookup labels of converted vertex indices in bnd image

    vtxLabels = bndImage.get_fdata()[vtxVxIdx[:, 0], vtxVxIdx[:, 1], vtxVxIdx[:, 2]]

    # write out psol file

    with open(os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".tetra.psol"), "w") as f:
        print("Solution:", file=f)
        [print(str(int(x)), file=f) for x in vtxLabels]

    # write out visualization files

    tetMeshBnd = tetMesh.boundary_tria()
    vRmBndKeep, vBndRmDel = tetMeshBnd.rm_free_vertices_()
    tetMeshBnd.orient_()
    TriaMesh.write_vtk(tetMeshBnd, os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".rm.bnd.tetra.vtk"))
    io.write_vfunc(
        os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".rm.bnd.tetra.psol"),
        np.array([float(x) for x in vtxLabels])[vRmBndKeep],
    )
    nb.freesurfer.save(
        nb.freesurfer.MGHImage(
            dataobj=np.array([float(x) for x in vtxLabels])[vRmBndKeep].astype("float32"), affine=None
        ),
        filename=os.path.join(params.OUTDIR, "tetra-labels", params.HEMI + ".rm.bnd.label.mgh"),
    )

    # return

    return params
