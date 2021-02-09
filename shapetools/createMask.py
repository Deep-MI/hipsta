"""
This module provides a function to fill holes and create binary masks

"""

# ------------------------------------------------------------------------------
# main function

def fillHoles(params):

    """

    """

    # imports

    import os
    import subprocess

    # message

    print()
    print("-------------------------------------------------------------------------")
    print()
    print("Fill holes")
    print()
    print("-------------------------------------------------------------------------")
    print()

    # dilate/erode by 1 voxel to fill holes; we use the FSL tool because
    # mri_binarize cannot be applied without binarization, which we don't want yet.
    # outputs ${HEMI}.de.${HSFLABEL_01}.mgz

    # convert forward

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_convert") + " " \
        + "--in_type mgz " \
        + "--out_type nii " \
        + os.path.join(params.OUTDIR, "merge-ml", params.HEMI + "." + params.internal.HSFLABEL_01 + "_merged.mgz") + " " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + "." + params.internal.HSFLABEL_01 + "_merged.nii")

    print(cmd)

    subprocess.run(cmd.split())

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_convert") + " " \
        + "--in_type mgz " \
        + "--out_type nii " \
        + os.path.join(params.OUTDIR, "merge-ml", params.HEMI + "." + params.internal.HSFLABEL_01 + "_assigned.mgz") + " " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + "." + params.internal.HSFLABEL_01 + "_assigned.nii")

    print(cmd)

    subprocess.run(cmd.split())

    # do the math

    os.environ["FSLOUTPUTTYPE"] = "NIFTI"

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "fslmaths.fsl") + " " \
        + "-dt input " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + "." + params.internal.HSFLABEL_01 + "_merged.nii") + " " \
        + "-kernel box 1 -dilD -ero " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_merged.nii")

    print(cmd)

    subprocess.run(cmd.split())

    # cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "fslmaths.fsl") + " " \
    #     + "-dt input " \
    #     + os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_merged.nii") + " " \
    #     + "-kernel box 1 -ero -dilD " \
    #     + os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_merged.nii")
    #
    # print(cmd)
    #
    # subprocess.run(cmd.split())

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "fslmaths.fsl") + " " \
        + "-dt input " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + "." + params.internal.HSFLABEL_01 + "_assigned.nii") + " " \
        + "-kernel box 1 -dilD -ero " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_assigned.nii")

    print(cmd)

    subprocess.run(cmd.split())

    # cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "fslmaths.fsl") + " " \
    #     + "-dt input " \
    #     + os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_assigned.nii") + " " \
    #     + "-kernel box 1 -ero -dilD " \
    #     + os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_assigned.nii")
    #
    # print(cmd)
    #
    # subprocess.run(cmd.split())

    # convert back

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_convert") + " " \
        + "--in_type nii " \
        + "--out_type mgz " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_merged.nii") + " " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_merged.mgz")

    print(cmd)

    subprocess.run(cmd.split())

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_convert") + " " \
        + "--in_type nii " \
        + "--out_type mgz " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_assigned.nii") + " " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_assigned.mgz")

    print(cmd)

    subprocess.run(cmd.split())

    # clean up

    if params.skipCLEANUP is False:

        os.remove(os.path.join(params.OUTDIR, "mask", params.HEMI + "." + params.internal.HSFLABEL_01 + "_merged.nii"))
        os.remove(os.path.join(params.OUTDIR, "mask", params.HEMI + "." + params.internal.HSFLABEL_01 + "_assigned.nii"))

        os.remove(os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_merged.nii"))
        os.remove(os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_assigned.nii"))

    # update params

    params.internal.HSFLABEL_02 = "de." + params.internal.HSFLABEL_01

    # return

    return(params)


def createMask(params):
    """

    """

    # imports

    import os
    import subprocess

    # message

    print()
    print("-------------------------------------------------------------------------")
    print()
    print("Create binary mask")
    print()
    print("-------------------------------------------------------------------------")
    print()

    # mri_binarize

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_binarize") + " " \
        + "--i " + os.path.join(params.OUTDIR, "mask", params.HEMI + "." + params.internal.HSFLABEL_02 + "_merged.mgz") + " " \
        + "--min 1 " \
        + "--o " + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_02 + ".mgz")

    print(cmd)

    subprocess.run(cmd.split())

    # return

    return(params)
