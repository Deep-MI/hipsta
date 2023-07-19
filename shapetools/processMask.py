"""
This module provides functions for image processing operations (such as binarization, filtering, closing operations)

"""

def gaussFilter(params):
    """

    """

    # imports

    import os
    import nibabel as nb

    from scipy import ndimage

    #

    if params.internal.GAUSSFILTER is True:

        # message

        print()
        print("-------------------------------------------------------------------------")
        print()
        print("Gaussian filtering")
        print()
        print("-------------------------------------------------------------------------")
        print()

        # gaussian filtering

        img = nb.load(params.FILENAME)
        dat = img.get_fdata()

        dat_filtered = ndimage.gaussian_filter((dat>0)*100, sigma=params.internal.GAUSSFILTER_SIZE[0], mode="constant", cval=0.0)>params.internal.GAUSSFILTER_SIZE[1] # TODO: need to set threshold after testing

        nb.freesurfer.save(nb.freesurfer.MGHImage(dataobj=dat_filtered.astype("float32"), affine=img.get_affine()), filename=os.path.join(params.OUTDIR, "mask", params.HEMI + ".gaussian_filter.mgz"))

        # update params

        params.FILENAME = os.path.join(params.OUTDIR, "mask", params.HEMI + ".gaussian_filter.mgz")

    # return

    return(params)


def longFilter(params):
    """

    """

    # imports

    import os
    from scipy import ndimage
    import nibabel as nb
    import numpy as np

    #

    if params.internal.LONGFILTER is True:

        # message

        print()
        print("-------------------------------------------------------------------------")
        print()
        print("Filtering along longitudinal axis")
        print()
        print("-------------------------------------------------------------------------")
        print()

        # longitudinal filtering

        img = nb.load(params.FILENAME)
        dat = img.get_fdata()

        k = np.zeros(params.internal.LONGFILTER_SIZE)
        k[2,2,:] = 1

        dat_filtered = ndimage.convolve(dat, k, mode="constant", cval=0.0)

        dat_filtered = dat_filtered>0

        nb.freesurfer.save(nb.freesurfer.MGHImage(dataobj=dat_filtered.astype("float32"), affine=img.get_affine()), filename=os.path.join(params.OUTDIR, "mask", params.HEMI + ".longitudinal_filter.mgz"))

        # update params

        params.FILENAME = os.path.join(params.OUTDIR, "mask", params.HEMI + ".longitudinal_filter.mgz")

    # return

    return(params)


def closeMask(params):
    """

    """

    # imports

    import os
    from scipy import ndimage
    import nibabel as nb

    #

    if params.internal.CLOSEMASK is True:

        # message

        print()
        print("-------------------------------------------------------------------------")
        print()
        print("Applying closing operation to mask")
        print()
        print("-------------------------------------------------------------------------")
        print()

        # get data

        img = nb.load(params.FILENAME)
        dat = img.get_fdata()

        dat_filtered = ndimage.binary_closing(dat, structure=ndimage.iterate_structure(ndimage.generate_binary_structure(rank=3, connectivity=1), iterations=1).astype(int), iterations=1)

        dat_filtered = dat_filtered>0

        nb.freesurfer.save(nb.freesurfer.MGHImage(dataobj=dat_filtered.astype("float32"), affine=img.get_affine()), filename=os.path.join(params.OUTDIR, "mask", params.HEMI + ".close_mask.mgz"))

        # update params

        params.FILENAME = os.path.join(params.OUTDIR, "mask", params.HEMI + ".close_mask.mgz")

    # return

    return(params)


def binarizeMask(params):

    # imports

    import os

    import nibabel as nb

    # binarize

    img = nb.load(params.FILENAME)
    dat = img.get_fdata()

    dat_filtered = dat!=0

    nb.freesurfer.save(nb.freesurfer.MGHImage(dataobj=dat_filtered.astype("float32"), affine=img.get_affine()), filename=os.path.join(params.OUTDIR, "mask", params.HEMI + ".initial_mask.mgz"))

    # update params

    params.FILENAME = os.path.join(params.OUTDIR, "mask", params.HEMI + ".initial_mask.mgz")

    # return

    return(params)


def copy_mask_to_main(params):
    """

    """

    #

    import os
    import shutil

    # copy to main directory

    shutil.copyfile(
        params.FILENAME,
        os.path.join(params.OUTDIR, params.HEMI + ".mask.mgz"))

    # update params

    params.FILENAME = os.path.join(params.OUTDIR, params.HEMI + ".mask.mgz")

    # return

    return(params)


def filterMaskRemove(params):
    """

    """

    # imports

    import os
    import subprocess

    #

    if params.internal.FILTERMASK is not None:

        # message

        print()
        print("-------------------------------------------------------------------------")
        print()
        print("Filter mask")
        print()
        print("-------------------------------------------------------------------------")
        print()

        #

        os.environ["FSLOUTPUTTYPE"] = "NIFTI"

        # convert

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_convert") + " " \
            + "--in_type mgz " \
            + "--out_type nii " \
            + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_02 + ".mgz") + " " \
            + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_02 + ".nii")

        print(cmd)

        subprocess.run(cmd.split())

        # filter

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "fslmaths.fsl") + " " \
            + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_02 + ".nii") + " " \
            + "-mul 100 -kernel gauss " + str(params.internal.FILTERMASK[0]) + " -fmean " \
            + os.path.join(params.OUTDIR, params.HEMI + ".filt." + params.internal.HSFLABEL_02 + ".nii")

        print(cmd)

        subprocess.run(cmd.split())

        # convert back

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_convert") + " " \
            + "--in_type nii " \
            + "--out_type mgz " \
            + os.path.join(params.OUTDIR, params.HEMI + ".filt." + params.internal.HSFLABEL_02 + ".nii") + " " \
            + os.path.join(params.OUTDIR, params.HEMI + ".filt." + params.internal.HSFLABEL_02 + ".mgz")

        print(cmd)

        subprocess.run(cmd.split())

        # binarize

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_binarize") + " " \
            + "--i " + os.path.join(params.OUTDIR, params.HEMI + ".filt." + params.internal.HSFLABEL_02 + ".mgz") + " " \
            + "--min " + str(params.internal.FILTERMASK[1]) + " --binval 1 " \
            + "--o " + os.path.join(params.OUTDIR, params.HEMI + ".filt." + params.internal.HSFLABEL_02 + ".mgz")

        print(cmd)

        subprocess.run(cmd.split())

        # clean up

        os.remove(os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_02 + ".nii"))

        os.remove(os.path.join(params.OUTDIR, params.HEMI + ".filt." + params.internal.HSFLABEL_02 + ".nii"))

        # update params

        params.internal.HSFLABEL_02 = "filt." + params.internal.HSFLABEL_02

    # return

    return(params)


def fillHolesRemove(params):

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

    # dilate/erode by n voxels to fill holes; we use the FSL tool because
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

    #

    os.environ["FSLOUTPUTTYPE"] = "NIFTI"

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "fslmaths.fsl") + " " \
        + "-dt int " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + "." + params.internal.HSFLABEL_01 + "_merged.nii") + " " \
        + "-kernel boxv " + str(params.internal.DIL) + " -dilD " \
        + "-kernel boxv " + str(params.internal.ERO) + " -ero " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_merged.nii")

    print(cmd)

    subprocess.run(cmd.split())

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "fslmaths.fsl") + " " \
        + "-dt int " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + "." + params.internal.HSFLABEL_01 + "_assigned.nii") + " " \
        + "-kernel boxv " + str(params.internal.DIL) + " -dilD " \
        + "-kernel boxv " + str(params.internal.ERO) + " -ero " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_assigned.nii")

    print(cmd)

    subprocess.run(cmd.split())

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


def createMaskRemove(params):
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