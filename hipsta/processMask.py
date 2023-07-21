"""
This module provides functions for image processing operations (such as binarization, filtering, closing operations)

"""

import os
import shutil

import nibabel as nb
import numpy as np

from scipy import ndimage


# ==============================================================================
# FUNCTIONS

def gaussFilter(params):
    """

    """

    #

    if params.internal.GAUSSFILTER is True:

        # message

        print()
        print("--------------------------------------------------------------------------------")
        print("Gaussian filtering")
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

    #

    if params.internal.LONGFILTER is True:

        # message

        print()
        print("--------------------------------------------------------------------------------")
        print("Filtering along longitudinal axis")
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

    #

    if params.internal.CLOSEMASK is True:

        # message

        print()
        print("--------------------------------------------------------------------------------")
        print("Applying closing operation to mask")
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

    # copy to main directory

    shutil.copyfile(
        params.FILENAME,
        os.path.join(params.OUTDIR, params.HEMI + ".mask.mgz"))

    # update params

    params.FILENAME = os.path.join(params.OUTDIR, params.HEMI + ".mask.mgz")

    # return

    return(params)
