"""
This module provides functions for basic image processing operations (conversion, cropping, upsampling)

"""

import os
import shutil
import subprocess

import nibabel as nb

from nilearn import image as nli

# ==============================================================================
# FUNCTIONS


def convertFormat(params):
    """

    """

    # message

    print()
    print("--------------------------------------------------------------------------------")
    print("Convert to mgz and copy to output directory")
    print()

    # convert and copy

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_convert") + " " \
        + params.FILENAME + " " \
        + os.path.join(params.OUTDIR, "image", params.HEMI + ".orig.mgz")

    print(cmd)

    subprocess.run(cmd.split())

    # update params

    params.FILENAME = os.path.join(params.OUTDIR, "image", params.HEMI + ".orig.mgz")

    # return

    return(params)


def cropImage(params):
    """

    """

    if params.internal.CROP is True:

        # message

        print()
        print("--------------------------------------------------------------------------------")
        print("Cropping")
        print()

        # crop

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_mask") + " " \
            + "-bb 5 " \
            + params.FILENAME + " " \
            + params.FILENAME + " " \
            + os.path.join(params.OUTDIR, "image", params.HEMI + ".cropped.mgz")

        print(cmd)

        subprocess.run(cmd.split())

        # update params

        params.FILENAME = os.path.join(params.OUTDIR, "image", params.HEMI + ".cropped.mgz")

    # return

    return(params)


def upsampleImage(params):
    """

    """

    if params.internal.UPSAMPLE is not None:

        # message

        print()
        print("--------------------------------------------------------------------------------")
        print("Upsampling")
        print()

        # upsample

        if len(params.internal.UPSAMPLE) == 3:

            logging.info("Upsampling with custom parameters ...")

            cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_convert") + " " \
                + " -ds " + str(params.internal.UPSAMPLE[0]) + " " + str(params.internal.UPSAMPLE[1]) + "  " + str(params.internal.UPSAMPLE[2]) + " " \
                + " -rt nearest " \
                + params.FILENAME + " " \
                + os.path.join(params.OUTDIR, "image", params.HEMI + ".upsampled.mgz")

            print(cmd)

            subprocess.run(cmd.split())

        elif len(params.internal.UPSAMPLE) == 0:

            #
            logging.info("Upsampling to min voxel size ...")

            # get image and info
            img = nb.load(params.FILENAME)
            vxsz = img.header["delta"]
            affn = img.affine

            # compute new affine
            target_affn = affn.copy()
            target_affn[0:4,0] *= min(vxsz)/vxsz[0]
            target_affn[0:4,1] *= min(vxsz)/vxsz[1]
            target_affn[0:4,2] *= min(vxsz)/vxsz[2]

            # resample
            img_int = nli.resample_img(img, target_affine=target_affn, interpolation="nearest")

            # save
            nb.save(img_int, os.path.join(params.OUTDIR, "image", params.HEMI + ".upsampled.mgz"))

        # update params
        params.FILENAME = os.path.join(params.OUTDIR, "image", params.HEMI + ".upsampled.mgz")

    # return

    return(params)


def copy_image_to_main(params):
    """
    
    """

    # copy to main directory

    shutil.copyfile(
        params.FILENAME,
        os.path.join(params.OUTDIR, params.HEMI + ".image.mgz"))

    # update params

    params.FILENAME = os.path.join(params.OUTDIR, params.HEMI + ".image.mgz")

    # return

    return(params)