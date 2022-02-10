"""
This module provides a function to convert images to mgz format

"""

# ------------------------------------------------------------------------------
# main functions

def convertFormat(params):
    """

    """

    # imports

    import os
    import sys
    import logging
    import subprocess

    # message

    print()
    print("-------------------------------------------------------------------------")
    print("Convert to mgz and copy to output directory")
    print()

    # convert and copy

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_convert") + " " \
        + params.FILENAME + " " \
        + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_00 + ".orig.mgz")

    print(cmd)

    subprocess.run(cmd.split())

    # update params

    params.FILENAME = os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_00 + ".orig.mgz")

    # return

    return(params)


def cropImage(params):
    """

    """

    # imports

    import os
    import sys
    import logging
    import shutil
    import subprocess

    if params.internal.CROP is True:

        # message

        print()
        print("-------------------------------------------------------------------------")
        print("Cropping")
        print()

        # crop

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_mask") + " " \
            + "-bb 5 " \
            + params.FILENAME + " " \
            + params.FILENAME + " " \
            + os.path.splitext(params.FILENAME)[0] + ".crop" + os.path.splitext(params.FILENAME)[1]

        print(cmd)

        subprocess.run(cmd.split())

        # update params

        params.FILENAME = os.path.splitext(params.FILENAME)[0] + ".crop" + os.path.splitext(params.FILENAME)[1]

    # return

    return(params)


def upsampleImage(params):
    """

    """

    # imports

    import os
    import sys
    import logging
    import shutil
    import subprocess

    #

    if params.internal.UPSAMPLE is not None:

        # message

        print()
        print("-------------------------------------------------------------------------")
        print("Upsampling")
        print()

        # crop

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_convert") + " " \
            + " -ds " + str(params.internal.UPSAMPLE[0]) + " " + str(params.internal.UPSAMPLE[1]) + "  " + str(params.internal.UPSAMPLE[2]) + " " \
            + " -rt nearest " \
            + params.FILENAME + " " \
            + os.path.splitext(params.FILENAME)[0] + ".ups" + os.path.splitext(params.FILENAME)[1]

        print(cmd)

        subprocess.run(cmd.split())

        # update params

        params.FILENAME = os.path.splitext(params.FILENAME)[0] + ".ups" + os.path.splitext(params.FILENAME)[1]

    # return

    return(params)



def autoMask(params):
    """

    """

    # imports

    import os
    import numpy as np
    import nibabel as nb
    import subprocess

    # message

    print()
    print("-------------------------------------------------------------------------")
    print()
    print("Creating auto-mask")
    print()
    print("-------------------------------------------------------------------------")
    print()

    def autoMaskDo(imageFilename, params):

        # get image
        img = nb.freesurfer.load(imageFilename)

        # get image data
        dat = img.get_fdata()

        # get orientation
        imgDims = nb.aff2axcodes(img.affine)

        # get dimension (reference is RAS)
        imgDimsAP = np.where(np.isin(imgDims, ["A", "P"]))[0][0]

        # get direction
        if imgDims[imgDimsAP] == "A":
            imgDimsAPDir = 1
        elif imgDims[imgDimsAP] == "P":
            imgDimsAPDir = -1

        # get LUT indices for head, tail and CA2/CA3
        labelHead = params.LUTDICT["head"]
        labelTail = params.LUTDICT["tail"]
        labelSbc = params.LUTDICT["subiculum"]
        labelCA1 = params.LUTDICT["ca1"]
        labelCA2 = params.LUTDICT["ca2"]
        labelCA3 = params.LUTDICT["ca3"]
        labelPHC = params.LUTDICT["parahippocampal"]

        # head
        if params.internal.BNDHEAD is True:
            #
            print("Using auto-mask for head")
            # remove any potentially existing head labels
            dat[dat==labelHead] = 0
            # find instance of CA2/CA3
            idxCA23 = np.argwhere(np.logical_or(dat==labelCA2, dat==labelCA3))
            # get min/max row/col/slice
            if imgDimsAPDir == -1:
                cutFrom = np.min(idxCA23[:, imgDimsAP])
            elif imgDimsAPDir == 1:
                cutFrom = np.max(idxCA23[:, imgDimsAP])
            # set to zero
            if imgDimsAP == 0 and imgDimsAPDir == -1:
                dat[:cutFrom, :, :] = 0
            elif imgDimsAP == 0 and imgDimsAPDir == 1:
                dat[cutFrom+1:, :, :] = 0
            elif imgDimsAP == 1 and imgDimsAPDir == -1:
                dat[:, :cutFrom, :] = 0
            elif imgDimsAP == 1 and imgDimsAPDir == 1:
                dat[:, cutFrom+1:,:] = 0
            elif imgDimsAP == 2 and imgDimsAPDir == -1:
                dat[:, :, :cutFrom] = 0
            elif imgDimsAP == 2 and imgDimsAPDir == 1:
                dat[:, :, cutFrom+1:] = 0
            # set to label
            if imgDimsAP == 0:
                dat[cutFrom, np.nonzero(dat[cutFrom, :, :])[0], np.nonzero(dat[cutFrom, :, :])[1]] = labelHead
            elif imgDimsAP == 1:
                dat[np.nonzero(dat[:, cutFrom, :])[0], cutFrom, np.nonzero(dat[:, cutFrom, :])[1]] = labelHead
            elif imgDimsAP == 2:
                dat[np.nonzero(dat[:, :, cutFrom])[0], np.nonzero(dat[:, :, cutFrom])[1], cutFrom] = labelHead

        # tail
        if params.internal.BNDTAIL is True:
            #
            print("Using auto-mask for tail")
            # remove any potentially existing tail labels
            dat[dat==labelTail] = 0
            # find instance of Tail
            if np.any(dat==labelPHC):
                idxTail = np.intersect1d(np.argwhere(dat==labelPHC)[:,imgDimsAP],
                    np.intersect1d(np.argwhere(dat==labelSbc)[:,imgDimsAP],
                    np.intersect1d(np.argwhere(dat==labelCA1)[:,imgDimsAP],
                    np.intersect1d(np.argwhere(dat==labelCA2)[:,imgDimsAP],
                    np.argwhere(dat==labelCA3)[:,imgDimsAP]))))
            else:
                idxTail = np.intersect1d(np.argwhere(dat==labelSbc)[:,imgDimsAP],
                    np.intersect1d(np.argwhere(dat==labelCA1)[:,imgDimsAP],
                    np.intersect1d(np.argwhere(dat==labelCA2)[:,imgDimsAP],
                    np.argwhere(dat==labelCA3)[:,imgDimsAP])))
            # get min/max row/col/slice
            if imgDimsAPDir == -1:
                cutFrom = np.max(idxTail)
            elif imgDimsAPDir == 1:
                cutFrom = np.min(idxTail)
            # set to zero
            if imgDimsAP == 0 and imgDimsAPDir == -1:
                dat[cutFrom+1:, :, :] = 0
            elif imgDimsAP == 0 and imgDimsAPDir == 1:
                dat[:cutFrom, :, :] = 0
            elif imgDimsAP == 1 and imgDimsAPDir == -1:
                dat[:, cutFrom+1:, :] = 0
            elif imgDimsAP == 1 and imgDimsAPDir == 1:
                dat[:, :cutFrom:,:] = 0
            elif imgDimsAP == 2 and imgDimsAPDir == -1:
                dat[:, :, cutFrom+1:] = 0
            elif imgDimsAP == 2 and imgDimsAPDir == 1:
                dat[:, :, :cutFrom] = 0
            # set to label
            if imgDimsAP == 0:
                dat[cutFrom, np.nonzero(dat[cutFrom, :, :])[0], np.nonzero(dat[cutFrom, :, :])[1]] = labelTail
            elif imgDimsAP == 1:
                dat[np.nonzero(dat[:, cutFrom, :])[0], cutFrom, np.nonzero(dat[:, cutFrom, :])[1]] = labelTail
            elif imgDimsAP == 2:
                dat[np.nonzero(dat[:, :, cutFrom])[0], np.nonzero(dat[:, :, cutFrom])[1], cutFrom] = labelTail

        # create image
        out = nb.MGHImage(dat, img.affine, header=img.header)

        # return image
        return out

    if params.internal.BNDHEAD is True or params.internal.BNDTAIL is True:

        #

        out = autoMaskDo(params.FILENAME, params)

        nb.save(out, os.path.splitext(params.FILENAME)[0] + ".am" + os.path.splitext(params.FILENAME)[1])

        # update params
        params.FILENAME = os.path.splitext(params.FILENAME)[0] + ".am" + os.path.splitext(params.FILENAME)[1]

    else:

        #
        print("Not using auto-mask")

    # return

    return(params)
