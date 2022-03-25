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

        # upsample

        if len(params.internal.UPSAMPLE) == 3:

            print("Upsampling with custom parameters ...")
            print()

            cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_convert") + " " \
                + " -ds " + str(params.internal.UPSAMPLE[0]) + " " + str(params.internal.UPSAMPLE[1]) + "  " + str(params.internal.UPSAMPLE[2]) + " " \
                + " -rt nearest " \
                + params.FILENAME + " " \
                + os.path.splitext(params.FILENAME)[0] + ".ups" + os.path.splitext(params.FILENAME)[1]

            print(cmd)

            subprocess.run(cmd.split())

        elif len(params.internal.UPSAMPLE) == 0:

            #
            import nibabel as nb
            import numpy as np
            import plotly.express as px
            from scipy import ndimage
            from nilearn import image as nli

            #
            print("Upsampling to min voxel size ...")
            print()

            # get image and info
            img = nb.load(params.FILENAME)
            dat = img.get_fdata()
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
            nb.save(img_int, os.path.splitext(params.FILENAME)[0] + ".ups" + os.path.splitext(params.FILENAME)[1])

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
        if "parahippocampal" in params.LUTDICT.keys():
            labelPHC = params.LUTDICT["parahippocampal"]
        else:
            labelPHC = None
        if "entorhinal" in params.LUTDICT.keys():
            labelENT = params.LUTDICT["entorhinal"]
        else:
            labelENT = None
        if "ba35" in params.LUTDICT.keys():
            labelBA35 = params.LUTDICT["ba35"]
        else:
            labelBA35 = None
        if "ba36" in params.LUTDICT.keys():
            labelBA36= params.LUTDICT["ba36"]
        else:
            labelBA36 = None

        # head
        if params.internal.BNDHEAD is True:
            #
            print("Using auto-mask for head")
            # remove any potentially existing head labels
            dat[dat==labelHead] = 0
            #
            if params.LUT == "ashs-ctx":
                if labelCA2 is not None and labelCA3 is not None and labelENT is not None and labelBA35 is not None:
                    # find instance of CA2/CA3
                    idxHead = np.argwhere(np.logical_or(dat==labelBA35, np.logical_or(dat==labelENT, np.logical_or(dat==labelCA2, dat==labelCA3))))
                else:
                    logging.info("Insufficient label information, exiting.")
                    sys.exit(1)
            elif params.LUT == "ashs-ctx-nohc":
                if labelENT is not None and labelBA35 is not None:
                    # find instance of ENT or BA35
                    idxHead = np.argwhere(np.logical_or(dat==labelENT, dat==labelBA35))
                else:
                    logging.info("Insufficient label information, exiting.")
                    sys.exit(1)
            elif params.LUT == "ashs-ent":
                if labelENT is not None:
                    # find instance of ENT
                    idxHead = np.argwhere(dat==labelENT)
                else:
                    logging.info("Insufficient label information, exiting.")
                    sys.exit(1)
            else:
                if labelCA2 is not None and labelCA3 is not None:
                    # find instance of CA2/CA3
                    idxHead = np.argwhere(np.logical_or(dat==labelCA2, dat==labelCA3))
                else:
                    logging.info("Insufficient label information, exiting.")
                    sys.exit(1)
            # get min/max row/col/slice
            if imgDimsAPDir == -1:
                cutFrom = np.min(idxHead[:, imgDimsAP]) + params.internal.BNDHEADMARGIN
            elif imgDimsAPDir == 1:
                cutFrom = np.max(idxHead[:, imgDimsAP]) - params.internal.BNDHEADMARGIN
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
            # TODO: look for consistent mask definitions, i.e. coords vs indices
            #
            print("Using auto-mask for tail")
            # remove any potentially existing tail labels
            dat[dat==labelTail] = 0
            #
            if params.LUT == "ashs-ctx":
                if labelPHC is not None and labelSbc is not None and labelCA1 is not None and labelCA2 is not None and labelCA3 is not None:
                    # find instance of Sbc/CA1/CA2/CA3 and PHC
                    idxTail = np.intersect1d(np.argwhere(dat==labelPHC)[:,imgDimsAP],
                        np.intersect1d(np.argwhere(dat==labelSbc)[:,imgDimsAP],
                        np.intersect1d(np.argwhere(dat==labelCA1)[:,imgDimsAP],
                        np.intersect1d(np.argwhere(dat==labelCA2)[:,imgDimsAP],
                        np.argwhere(dat==labelCA3)[:,imgDimsAP]))))
                else:
                    logging.info("Insufficient label information, exiting.")
                    sys.exit(1)
            elif params.LUT == "ashs-ctx-nohc":
                if labelPHC is not None:
                    # find instance of PHC
                    idxTail = np.argwhere(dat==labelPHC)[:,imgDimsAP]
                else:
                    logging.info("Insufficient label information, exiting.")
                    sys.exit(1)
            elif params.LUT == "ashs-ent":
                if labelENT is not None:
                    # find instance of ENT
                    idxTail = np.argwhere(dat==labelENT)[:,imgDimsAP]
                else:
                    logging.info("Insufficient label information, exiting.")
                    sys.exit(1)
            elif params.LUT == "ashs-noca3":
                if labelSbc is not None and labelCA1 is not None and labelCA2 is not None:
                    # find instance of Sbc/CA1/CA2/CA3
                    idxTail = np.intersect1d(np.argwhere(dat==labelSbc)[:,imgDimsAP],
                        np.intersect1d(np.argwhere(dat==labelCA1)[:,imgDimsAP],
                        np.argwhere(dat==labelCA2)[:,imgDimsAP]))
                else:
                    logging.info("Insufficient label information, exiting.")
                    sys.exit(1)
            elif params.LUT == "ashs-ca2ca3":
                if labelSbc is not None and labelCA1 is not None and labelCA2 is not None and labelCA3 is not None:
                    # find instance of Sbc/CA1/(CA2orCA3)
                    idxTail = np.intersect1d(np.argwhere(dat==labelSbc)[:,imgDimsAP],
                        np.intersect1d(np.argwhere(dat==labelCA1)[:,imgDimsAP],
                        np.argwhere(np.logical_or(dat==labelCA2, dat==labelCA3))[:,imgDimsAP]))
                else:
                    logging.info("Insufficient label information, exiting.")
                    sys.exit(1)
            else:
                if labelSbc is not None and labelCA1 is not None and labelCA2 is not None and labelCA3 is not None:
                    # find instance of Sbc/CA1/CA2/CA3
                    idxTail = np.intersect1d(np.argwhere(dat==labelSbc)[:,imgDimsAP],
                        np.intersect1d(np.argwhere(dat==labelCA1)[:,imgDimsAP],
                        np.intersect1d(np.argwhere(dat==labelCA2)[:,imgDimsAP],
                        np.argwhere(dat==labelCA3)[:,imgDimsAP])))
                else:
                    logging.info("Insufficient label information, exiting.")
                    sys.exit(1)
            # get min/max row/col/slice
            if imgDimsAPDir == -1:
                cutFrom = np.max(idxTail) - params.internal.BNDTAILMARGIN
            elif imgDimsAPDir == 1:
                cutFrom = np.min(idxTail) + params.internal.BNDTAILMARGIN
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
