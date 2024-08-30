"""
This module provides functions for processing label images

"""

import logging
import os
import shutil
import subprocess

import nibabel as nb
import numpy as np

# ==============================================================================
# LOGGING

LOGGER = logging.getLogger(__name__)

# ==============================================================================
# FUNCTIONS


def autoMask(params):
    """ """

    if params.internal.AUTOMASK_HEAD is True or params.internal.AUTOMASK_TAIL is True:
        # message

        print()
        print("--------------------------------------------------------------------------------")
        print("Creating auto-mask")
        print()

        # get image
        img = nb.freesurfer.load(params.FILENAME)

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
        # if "parahippocampal" in params.LUTDICT.keys():
        #     labelPHC = params.LUTDICT["parahippocampal"]
        # else:
        #     labelPHC = None
        # if "entorhinal" in params.LUTDICT.keys():
        #     labelENT = params.LUTDICT["entorhinal"]
        # else:
        #    labelENT = None
        # if "ba35" in params.LUTDICT.keys():
        #     labelBA35 = params.LUTDICT["ba35"]
        # else:
        #     labelBA35 = None
        # if "ba36" in params.LUTDICT.keys():
        #     labelBA36 = params.LUTDICT["ba36"]
        # else:
        #    labelBA36 = None

        # head
        if params.internal.AUTOMASK_HEAD is True:
            #
            LOGGER.info("Using auto-mask for head")
            # remove any potentially existing head labels
            dat[dat == labelHead] = 0
            #
            if labelCA2 is not None and labelCA3 is not None:
                # find instance of CA2/CA3
                idxHead = np.argwhere(np.logical_or(dat == labelCA2, dat == labelCA3))
            else:
                raise RuntimeError("Insufficient label information, exiting.")
            # get min/max row/col/slice
            if imgDimsAPDir == -1:
                cutFrom = np.min(idxHead[:, imgDimsAP]) + params.internal.AUTOMASK_HEAD_MARGIN
            elif imgDimsAPDir == 1:
                cutFrom = np.max(idxHead[:, imgDimsAP]) - params.internal.AUTOMASK_HEAD_MARGIN
            # set to zero
            if imgDimsAP == 0 and imgDimsAPDir == -1:
                dat[:cutFrom, :, :] = 0
            elif imgDimsAP == 0 and imgDimsAPDir == 1:
                dat[cutFrom + 1 :, :, :] = 0
            elif imgDimsAP == 1 and imgDimsAPDir == -1:
                dat[:, :cutFrom, :] = 0
            elif imgDimsAP == 1 and imgDimsAPDir == 1:
                dat[:, cutFrom + 1 :, :] = 0
            elif imgDimsAP == 2 and imgDimsAPDir == -1:
                dat[:, :, :cutFrom] = 0
            elif imgDimsAP == 2 and imgDimsAPDir == 1:
                dat[:, :, cutFrom + 1 :] = 0
            # set to label
            if imgDimsAP == 0:
                dat[cutFrom, np.nonzero(dat[cutFrom, :, :])[0], np.nonzero(dat[cutFrom, :, :])[1]] = labelHead
            elif imgDimsAP == 1:
                dat[np.nonzero(dat[:, cutFrom, :])[0], cutFrom, np.nonzero(dat[:, cutFrom, :])[1]] = labelHead
            elif imgDimsAP == 2:
                dat[np.nonzero(dat[:, :, cutFrom])[0], np.nonzero(dat[:, :, cutFrom])[1], cutFrom] = labelHead

        # tail
        if params.internal.AUTOMASK_TAIL is True:
            #
            LOGGER.info("Using auto-mask for tail")
            # remove any potentially existing tail labels
            dat[dat == labelTail] = 0
            #
            if labelSbc is not None and labelCA1 is not None and labelCA2 is not None and labelCA3 is not None:
                ## find instance of Sbc/CA1/CA2/CA3
                # idxTail = np.intersect1d(
                #    np.argwhere(dat == labelSbc)[:, imgDimsAP],
                #    np.intersect1d(
                #        np.argwhere(dat == labelCA1)[:, imgDimsAP],
                #        np.intersect1d(
                #            np.argwhere(dat == labelCA2)[:, imgDimsAP], np.argwhere(dat == labelCA3)[:, imgDimsAP]
                #        ),
                #    ),
                # )
                # find instance of Sbc/CA1
                idxTail = np.intersect1d(
                    np.argwhere(dat == labelSbc)[:, imgDimsAP], np.argwhere(dat == labelCA1)[:, imgDimsAP]
                )
            else:
                raise RuntimeError("Insufficient label information, exiting.")
            # get min/max row/col/slice
            if imgDimsAPDir == -1:
                cutFrom = np.max(idxTail) - params.internal.AUTOMASK_TAIL_MARGIN
            elif imgDimsAPDir == 1:
                cutFrom = np.min(idxTail) + params.internal.AUTOMASK_TAIL_MARGIN
            # set to zero
            if imgDimsAP == 0 and imgDimsAPDir == -1:
                dat[cutFrom + 1 :, :, :] = 0
            elif imgDimsAP == 0 and imgDimsAPDir == 1:
                dat[:cutFrom, :, :] = 0
            elif imgDimsAP == 1 and imgDimsAPDir == -1:
                dat[:, cutFrom + 1 :, :] = 0
            elif imgDimsAP == 1 and imgDimsAPDir == 1:
                dat[:, :cutFrom:, :] = 0
            elif imgDimsAP == 2 and imgDimsAPDir == -1:
                dat[:, :, cutFrom + 1 :] = 0
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

        # save
        nb.save(out, os.path.join(params.OUTDIR, "labels", params.HEMI + ".automask_image.mgz"))

        # update params
        params.FILENAME = os.path.join(params.OUTDIR, "labels", params.HEMI + ".automask_image.mgz")

    # return

    return params


def createLabels(params):
    """ """

    # message

    print()
    print("--------------------------------------------------------------------------------")
    print("Create labels")
    print()

    # create a mask
    cmd = (
        os.path.join(os.environ.get("FREESURFER_HOME"), "bin", "mri_binarize")
        + " "
        + "--i "
        + params.FILENAME
        + " "
        + "--match "
        + " ".join([str(x) for x in params.HSFLIST])
        + " "
        + "--o "
        + os.path.join(params.OUTDIR, "labels", params.HEMI + ".initial_labels.mgz")
    )

    print(cmd)

    subprocess.run(cmd.split(), capture_output=True)

    # multiply original image with mask
    cmd = (
        os.path.join(os.environ.get("FREESURFER_HOME"), "bin", "fscalc")
        + " "
        + params.FILENAME
        + " mul "
        + os.path.join(params.OUTDIR, "labels", params.HEMI + ".initial_labels.mgz")
        + " --o "
        + os.path.join(params.OUTDIR, "labels", params.HEMI + ".initial_labels.mgz")
    )

    print(cmd)

    subprocess.run(cmd.split(), capture_output=True)

    # update params
    params.FILENAME = os.path.join(params.OUTDIR, "labels", params.HEMI + ".initial_labels.mgz")

    #
    return params


def mergeMolecularLayer(params):
    """ """

    if params.internal.MERGE_MOLECULAR_LAYER is True:
        # message
        print()
        print("--------------------------------------------------------------------------------")
        print("Attaching the molecular layer")
        print()

        # settings
        offset = 0

        # load image
        im = nb.load(params.FILENAME)

        # get voxel data
        vx = im.get_fdata()

        # get all voxels that are not zero; sort by dimensions 1,2,3
        vxNZ = list(np.nonzero(vx))
        vxNZ.append(vx[np.nonzero(vx)])
        vxNZ = np.transpose(vxNZ)
        vxNZ = vxNZ[np.argsort(vxNZ[:, 0]), :]
        vxNZ = vxNZ[np.argsort(vxNZ[:, 1], kind="mergesort"), :]
        vxNZ = vxNZ[np.argsort(vxNZ[:, 2], kind="mergesort"), :]

        # get all ML voxels
        vxML = list(np.where((vx == 245) | (vx == 246) | (vx == 214)))
        vxML.append(vx[np.where((vx == 245) | (vx == 246) | (vx == 214))])
        vxML = np.transpose(vxML)
        vxML = vxML[np.argsort(vxML[:, 0]), :]
        vxML = vxML[np.argsort(vxML[:, 1], kind="mergesort"), :]
        vxML = vxML[np.argsort(vxML[:, 2], kind="mergesort"), :]

        #
        for i in range(np.shape(vxML)[0]):
            n = 1
            tmpML = np.transpose((np.empty([0]), np.empty([0]), np.empty([0])))
            while np.shape(tmpML)[0] == 0:
                # break if iterations are exceeded or if coordinate is below zero or equal / above image dims
                if (
                    n > 10
                    or (int(vxML[i, 0] - n) < 0)
                    or (int(vxML[i, 0] + n) > (vx.shape[0] - 1))
                    or (int(vxML[i, 1] - n) < 0)
                    or (int(vxML[i, 1] + n) > (vx.shape[1] - 1))
                    or (int(vxML[i, 2] - n) < 0)
                    or (int(vxML[i, 2] + n) > (vx.shape[2] - 1))
                ):
                    vxML[i, 3] = 0
                    break
                tmp = vx[
                    np.ix_(
                        range(int(vxML[i, 0] - n), int((vxML[i, 0] + n) + 1)),
                        range(int(vxML[i, 1] - n), int((vxML[i, 1] + n) + 1)),
                        range(int(vxML[i, 2] - n), int((vxML[i, 2] + n) + 1)),
                    )
                ]
                tmpML = np.transpose(np.where((tmp != 0) & (tmp != 246) & (tmp != 245) & (tmp != 214)))
                if np.shape(tmpML)[0] != 0:
                    tmp0 = np.asarray(tmpML[np.where(np.sum(tmpML, axis=1) == min(np.sum(tmpML, axis=1)))[0], :])
                    vxML[i, 3] = offset + min(tmp[tmp0[:, 0], tmp0[:, 1], tmp0[:, 2]])
                    tmp0 = None
                else:
                    n = n + 1
                tmp = None
            tmpML = None

        # write back to vx and im
        for i in range(np.shape(vxML)[0]):
            vx[int(vxML[i, 0]), int(vxML[i, 1]), int(vxML[i, 2])] = vxML[i, 3]

        om = nb.MGHImage(dataobj=vx, affine=im.affine, header=im.header)
        nb.save(img=om, filename=os.path.join(params.OUTDIR, "labels", params.HEMI + ".ml_labels.mgz"))

        # update params
        params.FILENAME = os.path.join(params.OUTDIR, "labels", params.HEMI + ".ml_labels.mgz")

    # return
    return params


def copy_labels_to_main(params):
    """ """

    # copy to main directory

    shutil.copyfile(params.FILENAME, os.path.join(params.OUTDIR, params.HEMI + ".labels.mgz"))

    # update params

    params.FILENAME = os.path.join(params.OUTDIR, params.HEMI + ".labels.mgz")

    # return

    return params
