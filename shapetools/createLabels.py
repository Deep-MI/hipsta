"""
This module provides a function to process image labels

"""

# ------------------------------------------------------------------------------
# main function

def createLabels(params):
    """

    """

    # imports

    import os
    import sys
    import shutil
    import subprocess

    # message

    print()
    print("-------------------------------------------------------------------")
    print()
    print("Process HSF label image")
    print()
    print("-------------------------------------------------------------------")
    print()

    # process HSF label image

    # This step extracts HSFs from the label image (note that dilate and erode may
    # create masks that are not disjunct; this will be taken care of). It outputs
    # several files that contain binary masks for subfields, i.e.
    # (XXX_XX_<label>_<original-mri>.mgz)

    for i in params.HSFLIST:

        print()
        print("---------------------------------------------------------------")
        print("Working on label "+i)
        print()

        # determine which anatomical label corresponds to the given index

        # check for errors
        if not i in params.LUTDICT.values():

            print("ERROR: Label "+i+" not found within look-up table, exiting.")
            print()
            sys.exit(1)

        # create a mask after erosion and dilation
        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_binarize") + " " \
            + "--i " +  params.FILENAME + " " \
            + "--match " + i + " " \
            + "--dilate " + str(params.internal.DIL) + " " \
            + "--erode " + str(params.internal.ERO) + " " \
            + "--binval " + i + " " \
            + "--o " + os.path.join(params.OUTDIR, "labels", i + "." + os.path.basename(params.FILENAME))

        print(cmd)

        subprocess.run(cmd.split())

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_binarize") + " " \
            + "--i " +  params.FILENAME + " " \
            + "--match " + i + " " \
            + "--erode 1 " \
            + "--dilate 1 " \
            + "--binval " + i + " " \
            + "--o " + os.path.join(params.OUTDIR, "labels", i + "." + os.path.basename(params.FILENAME))

        print(cmd)

        subprocess.run(cmd.split())

    # merge binary subfield masks

    # This step merges the binary subfield masks (but does not binarize yet) and
    # outputs one volume that contains all labels (${HEMI}.${HSFLABEL_00}.mgz)

    shutil.copy(
        os.path.join(params.OUTDIR, "labels", params.HSFLIST[0] + "." + os.path.basename(params.FILENAME)),
        os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_00 + ".mgz"))

    for i in params.HSFLIST:

        print()
        print("---------------------------------------------------------------")
        print("Merging label "+i)
        print()

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "fscalc") + " " \
            + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_00 + ".mgz") + " " \
            + "and " \
            + os.path.join(params.OUTDIR, "labels", i + "." + os.path.basename(params.FILENAME)) + " " \
            + "mul " + i + " " \
            + "--o " + os.path.join(params.OUTDIR,  "labels", params.HEMI + "." + params.internal.HSFLABEL_00 + ".tmp.mgz")

        print(cmd)

        subprocess.run(cmd.split())

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "fscalc") + " " \
            + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_00 + ".mgz") + " " \
            + "sub " \
            + os.path.join(params.OUTDIR,  "labels", params.HEMI + "." + params.internal.HSFLABEL_00 + ".tmp.mgz") + " " \
            + "--o " + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_00 + ".mgz")

        print(cmd)

        subprocess.run(cmd.split())

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "fscalc") + " " \
            + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_00 + ".mgz") + " " \
            + "add " \
            + os.path.join(params.OUTDIR, "labels", i + "." + os.path.basename(params.FILENAME)) + " " \
            + "--o " + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_00 + ".mgz")

        print(cmd)

        subprocess.run(cmd.split())

    # clean up

    if params.skipCLEANUP is False:

        os.remove(os.path.join(params.OUTDIR,  "labels", params.HEMI + "." + params.internal.HSFLABEL_00 + ".tmp.mgz"))

    # return

    return(params)
