"""
This module provides a function to convert images to mgz format

"""

# ------------------------------------------------------------------------------
# main function

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
