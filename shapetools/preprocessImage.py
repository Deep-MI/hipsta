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
    
    
    
