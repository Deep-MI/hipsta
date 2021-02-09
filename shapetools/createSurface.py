"""
This module provides a function to create surfaces, which includes the marching
cube algorithm, optional topology fixing, and remeshing and smoothing.

"""

# ------------------------------------------------------------------------------
# main function

def createSurface(params):

    """

    """

    # imports

    import os
    import uuid
    import shutil
    import subprocess

    import shapetools.triaUtils as triaUtils
    #import shapetools.sphericalProject as sphericalProject

    # message

    print()
    print("-------------------------------------------------------------------------")
    print()
    print("Creating surface via marching cube algorithm")
    print()
    print("-------------------------------------------------------------------------")
    print()

    # create surface via marching cube algorithm

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_mc") + " " \
        + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_02 + ".mgz") \
        + " 1 " \
        + os.path.join(params.OUTDIR, params.HEMI + ".mc." + params.internal.HSFLABEL_02 + ".vtk") + " " \
        + str(params.internal.MCC)

    print(cmd)

    subprocess.run(cmd.split())

    # update params

    params.internal.HSFLABEL_03 = "mc."+params.internal.HSFLABEL_02
    params.internal.HSFLABEL_04 = params.internal.HSFLABEL_03

    # fix topology (optional)

    if params.skipTFX is False:

        # message

        print()
        print("-------------------------------------------------------------------------")
        print()
        print("Topology fixing")
        print()
        print("-------------------------------------------------------------------------")
        print()

        # create temporary subjects dir and copy the required files

        OLD_SUBJECTS_DIR = os.environ.get('SUBJECTS_DIR')
        os.environ['SUBJECTS_DIR'] = params.OUTDIR

        TMPSUBJ = 'tmp.' + str(uuid.uuid4())

        os.makedirs(os.path.join(params.OUTDIR, TMPSUBJ, "mri"))
        os.makedirs(os.path.join(params.OUTDIR, TMPSUBJ, "surf"))

        shutil.copyfile(os.path.join(params.SUBJDIR, params.SUBJID, "mri", "brain.mgz"), os.path.join(params.OUTDIR, TMPSUBJ, "mri", "brain.mgz"))
        shutil.copyfile(os.path.join(params.SUBJDIR, params.SUBJID, "mri", "wm.mgz"), os.path.join(params.OUTDIR, TMPSUBJ, "mri", "wm.mgz"))

        # convert to freesurfer binary surface format

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mris_convert") + " " \
            + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_04 + ".vtk") + " " \
            + os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".fsmesh")

        print(cmd)

        subprocess.run(cmd.split())

        ## run spherical project
	#
        #sphericalProject.spherically_project_surface(
        #    os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".fsmesh"),
        #    os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".sphere"),
        #    os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".eig3d"))
	#
        #shutil.copyfile(
        #    os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".eig3d"),
        #    os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + ".inflated"))

        # run mris_inflate etc.

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mris_inflate") + " " \
            + os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".fsmesh") + " " \
            + os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".inflated")

        print(cmd)

        subprocess.run(cmd.split())

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mris_sphere") + " " \
            + "-q " + os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".inflated") + " " \
            + os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".sphere")

        print(cmd)

        subprocess.run(cmd.split())

        shutil.copyfile(
            os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".inflated"),
            os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + ".inflated"))

        shutil.copyfile(
            os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".fsmesh"),
            os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".fsmesh.nofix"))

        # run topology fixer (will overwrite ${HEMI}.${HSFLABEL_04}.fsmesh - or leave it unchanged)

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mris_fix_topology") + " " \
            + "-orig "  + params.internal.HSFLABEL_04 + ".fsmesh" + " " \
            + "-sphere " + params.internal.HSFLABEL_04 + ".sphere" + " " \
            + TMPSUBJ + " " \
            + params.HEMI

        print(cmd)

        subprocess.run(cmd.split())

        # copy files from temporary subjects dir to output dir

        shutil.move(
            os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".fsmesh"),
            os.path.join(params.OUTDIR, "fixed-surface", params.HEMI + ".tf." + params.internal.HSFLABEL_04 + ".fsmesh"))

        shutil.move(
            os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".sphere"),
            os.path.join(params.OUTDIR, "fixed-surface", params.HEMI + ".tf." + params.internal.HSFLABEL_04 + ".sphere"))

        #shutil.move(
        #    os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".eig3d"),
        #    os.path.join(params.OUTDIR, "fixed-surface", params.HEMI + ".tf." + params.internal.HSFLABEL_04 + ".eig3d"))

        shutil.move(
            os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".inflated"),
            os.path.join(params.OUTDIR, "fixed-surface", params.HEMI + ".tf." + params.internal.HSFLABEL_04 + ".inflated"))

        shutil.move(
            os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".fsmesh.nofix"),
            os.path.join(params.OUTDIR, "fixed-surface", params.HEMI + ".tf." + params.internal.HSFLABEL_04 + ".fsmesh.nofix"))

        # convert back to vtk

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mris_convert") + " " \
            + os.path.join(params.OUTDIR, "fixed-surface", params.HEMI + ".tf." + params.internal.HSFLABEL_04 + ".fsmesh") + " " \
            + os.path.join(params.OUTDIR, params.HEMI + ".tf." + params.internal.HSFLABEL_04 + ".vtk")

        print(cmd)

        subprocess.run(cmd.split())

        # delete temporary subjects dir

        if params.skipCLEANUP is False:

            os.remove(os.path.join(params.OUTDIR, TMPSUBJ, "mri", "brain.mgz"))
            os.remove(os.path.join(params.OUTDIR, TMPSUBJ, "mri", "wm.mgz"))
            os.remove(os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + ".inflated"))

            if os.path.exists(os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + ".defect_borders")):
                os.remove(os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + ".defect_borders"))

            if os.path.exists(os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + ".defect_chull")):
                os.remove(os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + ".defect_chull"))

            if os.path.exists(os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + ".defect_labels")):
                os.remove(os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + ".defect_labels"))

            if os.path.exists(os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + ".sulc")):
                os.remove(os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + ".sulc"))

            os.rmdir(os.path.join(params.OUTDIR, TMPSUBJ, "mri"))
            os.rmdir(os.path.join(params.OUTDIR, TMPSUBJ, "surf"))
            os.rmdir(os.path.join(params.OUTDIR, TMPSUBJ))

        # revert changes

        os.environ['SUBJECTS_DIR'] = OLD_SUBJECTS_DIR

        # update HSFLABEL 5

        params.internal.HSFLABEL_05 = "tf." + params.internal.HSFLABEL_04

    else:

        # update HSFLABEL 5

        params.internal.HSFLABEL_05 = params.internal.HSFLABEL_04

    # message

    print()
    print("-------------------------------------------------------------------------")
    print()
    print("Smooth surface")
    print()
    print("-------------------------------------------------------------------------")
    print()

    # smooth surface

    v, t = triaUtils.readVTK(os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_05 + ".vtk"))

    vs = triaUtils.smooth(v, t, n=params.internal.SMO)

    triaUtils.writeVTK(os.path.join(params.OUTDIR, params.HEMI + ".rs." + params.internal.HSFLABEL_05 + ".vtk"), vs, t)

    # remesh surface

    # ${SHAPEDNA_HOME}/triaIO --infile ${OUTDIR}/${HEMI}.${HSFLABEL_05}.vtk --outfile ${OUTDIR}/${HEMI}.rs.${HSFLABEL_05}.vtk --smooth $SMO --remeshbk 5 --info --check

    # update HSFLABEL 6

    params.internal.HSFLABEL_06 = "rs." + params.internal.HSFLABEL_05

    # return

    return(params)
