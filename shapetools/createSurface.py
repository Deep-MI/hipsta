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
    import sys
    import uuid
    import shutil
    import subprocess

    import pyacvd
    import pyvista as pv
    import numpy as np
    import lapy as lp

    from lapy import TriaIO as lpio
    from scipy import sparse as sp

    # message

    print()
    print("-------------------------------------------------------------------------")
    print()
    print("Creating surface via marching cube algorithm")
    print()
    print("-------------------------------------------------------------------------")
    print()

    # create surface via marching cube algorithm
    
    print(params.internal.MCA)

    if params.internal.MCA == "mri_mc":

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_mc") + " " \
            + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_02 + ".mgz") \
            + " 1 " \
            + os.path.join(params.OUTDIR, params.HEMI + ".mc." + params.internal.HSFLABEL_02 + ".vtk") + " " \
            + str(params.internal.MCC)

        print(cmd)

        subprocess.run(cmd.split())

    elif params.internal.MCA == "mri_tessellate":

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_pretess") + " " \
            + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_02 + ".mgz") \
            + " xyz " \
            + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_02 + ".mgz") \
            + " " \
            + os.path.join(params.OUTDIR, params.HEMI + ".pt." + params.internal.HSFLABEL_02 + ".mgz")

        print(cmd)

        subprocess.run(cmd.split())

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_tessellate") + " " \
            + os.path.join(params.OUTDIR, params.HEMI + ".pt." + params.internal.HSFLABEL_02 + ".mgz") \
            + " 1 " \
            + os.path.join(params.OUTDIR, params.HEMI + ".mc." + params.internal.HSFLABEL_02 + ".fsmesh")

        print(cmd)

        subprocess.run(cmd.split())

        # convert from freesurfer binary surface format to vtk

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mris_convert") + " " \
            + os.path.join(params.OUTDIR, params.HEMI + ".mc." + params.internal.HSFLABEL_02 + ".fsmesh") + " " \
            + os.path.join(params.OUTDIR, params.HEMI + ".mc." + params.internal.HSFLABEL_02 + ".vtk")

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

        # set preliminary SUBJECTS_DIR

        OLD_SUBJ_DIR = os.environ["SUBJECTS_DIR"]
        os.environ["SUBJECTS_DIR"] = params.OUTDIR

        # create temporary subjects dir and copy the required files

        TMPSUBJ = 'tmp.' + str(uuid.uuid4())

        os.makedirs(os.path.join(params.OUTDIR, TMPSUBJ, "mri"))
        os.makedirs(os.path.join(params.OUTDIR, TMPSUBJ, "surf"))

        if params.internal.tfx[0] != "custom_brain":
            shutil.copyfile(params.internal.tfx[0], os.path.join(params.OUTDIR, TMPSUBJ, "mri", "brain.mgz"))
        else:
            shutil.copyfile(
                os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_02 + ".mgz"),
                os.path.join(params.OUTDIR, TMPSUBJ, "mri", "brain.mgz"))

        if params.internal.tfx[1] != "custom_wm":
            shutil.copyfile(params.internal.tfx[1], os.path.join(params.OUTDIR, TMPSUBJ, "mri", "wm.mgz"))
        else:
            shutil.copyfile(
                os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_02 + ".mgz"),
                os.path.join(params.OUTDIR, TMPSUBJ, "mri", "wm.mgz"))

            cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "fscalc") + " " \
                + os.path.join(params.OUTDIR, TMPSUBJ, "mri", "wm.mgz") + " " \
                + "mul 0" + " " \
                + "--o " + os.path.join(params.OUTDIR, TMPSUBJ, "mri", "wm.mgz")

            print(cmd)

            subprocess.run(cmd.split())

        # convert to freesurfer binary surface format

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mris_convert") + " " \
            + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_04 + ".vtk") + " " \
            + os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".fsmesh")

        print(cmd)

        subprocess.run(cmd.split())

        # remove intersections

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mris_remove_intersection") + " " \
            + os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".fsmesh") + " " \
            + os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".fsmesh")

        print(cmd)

        subprocess.run(cmd.split())

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
            os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + ".inflated"),
            os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + ".inflated.nofix"))

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

        shutil.move(
            os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".inflated"),
            os.path.join(params.OUTDIR, "fixed-surface", params.HEMI + ".tf." + params.internal.HSFLABEL_04 + ".inflated"))

        shutil.move(
            os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + "." + params.internal.HSFLABEL_04 + ".fsmesh.nofix"),
            os.path.join(params.OUTDIR, "fixed-surface", params.HEMI + ".tf." + params.internal.HSFLABEL_04 + ".fsmesh.nofix"))

        # convert back to vtk

#        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mris_convert") + " " \
#            + os.path.join(params.OUTDIR, "fixed-surface", params.HEMI + ".tf." + params.internal.HSFLABEL_04 + ".fsmesh") + " " \
#            + os.path.join(params.OUTDIR, params.HEMI + ".tf." + params.internal.HSFLABEL_04 + ".vtk")

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mris_convert") + " " \
            + os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + ".orig") + " " \
            + os.path.join(params.OUTDIR, params.HEMI + ".tf." + params.internal.HSFLABEL_04 + ".vtk")

        print(cmd)

        subprocess.run(cmd.split())

        # delete temporary subjects dir

        if params.skipCLEANUP is False:

            os.remove(os.path.join(params.OUTDIR, TMPSUBJ, "mri", "brain.mgz"))
            os.remove(os.path.join(params.OUTDIR, TMPSUBJ, "mri", "wm.mgz"))
            os.remove(os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + ".inflated"))
            os.remove(os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + ".inflated.nofix"))
            os.remove(os.path.join(params.OUTDIR, TMPSUBJ, "surf", params.HEMI + ".orig"))

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

        # revert

        os.environ["SUBJECTS_DIR"] = OLD_SUBJ_DIR

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

    # remesh

    if params.internal.REMESH is True:

        Mesh = pv.PolyData(os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_05 + ".vtk"))

        clustered = pyacvd.Clustering(Mesh)
        clustered.subdivide(4)
        clustered.cluster(Mesh.n_points)

        remeshed = clustered.create_mesh()

        vr = remeshed.points
        tr = np.reshape(remeshed.faces, (int(len(remeshed.faces)/4), 4))[:,1:]

        # check remesh results

        # creates list of edges
        trSort = np.sort(tr, axis=1)
        trSortEdges = np.concatenate((trSort[:,[0,1]],  trSort[:,[1,2]], trSort[:,[0,2]]), axis=0)

        # remove trias that have an edge that only occurs once (boundary trias)
        countEdges = np.unique(trSortEdges, axis=0, return_counts=True)
        removeEdges = countEdges[0][np.where(countEdges[1]==1)[0]]
        if len(removeEdges)>0:
            removeTrias = np.unique([ np.where(np.sum(np.logical_or(trSort==i[0], trSort==i[1]), axis=1)==2)[0] for i in removeEdges ])
            trSortRmBnd = np.delete(trSort, removeTrias, axis=0)
        else:
            trSortRmBnd = trSort.copy()

        # assure that any edge occurs exactly two times (duplicates)
        trSortRmBndEdges = np.concatenate((trSortRmBnd[:,[0,1]],  trSortRmBnd[:,[1,2]], trSortRmBnd[:,[0,2]]), axis=0)
        if len(np.where(np.unique(trSortRmBndEdges, axis=0, return_counts=True)[1]!=2)[0])!=0:
            print("Duplicate edges in mesh, exiting.")
            sys.exit(1)

        # assure that every edge must be part of exactly two different triangles (no boundary edges, no duplicates)
        countEdgesInTrias = np.array([ np.sum(np.sum(np.logical_or(trSortRmBnd==trSortRmBndEdges[i,0], trSortRmBnd==trSortRmBndEdges[i,1]), axis=1)==2) for i in range(0, len(trSortRmBndEdges)) ])
        if (countEdgesInTrias!=2).any():
            print("Boundary or duplicate edges in mesh, exiting.")
            sys.exit(1)

        # restrict to largest connected component

        triaMesh = lp.TriaMesh(v=vr, t=trSortRmBnd)
        comps = sp.csgraph.connected_components(triaMesh.adj_sym, directed=False)
        if comps[0]>1:
            compsLargest = np.argmax(np.unique(comps[1], return_counts=True)[1])
            vtcsRemove = np.where(comps[1]!=compsLargest)
            triaKeep = np.sum(np.isin(trSortRmBnd, vtcsRemove), axis=1)==0
            trSortRmBndRmComps = trSortRmBnd[triaKeep,:]
        else:
            trSortRmBndRmComps = trSortRmBnd

        # remove free vertices and re-orient mesh

        triaMesh = lp.TriaMesh(v=vr, t=trSortRmBndRmComps)
        triaMesh.rm_free_vertices_()
        triaMesh.orient_()

    else:

        # remove free vertices and re-orient mesh

        triaMesh = lpio.import_vtk(os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_05 + ".vtk"))
        triaMesh.rm_free_vertices_()
        triaMesh.orient_()

    # save

    lpio.export_vtk(triaMesh, os.path.join(params.OUTDIR, params.HEMI + ".rm." + params.internal.HSFLABEL_05 + ".vtk"))

    # smoothing

    triaMesh = lpio.import_vtk(os.path.join(params.OUTDIR, params.HEMI + ".rm." + params.internal.HSFLABEL_05 + ".vtk"))

    triaMesh.smooth_(n=params.internal.SMO)

    lpio.export_vtk(triaMesh, os.path.join(params.OUTDIR, params.HEMI + ".rs." + params.internal.HSFLABEL_05 + ".vtk"))

    # update HSFLABEL 6

    params.internal.HSFLABEL_06 = "rs." + params.internal.HSFLABEL_05

    # return

    return(params)
