"""
This module provides a function to create surfaces, which includes the marching
cube algorithm, remeshing and smoothing.

"""

import os
import shutil
import subprocess

import nibabel as nb
import numpy as np
import pyacvd
import pyvista as pv
from lapy import TriaMesh
from scipy import sparse as sp
from skimage import measure as skm

# ==============================================================================
# FUNCTIONS

def extractSurface(params):

    """

    """

    # message

    print()
    print("--------------------------------------------------------------------------------")
    print("Creating surface via marching cube algorithm")
    print()

    # create surface via marching cube algorithm

    if params.internal.MCA == "mri_mc":

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_mc") + " " \
            + os.path.join(params.OUTDIR, params.HEMI + ".mask.mgz") \
            + " 1 " \
            + os.path.join(params.OUTDIR, "surface", params.HEMI + ".initial_surf.vtk") + " " \
            + str(params.internal.MCC)

        print(cmd)

        subprocess.run(cmd.split(), capture_output=True)

    elif params.internal.MCA == "mri_tessellate":

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_pretess") + " " \
            + os.path.join(params.OUTDIR, params.HEMI + ".mask.mgz") \
            + " xyz " \
            + os.path.join(params.OUTDIR, params.HEMI + ".mask.mgz") \
            + " " \
            + os.path.join(params.OUTDIR, "surface", params.HEMI + ".pretess.mgz")

        print(cmd)

        subprocess.run(cmd.split(), capture_output=True)

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_tessellate") + " " \
            + os.path.join(params.OUTDIR, "surface", params.HEMI + ".pretess.mgz") \
            + " 1 " \
            + os.path.join(params.OUTDIR, "surface", params.HEMI + ".surf.vtk")

        print(cmd)

        subprocess.run(cmd.split(), capture_output=True)

        # convert from freesurfer binary surface format to vtk

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mris_convert") + " " \
            + os.path.join(params.OUTDIR, "surface", params.HEMI + ".surf.vtk") + " " \
            + os.path.join(params.OUTDIR, "surface", params.HEMI + ".initial_surf.vtk")

        print(cmd)

        subprocess.run(cmd.split(), capture_output=True)

    elif params.internal.MCA == "skimage":

        img = nb.load(os.path.join(params.OUTDIR, params.HEMI + ".mask.mgz"))
        dat = img.get_fdata()

        msh = skm.marching_cubes(dat)

        v = np.matmul(img.header.get_vox2ras_tkr(), np.concatenate((msh[0], np.ones((msh[0].shape[0],1))), axis=1).T).T[:, 0:3]
        t = msh[1]

        TriaMesh(v, t).write_vtk(os.path.join(params.OUTDIR, "surface", params.HEMI + ".initial_surf.vtk"))

    # update params

    params.SURFNAME = os.path.join(params.OUTDIR, "surface", params.HEMI + ".initial_surf.vtk")

    # return

    return(params)


def remeshSurface(params):

    if params.internal.REMESH is True:

        # message

        print()
        print("--------------------------------------------------------------------------------")
        print("Remesh surface")
        print()

        if shutil.which("mris_remesh") is None:

            Mesh = pv.PolyData(params.SURFNAME)

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
                raise RuntimeError("Duplicate edges in mesh, exiting.")

            # assure that every edge must be part of exactly two different triangles (no boundary edges, no duplicates)
            countEdgesInTrias = np.array([ np.sum(np.sum(np.logical_or(trSortRmBnd==trSortRmBndEdges[i,0], trSortRmBnd==trSortRmBndEdges[i,1]), axis=1)==2) for i in range(0, len(trSortRmBndEdges)) ])
            if (countEdgesInTrias!=2).any():
                raise RuntimeError("Boundary or duplicate edges in mesh, exiting.")
                
            # restrict to largest connected component

            triaMesh = TriaMesh(v=vr, t=trSortRmBnd)
            comps = sp.csgraph.connected_components(triaMesh.adj_sym, directed=False)
            if comps[0]>1:
                compsLargest = np.argmax(np.unique(comps[1], return_counts=True)[1])
                vtcsRemove = np.where(comps[1]!=compsLargest)
                triaKeep = np.sum(np.isin(trSortRmBnd, vtcsRemove), axis=1)==0
                trSortRmBndRmComps = trSortRmBnd[triaKeep,:]
            else:
                trSortRmBndRmComps = trSortRmBnd

            # remove free vertices and re-orient mesh

            triaMesh = TriaMesh(v=vr, t=trSortRmBndRmComps)
            triaMesh.rm_free_vertices_()
            triaMesh.orient_()

        else:

            if  params.internal.REMESH_SIZE == 0:
                cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mris_remesh") + " " \
                    + "--remesh  -i " + params.SURFNAME + " " \
                    + "-o " + os.path.join(params.OUTDIR, "surface", params.HEMI + ".remeshed_surf.vtk")
            elif params.internal.REMESH_SIZE > 0:
                cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mris_remesh") + " " \
                    + "--nvert " + str(params.internal.REMESH_SIZE) + "  -i " + params.SURFNAME + " " \
                    + "-o " + os.path.join(params.OUTDIR, "surface", params.HEMI + ".remeshed_surf.vtk")

            print(cmd)

            subprocess.run(cmd.split(), capture_output=True)

            # remove free vertices and re-orient mesh

            triaMesh = TriaMesh.read_vtk(os.path.join(params.OUTDIR, "surface", params.HEMI + ".remeshed_surf.vtk"))
            triaMesh.rm_free_vertices_()
            triaMesh.orient_()

        # save

        TriaMesh.write_vtk(triaMesh, os.path.join(params.OUTDIR, "surface", params.HEMI + ".remeshed_surf.vtk"))

        # update params

        params.SURFNAME = os.path.join(params.OUTDIR, "surface", params.HEMI + ".remeshed_surf.vtk")

    # return

    return(params)


def smoothSurface(params):

    # message

    print()
    print("--------------------------------------------------------------------------------")
    print("Smooth surface")
    print()

    # read mesh

    triaMesh = TriaMesh.read_vtk(params.SURFNAME)

    # remove free vertices and re-orient mesh

    triaMesh.rm_free_vertices_()
    triaMesh.orient_()

    # smoothing

    triaMesh.smooth_(n=params.internal.SMO)

    # write mesh

    TriaMesh.write_vtk(triaMesh, os.path.join(params.OUTDIR, params.HEMI + ".surf.vtk"))

    # update params

    params.SURFNAME = os.path.join(params.OUTDIR, params.HEMI + ".surf.vtk")

    # return

    return(params)
