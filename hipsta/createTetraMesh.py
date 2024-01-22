"""
This module provides a function to create a tetrahedral mesh

"""

import os
import shutil
import subprocess

from lapy import TetMesh, TriaMesh

# ==============================================================================
# FUNCTIONS

# ------------------------------------------------------------------------------
# _createSTL


def _createSTL(filename, v, t):
    """
    createSTL(filename, v, t)

    A function to write STL files

    """

    import numpy as np

    with open(filename, "w") as f:
        print("solid lap-data", file=f)

        for i in range(0, len(t)):
            v1mv0 = v[t[i, 1], :] - v[t[i, 0], :]
            v2mv0 = v[t[i, 2], :] - v[t[i, 0], :]

            normal = np.cross(v1mv0, v2mv0) / np.linalg.norm(np.cross(v1mv0, v2mv0))

            print("  facet normal", file=f)

            print("    outer loop", file=f)

            print("      vertex %f %f %f " % (v[t[i, 0], 0], v[t[i, 0], 1], v[t[i, 0], 2]), file=f)
            print("      vertex %f %f %f " % (v[t[i, 1], 0], v[t[i, 1], 1], v[t[i, 1], 2]), file=f)
            print("      vertex %f %f %f " % (v[t[i, 2], 0], v[t[i, 2], 1], v[t[i, 2], 2]), file=f)

            print("    endloop", file=f)

            print("  endfacet", file=f)

        print("end solid lap-data", file=f)


# ------------------------------------------------------------------------------
# createTetraMesh


def createTetraMesh(params):
    """ """

    # message

    print()
    print("--------------------------------------------------------------------------------")
    print("Create tetrahedral mesh")
    print()

    # export mesh as STL

    triaMesh = TriaMesh.read_vtk(os.path.join(params.OUTDIR, params.HEMI + ".surf.vtk"))

    _createSTL(os.path.join(params.OUTDIR, "tetra-mesh", params.HEMI + ".tetra.stl"), v=triaMesh.v, t=triaMesh.t)

    # create geofile (test.geo)

    with open(os.path.join(params.OUTDIR, "tetra-mesh", params.HEMI + ".tetra.geo"), "w") as f:
        print("Mesh.Algorithm3D=4;", file=f)
        print("Mesh.Optimize=1;", file=f)
        print("Mesh.OptimizeNetgen=1;", file=f)
        print('Merge "' + params.HEMI + ".tetra.stl" + '";', file=f)
        print("Surface Loop(1) = {1};", file=f)
        print("Volume(1) = {1};", file=f)
        print("Physical Volume(1) = {1};", file=f)

    # use gmsh to get tetras

    cmd = (
        shutil.which("gmsh")
        + " "
        + "-3 "
        + "-o "
        + os.path.join(params.OUTDIR, "tetra-mesh", params.HEMI + ".tetra-gmsh.vtk")
        + " "
        + os.path.join(params.OUTDIR, "tetra-mesh", params.HEMI + ".tetra.geo")
    )

    print(cmd)

    subprocess.run(cmd.split(), capture_output=True)

    # convert vtk2 to vtk1

    tetMesh = TetMesh.read_vtk(os.path.join(params.OUTDIR, "tetra-mesh", params.HEMI + ".tetra-gmsh.vtk"))

    TetMesh.write_vtk(tetMesh, os.path.join(params.OUTDIR, params.HEMI + ".tetra.vtk"))

    # return

    return params
