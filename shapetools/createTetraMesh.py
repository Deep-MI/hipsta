"""
This module provides a function to create a tetrahedral mesh

"""

# ------------------------------------------------------------------------------
# main function

def createTetraMesh(params):
    """

    """

    # imports

    import os
    import shutil
    import subprocess

    from lapy import TriaMesh, TetMesh

    from shapetools.triaUtils import createSTL

    # message

    print()
    print("-------------------------------------------------------------------")
    print()
    print("Create tetrahedral mesh")
    print()
    print("-------------------------------------------------------------------")
    print()

    # export mesh as STL

    triaMesh = TriaMesh.read_vtk(os.path.join(params.OUTDIR, params.HEMI + ".surf.vtk"))

    createSTL(os.path.join(params.OUTDIR, "tetra-mesh", params.HEMI + ".tetra.stl"), v=triaMesh.v, t=triaMesh.t)

    # create geofile (test.geo)

    with open(os.path.join(params.OUTDIR, "tetra-mesh", params.HEMI + ".tetra.geo"), 'w') as f:

        print("Mesh.Algorithm3D=4;", file=f)
        print("Mesh.Optimize=1;", file=f)
        print("Mesh.OptimizeNetgen=1;", file=f)
        print("Merge \"" + params.HEMI + ".tetra.stl" + "\";", file=f)
        print("Surface Loop(1) = {1};", file=f)
        print("Volume(1) = {1};", file=f)
        print("Physical Volume(1) = {1};", file=f)

    # use gmsh to get tetras

    cmd = shutil.which("gmsh") + " " \
        + "-3 " \
        + "-o " + os.path.join(params.OUTDIR, "tetra-mesh", params.HEMI + ".tetra-gmsh.vtk") + " " \
        + os.path.join(params.OUTDIR, "tetra-mesh", params.HEMI + ".tetra.geo")

    print(cmd)

    subprocess.run(cmd.split())

    # convert vtk2 to vtk1

    tetMesh = TetMesh.read_vtk(os.path.join(params.OUTDIR, "tetra-mesh", params.HEMI + ".tetra-gmsh.vtk"))

    TetMesh.write_vtk(tetMesh, os.path.join(params.OUTDIR, params.HEMI + ".tetra.vtk"))

    # return

    return(params)
