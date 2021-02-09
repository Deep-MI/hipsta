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

    from shapetools.triaUtils import createSTL
    from shapetools.triaUtils import readVTK, writeVTK

    # message

    print()
    print("-------------------------------------------------------------------")
    print()
    print("Create tetrahedral mesh")
    print()
    print("-------------------------------------------------------------------")
    print()

    # smooth, remesh BK, check, etc.
    # (this should have happened already)

    # create temporary directory

    os.mkdir(os.path.join(params.OUTDIR, "tetra-tmp"))

    # export mesh as STL

    v, t = readVTK(os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_06 + ".vtk"))

    createSTL(os.path.join(params.OUTDIR, "tetra-tmp", params.HEMI + "." + params.internal.HSFLABEL_06 + ".stl"), v=v, t=t)

    # create geofile (test.geo)

    with open(os.path.join(params.OUTDIR, "tetra-tmp", params.HEMI + "." + params.internal.HSFLABEL_06 + ".geo"), 'w') as f:

        print("Mesh.Algorithm3D=4;", file=f)
        print("Mesh.Optimize=1;", file=f)
        print("Mesh.OptimizeNetgen=1;", file=f)
        print("Merge \"" + params.HEMI + "." + params.internal.HSFLABEL_06 + ".stl" + "\";", file=f)
        print("Surface Loop(1) = {1};", file=f)
        print("Volume(1) = {1};", file=f)
        print("Physical Volume(1) = {1};", file=f)

    # use gmsh to get tetras

    cmd = shutil.which("gmsh") + " " \
        + "-3 " \
        + "-o " + os.path.join(params.OUTDIR, "tetra-tmp", params.HEMI + ".tet." + params.internal.HSFLABEL_06 + ".vtk2.vtk") + " " \
        + os.path.join(params.OUTDIR, "tetra-tmp", params.HEMI + "." + params.internal.HSFLABEL_06 + ".geo")

    print(cmd)

    subprocess.run(cmd.split())

    # convert vtk2 to vtk1

    v, t = readVTK(os.path.join(params.OUTDIR, "tetra-tmp", params.HEMI + ".tet." + params.internal.HSFLABEL_06 + ".vtk2.vtk"))

    writeVTK(os.path.join(params.OUTDIR, params.HEMI + ".tet." + params.internal.HSFLABEL_06 + ".vtk"), v, t)

    # clean up

    if params.skipCLEANUP is False:

        os.remove(os.path.join(params.OUTDIR, "tetra-tmp", params.HEMI + "." + params.internal.HSFLABEL_06 + ".geo"))
        os.remove(os.path.join(params.OUTDIR, "tetra-tmp", params.HEMI + "." + params.internal.HSFLABEL_06 + ".stl"))
        os.remove(os.path.join(params.OUTDIR, "tetra-tmp", params.HEMI + ".tet." + params.internal.HSFLABEL_06 + ".vtk2.vtk"))
        os.rmdir(os.path.join(params.OUTDIR, "tetra-tmp"))

    # update params

    params.internal.HSFLABEL_07 = "tet." + params.internal.HSFLABEL_06

    # return

    return(params)
