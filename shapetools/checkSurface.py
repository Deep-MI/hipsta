"""
This module provides a function to check surfaces

"""

def checkSurface(params, stage=None):

    """

    """

    # imports

    import os
    import sys
    import logging

    from lapy import TriaMesh

    # message

    print()
    print("-------------------------------------------------------------------------")
    print()
    print("Check surfaces")
    print()
    print("-------------------------------------------------------------------------")
    print()

    if params.internal.CHECKSURFACE is not None and stage=="check_surface":

        triaMesh = TriaMesh.read_vtk(os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_06 + ".vtk"))

        euler = triaMesh.euler()

        logging.info("Euler number for " + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_06 + ".vtk") + " is " + str(euler))

        if euler!=2:

            logging.info("Surface contains holes. Please edit the corresponding hippocampal segmentation and re-run.")
            continue_program = False

        else:

            continue_program = True

    elif params.internal.CHECKBOUNDARIES is not None and stage=="check_boundaries":

        triaMesh = TriaMesh.read_vtk(os.path.join(os.path.join(params.OUTDIR, "tetra-cut"), params.HEMI + '.rm.open.bnd.cut.tetra.vtk'))

        bnd_loops = triaMesh.boundary_loops()

        logging.info("There are " + str(len(bnd_loops)) + " boundary loops for " + os.path.join(os.path.join(params.OUTDIR, "tetra-cut"), params.HEMI + '.rm.open.bnd.cut.tetra.vtk'))

        if len(bnd_loops)!=2:

            logging.info("Surface contains does not contain 2 boundary loops. Please retry with different cutting parameters.")
            continue_program = False

        else:

            continue_program = True

    else:

        continue_program = True

    #

    return continue_program, params
