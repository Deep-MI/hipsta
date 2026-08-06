"""
This module provides a function to check surfaces

"""

import logging
import os

from lapy import TriaMesh

from ..cfg.config import get_defaults

# ==============================================================================
# LOGGING

LOGGER = logging.getLogger(__name__)

# ==============================================================================
# FUNCTIONS


def checkSurface(params, stage=None):
    """ """

    # message

    print()
    print("--------------------------------------------------------------------------------")
    print("Check surfaces")
    print()

    if params.internal.CHECKSURFACE is True and stage == "check_surface":
        triaMesh = TriaMesh.read_vtk(os.path.join(params.OUTDIR, params.HEMI + ".surf.vtk"))

        euler = triaMesh.euler()

        LOGGER.info("Euler number for " + os.path.join(params.OUTDIR, params.HEMI + ".surf.vtk") + " is " + str(euler))

        if euler != 2:
            LOGGER.info("Surface contains holes. Please edit the corresponding hippocampal segmentation and re-run.")

            voxel_size = getattr(params.internal, "VOXEL_SIZE", None)

            if voxel_size is not None and max(voxel_size) > get_defaults("voxel_size_threshold"):
                LOGGER.info(
                    "Note that the input image has a voxel size of %s mm, which is coarse for this method. "
                    "Holes are common at this resolution even for an otherwise correct segmentation, so "
                    "please check whether a higher-resolution version of the segmentation is available "
                    "before editing it.",
                    " x ".join(format(x, ".3f") for x in voxel_size),
                )

            continue_program = False

        else:
            continue_program = True

    elif params.internal.CHECKBOUNDARIES is True and stage == "check_boundaries":
        triaMesh = TriaMesh.read_vtk(
            os.path.join(
                os.path.join(params.OUTDIR, "tetra-cut"),
                params.HEMI + ".rm.open.bnd.cut.vtk",
            )
        )

        bnd_loops = triaMesh.boundary_loops()

        LOGGER.info(
            "There are "
            + str(len(bnd_loops))
            + " boundary loops for "
            + os.path.join(
                os.path.join(params.OUTDIR, "tetra-cut"),
                params.HEMI + ".rm.open.bnd.cut.vtk",
            )
        )

        if len(bnd_loops) != 2:
            LOGGER.info(
                "Surface contains does not contain 2 boundary loops. Please retry with different cutting parameters."
            )
            continue_program = False

        else:
            continue_program = True

    else:
        continue_program = True

    #

    params.internal.continue_program = continue_program

    #

    return params
