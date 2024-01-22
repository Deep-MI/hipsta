"""
This module provides the main functionality of the hippocampal shape and
thickness analysis package.

"""

import argparse
import logging
import os
import shutil
import sys
import time

from .cfg.atlases import get_atlases
from .cfg.config import get_defaults
from .cfg.logging import setup_logging
from .cfg.version import get_version
from .computeCubeParam import computeCubeParam
from .computeThickness import computeThickness
from .createSurface import extractSurface, remeshSurface, smoothSurface
from .createTetraLabels import createTetraLabels
from .createTetraMesh import createTetraMesh
from .cutTetra import cutTetra
from .doc.documentation import get_help_text
from .processImage import convertFormat, copy_image_to_main, cropImage, upsampleImage
from .processLabels import autoMask, copy_labels_to_main, createLabels, mergeMolecularLayer
from .processMask import binarizeMask, closeMask, copy_mask_to_main, gaussFilter, longFilter
from .removeBoundaryMask import removeBoundaryMask
from .utils.check_surface import checkSurface
from .utils.create_supplementary_files import createSupplementaryFiles
from .utils.map_values import mapValues
from .utils.qc_plots import qcPlots

# ==============================================================================
# LOGGING

LOGGER = logging.getLogger(__name__)

# ==============================================================================
# FUNCTIONS

# ------------------------------------------------------------------------------
# get_help()


def get_help(print_help=True, return_help=False):
    """
    a function to return a help message

    """

    HELPTEXT = get_help_text()

    if print_help:
        print(HELPTEXT)

    if return_help:
        return HELPTEXT


# ------------------------------------------------------------------------------
# check environment and packages


def _check_environment_and_packages():
    # check environment variables
    if os.environ.get("FREESURFER_HOME") is None:
        raise RuntimeError("Need to set the FreeSurfer_HOME environment variable")

    # check python version
    if sys.version_info <= (3, 5):
        raise RuntimeError("Python version must be 3.5 or greater")

    # check for gmsh
    if shutil.which("gmsh") is None:
        raise RuntimeError("Could not find a 'gmsh' executable")


# ------------------------------------------------------------------------------
# check directories


def _create_directories(args):
    # note that the main output directory has been created during
    # _start_logging() already

    # define list of directories

    list_of_directories = [
        os.path.join(args.outputdir, "image"),
        os.path.join(args.outputdir, "labels"),
        os.path.join(args.outputdir, "mask"),
        os.path.join(args.outputdir, "surface"),
        os.path.join(args.outputdir, "tetra-mesh"),
        os.path.join(args.outputdir, "tetra-labels"),
        os.path.join(args.outputdir, "tetra-cut"),
        os.path.join(args.outputdir, "tetra-cube"),
        os.path.join(args.outputdir, "thickness"),
        os.path.join(args.outputdir, "qc"),
    ]

    # loop over list of directories

    for directory in list_of_directories:
        if not os.path.isdir(directory):
            try:
                LOGGER.info("Creating output directory " + directory)
                os.mkdir(directory)
            except OSError as e:
                LOGGER.error("Cannot create output directory " + args.outputdir)
                LOGGER.error("Reason: " + str(e))
                raise


# ------------------------------------------------------------------------------
# parse_arguments


def _parse_arguments():
    """
    an internal function to parse input arguments

    """

    # setup parser
    parser = argparse.ArgumentParser(
        description="This program conducts a thickness analysis of the hippocampus, based on a FreeSurfer, ASHS, or custom hippocampal subfield segmentation.",
        add_help=False,
    )

    # required arguments (set to required=False, because not required when using --more-help)
    required = parser.add_argument_group("Required arguments")

    required.add_argument(
        "--filename",
        dest="filename",
        help="Filename of a segmentation file.",
        default=None,
        metavar="<filename>",
        required=False,
    )
    required.add_argument(
        "--hemi", dest="hemi", help="Hemisphere. Either 'lh' or 'rh'.", default=None, metavar="<lh|rh>", required=False
    )
    required.add_argument(
        "--lut",
        dest="lut",
        help="Look-up table: a text file with numeric and verbal segmentation labels. 'freesurfer', 'ashs-penn_abc_3t_t2', and 'ashs-umcutrecht_7t' are keywords for built-in tables.",
        default=None,
        metavar="<freesurfer|ashs-penn_abc_3t_t2|ashs-umcutrecht_7t|filename>",
        required=False,
    )
    required.add_argument(
        "--outputdir",
        dest="outputdir",
        help="Directory where the results will be written.",
        default=None,
        metavar="<directory>",
        required=False,
    )

    # optional arguments
    optional = parser.add_argument_group("Optional arguments")

    optional.add_argument(
        "--no-cleanup",
        dest="no_cleanup",
        help="Do not remove files that may be useful for diagnostic or debugging purposes, but are not necessary otherwise.",
        default=get_defaults("no_cleanup"),
        action="store_true",
        required=False,
    )
    optional.add_argument(
        "--no-crop",
        dest="no_crop",
        help="Do not crop image.",
        default=get_defaults("no_crop"),
        action="store_true",
        required=False,
    )
    optional.add_argument(
        "--upsample",
        dest="upsample",
        help="Upsample to the smallest voxel edge length.",
        default=get_defaults("upsample"),
        action="store_true",
        required=False,
    )
    optional.add_argument(
        "--upsample-size",
        dest="upsample_size",
        help="Upsampling factors. Should be between 0 and 1. If all zeros, upsample to the smallest voxel edge length. Default: 0 0 0",
        default=get_defaults("upsample_size"),
        metavar="<float>",
        nargs=3,
        required=False,
        type=float,
    )
    optional.add_argument(
        "--no-merge-molecular-layer",
        dest="no_merge_molecular_layer",
        help="Do not merge molecular layer (only applicable for FreeSurfer segmentations).",
        default=get_defaults("no_merge_molecular_layer"),
        action="store_true",
        required=False,
    )
    optional.add_argument(
        "--automask-head",
        dest="automask_head",
        help="Automated boundary detection for hippocampal head.",
        default=get_defaults("automask_head"),
        action="store_true",
        required=False,
    )
    optional.add_argument(
        "--automask-tail",
        dest="automask_tail",
        help="Automated boundary detection for hippocampal tail.",
        default=get_defaults("automask_tail"),
        action="store_true",
        required=False,
    )
    optional.add_argument(
        "--automask-head-margin",
        dest="automask_head_margin",
        help="Margin for automated boundary detection for hippocampal head. Default: 0",
        default=get_defaults("automask_head_margin"),
        metavar="<int>",
        required=False,
        type=int,
    )
    optional.add_argument(
        "--automask-tail-margin",
        dest="automask_tail_margin",
        help="Margin for automated boundary detection for hippocampal tail. Default: 0",
        default=get_defaults("automask_tail_margin"),
        metavar="<int>",
        required=False,
        type=int,
    )
    optional.add_argument(
        "--no-gauss-filter",
        dest="no_gauss_filter",
        help="Do not apply gaussian filter.",
        default=get_defaults("no_gauss_filter"),
        action="store_true",
        required=False,
    )
    optional.add_argument(
        "--gauss-filter-size",
        dest="gauss_filter_size",
        help="Filter width and threshold for gaussian filtering. Default: 1 50.",
        default=get_defaults("gauss_filter_size"),
        metavar="<float>",
        nargs=2,
        required=False,
        type=float,
    )
    optional.add_argument(
        "--long-filter",
        dest="long_filter",
        help="Apply filter along longitudinal axis, i.e. attempt to create smooth transitions between slices.",
        default=get_defaults("long_filter"),
        action="store_true",
        required=False,
    )
    optional.add_argument(
        "--long-filter-size",
        dest="long_filter_size",
        help="Size of longitudinal filter. Default: 5",
        default=get_defaults("long_filter_size"),
        metavar="<int>",
        required=False,
        type=int,
    )
    optional.add_argument(
        "--no-close-mask",
        dest="no_close_mask",
        help="Do not apply closing operation to mask, i.e. do not attempt to close small holes.",
        default=get_defaults("no_close_mask"),
        action="store_true",
        required=False,
    )
    optional.add_argument(
        "--mca",
        dest="mca",
        help="Type of marching-cube algorithm. Either 'mri_mc' or 'skimage'. Default: 'skimage'",
        default=get_defaults("mca"),
        metavar="<mri_mc|skimage>",
        required=False,
    )
    optional.add_argument(
        "--remesh",
        dest="remesh",
        help="Apply remeshing operation to surface, i.e. create a regular, evenly spaced surface grid.",
        default=get_defaults("remesh"),
        action="store_true",
        required=False,
    )
    optional.add_argument(
        "--smooth",
        dest="smooth",
        help="Mesh smoothing iterations. Default: 5",
        default=get_defaults("smooth"),
        metavar="<int>",
        required=False,
        type=int,
    )
    optional.add_argument(
        "--cut-range",
        dest="cut_range",
        help="Range for tetrahedral boundary cutting. Default: -0.975, 0.975",
        default=get_defaults("cut_range"),
        metavar="<float>",
        nargs=2,
        required=False,
        type=float,
    )
    optional.add_argument(
        "--aniso-alpha",
        dest="aniso_alpha",
        help="Anisotropy parameter(s). Can be one or two or two numbers. Default: 40",
        default=get_defaults("aniso_alpha"),
        metavar="<float>",
        nargs="+",
        required=False,
        type=float,
    )
    optional.add_argument(
        "--aniso-smooth",
        dest="aniso_smooth",
        help="Anisotropy smoothing iterations. Default: 3",
        default=get_defaults("aniso_smooth"),
        metavar="<int>",
        required=False,
        type=int,
    )
    optional.add_argument(
        "--thickness-grid",
        dest="thickness_grid",
        help="Extent and resolution of the grid used for thickness computation; three lists of three numbers: negative extent of x axis, positive extent of x axis, resolution on x axis. Repeat for the y and z axes. Default: -0.9 0.9 41 -0.975 0.975 21 -0.9 0.9 11",
        default=get_defaults("thickness_grid"),
        metavar="<float>",
        nargs=9,
        required=False,
        type=float,
    )

    # expert options
    expert = parser.add_argument_group("Expert options")

    expert.add_argument(
        "--mcc",
        dest="mcc",
        help=argparse.SUPPRESS,
        default=get_defaults("mcc"),
        metavar="<int>",
        required=False,
        type=int,
    )  # help="Marching-cube connectivity. Only used for \'mri_mc\' algorithm. Default: 1",
    expert.add_argument(
        "--remesh-size",
        dest="remesh_size",
        help=argparse.SUPPRESS,
        default=get_defaults("remesh_size"),
        metavar="<int>",
        required=False,
        type=int,
    )  # help="Target number of vertices for remeshing operation. If zero, keep original number of vertices. Default: 0",
    expert.add_argument(
        "--no-check-surface",
        dest="no_check_surface",
        help=argparse.SUPPRESS,
        default=get_defaults("no_check_surface"),
        action="store_true",
        required=False,
    )  # help="Do not check surface and proceed if there are holes.",
    expert.add_argument(
        "--no-check-boundaries",
        dest="no_check_boundaries",
        help=argparse.SUPPRESS,
        default=get_defaults("no_check_boundaries"),
        action="store_true",
        required=False,
    )  # help="Do not check boundaries and proceed if there are less / more than two continuous boundary loops.",
    expert.add_argument(
        "--no-qc",
        dest="no_qc",
        help=argparse.SUPPRESS,
        default=get_defaults("no_qc"),
        action="store_true",
        required=False,
    )  # help="Do not perform QC",
    expert.add_argument(
        "--allow-ragged-surfaces",
        dest="allow_ragged_surfaces",
        help=argparse.SUPPRESS,
        default=get_defaults("allow_ragged_surfaces"),
        action="store_true",
        required=False,
    )  # help="Allow ragged mid-surfaces.",
    expert.add_argument(
        "--allow-ragged-trias",
        dest="allow_ragged_trias",
        help=argparse.SUPPRESS,
        default=get_defaults("allow_ragged_trias"),
        action="store_true",
        required=False,
    )  # help="Allow triangles for ragged mid-surfaces.",

    # deprecated options - to be removed in the future
    deprecated = parser.add_argument_group("Deprecated options")

    deprecated.add_argument(
        "--no-orient",
        dest="no_orient",
        help=argparse.SUPPRESS,
        default=get_defaults("no_orient"),
        action="store_true",
        required=False,
    )  # help="Do not check surface.",

    # define help
    help = parser.add_argument_group("Getting help")

    help.add_argument("--help", help="Display this help message and exit", action="help")
    help.add_argument(
        "--more-help",
        dest="more_help",
        help="Display extensive help message and exit",
        default=False,
        action="store_true",
        required=False,
    )
    help.add_argument(
        "--version", help="Display version number and exit", action="version", version="%(prog)s " + get_version()
    )

    # check if there are any inputs; if not, print help
    if len(sys.argv) == 1:
        args = parser.parse_args(["--help"])
    else:
        args = parser.parse_args()

    #
    if args.more_help is False:
        # check for required arguments, print message and exit
        if args.filename is None:
            print("ERROR: the --filename argument is required, exiting. Use --help to see details.")
            print()
            sys.exit(1)

        if args.hemi is None:
            print("ERROR: the --hemi argument is required, exiting. Use --help to see details.")
            print()
            sys.exit(1)

        if args.lut is None:
            print("ERROR: the --lut argument is required, exiting. Use --help to see details.")
            print()
            sys.exit(1)

        if args.outputdir is None:
            print("ERROR: the --outputdir argument is required, exiting. Use --help to see details.")
            print()
            sys.exit(1)

    #
    return args


# ------------------------------------------------------------------------------
# evaluate arguments


def _evaluate_args(args):
    """
    an internal function to set and return internal parameters

    """

    # message

    LOGGER.info("Evaluating arguments ...")

    # initialize settings

    class settings:
        pass

    # clean up

    settings.CLEANUP = not (args.no_cleanup)

    # processImage

    settings.CROP = not (args.no_crop)
    settings.UPSAMPLE = args.upsample
    settings.UPSAMPLE_SIZE = args.upsample_size

    # processLabels

    ## mergeMolecularLayer

    if args.lut == "freesurfer" and args.no_merge_molecular_layer is False:
        settings.MERGE_MOLECULAR_LAYER = True
    elif args.lut == "freesurfer" and args.no_merge_molecular_layer is True:
        settings.MERGE_MOLECULAR_LAYER = False
    else:
        settings.MERGE_MOLECULAR_LAYER = False

    ## automask

    settings.AUTOMASK_HEAD = args.automask_head
    settings.AUTOMASK_TAIL = args.automask_tail
    settings.AUTOMASK_HEAD_MARGIN = args.automask_head_margin
    settings.AUTOMASK_TAIL_MARGIN = args.automask_tail_margin

    # processMask

    settings.GAUSSFILTER = not (args.no_gauss_filter)
    settings.GAUSSFILTER_SIZE = args.gauss_filter_size
    settings.LONGFILTER = args.long_filter
    settings.LONGFILTER_SIZE = args.long_filter_size
    settings.CLOSEMASK = not (args.no_close_mask)

    # createSurface

    ## marching cube algorithm

    settings.MCA = args.mca

    ## marching cube connectivity; 1=6+, 2=18, 3=6, 4=26

    settings.MCC = args.mcc

    ## remeshing

    settings.REMESH = args.remesh
    settings.REMESH_SIZE = args.remesh_size

    ## smoothing iterations

    settings.SMO = args.smooth

    # qc params

    settings.no_qc = args.no_qc

    # checkSurface

    settings.CHECKSURFACE = not (args.no_check_surface)
    settings.CHECKBOUNDARIES = not (args.no_check_boundaries)

    # cutTetra

    settings.cut_range = args.cut_range

    # computeCubeParam

    settings.aniso_alpha = args.aniso_alpha
    settings.aniso_smooth = args.aniso_smooth

    # computeThickness

    settings.THXn = args.thickness_grid[0]
    settings.THXp = args.thickness_grid[1]
    settings.THXk = int(args.thickness_grid[2])
    settings.THYn = args.thickness_grid[3]
    settings.THYp = args.thickness_grid[4]
    settings.THYk = int(args.thickness_grid[5])
    settings.THZn = args.thickness_grid[6]
    settings.THZp = args.thickness_grid[7]
    settings.THZk = int(args.thickness_grid[8])
    settings.allow_ragged_surfaces = args.allow_ragged_surfaces
    settings.allow_ragged_trias = args.allow_ragged_trias
    settings.no_orient = args.no_orient

    # mapValues

    settings.mapValuesIntegrate = get_defaults("map_values_integrate")
    settings.mapValuesSelect = get_defaults("map_values_select")
    settings.mapValuesInterp = get_defaults("map_values_interp")
    settings.mapValuesWritePSOL = get_defaults("map_values_write_psol")
    settings.mapValuesWriteMGH = get_defaults("map_values_write_mgh")
    settings.mapValuesWriteANNOT = get_defaults("map_values_write_annot")

    # create the LUT

    LUTDICT, HSFLIST = get_atlases(args.lut)

    # assemble and return params

    class params:
        pass

    params.FILENAME = os.path.abspath(args.filename)
    params.HEMI = args.hemi
    params.LUT = args.lut
    params.OUTDIR = args.outputdir

    params.LUTDICT = LUTDICT
    params.HSFLIST = HSFLIST

    params.internal = settings

    LOGGER.info("Using " + params.FILENAME + " as input file")
    LOGGER.info("Using " + params.OUTDIR + " as output directory")

    return params


# ------------------------------------------------------------------------------
# check params


def _check_params(params):
    """
    an internal function to check input arguments

    """

    # check for subfield segmentation file

    if os.path.isfile(params.FILENAME):
        LOGGER.info("Found " + params.FILENAME)
    else:
        raise RuntimeError("Could not find " + params.FILENAME)

    # check hemisphere

    if params.HEMI != "lh" and params.HEMI != "rh":
        raise RuntimeError("Hemisphere must be either lh or rh, but not " + params.HEMI + ", exiting.")

    # check MC algorithm

    if params.internal.MCA != "mri_mc" and params.internal.MCA != "mri_tessellate" and params.internal.MCA != "skimage":
        raise RuntimeError("Could not recognise marching cube algorithm " + params.internal.MCA + ", exiting.")

    # check upsampling

    if any(x < 0 for x in params.internal.UPSAMPLE_SIZE) or any(x > 1 for x in params.internal.UPSAMPLE_SIZE):
        raise RuntimeError("All upsampling parameters should be between 0 and 1.")

    # check ML

    if params.internal.MERGE_MOLECULAR_LAYER is True and params.LUT != "freesurfer":
        raise RuntimeError("Cannot use no-merge-molecular-layer with " + params.LUT + ".")

    if params.internal.MERGE_MOLECULAR_LAYER is True and ("ml" in params.LUTDICT.keys()) is False:
        LOGGER.info("Could not find molecular layer in lookup table.")
        raise AssertionError()

    # check longfilter size

    if params.internal.LONGFILTER_SIZE < 2 or (params.internal.LONGFILTER_SIZE % 2) == 0:
        raise RuntimeError("Longfilter size must be an uneven integer greater than 2.")

    # check aniso alpha

    if len(params.internal.aniso_alpha) > 2:
        raise RuntimeError("Length of aniso-alpha must be 1 or 2.")

    # check LUT

    if (
        params.LUT != "freesurfer"
        and params.LUT != "ashs-penn_abc_3t_t2"
        and params.LUT != "ashs-umcutrecht_7t"
        and not os.path.isfile(params.LUT)
    ):
        raise RuntimeError(
            "Look-up table can only be 'freesurfer', 'ashs-penn_abc_3t_t2', 'ashs-umcutrecht_7t', or an existing file, but not "
            + params.LUT
        )

    # return

    return params


# ------------------------------------------------------------------------------
# run analysis


def _run_analysis(params):
    """ """

    # process image (1)

    LOGGER.info("Starting convertFormat() ...")
    params = convertFormat(params)

    LOGGER.info("Starting cropImage() ...")
    params = cropImage(params)

    LOGGER.info("Starting upsampleImage() ...")
    params = upsampleImage(params)

    LOGGER.info("Starting copy_image_to_main() ...")
    params = copy_image_to_main(params)

    # process labels (2)

    LOGGER.info("Starting autoMask() ...")
    params = autoMask(params)

    LOGGER.info("Starting createLabels() ...")
    params = createLabels(params)

    LOGGER.info("Starting mergeMolecularLayer() ...")
    params = mergeMolecularLayer(params)

    LOGGER.info("Starting copy_labels_to_main() ...")
    params = copy_labels_to_main(params)

    # process mask (3)

    LOGGER.info("Starting binarizeMask() ...")
    params = binarizeMask(params)

    LOGGER.info("Starting gaussFilter() ...")
    params = gaussFilter(params)

    LOGGER.info("Starting longFilter() ...")
    params = longFilter(params)

    LOGGER.info("Starting closeMask() ...")
    params = closeMask(params)

    LOGGER.info("Starting copy_mask_to_main() ...")
    params = copy_mask_to_main(params)

    # create surface (4)

    LOGGER.info("Starting extractSurface() ...")
    params = extractSurface(params)

    LOGGER.info("Starting remeshSurface() ...")
    params = remeshSurface(params)

    LOGGER.info("Starting smoothSurface() ...")
    params = smoothSurface(params)

    LOGGER.info("Starting qcPlots() ...")
    params = qcPlots(params, stage="mesh")

    LOGGER.info("Starting checkSurface() ...")
    params = checkSurface(params, stage="check_surface")

    if params.internal.continue_program is False:
        LOGGER.info("Hipsta finished WITH ERRORS.")
        raise AssertionError("Check surface failed (stage: surface)")

    # create tetra mesh for whole hippocampal body (5)

    LOGGER.info("Starting createTetraMesh() ...")
    params = createTetraMesh(params)

    # create label files for tetra mesh (6)

    LOGGER.info("Starting createTetraLabels() ...")
    params = createTetraLabels(params)

    # cut hippocampal body (7)

    LOGGER.info("Starting removeBoundaryMask() ...")
    params = removeBoundaryMask(params)

    LOGGER.info("Starting cutTetra() ...")
    params = cutTetra(params)

    LOGGER.info("Starting checkSurface() ...")
    params = checkSurface(params, stage="check_boundaries")

    if params.internal.continue_program is False:
        LOGGER.info("Hipsta finished WITH ERRORS.")
        raise AssertionError("Check surface failed (stage: boundaries)")

    # compute cube parametrization (8)

    LOGGER.info("Starting computeCubeParam() ...")
    params = computeCubeParam(params)

    LOGGER.info("Starting qcPlots() ...")
    params = qcPlots(params, stage="profile")

    # compute thickness (9)

    LOGGER.info("Starting computeThickness() ...")
    params = computeThickness(params)

    LOGGER.info("Starting qcPlots() ...")
    params = qcPlots(params, stage="hull")

    # map subfield mapValues (10)

    LOGGER.info("Starting mapValues() ...")
    params = mapValues(params)

    # create supplementary files (11)

    LOGGER.info("Starting createSupplementaryFiles() ...")
    params = createSupplementaryFiles(params)

    # all done

    LOGGER.info("Date: %s", time.strftime("%d/%m/%Y %H:%M:%S"))
    LOGGER.info("Hipsta finished without errors.")


# ------------------------------------------------------------------------------
# _run_hipsta


def _run_hipsta(args):
    """
    a function to run the hipsta submodules

    """

    #
    try:
        #
        setup_logging(args)

        # check environment and packages
        _check_environment_and_packages()

        # create directories
        _create_directories(args)

        # convert arguments to params
        params = _evaluate_args(args)

        # check params
        _check_params(params)

        # run analysis
        _run_analysis(params)

    except Exception as e:
        LOGGER.error("Error Information:")
        LOGGER.error("Type: %s", type(e))
        LOGGER.error("Value: %s", e.args)

        raise


# ------------------------------------------------------------------------------
# run_hipsta


def run_hipsta(filename, hemi, lut, outputdir, **kwargs):
    """
    Run the hippocampal shape and thickness analysis

    Parameters
    ----------

    filename :
        Filename of a segmentation file.
    hemi :
        Hemisphere. Either \'lh\' or \'rh\'.
    lut :
        Look-up table: a text file with numeric and verbal segmentation labels. \'freesurfer\', 'ashs-penn_abc_3t_t2' and \'ashs-umcutrecht_7t\' are keywords for built-in tables.
    outputdir
        Directory where the results will be written.

    no_cleanup : optional
        Keep files that may be useful for diagnostic or debugging purposes, but are not necessary otherwise. Default: False
    no_crop : optional
        Do not crop image. Default: False.
    upsample : optional
        Upsample to the smallest voxel edge length. Default: False
    upsample_size : optional
        Upsampling factors. Should be between 0 and 1. If all zeros, upsample to the smallest voxel edge length. Default: [0, 0, 0]
    no_merge_molecular_layer : optional
        Do not merge molecular layer (only applicable for FreeSurfer segmentations). Default: False
    automask_head : optional
        Automated boundary detection for hippocampal head. Default: False
    automask_tail : optional
        Automated boundary detection for hippocampal tail. Default: False
    automask_head_margin : optional
        Margin for automated boundary detection for hippocampal head. Default: 0
    automask_tail_margin : optional
        Margin for automated boundary detection for hippocampal tail. Default: 0
    no_gauss_filter : optional
        Do not apply gaussian filter. Default: False
    gauss_filter_size : optional
        Filter width and threshold for gaussian filtering. Default: [1, 50]
    long_filter : optional
        Apply filter along longitudinal axis, i.e. attempt to create smooth transitions between slices. Default: False
    long_filter_size : optional
        Size of longitudinal filter. Default: 5
    no_close_mask : optional
        Do not apply closing operation to mask, i.e. do not attempt to close small holes. Default: False
    mca : optional
        Type of marching-cube algorithm. Either \'mri_mc\' or \'skimage\'. Default: \'mri_mc\'
    remesh : optional
        Apply remeshing operation to surface, i.e. create a regular, evenly spaced surface grid. Default: False
    smooth : optional
        Mesh smoothing iterations. Default: 5
    cut_range : optional
        Range for tetrahedral boundary cutting. Default: [-0.975, 0.975]
    aniso_alpha : optional
        Anisotropy parameter(s). Can be one or two numbers. Default: [40]
    aniso_smooth : optional
        Anisotropy smoothing iterations parameter. Default: 3
    thickness_grid : optional
        Extent and resolution of the grid used for thickness computation; three lists of three numbers: negative extent of x axis, positive extent of x axis, resolution on x axis. Repeat for the y and z axes. Default: [ -0.9, 0.9, 41, -0.975, 0.975, 21, -0.9, 0.9, 11]

    mcc : optional
        Marching-cube connectivity. Only used for \'mri_mc\' algorithm. Default: 1
    remesh_size : optional
        Target number of vertices for remeshing operation. If zero, keep original number of vertices. Default: 0
    no_check_surface : optional
        Do not check surface and proceed if there are holes. Default: False
    no_check_boundaries : optional
        Do not check boundaries and proceed if there are less / more than two continuous boundary loops. Default: False
    no_qc : optional
        Do not perform QC. Default: False
    allow_ragged_surfaces
        Allow ragged mid-surfaces. Default: False
    allow_ragged_trias
        Allow triangles for ragged mid-surfaces. Default: False

    no_orient : optional
        Do not orient surfaces. Default: False

    Returns
    -------
    dict
        a dictionary of input arguments
    """

    # define class

    class Args:
        def __init__(self, dct=None):
            # get defaults
            self.no_cleanup = get_defaults("no_cleanup")
            self.no_crop = get_defaults("no_crop")
            self.upsample = get_defaults("upsample")
            self.upsample_size = get_defaults("upsample_size")
            self.no_merge_molecular_layer = get_defaults("no_merge_molecular_layer")
            self.automask_head = get_defaults("automask_head")
            self.automask_tail = get_defaults("automask_tail")
            self.automask_head_margin = get_defaults("automask_head_margin")
            self.automask_tail_margin = get_defaults("automask_tail_margin")
            self.no_gauss_filter = get_defaults("no_gauss_filter")
            self.gauss_filter_size = get_defaults("gauss_filter_size")
            self.long_filter = get_defaults("long_filter")
            self.long_filter_size = get_defaults("long_filter_size")
            self.no_close_mask = get_defaults("no_close_mask")
            self.mca = get_defaults("mca")
            self.remesh = get_defaults("remesh")
            self.smooth = get_defaults("smooth")
            self.cut_range = get_defaults("cut_range")
            self.aniso_alpha = get_defaults("aniso_alpha")
            self.aniso_smooth = get_defaults("aniso_smooth")
            self.thickness_grid = get_defaults("thickness_grid")
            # expert options
            self.mcc = get_defaults("mcc")
            self.remesh_size = get_defaults("remesh_size")
            self.no_check_surface = get_defaults("no_check_surface")
            self.no_check_boundaries = get_defaults("no_check_boundaries")
            self.no_qc = get_defaults("no_qc")
            self.allow_ragged_surfaces = get_defaults("allow_ragged_surfaces")
            self.allow_ragged_trias = get_defaults("allow_ragged_trias")
            # deprecated options
            self.no_orient = get_defaults("no_orient")

            # parse kwargs
            if dct is not None:
                for key, value in dct.items():
                    if key in vars(self).keys():
                        setattr(self, key, value)
                    else:
                        raise AssertionError("Key " + key + " not included in list of default keys")

    #

    my_args = Args(kwargs)

    my_args.filename = filename
    my_args.hemi = hemi
    my_args.lut = lut
    my_args.outputdir = outputdir

    #

    _run_hipsta(my_args)

    #

    return vars(my_args)
