"""
This module provides the main functionality of the hippocampal shape and
thickness analysis package.

"""

import os
import sys
import shutil
import logging
import argparse
import pandas
import tempfile
import time
import traceback

from .processImage import convertFormat, cropImage, upsampleImage, copy_image_to_main
from .processLabels import createLabels, mergeMolecularLayer, autoMask, copy_labels_to_main
from .processMask import binarizeMask, gaussFilter, longFilter, closeMask, copy_mask_to_main
from .createSurface import extractSurface, remeshSurface, smoothSurface
from .createTetraMesh import createTetraMesh
from .createTetraLabels import createTetraLabels
from .removeBoundaryMask import removeBoundaryMask
from .cutTetra import cutTetra
from .computeCubeParam import computeCubeParam
from .computeThickness import computeThickness
from .utils.checkSurface import checkSurface
from .utils.mapValues import mapValues
from .utils.createSupplementaryFiles import createSupplementaryFiles
from .utils.qcPlots import qcPlots
from .cfg.config import get_defaults

# ==============================================================================
# FUNCTIONS

# ------------------------------------------------------------------------------
# get_version()

def get_version():

    from importlib import metadata

    try:
        # requires existing installation
        version = metadata.version("hipsta")
    except Exception:
        # fall-back if package is not installed, but run directly
        version = "unknown"

    return version


# ------------------------------------------------------------------------------
# get_help()

def get_help(print_help=True, return_help=False):
    """
    a function to return a help message

    """

    HELPTEXT = """

    Hippocampal shape and thickness analysis


    Purpose:

    This is a collection of scripts for hippocampal shape and thickness analysis.


    Description:

    This script performs the following major processing steps:

      1. process image
         1. conversion
         2. upsampling (optional)
         3. cropping
      2. process labels
         1. extract and merge subfield labels
         2. run automask (optional)
         3. merge molecular layer (optional)
      3. process mask
         1. binarize
         2. apply filter operations (optional)
         3. fill holes using dilation and erosion (optional)
      4. create surface
         1. extract initial surface using marching cube algorithm
         2. remesh surface (optional)
         3. smooth surface
      5. create tetrahedral mesh from triangular surface
      6. create label files for tetra mesh
      7. cut open tetrahedral mesh at anterior and posterior ends
         1. remove boundary mask
         2. mesh cutting
      8. create cube parametrization
      9. compute thickness and curvature
      10. map subfield values (and other volume-based data, optional)
      11. create supplementary files for visualization


    Usage:

      python3 shapetools.py --filename <filename> --hemi <hemi>
          --lut <lookup-table> [--outputdir <OUTDIR>] [further options]


    Input arguments:

      --filename  A segmentation file.
      --hemi      Hemisphere. Either 'lh' or 'rh'.
      --lut       A valid look-up table for hippocampal subfield segmentation.
                  Either 'freesurfer' or 'ashs' or a valid filename.
      --outputdir Directory where the results will be written. 

    Optional input arguments:

      --no-cleanup
                  Do not delete files that may be useful for diagnostic or
                  debugging purposes, but are not strictly necessary otherwise.

    Getting help:

      --help      Display this help and exit.

      --more-help Display extensive help and exit.

      --version   Display version number and exit.


    Example:

      python3 shapetools.py --filename <filename> --hemi lh --lut freesurfer
          --outputdir /my/output/directory


    Outputs:

      1. The primary output are thickness values that will be computed and stored as
         PSOL overlay files within the 'tickness' folder. These can be overlaid onto
         the mid-surface vtk file that is also found within the 'thickness' folder.

      2. Intermediate volume and surface files.

      3. An 'image' folder with intermediate results, primarily from basic image
         procesing.

      4. A 'labels' folder with intermediate results, primarily label files.

      5. A 'mask' folder with intermediate results, primarily files created
         during creation and post-processing of binary masks.

      6. A 'surface' folder with intermediate results, primarily files
         created during surface construction

      7. A 'tetra-labels' folder with intermediate results, primarily files
         created during tetrahedral mesh construction.

      8. A 'tetra-cut' folder with intermediate results, primarily files
         created during mesh cutting.

      9. A 'tetra-cube' folder with intermediate results, primarily files
         created during cube parametrization.

      10. A 'thickness' folder with thickness overlays and mid-surfaces.


    Custom segmentations:

    - If using the `ashs` segmentation, additional labels for the hippocampal head
      (label value: 20) and tail (label value: 5) labels are required. Use `--lut ashs`.
    - If using a custom look-up table (`--lut <filename>`), the expected format for
      the file is: `<numeric ID> <name> <R> <G> <B> <A>`. `R`, `G`, `B`, `A` are
      numerical values for RGB colors and transparency (alpha) and will be ignored.
      For example, `236 Subiculum 255 0 0 0`. Each line may only contain a single
      anatomical structure. Do not use a header line. The following labels need to
      be present in the look-up table: `presubiculum`,  `subiculum`, `ca1`, `ca2`,
      `ca3`, `ca4`, `head`, and `tail`. Additional labels may be present, but will
      be ignored. Lines starting with `#` will be ignored (and can be used for
      comments).
    - Multiple substructures can have the same numeric ID, e.g. `presubiculum` and
      `subiculum` can have the same numeric ID if these substructures are not
      distinguished in the segmentation. `head` and `tail` can have multiple labels
      if these are distinguished in the segmentation, but should be combined for
      the processing within the hippocampal shapetools.


    Requirements:

    1. A FreeSurfer version (6.x or 7.x) must be sourced, i.e. FREESURFER_HOME must
       exist as an environment variable and point to a valid FreeSurfer installation.

    2. A hippocampal subfield segmentation created by FreeSurfer 7.11 or later
       or the ASHS software. A custom segmentation is also permissible (some restrictions
       and settings apply; see `Custom Segmentations`).

    3. Python 3.8 or higher including the lapy, numpy, scipy, nibabel, pyvista, and
       pyacvd libraries. See `requirements.txt` for a full list.

    4. The gmsh package (verson 2.x; http://gmsh.info) must be installed. Can be
       downloaded from e.g. https://gmsh.info/bin/Linux/gmsh-2.16.0-Linux64.tgz
       or https://gmsh.info/bin/MacOSX/gmsh-2.16.0-MacOSX.dmg. The 'gmsh' binary
       must also be on the $PATH, i.e `export PATH=${PATH}:/path/to/my/gmsh`

    5. The PYTHONPATH environment variable should include the toolbox directory,
       e.g. `export PYTHONPATH=${PYTHONPATH}:/path/to/hipsta-package`.


    References:

    Please cite the following publications if you use these scripts in your work:

    - Diers, K., Baumeister, H., Jessen, F., Düzel, E., Berron, D., & Reuter, M. (2023).
      An automated, geometry-based method for hippocampal shape and thickness analysis.
      Neuroimage, 276:120182. doi: 10.1016/j.neuroimage.2023.120182.

    Please also consider citing the these publications:

    - Geuzaine, C., & Remacle, J.-F. (2009). Gmsh: a three-dimensional finite element mesh
      generator with built-in pre- and post-processing facilities. International Journal
      for Numerical Methods in Engineering, 79, 1309-1331.

    - Andreux, M., Rodola, E., Aubry, M., & Cremers, D. (2014). Anisotropic Laplace-Beltrami
      operators for shape analysis. In European Conference on Computer Vision (pp. 299-312).
      Springer, Cham.

    - Iglesias, J. E., Augustinack, J. C., Nguyen, K., Player, C. M., Player, A., Wright, M.,
      ... & Fischl, B. (2015). A computational atlas of the hippocampal formation using ex vivo,
      ultra-high resolution MRI: application to adaptive segmentation of in vivo MRI. Neuroimage,
      115, 117-137.

    - Yushkevich, P. A., Pluta, J., Wang, H., Ding, S.L., Xie, L., Gertje, E., Mancuso, L.,
      Kliot, D., Das, S. R., & Wolk, D.A. (2015). Automated Volumetry and Regional Thickness
      Analysis of Hippocampal Subfields and Medial Temporal Cortical Structures in Mild Cognitive
      Impairment. Human Brain Mapping, 36, 258-287.


    """

    if print_help:
        print(HELPTEXT)

    if return_help:
        return HELPTEXT


# ------------------------------------------------------------------------------
# start logging

def _start_logging(args):
    """
    start logging

    """

    # setup function to log uncaught exceptions
    def foo(exctype, value, tb):
        # log
        logging.error('Error Information:')
        logging.error('Type: %s', exctype)
        logging.error('Value: %s', value)
        for i in traceback.format_list(traceback.extract_tb(tb)):
            logging.error('Traceback: %s', i)
        # message
        logging.error('Program exited with ERRORS.')
    sys.excepthook = foo

    # check if output directory exists or can be created
    if os.path.isdir(args.outputdir):
        logging.info("Found output directory " + args.outputdir)
    else:
        try:
            os.mkdir(args.outputdir)
        except OSError as e:
            logging.error("Cannot create output directory " + args.outputdir)
            logging.error("Reason: " + str(e))
            raise

    # check if logfile can be written in output directory
    try:
        testfile = tempfile.TemporaryFile(dir=args.outputdir)
        testfile.close()
    except OSError as e:
        logging.error(args.outputdir + " not writeable")
        logging.error("Reason: " + str(e))
        raise

    # set up logging
    logfile_format = "[%(levelname)s: %(filename)s] %(message)s"
    logfile_handlers = [logging.StreamHandler(sys.stdout)]
    logging.basicConfig(level=logging.INFO, format=logfile_format, handlers=logfile_handlers)    

    # start logging
    logfile =  os.path.join(args.outputdir, 'logfile.txt')
    logging.getLogger().addHandler(logging.FileHandler(filename=logfile, mode="w"))

    # intial messages
    logging.info("Starting logging for hippocampal shapetools ...")
    logging.info("Logfile: %s", logfile)    
    logging.info("Version: %s", get_version())
    logging.info("Date: %s", time.strftime('%d/%m/%Y %H:%M:%S'))

    # log args
    logging.info("Command: " + " ".join(sys.argv))

    # update args
    args.logfile = logfile

    # return
    return args


# ------------------------------------------------------------------------------
# check environment and packages

def _check_environment_and_packages():

    # check environment variables
    if os.environ.get('FREESURFER_HOME') is None:
        raise RuntimeError('Need to set the FreeSurfer_HOME environment variable')

    # check python version
    if sys.version_info <= (3, 5):
        raise RuntimeError('Python version must be 3.5 or greater')

    # check for gmsh
    if shutil.which("gmsh") is None:
        raise RuntimeError('Could not find a \'gmsh\' executable')


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
        os.path.join(args.outputdir, "qc")
        ]

    # loop over list of directories

    for directory in list_of_directories:

        if not os.path.isdir(directory):

            try:
                logging.info('Creating output directory ' + directory)
                os.mkdir(directory)
            except OSError as e:
                logging.error("Cannot create output directory " + args.outputdir)
                logging.error("Reason: " + str(e))
                raise


# ------------------------------------------------------------------------------
# parse_arguments

def _parse_arguments():
    """
    an internal function to parse input arguments

    """

    # setup parser
    parser = argparse.ArgumentParser(description="This program conducts a thickness analysis of the hippocampus, based on a FreeSurfer, ASHS, or custom hippocampal subfield segmentation.",
        add_help=False)

    # required arguments (set to required=False, because not required when using --more-help)
    required = parser.add_argument_group('Required arguments')

    required.add_argument('--filename', dest="filename", help="Filename of a segmentation file.",
        default=None, metavar="<filename>", required=False)
    required.add_argument('--hemi', dest="hemi", help="Specify hemisphere, either \'lh\' or \'rh\'.",
        default=None, metavar="<lh|rh>", required=False)
    required.add_argument('--lut', dest="lut", help="Look-up table: a text file with numeric and verbal segmentation labels. \'freesurfer\' and \'ashs\' are keywords for built-in tables.",
        default=None, metavar="<freesurfer|ashs|filename>", required=False)
    required.add_argument('--outputdir', dest="outputdir", help="Directory where the results will be written.",
        default=None, metavar="<directory>", required=False)    

    # optional arguments
    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument('--no-cleanup', dest='cleanup', help="Keep files that may be useful for diagnostic or debugging purposes, but are not strictly necessary otherwise.",
        default=get_defaults("cleanup"), action="store_true", required=False)

    # expert options
    expert = parser.add_argument_group('Expert options')

    expert.add_argument('--no-crop', dest='nocrop', help=argparse.SUPPRESS,
        default=get_defaults("nocrop"), action="store_true", required=False) # help="Do not crop image.",
    expert.add_argument('--upsample', dest='upsample', help=argparse.SUPPRESS,
        default=get_defaults("upsample"), metavar="<factor x> <factor y> <factor z>", nargs='*', required=False, type=float) #  help="A list of parameters to fine-tune image upsampling.",
    expert.add_argument('--no-merge-molecular-layer', dest='noml', help=argparse.SUPPRESS, 
        default=get_defaults("noml"), action="store_true", required=False) # help="Do not merge molecular layer."
    expert.add_argument('--automated-head', dest='automaskhead', help=argparse.SUPPRESS,
        default=get_defaults("automaskhead"), action="store_true", required=False) # help="",
    expert.add_argument('--automated-tail', dest='automasktail', help=argparse.SUPPRESS,
        default=get_defaults("automasktail"), action="store_true", required=False) # help="",
    expert.add_argument('--margin-head', dest="automaskheadmargin", help=argparse.SUPPRESS,
        default=get_defaults("automaskheadmargin"), metavar="<params>", nargs=1, required=False) # help="",
    expert.add_argument('--margin-tail', dest="automasktailmargin", help=argparse.SUPPRESS,
        default=get_defaults("automasktailmargin"), metavar="<params>", nargs=1, required=False) # help="",
    expert.add_argument('--no-filter', dest='nofilter', help=argparse.SUPPRESS,
        default=get_defaults("nofilter"), action="store_true", required=False) # help="Do not filter image.",
    expert.add_argument('--no-gauss-filter', dest='nogaussfilter', help=argparse.SUPPRESS,
        default=get_defaults("nogaussfilter"), action="store_true", required=False) # help="Filter image.",
    expert.add_argument('--gauss-filter-size', dest='gaussfilter_size', help=argparse.SUPPRESS,
        default=get_defaults("gaussfilter_size"), metavar="<params>", nargs=2, required=False)
    expert.add_argument('--long-filter', dest='longfilter', help=argparse.SUPPRESS,
        default=get_defaults("longfilter"), action="store_true", required=False) # help="Filter image along longitudinal axis.",
    expert.add_argument('--long-filter-size', dest='longfilter_size', help=argparse.SUPPRESS,
        default=get_defaults("longfilter_size"), metavar="<params>", nargs=1, required=False)
    expert.add_argument('--no-close-mask', dest='noclosemask', help=argparse.SUPPRESS,
        default=get_defaults("noclosemask"), action="store_true", required=False) # help="Apply closing operation to mask.",
    expert.add_argument('--mca', dest="mca", help=argparse.SUPPRESS,
        default=get_defaults("mca"), metavar="<params>", required=False) # help="Marching-cube algorithm.",
    expert.add_argument('--mcc', dest="mcc", help=argparse.SUPPRESS,
        default=get_defaults("mcc"), metavar="<params>", required=False) # help="Marching-cube connectivity.",
    expert.add_argument('--smooth', dest="smooth", help=argparse.SUPPRESS,
        default=get_defaults("smooth"), metavar="<params>", required=False) # help="Mesh smoothing iterations.",
    expert.add_argument('--remesh', dest='remesh', help=argparse.SUPPRESS,
        default=get_defaults("remesh"), metavar="<params>", nargs='*', required=False) # help="Switch on remeshing.",
    expert.add_argument('--no-check-surface', dest='nochecksurface', help=argparse.SUPPRESS,
        default=get_defaults("nochecksurface"), action="store_true", required=False) # help="Do not check surface.",
    expert.add_argument('--no-check-boundaries', dest='nocheckboundaries', help=argparse.SUPPRESS,
        default=get_defaults("nocheckboundaries"), action="store_true", required=False) # help="Do not check boundaries.",
    expert.add_argument('--no-qc', dest='noqc', help=argparse.SUPPRESS,
        default=get_defaults("noqc"), action="store_true", required=False) # help="Skip QC",
    expert.add_argument('--cut-params', dest="cutrange", help=argparse.SUPPRESS,
        default=get_defaults("cutrange"), metavar="<params>", nargs='+', required=False, type=float) #  help="A list of parameters to fine-tune the cut operation.",
    expert.add_argument('--aniso-alpha', dest="anisoAlpha", help=argparse.SUPPRESS,
        default=get_defaults("anisoAlpha"), metavar="<params>", nargs='+', required=False) # help="Anisotropy parameter.",
    expert.add_argument('--aniso-smooth', dest="anisoSmooth", help=argparse.SUPPRESS,
        default=get_defaults("anisoSmooth"), metavar="<params>", required=False) # help="Anisotropy smoothing.",
    expert.add_argument('--thickness-params', dest="thxyz", help=argparse.SUPPRESS,
        default=get_defaults("thxyz"), metavar="<params>", nargs='+', required=False, type=float) #  help="A list of parameters to fine-tune the thickness computation.",
    expert.add_argument('--allow-ragged', dest='allowRagged', help=argparse.SUPPRESS,
        default=get_defaults("allowRagged"), action="store_true", required=False) # help="Allow ragged mid-surfaces.",
    expert.add_argument('--allow-ragged-triangles', dest='allowRaggedTriangles', help=argparse.SUPPRESS,
        default=get_defaults("allowRaggedTriangles"), action="store_true", required=False) # help="Allow triangles for ragged mid-surfaces.",
    expert.add_argument('--logfiledir', dest='logfiledir', help=argparse.SUPPRESS,
        default=get_defaults("logfiledir"), metavar="<directory>", required=False) # help="Where to store temporary logfile. Default: current working directory.",

    # Deprecated options
    expert.add_argument('--skip-orient', dest='skiporient', help=argparse.SUPPRESS,
        default=get_defaults("skiporient"), action="store_true", required=False) # help="Do not check surface.",

    # define help
    help = parser.add_argument_group('Getting help')

    help.add_argument('--help', help="Display this help message and exit", action='help')
    help.add_argument('--more-help', dest='more_help', help="Display extensive help message and exit", default=False, action="store_true", required=False)
    help.add_argument('--version', help="Display version number and exit", action='version', version='%(prog)s ' + get_version())

    # check if there are any inputs; if not, print help
    if len(sys.argv) == 1:
        args = parser.parse_args(['--help'])
    else:
        args = parser.parse_args()

    #
    if args.more_help is False:
        # check for required arguments print help and exit
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

    logging.info("Evaluating arguments ...")

    # initialize settings

    class settings:
        pass

    # filename

    FILENAME = os.path.abspath(args.filename)

    logging.info("Using " + FILENAME + " as input file.")

    # output directory

    OUTDIR = args.outputdir

    logging.info("Using " + OUTDIR + " as output directory")

    # processImage

    settings.CROP = not(args.nocrop)

    if args.upsample is None:
        settings.UPSAMPLE = None
    else:
        settings.UPSAMPLE = args.upsample

    # processLabels

    ## mergeMolecularLayer

    if args.lut == "freesurfer" and args.noml is False:
        settings.MERGE_MOLECULAR_LAYER = True
    elif args.lut == "freesurfer" and args.noml is True:
        settings.MERGE_MOLECULAR_LAYER = False
    else:
        settings.MERGE_MOLECULAR_LAYER = False

    ## automask

    settings.AUTOMASK_HEAD = args.automaskhead
    settings.AUTOMASK_TAIL = args.automasktail
    settings.AUTOMASK_HEAD_MARGIN = int(args.automaskheadmargin[0])
    settings.AUTOMASK_TAIL_MARGIN = int(args.automasktailmargin[0])

    # processMask

    settings.GAUSSFILTER = not(args.nofilter)
    settings.GAUSSFILTER_SIZE = [ float(x) for x in args.gaussfilter_size ]
    settings.LONGFILTER = args.longfilter
    settings.LONGFILTER_SIZE = (int(args.longfilter_size[0]), int(args.longfilter_size[0]), int(args.longfilter_size[0]))
    settings.CLOSEMASK = not(args.noclosemask)

    # createSurface

    ## marching cube algorithm

    if args.mca is None:
        settings.MCA = "mri_mc"
    else:
        settings.MCA = args.mca

    ## marching cube connectivity; 1=6+, 2=18, 3=6, 4=26

    if args.mcc is None:
        settings.MCC = 1
    else:
        settings.MCC = int(args.mcc)

    ## remeshing

    if args.remesh is None:
        settings.REMESH = None
    elif not args.remesh:
        settings.REMESH = 0
    else:
        settings.REMESH = int(args.remesh[0])

    ## smoothing iterations

    if args.smooth is None:
        settings.SMO = 5
    else:
        settings.SMO = int(args.smooth)

    # checkSurface

    if args.nochecksurface is True:
        settings.CHECKSURFACE = None
    else:
        settings.CHECKSURFACE = args.nochecksurface

    if args.nocheckboundaries is True:
        settings.CHECKBOUNDARIES = None
    else:
        settings.CHECKBOUNDARIES = args.nocheckboundaries

    # qc

    settings.noqc = args.noqc

    # anisotropy parameters

    if args.anisoAlpha is None:
        settings.anisoAlpha = None
    else:
        settings.anisoAlpha = [ float(x) for x in args.anisoAlpha ]

    if args.anisoSmooth is None:
        settings.anisoSmooth = 3
    else:
        settings.anisoSmooth = int(args.anisoSmooth)

    # legacy parameters

    settings.cubeWriteLegacyVTK = True

    # cutting params

    if args.cutrange is None:
        settings.cutrange = [-0.975, 0.975]
    else:
        settings.cutrange = [float(args.cutrange[0]), float(args.cutrange[1])]

    # thickness params

    if args.thxyz is None:
        settings.THXn = -0.9
        settings.THXp =  0.9
        settings.THXk =   41
        settings.THYn = -0.975
        settings.THYp =  0.975
        settings.THYk =   21
        settings.THZn = -0.9
        settings.THZp =  0.9
        settings.THZk =   11
    else:
        settings.THXn = float(args.thxyz[0])
        settings.THXp = float(args.thxyz[1])
        settings.THXk = int(args.thxyz[2])
        settings.THYn = float(args.thxyz[3])
        settings.THYp = float(args.thxyz[4])
        settings.THYk = int(args.thxyz[5])
        settings.THZn = float(args.thxyz[6])
        settings.THZp = float(args.thxyz[7])
        settings.THZk = int(args.thxyz[8])

    settings.allowRagged = args.allowRagged
    settings.allowRaggedTriangles = args.allowRaggedTriangles
    settings.skipOrient = args.skiporient

    # mapValues params

    settings.mapValuesIntegrate = "none"
    settings.mapValuesSelect = None
    settings.mapValuesInterp = "nearest"
    settings.mapValuesWritePSOL = True
    settings.mapValuesWriteMGH = True
    settings.mapValuesWriteANNOT = False

    # create the LUT

    if args.lut == "freesurfer":

        logging.info("Found internal, modified look-up table for FreeSurfer.")

        LUTLABEL = [ "tail", "head", "presubiculum", "subiculum", "ca1", "ca2", "ca3", "ca4", "dg", "ml" ]

        LUTINDEX = [ 226, [233, 235, 237, 239, 245], 234, 236, 238, 240, 240, 242, 244, 246 ]

        LUTDICT = dict(zip(LUTLABEL, LUTINDEX))

        HSFLIST = [ 234, 236, 238, 240, 246 ]

    elif args.lut == "freesurfer-noml":

        logging.info("Found internal, modified look-up table for FreeSurfer.")

        LUTLABEL = [ "tail", "head", "presubiculum", "subiculum", "ca1", "ca2", "ca3", "ca4", "dg" ]

        LUTINDEX = [ 226, [233, 235, 237, 239, 245], 234, 236, 238, 240, 240, 242, 244 ]

        LUTDICT = dict(zip(LUTLABEL, LUTINDEX))

        HSFLIST = [ 234, 236, 238, 240]

    elif args.lut == "ashs":

        logging.info("Found internal, modified look-up table for ASHS atlas.")

        LUTLABEL = [ "ca1", "ca2", "ca3", "ca4", "dg", "tail_orig", "subiculum",
            "presubiculum", "entorhinal", "ba35", "ba36", "parahippocampal",
            "head", "tail" ]

        LUTINDEX = [ 1, 2, 4, 3, 3, 5, 8, 8, 9, 10, 11, 12, 20, 5 ]

        LUTDICT = dict(zip(LUTLABEL, LUTINDEX))

        HSFLIST = [ 8, 1, 2, 4 ]

    elif os.path.isfile(args.lut):

        logging.info("Found look-up table " + args.lut)

        lut = pandas.read_csv(args.lut, sep=' ', comment='#',
            header=None, skipinitialspace=True, skip_blank_lines=True,
            error_bad_lines=False, warn_bad_lines=True)

        LUTDICT = dict(zip(lut[0], lut[1]))
        HSFLIST = list(lut[1])

    else:

        LUTDICT = dict()
        HSFLIST = []

    # add entries for tetra-labels

    LUTDICT['jointtail'] = 226
    LUTDICT['jointhead'] = 232
    LUTDICT['bndtail'] = 2260
    LUTDICT['bndhead'] = 2320
    LUTDICT['bndca4'] = 2420

    # assemble and return params

    class params:
        pass

    params.FILENAME = FILENAME
    params.OUTDIR = OUTDIR
    params.LUTDICT = LUTDICT
    params.HSFLIST = HSFLIST

    params.HEMI = args.hemi
    params.LUT = args.lut
    params.skipCLEANUP = args.cleanup
    params.LOGFILE = args.logfile

    params.internal = settings

    return params


# ------------------------------------------------------------------------------
# check params

def _check_params(params):
    """
    an internal function to check input arguments

    """

    # check for subfield segmentation file

    if os.path.isfile(params.FILENAME):
        logging.info("Found " + params.FILENAME)
    else:
        raise RuntimeError("Could not find " + params.FILENAME)
    
    # check hemisphere

    if params.HEMI != "lh" and params.HEMI != "lr":
        raise RuntimeError("Hemisphere must be either lh or rh, but not " + params.HEMI + ", exiting.")

    # check MC algorithm

    if params.internal.MCA != "mri_mc" and params.internal.MCA != "mri_tessellate" and params.internal.MCA != "skimage":
        raise RuntimeError("Could not recognise algorithm " + params.internal.MCA + ", exiting.")

    # check upsampling

    if params.internal.UPSAMPLE is not None:
        if len(params.internal.UPSAMPLE ) != 0 and len(params.internal.UPSAMPLE ) != 3:
            raise RuntimeError("Incorrect number of --upsampling parameters, must be 0 or 3.")

    # check ML

    if params.internal.MERGE_MOLECULAR_LAYER is True and params.LUT != "freesurfer":
        raise RuntimeError("Cannot use --no-merge-molecular-layer with " + params.LUT + ".")

    # check aniso alpha

    if params.internal.anisoAlpha is not None:
        if len(params.internal.anisoAlpha)>2:
            raise RuntimeError("Length of --aniso-alpha must be 1 or 2.")
        
    # check LUT

    if params.LUT != "freesurfer" and params.LUT != "ashs" and not os.path.isfile(params.LUT):
        raise RuntimeError("Look-up table can only be \'fs711\', \'ashs\', or an existing file, but not " + params.LUT)

    # return

    return params


# ------------------------------------------------------------------------------
# run analysis

def _run_analysis(params):
    """

    """

    # process image (1)

    logging.info("Starting convertFormat() ...")
    params = convertFormat(params)

    logging.info("Starting cropImage() ...")
    params = cropImage(params)

    logging.info("Starting upsampleImage() ...")
    params = upsampleImage(params)

    logging.info("Starting copy_image_to_main() ...")
    params = copy_image_to_main(params)

    # process labels (2)

    logging.info("Starting autoMask() ...")
    params = autoMask(params)

    logging.info("Starting createLabels() ...")
    params = createLabels(params)

    logging.info("Starting mergeMolecularLayer() ...")
    params = mergeMolecularLayer(params)

    logging.info("Starting copy_labels_to_main() ...")
    params = copy_labels_to_main(params)

    # process mask (3)

    logging.info("Starting binarizeMask() ...")
    params = binarizeMask(params)

    logging.info("Starting gaussFilter() ...")
    params = gaussFilter(params)

    logging.info("Starting longFilter() ...")
    params = longFilter(params)

    logging.info("Starting closeMask() ...")
    params = closeMask(params)

    logging.info("Starting copy_mask_to_main() ...")
    params = copy_mask_to_main(params)

    # create surface (4)

    logging.info("Starting extractSurface() ...")
    params = extractSurface(params)

    logging.info("Starting remeshSurface() ...")
    params = remeshSurface(params)

    logging.info("Starting smoothSurface() ...")
    params = smoothSurface(params)

    logging.info("Starting qcPlots() ...")
    params = qcPlots(params, stage="mesh")

    logging.info("Starting checkSurface() ...")
    params = checkSurface(params, stage="check_surface")

    if params.internal.continue_program is False:
        logging.info("Hippocampal shapetools finished WITH ERRORS.")
        raise AssertionError()

    # create tetra mesh for whole hippocampal body (5)

    logging.info("Starting createTetraMesh() ...")
    params = createTetraMesh(params)

    # create label files for tetra mesh (6)

    logging.info("Starting createTetraLabels() ...")
    params = createTetraLabels(params)

    # cut hippocampal body (7)

    logging.info("Starting removeBoundaryMask() ...")
    params = removeBoundaryMask(params)

    logging.info("Starting cutTetra() ...")
    params = cutTetra(params)

    logging.info("Starting checkSurface() ...")
    params = checkSurface(params, stage="check_boundaries")

    if params.internal.continue_program is False:
        logging.info("Hippocampal shapetools finished WITH ERRORS.")
        raise AssertionError()

    # compute cube parametrization (8)

    logging.info("Starting computeCubeParam() ...")
    params = computeCubeParam(params)

    logging.info("Starting qcPlots() ...")
    params = qcPlots(params, stage="profile")

    # compute thickness (9)

    logging.info("Starting computeThickness() ...")
    params = computeThickness(params)

    logging.info("Starting qcPlots() ...")
    params = qcPlots(params, stage="hull")

    # map subfield mapValues (10)

    logging.info("Starting mapValues() ...")
    params = mapValues(params)

    # create supplementary files (11)

    logging.info("Starting createSupplementaryFiles() ...")
    params = createSupplementaryFiles(params)

    # all done

    logging.info("Date: %s", time.strftime('%d/%m/%Y %H:%M:%S'))
    logging.info("Hippocampal shapetools finished without errors.")


# ------------------------------------------------------------------------------
# _run_hipsta

def _run_hipsta(args):
    """
    a function to run the shapetools submodules

    """

    # start logging
    args = _start_logging(args)

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

# ------------------------------------------------------------------------------
# run_hipsta

def run_hipsta(filename, hemisphere, lut, output_directory, **kwargs):
    """
    run the hippocampal shape and thickness analysis

    """

    # define class

    class Args:
        def __init__(self, dct=None):

            # get defaults
            self.cleanup = get_defaults("cleanup")
            self.nocrop = get_defaults("nocrop")
            self.upsample = get_defaults("upsample")
            self.noml = get_defaults("noml")
            self.automaskhead = get_defaults("automaskhead")
            self.automasktail = get_defaults("automasktail")
            self.automaskheadmargin = get_defaults("automaskheadmargin")
            self.automasktailmargin = get_defaults("automasktailmargin")
            self.nofilter = get_defaults("nofilter")
            self.nogaussfilter = get_defaults("nogaussfilter")
            self.gaussfilter_size = get_defaults("gaussfilter_size")
            self.longfilter = get_defaults("longfilter")
            self.longfilter_size = get_defaults("longfilter_size")
            self.noclosemask = get_defaults("noclosemask")
            self.mca = get_defaults("mca")
            self.mcc = get_defaults("mcc")
            self.smooth = get_defaults("smooth")
            self.remesh = get_defaults("remesh")
            self.nochecksurface = get_defaults("nochecksurface")
            self.nocheckboundaries = get_defaults("nocheckboundaries")
            self.noqc = get_defaults("noqc")
            self.cutrange = get_defaults("cutrange")
            self.anisoAlpha = get_defaults("anisoAlpha")
            self.anisoSmooth = get_defaults("anisoSmooth")
            self.thxyz = get_defaults("thxyz")
            self.allowRagged = get_defaults("allowRagged")
            self.allowRaggedTriangles = get_defaults("allowRaggedTriangles")
            self.logfiledir = get_defaults("logfiledir")
            self.skiporient = get_defaults("skiporient")

            # parse kwargs
            if dct is not None:
               for key, value in dct.items():
                   setattr(self, key, value)

    #

    my_args = Args(kwargs)

    my_args.filename = filename
    my_args.hemi = hemisphere
    my_args.lut = lut
    my_args.outputdir = output_directory

    #

    _run_hipsta(my_args)