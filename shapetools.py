"""
This module provides the main functionality of the hippocampal shape tools
package.

"""


# ==============================================================================
# FUNCTIONS

# ------------------------------------------------------------------------------
# get_version()

def get_version():

    VERSION = "v0.7.0-beta"

    return VERSION

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

      1. create labels
         1. extract subfield labels from label image
         2. merge individual labels
      2. attach molecular layer to neighboring regions
      3. create mask
         1. fill holes using dilation and erosion
         2. create final mask
      4. create surface
         1. create initial surface using marching cube algorithm
         2. topology fixing (optional)
         3. smooth (fixed) surface
      5. create tetrahedral mesh from triangular surface
      6. create label files for tetra mesh
      7. remove boundary mask (preprocessing for mesh cutting)
      8. cut open tetrahedral mesh at anterior and posterior ends
      9. create cube parametrization
      10. compute thickness and curvature
      11. map subfield values (and other volume-based data, optional)
      12. create supplementary files for visualization


    Usage:

      python3 shapetools.py --filename <filename> --hemi <hemi>
          --lut <lookup-table> [--outputdir <OUTDIR>] [further options]


    Input arguments:

      --filename  A segmentation file.
      --hemi      Hemisphere. Either 'lh' or 'rh'.
      --lut       A valid look-up table for hippocampal subfield segmentation.
                  Either 'fs711' or 'ashs' or a valid filename.

    Optional input arguments:

      --outputdir Directory where the results will be written. If not given, a
                  subfolder within each subject's directory will be created.

      --no-cleanup
                  Do not delete files that may be useful for diagnostic or
                  debugging purposes, but are not strictly necessary otherwise.

    Getting help:

      --help      Display this help and exit.

      --more-help Display extensive help and exit.

      --version   Display version number and exit.


    Example:

      python3 shapetools.py --filename <filename> --hemi lh --lut fs711
          --outputdir /my/output/directory


    Outputs:

      1. The primary output are thickness values that will be computed and stored as
         PSOL overlay files within the 'tickness' folder. These can be overlaid onto
         the mid-surface vtk file that is also found within the 'thickness' folder.

      2. Intermediate volume and surface files, each prefixed with the particular
         operation performed. The basic pattern for the filenames is:

         <hemisphere>.<prefixes>.<hsf-labels>.<suffixes>

         Hemisphere:
         - lh : left hemisphere
         - rh : right hemisphere

         Prefixes:
         - ml : merged molecular layer
         - de : dilation and erosion
         - mc : marching cube algorithm
         - tf : topological fixing
         - sm : smoothing
         - rs : remeshing

         HSF-Labels:
         - several numerics corresponding to the hippocampal subfield
           segmentation look-up table

         Suffixes:
         - assigned : molecular layer split into subregions
         - merged   : molecular layer subregions merged with other regions

      3. A 'labels' folder with intermediate results, primarily label files.

      4. A 'merge-ml' folder with intermediate results, primarily files created
         during the merging of the molecular layer.

      5. A 'mask' folder with intermediate results, primarily files created
         during creation and post-processing of binary masks.

      6. A 'fixed-surface' folder with intermediate results, primarily files
         created during topological fixing (optional).

      7. A 'tetra-labels' folder with intermediate results, primarily files
         created during tetrahedral mesh construction.

      8. A 'tetra-cut' folder with intermediate results, primarily files
         created during mesh cutting.

      8. A 'tetra-cube' folder with intermediate results, primarily files
         created during cube parametrization.

      10. A 'thickness' folder with thickness overlays and mid-surfaces.


    Custom segmentations:

    - If using the `ashs` segmentation, additional labels for the hippocampal head
      (label value: 20) and tail (label value: 5) labels are required. Use `--lut ashs`.
    - Topological fixing (`--topological-fixing <filename1> <filename2>`) should
      not be used unless working with FreeSurfer-processed data. `<filename1>` is
      the `brain.mgz` file and `<filename2>` is the `wm.mgz` file, both found within
      the `mri` subdirectory of an individual FreeSurfer output directory.
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
# check environment and packages

def _check_environment_and_packages():

    # imports
    import os
    import sys
    import shutil
    import logging
    import importlib.util
    import packaging.version

    # check environment variables
    if os.environ.get('FREESURFER_HOME') is None:
        logging.error('Need to set the FreeSurfer_HOME environment variable')
        print("Program exited with ERRORS.")
        sys.exit(1)

    # check python version
    if sys.version_info <= (3, 5):
        logging.error('Python version must be 3.5 or greater')
        print("Program exited with ERRORS.")
        sys.exit(1)

    # check for gmsh
    if shutil.which("gmsh") is None:
        logging.error('Could not find a \'gmsh\' executable')
        print("Program exited with ERRORS.")
        sys.exit(1)

    # check for LaPy
    if importlib.util.find_spec("lapy") is not None:
        import lapy as lp
        if packaging.version.parse(lp.__version__) < packaging.version.parse("1.0"):
            logging.error('A version >=1.0 is required for the \'lapy\' package (see README.md for details on installation)')
            print("Program exited with ERRORS.")
            sys.exit(1)
    else:
        logging.error('Could not find the \'lapy\' package (see README.md for details on installation)')
        print("Program exited with ERRORS.")
        sys.exit(1)


# ------------------------------------------------------------------------------
# parse_arguments

def _parse_arguments():
    """
    an internal function to parse input arguments

    """

    # imports
    import sys
    import argparse

    # setup parser
    parser = argparse.ArgumentParser(description="This program conducts a thickness analysis of the hippocampus, based on a FreeSurfer, ASHS, or custom hippocampal subfield segmentation.",
        add_help=False)

    # required arguments
    required = parser.add_argument_group('Required arguments')

    required.add_argument('--filename', dest="filename", help="Filename of a segmentation file.",
        default=None, metavar="<filename>", required=False)
    required.add_argument('--hemi', dest="hemi", help="Specify hemisphere, either \'lh\' or \'rh\'.",
        default=None, metavar="<lh|rh>", required=False)
    required.add_argument('--lut', dest="lut", help="Look-up table: a text file with numeric and verbal segmentation labels. \'fs711\' and \'ashs\' are keywords for built-in tables.",
        default=None, metavar="<fs711|ashs|filename>", required=False)

    # optional arguments
    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument('--outputdir', dest="outputdir", help="Directory where the results will be written. If not given, a subfolder within each subject's directory will be created.",
        default=None, metavar="<directory>", required=False)
    optional.add_argument('--no-cleanup', dest='cleanup', help="Keep files that may be useful for diagnostic or debugging purposes, but are not strictly necessary otherwise.",
        default=False, action="store_true", required=False)

    # expert options
    expert = parser.add_argument_group('Expert options')

    expert.add_argument('--no-crop', dest='nocrop', help=argparse.SUPPRESS,
        default=False, action="store_true", required=False) # help="Do not crop image.",
    expert.add_argument('--upsample', dest='upsample', help=argparse.SUPPRESS,
        default=None, metavar="<factor x> <factor y> <factor z>", nargs='*', required=False, type=float) #  help="A list of parameters to fine-tune image upsampling.",
    expert.add_argument('--automated-head', dest='bndHead', help=argparse.SUPPRESS,
        default=False, action="store_true", required=False) # help="",
    expert.add_argument('--automated-tail', dest='bndTail', help=argparse.SUPPRESS,
        default=False, action="store_true", required=False) # help="",
    expert.add_argument('--margin-head', dest="bndHeadMargin", help=argparse.SUPPRESS,
        default=[0], metavar="<params>", nargs=1, required=False)
    expert.add_argument('--margin-tail', dest="bndTailMargin", help=argparse.SUPPRESS,
        default=[0], metavar="<params>", nargs=1, required=False)
    expert.add_argument('--dilation-erosion', dest="dilero", help=argparse.SUPPRESS,
        default=None, metavar="<params>", nargs='+', required=False) # help="A list of parameters to fine-tune the surface creation.",
    expert.add_argument('--no-filter', dest='nofilter', help=argparse.SUPPRESS,
        default=False, action="store_true", required=False) # help="Do not filter image.",
    expert.add_argument('--filter', dest='filter', help=argparse.SUPPRESS,
        default=[0.5, 50], metavar="<kernel width> <threshold>", nargs=2, required=False, type=float) #  help="A list of parameters to fine-tune image filtering.",
    expert.add_argument('--long-filter', dest='longfilter', help=argparse.SUPPRESS,
        default=False, action="store_true", required=False) # help="Filter image along longitudinal axis.",
    expert.add_argument('--long-filter-size', dest='longfilter_size', help=argparse.SUPPRESS,
        default=[5], metavar="<params>", nargs=1, required=False)
    expert.add_argument('--close-mask', dest='closemask', help=argparse.SUPPRESS,
        default=None, metavar="regular|experimental", nargs="*", required=False) # help="Apply closing opreation to mask.",
    expert.add_argument('--mca', dest="mca", help=argparse.SUPPRESS,
        default=None, metavar="<params>", required=False) # help="Marching-cube algorithm.",
    expert.add_argument('--mcc', dest="mcc", help=argparse.SUPPRESS,
        default=None, metavar="<params>", required=False) # help="Marching-cube connectivity.",
    optional.add_argument('--topological-fixing', dest='tfx', help=argparse.SUPPRESS,
        default=None, metavar="<filename>", nargs=2, required=False) # help="Enable topological fixing algorithm for surfaces. Can only be used with FreeSurfer-processed data. Expects two files as input, brain.mgz and wm.mgz, from FreeSurfer\'s `mri` subdirectory.",
    expert.add_argument('--smooth', dest="smooth", help=argparse.SUPPRESS,
        default=None, metavar="<params>", required=False) # help="Mesh smoothing iterations.",
    expert.add_argument('--remesh', dest='remesh', help=argparse.SUPPRESS,
        default=None, metavar="<params>", nargs='*', required=False) # help="Switch on remeshing.",
    expert.add_argument('--no-check-surface', dest='nochecksurface', help=argparse.SUPPRESS,
        default=False, action="store_true", required=False) # help="Do not check surface.",
    expert.add_argument('--no-check-boundaries', dest='nocheckboundaries', help=argparse.SUPPRESS,
        default=False, action="store_true", required=False) # help="Do not boundaries.",
    expert.add_argument('--no-qc', dest='noqc', help=argparse.SUPPRESS,
        default=False, action="store_true", required=False) # help="Skip QC",
    expert.add_argument('--cut-params', dest="cutrange", help=argparse.SUPPRESS,
        default=None, metavar="<params>", nargs='+', required=False, type=float) #  help="A list of parameters to fine-tune the cut operation.",
    expert.add_argument('--aniso-alpha', dest="anisoAlpha", help=argparse.SUPPRESS,
        default=[40], metavar="<params>", nargs='+', required=False) # help="Anisotropy parameter.",
    expert.add_argument('--aniso-smooth', dest="anisoSmooth", help=argparse.SUPPRESS,
        default=None, metavar="<params>", required=False) # help="Anisotropy smoothing.",
    expert.add_argument('--thickness-params', dest="thxyz", help=argparse.SUPPRESS,
        default=None, metavar="<params>", nargs='+', required=False, type=float) #  help="A list of parameters to fine-tune the thickness computation.",
    expert.add_argument('--allow-ragged', dest='allowRagged', help=argparse.SUPPRESS,
        default=False, action="store_true", required=False) # help="Allow ragged mid-surfaces.",
    expert.add_argument('--allow-ragged-triangles', dest='allowRaggedTriangles', help=argparse.SUPPRESS,
        default=False, action="store_true", required=False) # help="Allow triangles for ragged mid-surfaces.",
    expert.add_argument('--logfiledir', dest='logfiledir', help=argparse.SUPPRESS,
        default=None, metavar="<directory>", required=False) # help="Where to store temporary logfile. Default: current working directory.",

    # Deprecated options
    expert.add_argument('--spherically-project', dest='spherically_project',  help=argparse.SUPPRESS,
        default=False, action="store_true", required=False) # help="Use eigenvalue-based spherical projection for topology fixing."
    expert.add_argument('--skip-orient', dest='skiporient', help=argparse.SUPPRESS,
        default=False, action="store_true", required=False) # help="Do not check surface.",

    # define help
    help = parser.add_argument_group('Getting help')

    help.add_argument('--help', help="Display this help message and exit", action='help')
    help.add_argument('--more-help', dest='more_help', help="Display extensive help message and exit",
        default=False, action="store_true", required=False)
    help.add_argument('--version', help="Display version number and exit", action='version', version='%(prog)s '+get_version())

    # check if there are any inputs; if not, print help and exit
    if len(sys.argv) == 1:
        args = parser.parse_args(['--help'])
    else:
        args = parser.parse_args()

    # return extensive helptext
    if args.more_help is True:
        get_help()
        sys.exit(0)

    # check for required arguments print help and exit
    if args.filename is None:
        print("ERROR: the --filename argument is required, exiting. Use --help to see details.")
        sys.exit(1)

    if args.hemi is None:
        print("ERROR: the --hemi argument is required, exiting. Use --help to see details.")
        sys.exit(1)

    if args.lut is None:
        print("ERROR: the --lut argument is required, exiting. Use --help to see details.")
        sys.exit(1)

    #
    return args


# ------------------------------------------------------------------------------
# evaluate arguments

def _evaluate_args(args):
    """
    an internal function to set and return internal parameters

    """

    # import

    import os
    import sys
    import logging

    # message

    logging.info("Evaluating arguments ...")

    # initialize settings

    class settings:
        pass

    settings.HSFLABEL_00 = 'hsf'

    # filename

    FILENAME = os.path.abspath(args.filename)

    logging.info("Using " + FILENAME + " as input file.")

    # output directory

    if args.outputdir is None:
        OUTDIR = os.path.join(os.path.dirname(os.path.abspath(args.filename)), "hsf-" + args.hemi)
    else:
        OUTDIR = args.outputdir

    logging.info("Using " + OUTDIR + " as output directory")

    # preprocessImage

    settings.CROP = not(args.nocrop)

    # processMask

    if args.upsample is None:
        settings.UPSAMPLE = None
    else:
        settings.UPSAMPLE = args.upsample

    if args.dilero is None:
        settings.DIL = 0 # was 2
        settings.ERO = 0 # was 2
    else:
        settings.DIL = int(args.dilero[0])
        settings.ERO = int(args.dilero[1])

    if args.nofilter is True:
        settings.FILTERMASK = None
    else:
        settings.FILTERMASK = args.filter

    if args.longfilter is False:
        settings.LONGFILTER = None
    else:
        settings.LONGFILTER = args.longfilter

    settings.LONGFILTER_SIZE = (int(args.longfilter_size[0]), int(args.longfilter_size[0]), int(args.longfilter_size[0]))

    if args.closemask is None:
        settings.CLOSEMASK = args.closemask
    else:
        if len(args.closemask) == 0:
            settings.CLOSEMASK = "regular"
        elif len(args.closemask) == 1:
            if args.closemask[0] != "regular" and args.closemask[0] != "experimental" and args.closemask[0] != "experimental_v2":
                logging.error("The --close-mask argument can be either 'regular' or 'experimental'")
                print("Program exited with ERRORS.")
                sys.exit(1)
            else:
                settings.CLOSEMASK = args.closemask[0]
        else:
            logging.error("The --close-mask argument can take only 0 or 1 options, but " + str(len(args.closemask)) + " were specified")
            print("Program exited with ERRORS.")
            sys.exit(1)

    if args.nochecksurface is True:
        settings.CHECKSURFACE = None
    else:
        settings.CHECKSURFACE = args.nochecksurface

    if args.nocheckboundaries is True:
        settings.CHECKBOUNDARIES = None
    else:
        settings.CHECKBOUNDARIES = args.nocheckboundaries

    settings.BNDHEAD = args.bndHead
    settings.BNDTAIL = args.bndTail
    settings.BNDHEADMARGIN = int(args.bndHeadMargin[0])
    settings.BNDTAILMARGIN = int(args.bndTailMargin[0])

    # marching cube algorithm

    if args.mca is None:
        settings.MCA = "mri_mc" # was mri_tesselate
    else:
        settings.MCA = args.mca

    # marching cube connectivity; 1=6+, 2=18, 3=6, 4=26

    if args.mcc is None:
        settings.MCC = 1 # was 3
    else:
        settings.MCC = int(args.mcc)

    # smoothing iterations

    if args.smooth is None:
        settings.SMO = 5
    else:
        settings.SMO = int(args.smooth)

    # remeshing

    if args.remesh is None:
        settings.REMESH = None
    elif not args.remesh:
        settings.REMESH = 0
    else:
        settings.REMESH = int(args.remesh[0])

    # topology fixing

    settings.tfx = args.tfx

    settings.spherically_project = args.spherically_project

    if args.tfx is None:
        args.skipTFX = True
    else:
        args.skipTFX = False

    # qc

    settings.noqc = args.noqc

    # anisotropy parameters

    if args.anisoAlpha is None:
        settings.anisoAlpha = None
    else:
        if len(args.anisoAlpha)==1:
            settings.anisoAlpha = [ float(x) for x in args.anisoAlpha ][0]
        elif len(args.anisoAlpha)==2:
            settings.anisoAlpha = [ float(x) for x in args.anisoAlpha ]
        else:
            logging.error("Length of --aniso-alpha must be 1 or 2.")
            print("Program exited with ERRORS.")
            sys.exit(1)

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

    # assemble and return params

    class params:
        pass

    params.FILENAME = FILENAME
    params.OUTDIR = OUTDIR

    params.HEMI = args.hemi
    params.LUT = args.lut
    params.skipTFX = args.skipTFX
    params.skipCLEANUP = args.cleanup
    params.LOGFILE = args.logfile

    params.internal = settings

    return params


# ------------------------------------------------------------------------------
# check arguments

def _check_arguments(params):
    """
    an internal function to check input arguments

    """

    # --------------------------------------------------------------------------
    # imports

    import os
    import sys
    import errno
    import pandas
    import logging
    import tempfile
    import importlib

    import  numpy as np

    # --------------------------------------------------------------------------
    # check arguments

    # check if output directory exists or can be created and is writable

    if os.path.isdir(params.OUTDIR):
        logging.info("Found output directory " + params.OUTDIR)
    else:
        try:
            os.mkdir(params.OUTDIR)
        except:
            logging.error('Cannot create output directory ' + params.OUTDIR)
            print("Program exited with ERRORS.")
            sys.exit(1)
        try:
            testfile = tempfile.TemporaryFile(dir=params.OUTDIR)
            testfile.close()
        except OSError as e:
            if e.errno != errno.EACCES:  # 13
                e.filename = params.OUTDIR
                raise
            logging.error('Directory ' + params.OUTDIR + ' not writeable')
            print("Program exited with ERRORS.")
            sys.exit(1)

    # check for subfield segmentation file

    if os.path.isfile(params.FILENAME):
        logging.info("Found " + params.FILENAME)
    else:
        logging.error("Could not find " + params.FILENAME)
        print("Program exited with ERRORS.")
        sys.exit(1)

    # check for topological fixing files if given

    if params.internal.tfx is not None:
        if len(params.internal.tfx) != 2:
            logging.error('Topology-fixing requires two input files.')
            print("Program exited with ERRORS.")
            sys.exit(1)
        else:
            if params.internal.tfx[0] == "custom_brain":
                logging.info("Using custom brain file")
            else:
                if os.path.isfile(params.internal.tfx[0]):
                    logging.info("Found " + params.internal.tfx[0])
                else:
                    logging.error("Could not find " + params.internal.tfx[0])
                    print("Program exited with ERRORS.")
                    sys.exit(1)
            if params.internal.tfx[1] == "custom_wm":
                logging.info("Using custom wm file")
            else:
                if os.path.isfile(params.internal.tfx[1]):
                    logging.info("Found " + params.internal.tfx[1])
                else:
                    logging.error("Could not find " + params.internal.tfx[1])
                    print("Program exited with ERRORS.")
                    sys.exit(1)

    # check MC algorithm

    if params.internal.MCA != "mri_mc" and params.internal.MCA != "mri_tessellate":
        print("Could not recognise algorithm " + params.internal.MCA + ", exiting.")
        sys.exit(1)

    # check upsampling

    if params.internal.UPSAMPLE is not None:
        if len(params.internal.UPSAMPLE ) != 0 and len(params.internal.UPSAMPLE ) != 3:
            logging.error("Incorrect number of --upsampling parameters, must be 0 or 3.")
            print("Program exited with ERRORS.")
            sys.exit(1)

    # create the LUT

    if params.LUT == "fs711":

        logging.info("Found internal, modified look-up table for FreeSurfer 7.11.")

        LUTLABEL = [ "tail", "head", "presubiculum", "subiculum", "ca1", "ca2",
            "ca3", "ca4", "dg", "ml" ]

        LUTINDEX = [ 226, [233, 235, 237, 239, 245], 234, 236, 238, 240, 240,
            242, 244, 246 ]

        lutDict = dict(zip(LUTLABEL, LUTINDEX))

        hsflist = [ 234, 236, 238, 240, 246 ]

    elif params.LUT == "fs711-noml":

        logging.info("Found internal, modified look-up table for FreeSurfer 7.11.")

        LUTLABEL = [ "tail", "head", "presubiculum", "subiculum", "ca1", "ca2",
            "ca3", "ca4", "dg" ]

        LUTINDEX = [ 226, [233, 235, 237, 239, 245], 234, 236, 238, 240, 240,
            242, 244 ]

        lutDict = dict(zip(LUTLABEL, LUTINDEX))

        hsflist = [ 234, 236, 238, 240]

    elif params.LUT == "ashs":

        logging.info("Found internal, modified look-up table for ASHS atlas.")

        LUTLABEL = [ "ca1", "ca2", "ca3", "ca4", "dg", "tail_orig", "subiculum",
            "presubiculum", "entorhinal", "ba35", "ba36", "parahippocampal",
            "head", "tail" ]

        LUTINDEX = [ 1, 2, 4, 3, 3, 5, 8, 8, 9, 10, 11, 12, 20, 5 ]

        lutDict = dict(zip(LUTLABEL, LUTINDEX))

        hsflist = [ 8, 1, 2, 4 ]

    elif os.path.isfile(params.LUT):

        try:

            logging.info("Found look-up table " + params.LUT)

            lut = pandas.read_csv(params.LUT, sep=' ', comment='#',
                header=None, skipinitialspace=True, skip_blank_lines=True,
                error_bad_lines=False, warn_bad_lines=True)

            lutDict = dict(zip(lut[0], lut[1]))
            hsflist = list(lut[1])

        except:

            logging.error("Could not read look-up table " + params.LUT)
            print("Program exited with ERRORS.")
            sys.exit(1)

    else:

        logging.error("Look-up table can only be \'fs711\', \'ashs\', or an existing file, but not " + params.LUT)
        print("Program exited with ERRORS.")
        sys.exit(1)

    params.LUTDICT = lutDict
    params.HSFLIST = hsflist

    # --------------------------------------------------------------------------
    # create directories if they don't exist

    # define list of directories

    list_of_directories = [
        os.path.join(params.OUTDIR, "labels"),
        os.path.join(params.OUTDIR, "merge-ml"),
        os.path.join(params.OUTDIR, "mask"),
        os.path.join(params.OUTDIR, "fixed-surface"),
        os.path.join(params.OUTDIR, "tetra-labels"),
        os.path.join(params.OUTDIR, "tetra-cut"),
        os.path.join(params.OUTDIR, "tetra-cube"),
        os.path.join(params.OUTDIR, "thickness"),
        os.path.join(params.OUTDIR, "qc")
        ]

    # loop over list of directories

    for directory in list_of_directories:

        if not os.path.isdir(directory):

            try:
                logging.info('Creating output directory ' + directory)
                os.mkdir(directory)
            except:
                logging.error('Cannot create output directory ' + directory)
                print("Program exited with ERRORS.")
                sys.exit(1)

    # return

    return params


# ------------------------------------------------------------------------------
# run analysis

def _run_analysis(params):
    """

    """

    # imports

    import os
    import time
    import logging

    from shapetools.preprocessImage import convertFormat
    from shapetools.preprocessImage import cropImage
    from shapetools.preprocessImage import upsampleImage
    from shapetools.preprocessImage import autoMask
    from shapetools.createLabels import createLabels
    from shapetools.mergeMolecularLayer import mergeMolecularLayer
    from shapetools.processMask import fillHoles
    from shapetools.processMask import createMask
    from shapetools.processMask import filterMask
    from shapetools.processMask import longFilter
    from shapetools.processMask import closeMask
    from shapetools.checkSurface import checkSurface
    from shapetools.createSurface import createSurface
    from shapetools.createTetraMesh import createTetraMesh
    from shapetools.createTetraLabels import createTetraLabels
    from shapetools.removeBoundaryMask import removeBoundaryMask
    from shapetools.cutTetra import cutTetra
    from shapetools.computeCubeParam import computeCubeParam
    from shapetools.computeThickness import computeThickness
    from shapetools.mapValues import mapValues
    from shapetools.createSupplementaryFiles import createSupplementaryFiles
    from shapetools.qcPlots import qcPlots

    # convert format (0)

    logging.info("Starting convertFormat() ...")
    params = convertFormat(params)

    logging.info("Starting cropImage() ...")
    params = cropImage(params)

    logging.info("Starting upsampleImage() ...")
    params = upsampleImage(params)

    logging.info("Starting autoMask() ...")
    params = autoMask(params)

    # create labels (1)

    logging.info("Starting createLabels() ...")
    params = createLabels(params)

    # merge molecular layer (2)

    logging.info("Starting mergeMolecularLayer() ...")
    params = mergeMolecularLayer(params)

    # create mask (3)

    logging.info("Starting fillHoles() ...")
    params = fillHoles(params)

    logging.info("Starting createMask() ...")
    params = createMask(params)

    logging.info("Starting longFilter() ...")
    params = longFilter(params)

    logging.info("Starting closeMask() ...")
    params = closeMask(params)

    logging.info("Starting filterMask() ...")
    params = filterMask(params)

    # create surface (4)

    logging.info("Starting createSurface() ...")
    params = createSurface(params)

    logging.info("Starting qcPlots() ...")
    params = qcPlots(params, stage="mesh")

    logging.info("Starting checkSurface() ...")
    cont, params = checkSurface(params, stage="check_surface")

    #

    if cont is True:

        # create tetra mesh for whole hippocampal body (5)

        logging.info("Starting createTetraMesh() ...")
        params = createTetraMesh(params)

        # create label files for tetra mesh (6)

        logging.info("Starting createTetraLabels() ...")
        params = createTetraLabels(params)

        # remove boundary mask (7)

        logging.info("Starting removeBoundaryMask() ...")
        params = removeBoundaryMask(params)

        # create tetra mesh for cut hippocampal body (8)

        logging.info("Starting cutTetra() ...")
        params = cutTetra(params)


        logging.info("Starting checkSurface() ...")
        cont, params = checkSurface(params, stage="check_boundaries")

        #

        if cont is True:

            # compute cube parametrization (9)

            logging.info("Starting computeCubeParam() ...")
            params = computeCubeParam(params)

            logging.info("Starting qcPlots() ...")
            params = qcPlots(params, stage="profile")

            # compute thickness (10)

            logging.info("Starting computeThickness() ...")
            params = computeThickness(params)

            logging.info("Starting qcPlots() ...")
            params = qcPlots(params, stage="hull")

            # map subfield mapValues (11)

            logging.info("Starting mapValues() ...")
            params = mapValues(params)

            # map subfield mapValues (12)

            logging.info("Starting createSupplementaryFiles() ...")
            params = createSupplementaryFiles(params)

    # all done

    logging.info("Date: %s", time.strftime('%d/%m/%Y %H:%M:%S'))
    if cont is True:
        logging.info("Hippocampal shapetools finished without errors.")
    else:
        logging.info("Hippocampal shapetools finished WITH ERRORS.")

# ------------------------------------------------------------------------------
# start logging

def _start_logging(args):
    """
    start logging

    """

    # imports
    import os
    import sys
    import uuid
    import time
    import logging
    import tempfile
    import traceback

    # setup function to log uncaught exceptions
    def foo(exctype, value, tb):
        # log
        logging.error('Error Information:')
        logging.error('Type: %s', exctype)
        logging.error('Value: %s', value)
        for i in traceback.format_list(traceback.extract_tb(tb)):
            logging.error('Traceback: %s', i)
        # message
        print('Program exited with ERRORS.')
        print('')
        sys.exit(1)
    sys.excepthook = foo

    # check if logfile can be written in logfiledir
    if args.logfiledir is not None:
        logfiledir = args.logfiledir
    else:
        logfiledir = os.getcwd()

    try:
        testfile = tempfile.TemporaryFile(dir=logfiledir)
        testfile.close()
    except OSError as e:
        print('Directory ' + logfiledir + ' not writeable for temporary logfile.')
        print("Program exited with ERRORS.")
        sys.exit(1)

    # start logging
    logfile =  os.path.join(logfiledir, 'logfile-' + str(uuid.uuid4()) + '.log')
    logging.basicConfig(filename=logfile, filemode='w', level=logging.INFO)
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout)) # this adds output to stdout

    # intial messages
    logging.info("Starting logging for hippocampal shapetools ...")
    logging.info("Version: %s", get_version())
    logging.info("Date: %s", time.strftime('%d/%m/%Y %H:%M:%S'))

    # log args
    logging.info("Command: " + " ".join(sys.argv))

    # update args
    args.logfile = logfile

    # return
    return args


# ------------------------------------------------------------------------------
# run analysis

def _stop_logging(params):
    """
    stop logging

    """

    # imports
    import os

    # move logfile to output dir
    os.replace(params.LOGFILE, os.path.join(params.OUTDIR, "logfile.txt"))


# ------------------------------------------------------------------------------
# run analysis

def run_analysis(args):
    """
    a function to run the shapetools submodules

    """

    # start logging
    args = _start_logging(args)

    # check environment and packages
    _check_environment_and_packages()

    # evaluate arguments
    params = _evaluate_args(args)

    # check arguments
    _check_arguments(params)

    # run analysis
    _run_analysis(params)

    # stop logging
    _stop_logging(params)


# ------------------------------------------------------------------------------
# main

if __name__ == "__main__":

    # message

    print('')
    print('----------------------------------------')
    print('Hippocampal shape and thickness analysis')
    print('----------------------------------------')
    print('')

    # parse arguments

    args = _parse_arguments()

    # run analysis

    run_analysis(args)
