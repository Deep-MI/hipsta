"""
This module provides the main functionality of the hippocampal shape tools
package.

"""


# ==============================================================================
# FUNCTIONS

# ------------------------------------------------------------------------------
# get_version()

def get_version():

    VERSION = "v0.1.0-beta"

    return VERSION

# ------------------------------------------------------------------------------
# get_help()

def get_help(print_help=True, return_help=False):
    """
    a function to return a help message

    """

    HELPTEXT = """

    Hippocampal shape tools


    Purpose:

    This is a script for the creation and analysis of hippocampal surfaces.


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


    Usage:

    python3 shapetools.py --subjid <SUBJID> --hemi <hemi> --lut <lookup-table>
        [--subjdir <SUBJDIR>] [--outputdir <OUTDIR>] [further options]


    Input arguments:

      --subjid    Subject ID. Must match a subdirectory in Freesurfer's
                  SUBJECTS_DIR or in a manually specified --subjdir. The
                  expected filename is <subjid>/mri/<lh|rh>.<suffix>.mgz.
      --hemi      Hemisphere. Either 'lh' or 'rh'.
      --lut       A valid look-up table for hippocampal subfield segmentation.
                  Either 'fs711' or 'ashs' or a valid filename.


    Optional input arguments:

      --subjdir   Subjects directory. This directory must contain a folder that
                  matches the subject ID and contains Freesurfer output
                  including the hippocampal subfield segmentation. If not
                  given, Freesurfer's SUBJECTS_DIR environment variable will be
                  used.

      --outputdir Directory where the results will be written. If not given, a
                  subfolder within each subject's directory will be created.

      --suffix    Subfield segmentation suffix. These are different for the
                  label files created with different versions and settings of
                  the subfield segmentation (default: hippoAmygLabels-T1.v20).

      --no-topological-fixing
                  Disable topological fixing. This is a simple switch and no
                  further parameter is required; if it is given, use
                  Freesurfer's topology fixing program to refine and fix
                  initial surfaces.

      --no-cleanup
                  Do not delete files that may be useful for diagnostic or
                  debugging purposes, but are not strictly necessary otherwise.

      --help      Display this help and exit.

      --version   Display version number and exit.


    Example:

    python3 shapetools.py --subjid my_subject --hemi lh --lut fs711
        --subjdir /my/directory --outputdir /my/output/directory
        --suffix hippoAmygLabels-T1.v20


    Outputs:

      1. The primary output are thickness values that will be computed and
         stored as PSOL overlay files within the 'tickness' folder. These can
         id onto the mid-surface vtk file that is also found within the
         'thickness' folder.

      2. Intermediate volume and surface files, each prefixed with the
         particular operation performed. The basic pattern for the filenames is:

         <hemisphere>.<prefixes>.<hsf-labels>.<suffixes>

         Hemisphere:
         - lh : left hemisphere
         - rh : right hemisphere

         Prefixes:
         - ml : merged molecular layer
         - de : dilation and erosion
         - mc : marching cube algorithm
         - tf : topological fixing
         - rs : remeshing and smoothing

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

      9. A 'thickness' folder with thickness overlays and mid-surfaces.


    Requirements:

      1. A Freesurfer version (6.x or 7.x) must be sourced, i.e.
         FREESURFER_HOME must exist as an environment variable and point to a
         valid Freesurfer installation.

      2. A hippocampal subfield segmentation created by Freesurfer 7.11
         (development versions after 6.0 will also work) must be present in the
         'mri' subdirectory of the Freesurfer 6.0 or newer processed
         SUBJECTS_DIR/<subjid> directory.

      3. Python 3.5 or higher including the lapy, numpy, scipy, and nibabel
         libraries.

      4. An installation of the gmsh package (verson 2.x; http://gmsh.info)
         must be on the PATH. E.g., http://gmsh.info/bin/Linux/older/gmsh-
         2.16.0-Linux64.tgz or http://gmsh.info/bin/MacOSX/older/gmsh-2.16.0-
         MacOSX.dmg


    References:

    Please cite the following publications if using this script in your work:

      M. Reuter, F.-E. Wolter, M. Shenton, M. Niethammer. Laplace-Beltrami
      Eigenvalues and Topological Features of Eigenfunctions for Statistical
      Shape Analysis. Computer-Aided Design 41 (10), pp.739-755, 2009.
      http://dx.doi.org/10.1016/j.cad.2009.02.007

      C. Geuzaine and J.-F. Remacle. Gmsh: a three-dimensional finite element
      mesh generator with built-in pre- and post-processing facilities.
      International Journal for Numerical Methods in Engineering 79(11), pp.
      1309-1331, 2009.

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
    import importlib
    import subprocess

    # check environment variables
    if os.environ.get('FREESURFER_HOME') is None:
        print('\nERROR: need to set the FREESURFER_HOME environment variable\n')
        sys.exit(1)

    if os.environ.get('SUBJECTS_DIR') is None:
        print('\nERROR: need to set the SUBJECTS_DIR environment variable\n')
        sys.exit(1)

    # check python version
    if sys.version_info <= (3, 5):
        print('\nERROR: Python version must be 3.5 or greater\n')
        sys.exit(1)

    # check for gmsh
    if shutil.which("gmsh") is None:
        print('\nERROR: could not find a \'gmsh\' executable\n')
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
    parser = argparse.ArgumentParser(description="This program takes returns a thickness analysis of the hippocampus, based on a Freesurfer 6.0 hippocampal subfield segmentation.",
        add_help=False)

    # required arguments
    required = parser.add_argument_group('required arguments')

    required.add_argument('--subjid', dest="subjid", help="Subject ID. Must match a subdirectory in Freesurfer's SUBJECTS_DIR or in a manually specified --subjdir.",
        metavar="<directory>", required=True)
    required.add_argument('--hemi', dest="hemi", help="Hemisphere.",
        metavar="<lh|rh>", required=True)
    required.add_argument('--lut', dest="lut", help="Look-up table.",
        metavar="<fs711|ashs|filename>", required=True)

    # optional arguments
    optional = parser.add_argument_group('optional arguments')

    optional.add_argument('--subjdir', dest="subjdir", help="Subjects directory. This directory must contain a folder that matches the subject ID and contains Freesurfer output including the hippocampal subfield segmentation. If not given, Freesurfer's SUBJECTS_DIR environment variable will be used.",
        default=None, metavar="<directory>", required=False)
    optional.add_argument('--outputdir', dest="outputdir", help="Directory where the results will be written. If not given, a subfolder within each subject's directory will be created.",
        default=None, metavar="<directory>", required=False)
    optional.add_argument('--suffix', dest="sfx", help="Subfield segmentation suffix. These are different for the label files created with different versions and settings of the subfield segmentation (default: hippoAmygLabels-T1.v20).",
        default="hippoAmygLabels-T1.v20", metavar="<string>", required=False)
    optional.add_argument('--no-topological-fixing', dest='tfx', help="Disable topological fixing. This is a simple switch and no further parameter is required; if it is given, don't use Freesurfer's topology fixing program to refine and fix initial surfaces.",
        default=False, action="store_true", required=False)
    optional.add_argument('--no-cleanup', dest='cleanup', help="Keep files that may be useful for diagnostic or debugging purposes, but are not strictly necessary otherwise.",
        default=False, action="store_true", required=False)

    # expert options
    expert = parser.add_argument_group('expert options')

    expert.add_argument('--thickness-params', dest="thxyz", help="A list of parameters to fine-tune the thickness computation.",
        default=None, metavar="<params>", nargs='+', required=False, type=float)
    expert.add_argument('--dilation-erosion', dest="dilero", help="A list of parameters to fine-tune the surface creation.",
        default=None, metavar="<params>", nargs='+', required=False)
    expert.add_argument('--mcc', dest="mcc", help="Marching-cube connectivity.",
        default=None, metavar="<params>", required=False)
    expert.add_argument('--smooth', dest="smooth", help="Mesh smoothing iterations.",
        default=None, metavar="<params>", required=False)
    expert.add_argument('--aniso-alpha', dest="anisoAlpha", help="Anisotropy parameter.",
        default=40, metavar="<params>", required=False)
    expert.add_argument('--aniso-smooth', dest="anisoSmooth", help="Anisotropy smoothing.",
        default=None, metavar="<params>", required=False)
    expert.add_argument('--mean-curv', dest="weightMeanCurv", help="Mean curvature weighting.",
        default=None, metavar="<params>", required=False)
    expert.add_argument('--allow-ragged', dest='allowRagged', help="Allow ragged mid-surfaces.",
        default=False, action="store_true", required=False)
    expert.add_argument('--allow-ragged-triangles', dest='allowRaggedTriangles', help="Allow triangles for ragged mid-surfaces.",
        default=False, action="store_true", required=False)

    # define help
    help = parser.add_argument_group('getting help')

    help.add_argument('--help', help="display this help message and exit", action='help')
    help.add_argument('--version', help="display version number and exit", action='version', version='%(prog)s '+get_version())

    # check if there are any inputs; if not, print help and exit
    if len(sys.argv) == 1:
        args = parser.parse_args(['--help', '--version'])
    else:
        args = parser.parse_args()

    return args


# ------------------------------------------------------------------------------
# evaluate arguments

def _evaluate_args(args):
    """
    an internal function to set and return internal parameters

    """

    # import

    import os

    # initialize settings

    class settings:
        pass

    settings.HSFLABEL_00 = 'hsf'

    # create the internal LUT

    if args.lut == "fs711":

        print("Using internal, modified look-up table for Freesurfer 7.11")

        LUTLABEL = [ "tail", "head", "presubiculum", "subiculum", "CA1", "CA3",
            "CA2", "CA4", "DG", "ML" ]

        LUTINDEX = [ "226", "232", "234", "236", "238", "240", "240", "242",
            "244", "246" ]

        lutDict = dict(zip(LUTLABEL, LUTINDEX))

        hsflist = [ "234", "236", "238", "240", "246" ]

    elif args.lut == "ashs":

        print("Using internal, modified look-up table for ASHS IKND Magdeburg Young Adult 7T Atlas:")

        LUTLABEL = [ "CA1", "CA2", "DG", "CA3", "tail_orig", "subiculum",
            "presubiculum", "entorhinal", "BA35", "BA36", "parahippocampal",
            "head", "tail" ]

        LUTINDEX = [ "1", "2", "3", "4", "5", "8", "8", "9", "10", "11", "12",
            "255", "254" ]

        lutDict = dict(zip(LUTLABEL, LUTINDEX))

        hsflist = [ "8", "4", "2", "1" ]

    else:

        lutDict = None
        hsflist = None

    # subjects dir

    if args.subjdir is None:
        SUBJDIR = os.environ.get('SUBJECTS_DIR')
    else:
        SUBJDIR = args.subjdir
    print("Using", SUBJDIR, "as subjects directory")

    # output directory

    if args.outputdir is None:
        OUTDIR = os.path.join(SUBJDIR, args.subjid, "hsf-" + args.hemi)
    else:
        OUTDIR = args.outputdir
    print("Using", OUTDIR, "as output directory")

    # filename

    FILENAME = os.path.join(SUBJDIR, args.subjid, "mri", args.hemi + "." + args.sfx + ".mgz")
    print("Using", FILENAME, "as input file.")

    # number of voxels for dilation and erosion
    if args.dilero is None:
        settings.DIL = 2
        settings.ERO = 2
    else:
        settings.DIL = int(args.dilero[0])
        settings.ERO = int(args.dilero[1])

    # marching cube connectivity; 1=6+, 2=18, 3=6, 4=26
    if args.mcc is None:
        settings.MCC = 4
    else:
        settings.MCC = int(args.mcc)

    # smoothing iterations
    if args.smooth is None:
        settings.SMO = 5
    else:
        settings.SMO = int(args.smooth)

    # anisotropy parameters
    if args.anisoAlpha is None:
        settings.anisoAlpha = None
    else:
        settings.anisoAlpha = float(args.anisoAlpha)

    if args.anisoSmooth is None:
        settings.anisoSmooth = 3
    else:
        settings.anisoSmooth = int(args.anisoSmooth)

    if args.weightMeanCurv is None:
        settings.weightMeanCurv = 0.0
    else:
        settings.weightMeanCurv = float(args.weightMeanCurv)

    settings.cubeAddNoise = False

    settings.cubeWriteLegacyVTK = True

    # thickness params
    if args.thxyz is None:
        settings.THXn = -0.9
        settings.THXp =  0.9
        settings.THXk =   41
        settings.THYn = -0.9
        settings.THYp =  0.9
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

    # mapValues params

    settings.mapValuesIntegrate = "none"
    settings.mapValuesSelect = None
    settings.mapValuesWritePSOL = True
    settings.mapValuesWriteMGH = True
    settings.mapValuesWriteANNOT = False

    # assemble and return params

    class params:
        pass

    params.SUBJDIR = SUBJDIR
    params.OUTDIR = OUTDIR
    params.FILENAME = FILENAME

    params.SUBJID = args.subjid
    params.HEMI = args.hemi
    params.LUT = args.lut
    params.SFX = args.sfx
    params.skipTFX = args.tfx
    params.skipCLEANUP = args.cleanup

    params.LUTDICT = lutDict
    params.HSFLIST = hsflist

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

    import tempfile
    import importlib

    # --------------------------------------------------------------------------
    # check arguments

    # check if subject directory exists

    if os.path.isdir(params.SUBJDIR):
        print("Found subjects directory", params.SUBJDIR)
    else:
        print('ERROR: subjects directory ' + params.SUBJDIR + ' is not an existing directory\n')
        sys.exit(1)

    # check if output directory exists or can be created and is writable

    if os.path.isdir(params.OUTDIR):
        print("Found output directory", params.OUTDIR)
    else:
        try:
            os.mkdir(params.OUTDIR)
        except:
            print('ERROR: cannot create output directory ' + params.OUTDIR + '\n')
            sys.exit(1)
        try:
            testfile = tempfile.TemporaryFile(dir=params.OUTDIR)
            testfile.close()
        except OSError as e:
            if e.errno != errno.EACCES:  # 13
                e.filename = params.OUTDIR
                raise
            print('\nERROR: ' + params.OUTDIR + ' not writeable (check access)!\n')
            sys.exit(1)

    # check for subfield segmentation file

    if os.path.isfile(params.FILENAME):
        print("Found", params.FILENAME)
    else:
        print("Could not find", params.FILENAME)
        sys.exit(1)

    # check if lookup-table exists if it was given, otherwise exit

    if params.LUTDICT is None:

        print("Lookup table must currently be fs711 or ashs.")
        sys.exit(1)

        # TODO: specify format for input lut and change the following code
        # accordingly; specifically, do not use range() for indices, but actual
        # values!

        if os.path.isfile(params.LUT):

            print("Found look-up table", params.LUT)

            lut = pandas.read_csv(params.LUT, sep=' ', comment='#',
                header=None, skipinitialspace=True, skip_blank_lines=True,
                error_bad_lines=False, warn_bad_lines=True)

            lut = np.array(lut)

            params.LUTDICT = dict(zip(range(len(lut[:, 0])), lut[:, 0]))
            params.HSFLIST = params.LUTDICT.keys()


        else:

            print("ERROR: Could not find look-up table ", params.LUT)
            sys.exit(1)

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
        os.path.join(params.OUTDIR, "thickness")
        ]

    # loop over list of directories

    for directory in list_of_directories:

        if not os.path.isdir(directory):

            try:
                print('Creating output directory ' + directory)
                os.mkdir(directory)
            except:
                print('ERROR: cannot create output directory ' + directory)
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
    import sys

    from shapetools.createLabels import createLabels
    from shapetools.mergeMolecularLayer import mergeMolecularLayer
    from shapetools.createMask import fillHoles
    from shapetools.createMask import createMask
    from shapetools.createSurface import createSurface
    from shapetools.createTetraMesh import createTetraMesh
    from shapetools.createTetraLabels import createTetraLabels
    from shapetools.removeBoundaryMask import removeBoundaryMask
    from shapetools.cutTetra import cutTetra
    from shapetools.computeCubeParam import computeCubeParam
    from shapetools.computeThickness import computeThickness
    from shapetools.mapValues import mapValues

    # create labels (1)

    params = createLabels(params)

    # merge molecular layer (2)

    params = mergeMolecularLayer(params)

    # create mask (3)

    params = fillHoles(params)

    params = createMask(params)

    # create surface (4)

    params = createSurface(params)

    # create tetra mesh for whole hippocampal body (5)

    params = createTetraMesh(params)

    # create label files for tetra mesh (6)

    params = createTetraLabels(params)

    # remove boundary mask (7)

    params = removeBoundaryMask(params)

    # create tetra mesh for cut hippocampal body (8)

    params = cutTetra(params)

    # compute cube parametrization (9)

    params = computeCubeParam(params)

    # compute thickness (10)

    params = computeThickness(params)

    # map subfield mapValues (11)

    params = mapValues(params)


# ------------------------------------------------------------------------------
# run analysis

def run_analysis(args):
    """
    a function to run the shapetools submodules

    """

    # check environment and packages
    _check_environment_and_packages()

    # evaluate arguments
    params = _evaluate_args(args)

    # check arguments
    _check_arguments(params)

    # run analysis
    _run_analysis(params)


# ------------------------------------------------------------------------------
#

if __name__ == "__main__":

    # message

    print('Hippocampal shape tools')

    # parse arguments

    args = _parse_arguments()

    # run analysis

    run_analysis(args)
