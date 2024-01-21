"""
This module provides a function to map values from a volume to the midsurface.

"""

import argparse
import logging
import os
import subprocess
import sys

import nibabel as nb
import numpy as np
import pandas as pd
from lapy import TriaMesh, io
from scipy import stats as st

# ==============================================================================
# LOGGING

LOGGER = logging.getLogger(__name__)

# ==============================================================================
# FUNCTIONS

# ------------------------------------------------------------------------------
# AUXILIARY FUNCTIONS
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# get_help()


def get_help(print_help=True, return_help=False):
    """
    a function to return a help message

    """

    # define helptext
    HELPTEXT = """
    SUMMARY

    This is a script for mapping volume-based data onto the hippocampal
    midsurface.

    USAGE

    The script requires four arguments:

    --volume <input volume>
    --surface <input surface>
    --label <input label file>
    --table <input csv table>
    --suffix <suffix>

    The following arguments are optional and determine which output files will
    be created:

    --integrate <none|mode|median|mean|min|max> (default: none)
    --select <0, 1, ..., n> (default: all)
    --interp <nearest|cubic> (default: nearest)
    --writePSOL
    --writeMGH

    EXAMPLES

    Example for mapping at the vertices of the midsurface:

    python3 mapValues.py --volume /my/path/to/<MY_VOLUME.mgz> \\
        --surface /my/path/to/<SUBJECT_ID>/thickness/lh.mid-surface.vtk \\
        --label /my/path/to/<MY_LABEL_FILE.mgz> \\
        --table /my/path/to/<SUBJECT_ID>/thickness/lh.mid-surface.csv \\
        --suffix MY_SUFFIX --writePSOL --writeMGH

    This will sample the vertices of the midsurface at the corresponding
    positions in the given volume. This mapping requires the '[lr]h.mid-
    surface.vtk' surface and the '[lr]h.mid-surface.csv' index table that are
    created by the hippocampal thickness tools. It also requires the original
    label file that contains the subfield labels. In some cases it will be
    necessary to reslice the input volume to match the input label file: in
    this case, there are two interpolation options, nearest neighbor and cubic.
    Nearest neighbor should be used for labels, and cubic for continuous
    volumetric data. The resliced file will be saved in the same directory as
    the input file. The suffix argument will be used to identify the output
    files and can be chosen arbitrarily. PSOL and MGH output files will be
    written to the same directory as the input surface.

    Example for mapping at the vertices of the streamlines:

    python3 mapValues.py --volume /my/path/to/<MY_VOLUME.mgz> \\
        --surface /my/path/to/<SUBJECT_ID>/thickness/lh.grid-lines-z.vtk \\
        --label /my/path/to/<SUBJECT_ID>/mri/lh.hippoAmygLabels-T1.v20.mgz \\
        --table /my/path/to/<SUBJECT_ID>/thickness/lh.grid-lines.csv \\
        --suffix MY_SUFFIX --integrate mean --select 0 1 2 3 4 5 \\
        --interp cubic --writePSOL --writeMGH

    This will sample the vertices of the interior/exterior streamlines at the
    corresponding positions in the given volume. That is, multiple vertices per
    streamlines will be integrated and projected onto the corresponding point
    of the midsurface. This mapping requires the '[lr]h.grid-lines-z.vtk' file
    and the '[lr]h.grid-lines.csv' index table that is created by the hippo-
    campal thickness tools. Also, the 'integrate' argument is required. In this
    example, the sampled values will be averaged ('mean') per streamline. The
    select argument can be used to restrict the integration to a subset of the
    integration points (0 ... 10 for the default z direction, 5 being the
    midsurface). We assume that our input volume contains continuous volumetric
    data, and hence use cubic interpolation. The suffix will be used to
    identify the output files and can be chosen arbitrarily. PSOL and MGH
    output files will be written to the same directory as the input surface.

    """

    if print_help:
        print(HELPTEXT)

    if return_help:
        return HELPTEXT


# ------------------------------------------------------------------------------
# parse_arguments


def _parse_arguments():
    """
    Command Line Options Parser:
    initiate the option parser and return the parsed object
    """

    # message
    print("\nReading input options ...")

    # setup parser
    parser = argparse.ArgumentParser(description="", add_help=False)

    # help text
    h_IN_VOL = "input volume"
    h_IN_SURF = "input surface"
    h_IN_LABEL = "input label file"
    h_IN_INDICES = "table with indices"
    h_IN_SUFFIX = "suffix for output files"
    h_writePSOL = "write out PSOL files"
    h_writeMGH = "write out MGH files"
    h_writeANNOT = "write out ANNOT files"
    h_integrate = "write out integrated values (default: none)"
    h_select = "select inegration points (default: all)"
    h_interp = "type of interpolation (cubic or nearest; default: nearest)"

    # required arguments
    required = parser.add_argument_group("required arguments")

    required.add_argument(
        "--volume",
        dest="IN_VOL",
        help=h_IN_VOL,
        default=None,
        required=False,
        metavar="<file>",
    )
    required.add_argument(
        "--surface",
        dest="IN_SURF",
        help=h_IN_SURF,
        default=None,
        required=False,
        metavar="<file>",
    )
    required.add_argument(
        "--label",
        dest="IN_LABEL",
        help=h_IN_LABEL,
        default=None,
        required=False,
        metavar="<file>",
    )
    required.add_argument(
        "--table",
        dest="IN_INDICES",
        help=h_IN_INDICES,
        default=None,
        required=False,
        metavar="<file>",
    )
    required.add_argument(
        "--suffix",
        dest="IN_SUFFIX",
        help=h_IN_SUFFIX,
        default=None,
        required=False,
        metavar="<string>",
    )

    # optional arguments
    optional = parser.add_argument_group("optional arguments")

    optional.add_argument(
        "--integrate",
        dest="INTEGRATE",
        help=h_integrate,
        default="none",
        required=False,
        metavar="<none|mode|median|mean|min|max>",
    )
    optional.add_argument(
        "--select",
        dest="SELECT",
        help=h_select,
        default=None,
        nargs="+",
        required=False,
        metavar="<list>",
    )
    optional.add_argument(
        "--interp",
        dest="INTERP",
        help=h_interp,
        default="nearest",
        required=False,
        metavar="<nearest|cubic>",
    )
    optional.add_argument(
        "--writePSOL",
        dest="writePSOL",
        help=h_writePSOL,
        default=False,
        action="store_true",
        required=False,
    )
    optional.add_argument(
        "--writeMGH",
        dest="writeMGH",
        help=h_writeMGH,
        default=False,
        action="store_true",
        required=False,
    )
    optional.add_argument(
        "--writeANNOT",
        dest="writeANNOT",
        help=argparse.SUPPRESS,
        default=False,
        action="store_true",
        required=False,
    )  # help=h_writeANNOT # currently hidden, might be added later

    # define help
    help = parser.add_argument_group("getting help")

    help.add_argument("--help", help="display this help message and exit", action="help")
    help.add_argument(
        "--more-help",
        dest="more_help",
        help="display extensive help message and exit",
        default=False,
        action="store_true",
        required=False,
    )

    # check if there are any inputs; if not, print help and exit
    if len(sys.argv) == 1:
        args = parser.parse_args(["--help"])
    else:
        args = parser.parse_args()

    return args


# ------------------------------------------------------------------------------
# check arguments


def _check_arguments(options):
    """
    an internal function to check input arguments

    """

    # more help

    if options.more_help is True:
        get_help()
        sys.exit(0)

    # check environment variables

    if os.environ.get("FREESURFER_HOME") is None:
        LOGGER.info("ERROR: need to set the FREESURFER_HOME environment variable and to source FREESURFER.")
        sys.exit(1)

    # required arguments

    if options.IN_VOL is None:
        LOGGER.info("The --volume option is required. See --help for details.")
        sys.exit(1)

    if options.IN_SURF is None:
        LOGGER.info("The --surface option is required. See --help for details.")
        sys.exit(1)

    if options.IN_LABEL is None:
        LOGGER.info("The --label option is required. See --help for details.")
        sys.exit(1)

    if options.IN_INDICES is None:
        LOGGER.info("The --table option is required. See --help for details.")
        sys.exit(1)

    if options.IN_SUFFIX is None:
        LOGGER.info("The --suffix option is required. See --help for details.")
        sys.exit(1)

    # change formats

    if options.SELECT is not None:
        options.SELECT = np.array(options.SELECT).astype(int)

    # add params

    options.params = None

    # return

    return options


# ------------------------------------------------------------------------------
# MAIN FUNCTION
# ------------------------------------------------------------------------------


def mapValues(
    params,
    IN_VOL=None,
    IN_SURF=None,
    IN_LABEL=None,
    IN_INDICES=None,
    IN_SUFFIX="hsf",
    INTEGRATE="mode",
    SELECT=None,
    INTERP="nearest",
    writePSOL=False,
    writeMGH=False,
    writeANNOT=False,
):
    # message

    print()
    print("--------------------------------------------------------------------------------")
    print("Mapping values")
    print()

    # -------------------------------------------------------
    # get params

    if params is not None:
        # these are the default settings for processing within the hippocampal
        # thickness toolbox; will override any other settings; set 'params=None'
        # for custom processing.
        IN_VOL = os.path.join(params.OUTDIR, params.HEMI + ".labels.mgz")
        IN_SURF = os.path.join(params.OUTDIR, "thickness", params.HEMI + ".mid-surface.vtk")  # or: grid-lines-z.vtk
        IN_LABEL = params.FILENAME
        IN_INDICES = os.path.join(params.OUTDIR, "thickness", params.HEMI + ".mid-surface.csv")  # or: grid-lines.csv
        IN_SUFFIX = "hsf"

        writePSOL = params.internal.mapValuesWritePSOL
        writeMGH = params.internal.mapValuesWriteMGH
        writeANNOT = params.internal.mapValuesWriteANNOT
        INTEGRATE = params.internal.mapValuesIntegrate
        SELECT = params.internal.mapValuesSelect
        INTERP = params.internal.mapValuesInterp

    # -------------------------------------------------------
    # check for differences

    cmd = (
        os.path.join(os.environ.get("FREESURFER_HOME"), "bin", "mri_diff")
        + " "
        + "--notallow-acq"
        + " "
        + "--notallow-prec"
        + " "
        + "--notallow-pix"
        + " "
        + IN_VOL
        + " "
        + IN_LABEL
    )
    # + " " + "--log" + " " + os.path.splitext(IN_VOL)[0] + ".diff"

    print(cmd)

    subproc = subprocess.run(cmd.split(), capture_output=True)

    print(subproc.returncode)

    # reslice input if necessary
    if subproc.returncode == 1:
        sys.exit(1)
    elif subproc.returncode > 100:
        if os.path.isdir(os.path.join(os.path.dirname(IN_VOL), "labels")):
            RESLICED_VOL = os.path.join(
                os.path.dirname(IN_VOL),
                "labels",
                os.path.splitext(os.path.basename(IN_VOL))[0] + "_reslicedToHighresHC" + os.path.splitext(IN_VOL)[1],
            )
        else:
            RESLICED_VOL = os.path.splitext(IN_VOL)[0] + "_reslicedToHighresHC" + os.path.splitext(IN_VOL)[1]

        cmd = (
            os.path.join(os.environ.get("FREESURFER_HOME"), "bin", "mri_convert")
            + " -rl "
            + IN_LABEL
            + " -rt "
            + INTERP
            + " "
            + IN_VOL
            + " "
            + RESLICED_VOL
        )

        print(cmd)

        subprocess.run(cmd.split(), capture_output=True)

    # read files
    if subproc.returncode > 100:
        vol = nb.freesurfer.load(RESLICED_VOL)
    else:
        vol = nb.freesurfer.load(IN_VOL)

    surf = TriaMesh.read_vtk(IN_SURF)

    # get data
    mat = vol.header.get_vox2ras_tkr()
    dat = vol.get_fdata()
    ind = np.array(np.nonzero(dat)).transpose()

    # variant 1: do it in surface RAS space (disadvantage: need to do distance
    # computation for many unnecessary values)
    # from scipy import spatial as sp
    # coordXYZ = np.concatenate((ind, np.ones((len(ind), 1))), axis=1)
    # coordSurfRAS = np.matmul(coordXYZ, mat.transpose())[:, 0:3]
    # dst = sp.distance_matrix(surf.v, coordSurfRAS)

    # variant 2: do it in XYZ space (easier lookup for relevant mid-surface
    # vertices)
    vtcs = np.concatenate((surf.v, np.ones((len(surf.v), 1))), axis=1)
    vtcsXYZ = np.matmul(vtcs, np.linalg.inv(mat).transpose())[:, 0:3]
    vtcsXYZ = np.round(vtcsXYZ).astype("int")
    lookup = dat[vtcsXYZ[:, 0], vtcsXYZ[:, 1], vtcsXYZ[:, 2]]

    # integrate along z-axis
    indices = np.array(pd.read_csv(IN_INDICES, header=None, index_col=None))
    indices3D = indices[:, 0:3].astype(int)
    lookup3D = np.full((np.max(indices3D, axis=0) + 1), np.nan)
    for i in range(0, len(lookup)):
        lookup3D[indices3D[i, 0], indices3D[i, 1], indices3D[i, 2]] = lookup[i]

    # if no selection is specified, set SELECT to all available indices
    if SELECT is None:
        SELECT = np.unique(indices3D[:, 2])

    if INTEGRATE == "mean":
        integr = np.nanmean(lookup3D[:, :, SELECT], axis=2)
    elif INTEGRATE == "median":
        integr = np.nanmedian(lookup3D[:, :, SELECT], axis=2)
    elif INTEGRATE == "mode":
        integr = st.mode(lookup3D[:, :, SELECT], axis=2)[0]
    elif INTEGRATE == "max":
        integr = np.nanmax(lookup3D[:, :, SELECT], axis=2)
    elif INTEGRATE == "min":
        integr = np.nanmin(lookup3D[:, :, SELECT], axis=2)
    elif INTEGRATE == "nonzeromin":
        lookup3D[lookup3D == 0] = np.Inf
        integr = np.nanmin(lookup3D[:, :, SELECT], axis=2)
    elif INTEGRATE == "none":
        if len(SELECT) > 1:
            LOGGER.info("Error: cannot use --integrate none with multiple sampling points, exiting.")
            sys.exit(1)
        else:
            integr = lookup3D[:, :, SELECT]

    integr = np.vstack(
        (
            np.indices(integr.shape)[0].flatten(),
            np.indices(integr.shape)[1].flatten(),
            integr.flatten(),
        )
    ).transpose()
    integr = pd.DataFrame(integr).sort_values([0, 1])

    # output
    OUT_CSV = IN_SURF.replace(".vtk", "." + IN_SUFFIX + ".csv")
    pd.DataFrame(integr).to_csv(OUT_CSV, header=False, index=False)

    if writePSOL is True:
        if (
            INTEGRATE == "mean"
            or INTEGRATE == "median"
            or INTEGRATE == "mode"
            or INTEGRATE == "max"
            or INTEGRATE == "min"
        ):
            OUT_PSOL_INTEGR = IN_SURF.replace(".vtk", "." + IN_SUFFIX + "-integrated.psol")
            io.write_vfunc(OUT_PSOL_INTEGR, np.asarray(integr)[:, 2])
        else:
            OUT_PSOL = IN_SURF.replace(".vtk", "." + IN_SUFFIX + ".psol")
            io.write_vfunc(OUT_PSOL, lookup)

    if writeMGH is True:
        if (
            INTEGRATE == "mean"
            or INTEGRATE == "median"
            or INTEGRATE == "mode"
            or INTEGRATE == "max"
            or INTEGRATE == "min"
        ):
            OUT_MGH_INTEGR = IN_SURF.replace(".vtk", "." + IN_SUFFIX + "-integrated.mgh")
            nb.freesurfer.save(
                nb.freesurfer.MGHImage(dataobj=np.asarray(integr)[:, 2].astype("float32"), affine=None),
                filename=OUT_MGH_INTEGR,
            )
        else:
            OUT_MGH = IN_SURF.replace(".vtk", "." + IN_SUFFIX + ".mgh")
            nb.freesurfer.save(
                nb.freesurfer.MGHImage(dataobj=lookup.astype("float32"), affine=None),
                filename=OUT_MGH,
            )

    if writeANNOT is True:
        OUT_ANNOT = IN_SURF.replace(".vtk", "." + IN_SUFFIX + ".annot")
        # ctab = np.array([(63,63,63,255,0), (255,0,0,255,234), (0,255,0,255,236), (0,0,255,255,238), (255,255,0,255,240), (255,255,0,255,246)])
        # names = ['Void', 'PrSbc', 'Sbc', 'CA1', 'CA2/3', 'ML']
        # labels = np.zeros(len(lookup))
        # labels[lookup==234] = 1
        # labels[lookup==236] = 2
        # labels[lookup==238] = 3
        # labels[lookup==240] = 4
        # labels[lookup==246] = 5
        # labels = labels.astype("int")
        # nb.freesurfer.write_annot(OUT_ANNOT, labels=labels, ctab=ctab, names=names, fill_ctab=False)

    # --------------------------------------------------------------------------
    # return

    return params


# ------------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# MAIN PART

if __name__ == "__main__":
    # Command line options and error checking

    options = _parse_arguments()

    options = _check_arguments(options)

    # Run analysis

    mapValues(
        params=options.params,
        IN_VOL=options.IN_VOL,
        IN_SURF=options.IN_SURF,
        IN_LABEL=options.IN_LABEL,
        IN_INDICES=options.IN_INDICES,
        IN_SUFFIX=options.IN_SUFFIX,
        writePSOL=options.writePSOL,
        writeMGH=options.writeMGH,
        writeANNOT=options.writeANNOT,
        INTEGRATE=options.INTEGRATE,
        SELECT=options.SELECT,
        INTERP=options.INTERP,
    )
