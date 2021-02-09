"""
This module provides a function to map values from a volume to the midsurface.

"""

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
    --table <input csv table>
    --suffix <suffix>

    The following arguments are optional and determine which output files will
    be created:

    --integrate <none|mode|median|mean> (default: none)
    --select <0, 1, ..., n> (default: all)
    --writePSOL
    --writeMGH

    EXAMPLES

    Example for mapping at the vertices of the midsurface:

    python3 mapValues.py --volume /my/path/to/<MY_VOLUME.mgz> \\
        --surface /my/path/to/<SUBJECT_ID>/thickness/lh.mid-surface.vtk \\
        --table /my/path/to/<SUBJECT_ID>/thickness/lh.mid-surface.csv \\
        --suffix MY_SUFFIX --writePSOL --writeMGH

    This will sample the vertices of the midsurface at the corresponding
    positions in the given volume. This mapping requires the '[lr]h.mid-
    surface.vtk' surface and the '[lr]h.mid-surface.csv' index table that are
    created by the hippocampal thickness tools. The suffix will be used to
    identify the output files and can be chosen arbitrarily. PSOL and MGH
    output files will be written to the same directory as the input surface.

    Example for mapping at the vertices of the streamlines:

    python3 mapValues.py --volume /my/path/to/<MY_VOLUME.mgz> \\
        --surface /my/path/to/<SUBJECT_ID>/thickness/lh.grid-lines-z.vtk \\
        --table /my/path/to/<SUBJECT_ID>/thickness/lh.grid-lines.csv \\
        --suffix MY_SUFFIX --integrate mean --select 0 1 2 3 4 5 --writePSOL \\
        --writeMGH

    This will sample the vertices of the interior/exterior streamlines at the
    corresponding positions in the given volume. That is, multiple vertices per
    streamlines will be integrated and projected onto the corresponding point
    of the midsurface. This mapping requires the '[lr]h.grid-lines-z.vtk' file
    and the '[lr]h.grid-lines.csv' index table that is created by the hippo-
    campal thickness tools. Also, the 'integrate' argument is required. In this
    example, the sampled values will be averaged ('mean') per streamline. The
    select argument can be used to restrict the integration to a subset of the
    integration points (0 ... 10 for the default z direction, 5 being the
    midsurface). The suffix will be used to identify the output files and can
    be chosen arbitrarily. PSOL and MGH output files will be written to the
    same directory as the input surface.

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

    # imports
    import sys
    import argparse

    # message
    print("\nReading input options ...")

    # setup parser
    parser = argparse.ArgumentParser(description="", add_help=False)

    # help text
    h_IN_VOL = 'input volume'
    h_IN_SURF = 'input surface'
    h_IN_INDICES = 'table with indices'
    h_IN_SUFFIX = 'suffix for output files'
    h_writePSOL = 'write out PSOL files'
    h_writeMGH = 'write out MGH files'
    h_writeANNOT = 'write out ANNOT files'
    h_integrate = 'write out integrated values (default: none)'
    h_select = 'select inegration points (default: all)'

    # required arguments
    required = parser.add_argument_group('required arguments')

    required.add_argument('--volume', dest='IN_VOL', help=h_IN_VOL, default=None, required=False, metavar='<file>')
    required.add_argument('--surface', dest='IN_SURF', help=h_IN_SURF, default=None, required=False, metavar='<file>')
    required.add_argument('--table', dest='IN_INDICES', help=h_IN_INDICES, default=None, required=False, metavar='<file>')
    required.add_argument('--suffix', dest='IN_SUFFIX', help=h_IN_SUFFIX, default=None, required=False, metavar='<string>')

    # optional argumentws
    optional = parser.add_argument_group('optional arguments')

    optional.add_argument('--integrate', dest='integrate', help=h_integrate, default='none', required=False, metavar= '<none|mode|median|mean>')
    optional.add_argument('--select', dest='select', help=h_select, default=None, nargs='+', required=False, metavar= '<list>')
    optional.add_argument('--writePSOL', dest='writePSOL', help=h_writePSOL, default=False, action="store_true", required=False)
    optional.add_argument('--writeMGH', dest='writeMGH', help=h_writeMGH, default=False, action="store_true", required=False)
    optional.add_argument('--writeANNOT', dest='writeANNOT', help=argparse.SUPPRESS, default=False, action="store_true", required=False) # help=h_writeANNOT, currently hidden might  be added later

    # define help
    help = parser.add_argument_group('getting help')

    help.add_argument('--help', help="display this help message and exit", action='help')
    help.add_argument('--more-help', dest='more_help', help="display extensive help message and exit", default=False, action="store_true", required=False)

    # check if there are any inputs; if not, print help and exit
    if len(sys.argv) == 1:
        args = parser.parse_args(['--help'])
    else:
        args = parser.parse_args()

    return args

# ------------------------------------------------------------------------------
# check arguments

def _check_arguments(options):
    """
    an internal function to check input arguments

    """

    # imports

    import sys
    import numpy as np

    # more help

    if options.more_help is True:
        get_help()
        sys.exit(0)

    # required arguments

    if options.IN_VOL is None:
        print("The --volume option is required. See --help for details.")
        sys.exit(1)

    if options.IN_SURF is None:
        print("The --surface option is required. See --help for details.")
        sys.exit(1)

    if options.IN_INDICES is None:
        print("The --table option is required. See --help for details.")
        sys.exit(1)

    if options.IN_SUFFIX is None:
        print("The --suffix option is required. See --help for details.")
        sys.exit(1)

    # change formats

    if options.select is not None:
        options.select = np.array(options.select).astype(int)

    # add params

    options.params = None

    # return

    return options

# ------------------------------------------------------------------------------
# MAIN FUNCTION
# ------------------------------------------------------------------------------

def mapValues(params, IN_VOL=None, IN_SURF=None, IN_INDICES=None, IN_SUFFIX='hsf', integrate='mode', select=['all'], writePSOL=False, writeMGH=False, writeANNOT=False):

    # imports

    import os
    import sys
    import nibabel as nb
    import numpy as np
    import pandas as pd
    import lapy as lp
    from lapy import TriaIO as lpio
    from lapy import FuncIO as lpfio
    from scipy import spatial as sp
    from scipy import stats as st
    from shapetools import triaUtils as tu

    # message

    print()
    print("-------------------------------------------------------------------")
    print()
    print("Mapping values")
    print()
    print("-------------------------------------------------------------------")
    print()

    #-------------------------------------------------------
    # get params

    if params is not None:

        # these are the default settings for processing within the hippocampal
        # thickness toolbox; will override any other settings; set 'params=None'
        # for custom processing.
        IN_VOL = os.path.join(params.OUTDIR, params.HEMI + ".ml." + params.internal.HSFLABEL_00 + ".mgz") # should be identical to the following line
        #IN_VOL = os.path.join(params.OUTDIR, 'merge-ml', params.HEMI + ".ml.234-236-238-240-246_assigned.mgz")
        IN_SURF = os.path.join(params.OUTDIR, 'thickness', params.HEMI + ".mid-surface.vtk") # or: grid-lines-z.vtk
        IN_INDICES = os.path.join(params.OUTDIR, 'thickness', params.HEMI + ".mid-surface.csv") # or: grid-lines.csv
        IN_SUFFIX = 'hsf'

        writePSOL = params.internal.mapValuesWritePSOL
        writeMGH = params.internal.mapValuesWriteMGH
        writeANNOT = params.internal.mapValuesWriteANNOT
        integrate = params.internal.mapValuesIntegrate
        select = params.internal.mapValuesSelect

    #-------------------------------------------------------
    # read files
    vol = nb.freesurfer.load(IN_VOL)
    surf = lpio.import_vtk(IN_SURF)

    # get data
    mat = vol.header.get_vox2ras_tkr()
    dat = vol.get_fdata()
    ind = np.array(np.nonzero(dat)).transpose()

    # variant 1: do it in surface RAS space (disadvantage: need to do distance
    # computation for many unnecessary values)
    #coordXYZ = np.concatenate((ind, np.ones((len(ind), 1))), axis=1)
    #coordSurfRAS = np.matmul(coordXYZ, mat.transpose())[:, 0:3]
    #dst = sp.distance_matrix(surf.v, coordSurfRAS)

    # variant 2: do it in XYZ space (easier lookup for relevant mid-surface
    # vertices)
    vtcs = np.concatenate((surf.v, np.ones((len(surf.v), 1))), axis=1)
    vtcsXYZ = np.matmul(vtcs, np.linalg.inv(mat).transpose())[:, 0:3]
    vtcsXYZ = np.round(vtcsXYZ).astype('int')
    lookup = dat[vtcsXYZ[:,0], vtcsXYZ[:,1], vtcsXYZ[:,2]]

    # integrate along z-axis
    integr = list()
    indices = np.array(pd.read_csv(IN_INDICES, header=None, index_col=None))
    labels = [ i + "_" + j for i, j in zip(indices[:,0].astype(int).astype(str).tolist(), indices[:,1].astype(int).astype(str).tolist()) ]
    uniqueLabels = np.unique(labels)
    for i in range(0, len(uniqueLabels)):
        if integrate == "mean":
            if options.select is not None:
               integr.append(np.mean(lookup[np.asarray(labels)==uniqueLabels[i]][options.select]))
            else:
               integr.append(np.mean(lookup[np.asarray(labels)==uniqueLabels[i]]))
        elif integrate == "median":
            if options.select is not None:
               integr.append(np.median(lookup[np.asarray(labels)==uniqueLabels[i]][options.select]))
            else:
               integr.append(np.median(lookup[np.asarray(labels)==uniqueLabels[i]]))
        elif integrate == "mode":
            if options.select is not None:
                integr.append(st.mode(lookup[np.asarray(labels)==uniqueLabels[i]][options.select])[0].item())
            else:
                integr.append(st.mode(lookup[np.asarray(labels)==uniqueLabels[i]])[0].item())
        elif integrate == "none":
            if options.select is not None:
                if len(lookup[np.asarray(labels)==uniqueLabels[i]][options.select])>1:
                    print("Error: cannot use --integrate none with multiple sampling points, exiting.")
                    sys.exit(1)
                else:
                    integr.append(lookup[np.asarray(labels)==uniqueLabels[i]][options.select].item())
            else:
                if len(lookup[np.asarray(labels)==uniqueLabels[i]])>1:
                    print("Error: cannot use --integrate none with multiple sampling points, exiting.")
                    sys.exit(1)
                else:
                    integr.append(lookup[np.asarray(labels)==uniqueLabels[i]].item())
    integr = np.concatenate((np.array([ x.split("_") for x in uniqueLabels ]).astype(int), np.array(integr, ndmin=2).transpose()), axis=1)
    integr = pd.DataFrame(integr).sort_values([0, 1])

    # output
    OUT_CSV = IN_SURF.replace('.vtk', '_'+IN_SUFFIX+'.csv')
    pd.DataFrame(integr).to_csv(OUT_CSV, header=False, index=False)

    if writePSOL is True:
        if integrate=="mean" or integrate=="median" or integrate=="mode":
            OUT_PSOL_INTEGR = IN_SURF.replace('.vtk', '_'+IN_SUFFIX+'-integrated.psol')
            lpfio.export_vfunc(OUT_PSOL_INTEGR, np.asarray(integr)[:,2])
        else:
            OUT_PSOL = IN_SURF.replace('.vtk', '_'+IN_SUFFIX+'.psol')
            lpfio.export_vfunc(OUT_PSOL, lookup)

    if writeMGH is True:
        if integrate=="mean" or integrate=="median" or integrate=="mode":
            OUT_MGH_INTEGR = IN_SURF.replace('.vtk', '_'+IN_SUFFIX+'-integrated.mgh')
            nb.freesurfer.save(nb.freesurfer.MGHImage(dataobj=np.asarray(integr)[:,2].astype("float32"), affine=None), filename=OUT_MGH_INTEGR)
        else:
            OUT_MGH = IN_SURF.replace('.vtk', '_'+IN_SUFFIX+'.mgh')
            nb.freesurfer.save(nb.freesurfer.MGHImage(dataobj=lookup.astype("float32"), affine=None), filename=OUT_MGH)

    if writeANNOT is True:
        OUT_ANNOT = IN_SURF.replace('.vtk', '_'+IN_SUFFIX+'.annot')
        #ctab = np.array([(63,63,63,255,0), (255,0,0,255,234), (0,255,0,255,236), (0,0,255,255,238), (255,255,0,255,240), (255,255,0,255,246)])
        #names = ['Void', 'PrSbc', 'Sbc', 'CA1', 'CA2/3', 'ML']
        #labels = np.zeros(len(lookup))
        #labels[lookup==234] = 1
        #labels[lookup==236] = 2
        #labels[lookup==238] = 3
        #labels[lookup==240] = 4
        #labels[lookup==246] = 5
        #labels = labels.astype("int")
        #nb.freesurfer.write_annot(OUT_ANNOT, labels=labels, ctab=ctab, names=names, fill_ctab=False)

    # --------------------------------------------------------------------------
    # return

    return(params)


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

    mapValues(params=options.params, IN_VOL=options.IN_VOL, IN_SURF=options.IN_SURF, IN_INDICES=options.IN_INDICES, IN_SUFFIX=options.IN_SUFFIX, writePSOL=options.writePSOL, writeMGH=options.writeMGH, writeANNOT=options.writeANNOT, integrate=options.integrate)
