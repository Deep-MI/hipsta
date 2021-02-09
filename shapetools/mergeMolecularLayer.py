"""
This module provides a function to merge voxels from the molecular layer

"""

# ------------------------------------------------------------------------------
# AUXILIARY FUNCTIONS

# options_parse()

def options_parse():
    """
    A function to parse command line options

    """

    # imports

    import optparse
    import sys
    import os

    # define help text
    HELPTEXT = """
    SUMMARY

    This is an auxiliary script for the shapetools.py script and is normally
    called from within that script. It attaches the molecular layer on a
    voxel-by-voxel basis to its neighboring structures.

    The script requires three arguments:

    --inputfile  <file>     input file. an mgh / mgz image file, preferably including the pathname
    --outputfile <filename> output filename. an mgh / mgz image filename, preferably including the pathname
    --offset     <int>      a numerical value indicating the molecular layer label offset

    """

    # initialize
    parser = optparse.OptionParser(usage=HELPTEXT)

    # help text
    h_inputfile  = 'input file. an mgh / mgz image file, preferably including the pathname'
    h_outputfile = 'output file. an mgh / mgz image filename, preferably including the pathname'
    h_offset     = 'a numerical value indicating the molecular layer label offset (0|1000)'
    h_suffix     = 'a suffix for the output file (assigned|merged)'

    # specify inputs
    group = optparse.OptionGroup(parser, "Required Options:", "...")
    group.add_option('--inputfile', dest='inputfile', help=h_inputfile)
    group.add_option('--outputfile', dest='outputfile', help=h_outputfile)
    group.add_option('--offset', dest='offset', help=h_offset, type='int' )
    group.add_option('--suffix', dest='suffix', help=h_suffix)
    parser.add_option_group(group)

    # parse arguments
    (options, args) = parser.parse_args()

    # check if there are any inputs
    if len(sys.argv)==1:
        print(HELPTEXT)
        sys.exit(0)

    # check if input file is given
    if options.inputfile is None:
        print('\nERROR: Specify --inputfile\n')
        sys.exit(1)
    else:
        print('... Found input file '+options.inputfile)

    # check if output file is given
    if options.outputfile is None:
        print('\nERROR: Specify --outputfile\n')
        sys.exit(1)
    else:
        print('... Found output filename '+options.outputfile)

    # check if label offset is given
    if options.offset is None:
        print('\nERROR: Specify --offset\n')
        sys.exit(1)
    else:
        print('... Found label offset '+str(options.offset))

    # check if suffix is given
    if options.suffix is None:
        print('\nERROR: Specify --suffix\n')
        sys.exit(1)
    else:
        print('... Found suffix '+str(options.suffix))

    # check if input file exists
    if not os.path.isfile(options.inputfile):
        print('ERROR: input list '+options.inputfile+' is not an existing regular file\n')
        sys.exit(1)

    # return
    return options


# -----------------------------------------------------------------------------
# MAIN FUNCTIONS

# mergeML()

def mergeML(params, inputfile=None, outputfile=None, offset=None, suffix=None):

    # imports
    import nibabel as nb
    import numpy as np
    import sys
    import os

    # evaluate params
    if params is not None:
        inputfile = os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_00 + ".mgz")
        outputfile = os.path.join(params.OUTDIR, "merge-ml", params.HEMI + ".ml." + params.internal.HSFLABEL_00 + "_" + suffix + ".mgz")

    # load image
    im=nb.load(inputfile)

    # get voxel data
    vx=im.get_data()

    # get all voxels that are not zero; sort by dimensions 1,2,3
    vxNZ=list(np.nonzero(vx))
    vxNZ.append(vx[np.nonzero(vx)])
    vxNZ=np.transpose(vxNZ)
    vxNZ=vxNZ[np.argsort(vxNZ[:,0]),:]
    vxNZ=vxNZ[np.argsort(vxNZ[:,1],kind='mergesort'),:]
    vxNZ=vxNZ[np.argsort(vxNZ[:,2],kind='mergesort'),:]


    # get all ML voxels
    vxML=list(np.where((vx==245)|(vx==246)|(vx==214)))
    vxML.append(vx[np.where((vx==245)|(vx==246)|(vx==214))])
    vxML=np.transpose(vxML)
    vxML=vxML[np.argsort(vxML[:,0]),:]
    vxML=vxML[np.argsort(vxML[:,1],kind='mergesort'),:]
    vxML=vxML[np.argsort(vxML[:,2],kind='mergesort'),:]

    # check if the following returns same results as R
    for i in range(np.shape(vxML)[0]):
        n=1
        tmpML=np.transpose((np.empty([0]),np.empty([0]),np.empty([0])))
        while np.shape(tmpML)[0]==0:
            tmp=vx[np.ix_(range(np.int((vxML[i,0]-n)),np.int((vxML[i,0]+n)+1)),range(np.int((vxML[i,1]-n)),np.int((vxML[i,1]+n)+1)),range(np.int((vxML[i,2]-n)),np.int((vxML[i,2]+n)+1)))]
            tmpML=np.transpose(np.where((tmp!=0)&(tmp!=246)&(tmp!=245)&(tmp!=214)))
            if np.shape(tmpML)[0]!=0:
                tmp0=np.asmatrix(tmpML[np.where(np.sum(tmpML,axis=1)==min(np.sum(tmpML,axis=1))),:])
                vxML[i,3]=offset+min(tmp[tmp0[:,0],tmp0[:,1],tmp0[:,2]])
                tmp0=None
            else:
                n=n+1
            tmp=None
        tmpML=None

    # write back to vx and im
    vxOld=vx
    for i in range(np.shape(vxML)[0]):
        vx[np.int(vxML[i,0]),np.int(vxML[i,1]),np.int(vxML[i,2])]=vxML[i,3]

    om=nb.MGHImage(dataobj=vx,affine=im.get_affine(),header=im.header)
    nb.save(img=om,filename=outputfile)

    # update params
    if params is not None:
        params.internal.HSFLABEL_01 = "ml."+params.internal.HSFLABEL_00

    # return
    return(params)


# mergeMolecularLayer()

def mergeMolecularLayer(params):
    """


    """

    # import
    import os
    import shutil

    # check if ML is present in lookup table and hsflist; if not, just copy
    # files

    if "ML" in params.LUTDICT.keys():

        #
        print("Attaching the molecular layer")

        # analysis 1
        params = mergeML(params=params, offset=0, suffix="assigned")

        # analyis 2
        params = mergeML(params=params, offset=1000, suffix="merged")

        # copy
        shutil.copyfile(
            os.path.join(params.OUTDIR, "merge-ml", params.HEMI + ".ml." + params.internal.HSFLABEL_00 + "_assigned.mgz"),
            os.path.join(params.OUTDIR, params.HEMI + ".ml." + params.internal.HSFLABEL_00 + ".mgz"))

    else:

            #
            print("Not attaching the molecular layer")

            # copy
            shutil.copyfile(
                os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_00 + ".mgz"),
                os.path.join(params.OUTDIR, "merge-ml", params.HEMI + ".ml." + params.internal.HSFLABEL_00 + "_merged.mgz"))

            shutil.copyfile(
                os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_00 + ".mgz"),
                os.path.join(params.OUTDIR, "merge-ml", params.HEMI + ".ml." + params.internal.HSFLABEL_00 + "_assigned.mgz"))

            shutil.copyfile(
                os.path.join(params.OUTDIR, "merge-ml", params.HEMI + ".ml." + params.internal.HSFLABEL_00 + "_assigned.mgz"),
                os.path.join(params.OUTDIR, params.HEMI + ".ml." + params.internal.HSFLABEL_00 + ".mgz"))

            # update params
            if params is not None:
                params.internal.HSFLABEL_01 = "ml."+params.internal.HSFLABEL_00

    # return
    return(params)

# -----------------------------------------------------------------------------
# MAIN PART

if __name__=="__main__":

    # message

    print("-----------------------------------------------------------------")
    print("Classification and assignment of the molecular layer")
    print("-----------------------------------------------------------------")

    print("\nReading input options ...")

    # command Line options and error checking

    options = options_parse()

    # run main function

    mergeML(inputfile=options.inputfile, outputfile=options.outputfile, offset=options.offset, suffix=options.suffix, params=None)
