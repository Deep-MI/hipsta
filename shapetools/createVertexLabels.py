"""


"""

# ------------------------------------------------------------------------------
# AUXILIARY FUNCTIONS

def options_parse():
    """
    Command Line Options Parser:
    initiate the option parser and return the parsed object
    """

    # imports
    import optparse
    import sys

    # define helptext
    HELPTEXT = """
    SUMMARY

    This is an auxiliary script for the shapetools.py script and is usually
    called from within that script.

    The script requires three arguments:

    --lst   <file>
    --mat   <file>
    --lab   <file>

    """

    # initialize parser
    parser = optparse.OptionParser(usage=HELPTEXT)

    # help text
    h_lst = 'lst file'
    h_mat = 'mat file'
    h_lab = 'lab file'

    # specify inputs
    group = optparse.OptionGroup(parser, "Required Options:", "...")
    group.add_option('--lst', dest='lst', help=h_lst)
    group.add_option('--mat', dest='mat', help=h_mat)
    group.add_option('--lab', dest='lab', help=h_lab)
    parser.add_option_group(group)

    # parse arguments
    options, args = parser.parse_args()

    # check if there are any inputs
    if len(sys.argv) == 1:
        print(HELPTEXT)
        sys.exit(0)

    # check if lst file is given
    if options.lst is None:
        print('\nERROR: Specify --lst\n')
        sys.exit(1)
    else:
        print('... Found lst file ' + options.lst)

    # check if mat file is given
    if options.mat is None:
        print('\nERROR: Specify --mat\n')
        sys.exit(1)
    else:
        print('... Found mat file ' + options.mat)

    # check if lab file is given
    if options.lab is None:
        print('\nERROR: Specify --lab\n')
        sys.exit(1)
    else:
        print('... Found lab file ' + options.lab)

    # return
    return options


# -----------------------------------------------------------------------------
# MAIN FUNCTION

def createVertexLabels(lstFile, matFile, labFile):
    """
    createVertexLabels(createVertexLabels(lstFile, matFile, labFile)

    """

    # imports

    import numpy as np

    #

    lst = np.loadtxt(lstFile, dtype="float")
    mat = np.loadtxt(matFile, dtype="float")
    lab = np.loadtxt(labFile, dtype="float", skiprows=2)

    crs = np.array(lst)  # to avoid mere referencing
    crs[:, 3] = 1

    ras = np.matmul(mat, crs.transpose()).transpose()
    ras = np.append(ras[:, 0:3], lst[:, 3:4], axis=1)

    vtx = np.zeros(shape=np.shape(lab)[0])
    for i in range(np.shape(vtx)[0]):
        tmp = np.linalg.norm(ras[:, 0:3] - np.repeat(lab[i:(i + 1), 1:4], np.shape(ras)[0], axis=0), ord=2, axis=1)
        vtx[i] = lst[np.where(tmp == tmp.min())[0][0], 3]  # note: we added [0][0] here

    # the following lines will produce zero-mean, step-one labels
    # key = np.array(range(0, len(np.unique(vtx)))) - np.mean(range(0, len(np.unique(vtx))))
    # index = np.digitize(vtx, np.unique(vtx), right=True)
    # np.savetxt(fname=labFile.replace(".label", ".asc"), X=key[index])

    # the following line will keep original labels
    np.savetxt(fname=labFile.replace(".label", ".asc"), X=vtx)


# -----------------------------------------------------------------------------
# MAIN PART

if __name__ == "__main__":

    # message

    print("-------------------------------------------------------------------")
    print("Mapping vertex label to surface")
    print("-------------------------------------------------------------------")

    #

    options = options_parse()

    #

    createVertexLabels(options.lst, options.mat, options.lab)
