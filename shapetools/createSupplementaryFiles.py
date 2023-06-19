"""
This module creates supplementary files, primarily for visualization.

"""

# ------------------------------------------------------------------------------
# MAIN FUNCTION
# ------------------------------------------------------------------------------

def createSupplementaryFiles(params):

    # imports

    import os
    import sys
    import logging

    import nibabel as nb
    import numpy as np

    from lapy import TriaMesh

    from shapetools import triaUtils

    # message

    print()
    print("-------------------------------------------------------------------")
    print()
    print("Creating supplementary files")
    print()
    print("-------------------------------------------------------------------")
    print()

    #-------------------------------------------------------
    # get params

    LUT = params.LUT
    HEMI = params.HEMI
    OUT_DIR = os.path.join(params.OUTDIR, 'thickness')

    # --------------------------------------------------------------------------
    # create overlay with subfield boundaries

    # list of boundary values:

    if LUT == "fs711":

        logging.info("Found internal, modified look-up table for FreeSurfer 7.11.")

        # get labels

        hsfBnd = nb.load(os.path.join(OUT_DIR, HEMI + '.mid-surface_hsf.mgh'))
        hsfBnd = np.array(hsfBnd.get_fdata()).flatten()

        # recode to consecutive values

        hsfBnd[hsfBnd==234] = 1
        hsfBnd[hsfBnd==236] = 2
        hsfBnd[hsfBnd==238] = 3
        hsfBnd[hsfBnd==240] = 4

        lstBnd = [ 1.5, 2.5, 3.5 ] # note that this implicitly also accounts for 246 labels if present

        #

        cont = True

    elif LUT == "ashs":

        logging.info("Found internal, modified look-up table for ASHS atlas.")

        # get labels

        hsfBnd = nb.load(os.path.join(OUT_DIR, HEMI + '.mid-surface_hsf.mgh'))
        hsfBnd = np.array(hsfBnd.get_fdata()).flatten()

        # recode to consecutive values

        hsfBnd[hsfBnd==8] = 11
        hsfBnd[hsfBnd==1] = 12
        hsfBnd[hsfBnd==2] = 13
        hsfBnd[hsfBnd==4] = 14

        lstBnd = [ 11.5, 12.5, 13.5 ]

        #

        cont = True

    else:

        logging.info("Creation of subfield boundary overlays is only supported for the fs711 and ashs look-up tables.")

        cont = False

    #

    if cont is True:

        # get mid-surface

        triaMidRm = TriaMesh.import_vtk(os.path.join(OUT_DIR, HEMI + '.mid-surface.vtk'))

        #

        vBnd = np.empty((0, 3))
        tBnd = np.empty((0, 3))

        for i in lstBnd:
            lvlBnd = triaUtils.levelsetsTria(triaMidRm.v, triaMidRm.t, hsfBnd, i)
            if len(lvlBnd[0][0]) > 0:
                vBnd_lvl = np.array(lvlBnd[0][0])
                tBnd_lvl = np.array(lvlBnd[1][0])[:,[0,1,1]] - 1
                tBnd = np.concatenate((tBnd, tBnd_lvl + len(vBnd)), axis=0)
                vBnd = np.concatenate((vBnd, vBnd_lvl), axis=0)

        # output

        if len(vBnd)>0:
            TriaMesh(v=vBnd, t=tBnd.astype('int')).write_vtk(outBndoutfile=os.path.join(OUT_DIR, HEMI + '.mid-surface_hsf-bnd.vtk'))
        else:
            logging.info("Could not create boundary overlays.")

