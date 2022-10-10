"""
This module provides a function to fill holes and create binary masks

"""

# ------------------------------------------------------------------------------
# main function

def fillHoles(params):

    """

    """

    # imports

    import os
    import subprocess

    # message

    print()
    print("-------------------------------------------------------------------------")
    print()
    print("Fill holes")
    print()
    print("-------------------------------------------------------------------------")
    print()

    # dilate/erode by n voxels to fill holes; we use the FSL tool because
    # mri_binarize cannot be applied without binarization, which we don't want yet.
    # outputs ${HEMI}.de.${HSFLABEL_01}.mgz

    # convert forward

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_convert") + " " \
        + "--in_type mgz " \
        + "--out_type nii " \
        + os.path.join(params.OUTDIR, "merge-ml", params.HEMI + "." + params.internal.HSFLABEL_01 + "_merged.mgz") + " " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + "." + params.internal.HSFLABEL_01 + "_merged.nii")

    print(cmd)

    subprocess.run(cmd.split())

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_convert") + " " \
        + "--in_type mgz " \
        + "--out_type nii " \
        + os.path.join(params.OUTDIR, "merge-ml", params.HEMI + "." + params.internal.HSFLABEL_01 + "_assigned.mgz") + " " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + "." + params.internal.HSFLABEL_01 + "_assigned.nii")

    print(cmd)

    subprocess.run(cmd.split())

    #

    os.environ["FSLOUTPUTTYPE"] = "NIFTI"

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "fslmaths.fsl") + " " \
        + "-dt int " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + "." + params.internal.HSFLABEL_01 + "_merged.nii") + " " \
        + "-kernel boxv " + str(params.internal.DIL) + " -dilD " \
        + "-kernel boxv " + str(params.internal.ERO) + " -ero " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_merged.nii")

    print(cmd)

    subprocess.run(cmd.split())

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "fslmaths.fsl") + " " \
        + "-dt int " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + "." + params.internal.HSFLABEL_01 + "_assigned.nii") + " " \
        + "-kernel boxv " + str(params.internal.DIL) + " -dilD " \
        + "-kernel boxv " + str(params.internal.ERO) + " -ero " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_assigned.nii")

    print(cmd)

    subprocess.run(cmd.split())

    # convert back

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_convert") + " " \
        + "--in_type nii " \
        + "--out_type mgz " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_merged.nii") + " " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_merged.mgz")

    print(cmd)

    subprocess.run(cmd.split())

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_convert") + " " \
        + "--in_type nii " \
        + "--out_type mgz " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_assigned.nii") + " " \
        + os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_assigned.mgz")

    print(cmd)

    subprocess.run(cmd.split())

    # clean up

    if params.skipCLEANUP is False:

        os.remove(os.path.join(params.OUTDIR, "mask", params.HEMI + "." + params.internal.HSFLABEL_01 + "_merged.nii"))
        os.remove(os.path.join(params.OUTDIR, "mask", params.HEMI + "." + params.internal.HSFLABEL_01 + "_assigned.nii"))

        os.remove(os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_merged.nii"))
        os.remove(os.path.join(params.OUTDIR, "mask", params.HEMI + ".de." + params.internal.HSFLABEL_01 + "_assigned.nii"))

    # update params

    params.internal.HSFLABEL_02 = "de." + params.internal.HSFLABEL_01

    # return

    return(params)


def createMask(params):
    """

    """

    # imports

    import os
    import subprocess

    # message

    print()
    print("-------------------------------------------------------------------------")
    print()
    print("Create binary mask")
    print()
    print("-------------------------------------------------------------------------")
    print()

    # mri_binarize

    cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_binarize") + " " \
        + "--i " + os.path.join(params.OUTDIR, "mask", params.HEMI + "." + params.internal.HSFLABEL_02 + "_merged.mgz") + " " \
        + "--min 1 " \
        + "--o " + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_02 + ".mgz")

    print(cmd)

    subprocess.run(cmd.split())

    # return

    return(params)


def filterMask(params):
    """

    """

    # imports

    import os
    import subprocess

    #

    if params.internal.FILTERMASK is not None:

        # message

        print()
        print("-------------------------------------------------------------------------")
        print()
        print("Filter mask")
        print()
        print("-------------------------------------------------------------------------")
        print()

        #

        os.environ["FSLOUTPUTTYPE"] = "NIFTI"

        # convert

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_convert") + " " \
            + "--in_type mgz " \
            + "--out_type nii " \
            + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_02 + ".mgz") + " " \
            + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_02 + ".nii")

        print(cmd)

        subprocess.run(cmd.split())

        # filter

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "fslmaths.fsl") + " " \
            + os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_02 + ".nii") + " " \
            + "-mul 100 -kernel gauss " + str(params.internal.FILTERMASK[0]) + " -fmean " \
            + os.path.join(params.OUTDIR, params.HEMI + ".filt." + params.internal.HSFLABEL_02 + ".nii")

        print(cmd)

        subprocess.run(cmd.split())

        # convert back

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_convert") + " " \
            + "--in_type nii " \
            + "--out_type mgz " \
            + os.path.join(params.OUTDIR, params.HEMI + ".filt." + params.internal.HSFLABEL_02 + ".nii") + " " \
            + os.path.join(params.OUTDIR, params.HEMI + ".filt." + params.internal.HSFLABEL_02 + ".mgz")

        print(cmd)

        subprocess.run(cmd.split())

        # binarize

        cmd = os.path.join(os.environ.get('FREESURFER_HOME'), "bin", "mri_binarize") + " " \
            + "--i " + os.path.join(params.OUTDIR, params.HEMI + ".filt." + params.internal.HSFLABEL_02 + ".mgz") + " " \
            + "--min " + str(params.internal.FILTERMASK[1]) + " --binval 1 " \
            + "--o " + os.path.join(params.OUTDIR, params.HEMI + ".filt." + params.internal.HSFLABEL_02 + ".mgz")

        print(cmd)

        subprocess.run(cmd.split())

        # clean up

        os.remove(os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_02 + ".nii"))

        os.remove(os.path.join(params.OUTDIR, params.HEMI + ".filt." + params.internal.HSFLABEL_02 + ".nii"))

        # update params

        params.internal.HSFLABEL_02 = "filt." + params.internal.HSFLABEL_02

    # return

    return(params)


def longFilter(params):
    """

    """

    # imports

    import os
    import subprocess
    from scipy import ndimage
    import nibabel as nb
    import numpy as np

    #

    if params.internal.LONGFILTER is not None:

        # message

        print()
        print("-------------------------------------------------------------------------")
        print()
        print("Filtering along longitudinal axis")
        print()
        print("-------------------------------------------------------------------------")
        print()

        # longitudinal filtering

        img = nb.load(os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_02 + ".mgz"))
        dat = img.get_fdata()

        #k = np.zeros((3,3,3))
        #k[1,1,:] = 1
        #k[:,1,1] = 1

        k = np.zeros(params.internal.LONGFILTER_SIZE)
        k[2,2,:] = 1

        dat_filtered = ndimage.convolve(dat, k, mode="constant", cval=0.0)

        dat_filtered = dat_filtered>0

        # opening
        #dat_filtered = ndimage.binary_opening(dat_filtered)

        # closing
        #dat_filtered = ndimage.binary_closing(dat_filtered)

        nb.freesurfer.save(nb.freesurfer.MGHImage(dataobj=dat_filtered.astype("float32"), affine=img.get_affine()), filename=os.path.join(params.OUTDIR, params.HEMI + ".lf." + params.internal.HSFLABEL_02 + ".mgz"))

        # update params

        params.internal.HSFLABEL_02 = "lf." + params.internal.HSFLABEL_02

    # return

    return(params)


def closeMask(params):
    """

    """

    # imports

    import os
    import subprocess
    from scipy import ndimage
    import nibabel as nb
    import numpy as np

    #

    if params.internal.CLOSEMASK is not None:

        # message

        print()
        print("-------------------------------------------------------------------------")
        print()
        print("Applying closing operation to mask")
        print()
        print("-------------------------------------------------------------------------")
        print()

        if params.internal.CLOSEMASK == "regular":

            # get data

            img = nb.load(os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_02 + ".mgz"))
            dat = img.get_fdata()

            #dat_filtered = ndimage.binary_closing(dat, iterations=1)
            dat_filtered = ndimage.binary_closing(dat, structure=ndimage.iterate_structure(ndimage.generate_binary_structure(rank=3, connectivity=1), iterations=1).astype(int), iterations=1)

            dat_filtered = dat_filtered>0

            nb.freesurfer.save(nb.freesurfer.MGHImage(dataobj=dat_filtered.astype("float32"), affine=img.get_affine()), filename=os.path.join(params.OUTDIR, params.HEMI + ".cm." + params.internal.HSFLABEL_02 + ".mgz"))

        elif params.internal.CLOSEMASK == "experimental":

            # set all voxels to .5 grid
            print("Experimental variant for closing operation")

            # get dataobj
            img = nb.load(os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_02 + ".mgz"))
            dat = img.get_fdata()

            # compute new grid
            dims = img.header["dims"]
            grd = np.array(np.meshgrid(
                np.arange(start=0.5, stop=dims[0]-1),
                np.arange(start=0.5, stop=dims[1]-1),
                np.arange(start=0.5, stop=dims[2]-1), indexing="ij"))

            # get values at new grid
            img_grd_dat = ndimage.map_coordinates(img.get_fdata(), grd)

            # binarize
            img_grd_dat = img_grd_dat > 0.1

            # opening
            #img_grd_dat = ndimage.binary_opening(img_grd_dat)

            # closing
            #img_grd_dat = ndimage.binary_closing(img_grd_dat)

            #  compute new affine
            off = np.matmul(img.affine[0:3,0:3], np.expand_dims(np.array([0.5,0.5,0.5]), axis=1))
            img_grd_affn = img.affine.copy()
            img_grd_affn[0:3,3] += off.flatten()

            # create new image and save
            img_grd = nb.MGHImage(dataobj=img_grd_dat.astype("float32"), affine=img_grd_affn)
            nb.save(img_grd, os.path.join(params.OUTDIR, params.HEMI + ".cm." + params.internal.HSFLABEL_02 + ".mgz"))

        elif params.internal.CLOSEMASK == "experimental_v2":

            # set voxels in longitudinal direction to .5 grid
            print("Experimental variant v2 for closing operation")

            # get dataobj
            img = nb.load(os.path.join(params.OUTDIR, params.HEMI + "." + params.internal.HSFLABEL_02 + ".mgz"))
            dat = img.get_fdata()

            # compute new grid
            dims = img.header["dims"]
            grd = np.array(np.meshgrid(
                np.arange(start=0, stop=dims[0]),
                np.arange(start=0, stop=dims[1]),
                np.arange(start=0.5, stop=dims[2]-1), indexing="ij"))

            # get values at new grid
            img_grd_dat = ndimage.map_coordinates(img.get_fdata(), grd)

            # binarize
            img_grd_dat = img_grd_dat > 0.1

            # opening
            img_grd_dat = ndimage.binary_opening(img_grd_dat)

            # closing
            img_grd_dat = ndimage.binary_closing(img_grd_dat)
            #img_grd_dat = ndimage.binary_closing(img_grd_dat, structure=ndimage.generate_binary_structure(rank=3, connectivity=3))

            # closing
            #img_grd_dat = ndimage.binary_closing(img_grd_dat)

            #  compute new affine
            off = np.matmul(img.affine[0:3,0:3], np.expand_dims(np.array([0,0,0.5]), axis=1))
            img_grd_affn = img.affine.copy()
            img_grd_affn[0:3,3] += off.flatten()

            # create new image and save
            img_grd = nb.MGHImage(dataobj=img_grd_dat.astype("float32"), affine=img_grd_affn)
            nb.save(img_grd, os.path.join(params.OUTDIR, params.HEMI + ".cm." + params.internal.HSFLABEL_02 + ".mgz"))

        # update params

        params.internal.HSFLABEL_02 = "cm." + params.internal.HSFLABEL_02

    # return

    return(params)
