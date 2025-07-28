"""
This module provides configuration defaults for the hippocampal shape and 
thickness analysis package.

"""


def get_defaults(x):
    defaults = dict(
        # note: if changing these values, remember to adjust the argparse help
        # messages as well as the run_hipsta docstring also
        # options
        start_with_edited_labels=False,
        no_cleanup=False,
        no_crop=False,
        upsample=False,
        upsample_size=[0, 0, 0],
        no_merge_molecular_layer=False,
        automask_head=False,
        automask_tail=False,
        automask_head_margin=0,
        automask_tail_margin=0,
        no_gauss_filter=False,
        gauss_filter_size=[1, 50],
        long_filter=False,
        long_filter_size=5,
        no_close_mask=False,
        mca="mri_mc",
        remesh=False,
        smooth=5,
        cut_range=[-0.975, 0.975],
        aniso_alpha=[40],
        aniso_smooth=3,
        thickness_grid=[-0.9, 0.9, 41, -0.975, 0.975, 21, -0.9, 0.9, 11],
        # expert options
        mcc=1,
        remesh_size=0,
        no_check_surface=False,
        no_check_boundaries=False,
        no_qc=False,
        allow_ragged_surfaces=False,
        allow_ragged_trias=False,
        # deprecated options
        no_orient=False,
        # internal options (not set during parsing or class definition, but during args evaluation)
        map_values_integrate="none",
        map_values_select=None,
        map_values_interp="nearest",
        map_values_write_psol=True,
        map_values_write_mgh=True,
        map_values_write_annot=False,
    )

    if x in defaults.keys():
        return defaults[x]
    else:
        raise ValueError(x + " is not a valid key in the defaults dict")
