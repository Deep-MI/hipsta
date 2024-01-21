"""
This module provides atlas information for the hippocampal shape and thickness 
analysis package.

"""

import logging
import os

import pandas

# ==============================================================================
# LOGGING

LOGGER = logging.getLogger(__name__)

# ==============================================================================
# FUNCTIONS


def get_atlases(lut):
    if lut == "freesurfer":
        LOGGER.info("Found internal, modified look-up table for FreeSurfer.")

        LUTLABEL = ["tail", "head", "presubiculum", "subiculum", "ca1", "ca2", "ca3", "ca4", "dg", "ml"]

        LUTINDEX = [226, [233, 235, 237, 239, 245], 234, 236, 238, 240, 240, 242, 244, 246]

        LUTDICT = dict(zip(LUTLABEL, LUTINDEX))

        HSFLIST = [234, 236, 238, 240, 246]

    elif lut == "freesurfer-no_ml":
        LOGGER.info("Found internal, modified look-up table for FreeSurfer.")

        LUTLABEL = ["tail", "head", "presubiculum", "subiculum", "ca1", "ca2", "ca3", "ca4", "dg"]

        LUTINDEX = [226, [233, 235, 237, 239, 245], 234, 236, 238, 240, 240, 242, 244]

        LUTDICT = dict(zip(LUTLABEL, LUTINDEX))

        HSFLIST = [234, 236, 238, 240]

    elif lut == "ashs-penn_abc_3t_t2":
        LOGGER.info("Found internal, modified look-up table the Penn ABC-3T ASHS Atlas for T2-weighted MRI.")

        LUTLABEL = [
            "ca1",
            "ca2",
            "ca3",
            "ca4",
            "dg",
            "tail_orig",
            "subiculum",
            "presubiculum",
            "entorhinal",
            "ba35",
            "ba36",
            "parahippocampal",
            "head",
            "tail",
        ]

        LUTINDEX = [1, 2, 4, 3, 3, 5, 8, 8, 9, 10, 11, 12, 20, 5]

        LUTDICT = dict(zip(LUTLABEL, LUTINDEX))

        HSFLIST = [8, 1, 2, 4]

    elif lut == "ashs-umcutrecht_7t":
        LOGGER.info("Found internal, modified look-up table for ASHS UMC Utrecht 7T atlas.")

        LUTLABEL = [
            "entorhinal",
            "subiculum",
            "presubiculum",
            "ca1",
            "ca2",
            "dg",
            "ca3",
            "ca4",
            "cyst",
            "tail",
            "head",
        ]

        LUTINDEX = [1, 2, 2, 3, 4, 5, 6, 5, 7, 8, 20]

        LUTDICT = dict(zip(LUTLABEL, LUTINDEX))

        HSFLIST = [2, 3, 4, 6]

    elif os.path.isfile(lut):
        LOGGER.info("Found look-up table " + lut)

        lut = pandas.read_csv(
            lut,
            sep=" ",
            comment="#",
            header=None,
            skipinitialspace=True,
            skip_blank_lines=True,
            error_bad_lines=False,
            warn_bad_lines=True,
        )

        LUTDICT = dict(zip(lut[0], lut[1]))

        HSFLIST = list(lut[1])

    else:
        LUTDICT = dict()

        HSFLIST = []

    # add entries for tetra-labels

    LUTDICT["jointtail"] = 226
    LUTDICT["jointhead"] = 232
    LUTDICT["bndtail"] = 2260
    LUTDICT["bndhead"] = 2320
    LUTDICT["bndca4"] = 2420

    # return

    return LUTDICT, HSFLIST
