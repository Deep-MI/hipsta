"""
This module provides configuration defaults for the hippocampal shape and 
thickness analysis package.

"""

def get_defaults(x):

    defaults = dict(
        cleanup = False,
        nocrop = False,
        upsample = None,
        noml = False,
        automaskhead = False,
        automasktail = False,
        automaskheadmargin = [0],
        automasktailmargin = [0],
        nofilter = False,
        nogaussfilter = False,
        gaussfilter_size = [1, 50],
        longfilter = False,
        longfilter_size = [5],
        noclosemask = False,
        mca = None,
        mcc = None,
        smooth = None,
        remesh = None,
        nochecksurface = False,
        nocheckboundaries = False,
        noqc = False,
        cutrange = None,
        anisoAlpha = [40],
        anisoSmooth = None,
        thxyz = None,
        allowRagged = False,
        allowRaggedTriangles = False,
        logfiledir = None,
        skiporient = False,
    )

    if x in defaults.keys():
        return defaults[x]
    else:
        raise ValueError(x + " is not a valid key in the defaults dict")