# Hippocampal Shape and Thickness Analysis


## Purpose:

This repository contains the Hipsta package, a collection of scripts for
hippocampal shape and thickness analysis as described in our recent [publication](https://doi.org/10.1016/j.neuroimage.2023.120182).


## Documentation:

Please see the documentation pages for a general overview, usage information and
examples, and output description. Brief usage information is also available [here](hipsta/doc/DOCUMENTATION.md).
Some suggestions for running the script can be found in the [tutorial](TUTORIAL.md).


## Current status:

The hippocampal shape and thickness analysis package is currently in its beta
stage, which means that it's open for testing, but may still contain unresolved
issues. Future changes with regard to the algorithms, interfaces, and package
structure are to be expected.


## Feedback:

Questions, suggestions, and feedback are welcome, and should preferably be
submitted as an [issue](https://github.com/Deep-MI/Hipsta/issues).


## Installation:

It is recommended to run this pipeline within its own virtual environment. A
virtual environment can, for example, be created using Python's `virtualenv`
command:

`virtualenv /path/to/a/new/directory/of/your/choice`

Activate the virtual environment as follows:

`source /path/to/a/new/directory/of/your/choice/bin/activate`

The package is available on `pypi.org`, and can be installed as follows 
(including all required dependencies):

`pip install hipsta`

Alternatively, the following code can be used to download this package from its
GitHub repository (this will create a 'Hipsta' directory within the current
working directory):

`git clone https://github.com/Deep-MI/Hipsta.git`


Use the following code to install the downloaded files as a Python packge (after
changing into the 'Hipsta' directory). It will also install all required
dependencies:

`pip install .`

The above steps are not necessary when running the [Docker](docker/Docker.md) or
[Singularity](singularity/Singularity.md) versions of the package.


## Requirements:

Unless using the [Docker](docker/Docker.md) or [Singularity](singularity/Singularity.md)
versions of the package, the following conditions need to be met for running an
analysis:

1. A FreeSurfer version (6.x or 7.x) must be sourced, i.e. FREESURFER_HOME must
exist as an environment variable and point to a valid FreeSurfer installation.

2. A hippocampal subfield segmentation created by FreeSurfer 7.11 or later
or the ASHS software. A custom segmentation is also permissible (some restrictions
and settings apply; see [Supported Segmentations](https://github.com/Deep-MI/Hipsta#supported-segmentations)).

3. Python 3.8 or higher including the lapy, numpy, scipy, nibabel, pyvista, and
pyacvd libraries, among others. See `requirements.txt` for a full list, and use
`pip install -r requirements.txt` to install.

4. The gmsh package (version 2.x; http://gmsh.info) must be installed. Can be
downloaded e.g. as binaries for [linux](https://gmsh.info/bin/Linux/gmsh-2.16.0-Linux64.tgz) or
[MacOSX](https://gmsh.info/bin/MacOSX/gmsh-2.16.0-MacOSX.dmg) . The 'gmsh' binary must
be on the $PATH:

    `export PATH=${PATH}:/path/to/gmsh-directory/bin`


## References:

Please cite the following publications if you use these scripts in your work:

- Diers, K., Baumeister, H., Jessen, F., Düzel, E., Berron, D., & Reuter, M. (2023). An automated, geometry-based method for hippocampal shape and thickness analysis. Neuroimage, 276:120182. doi: [10.1016/j.neuroimage.2023.120182](https://doi.org/10.1016/j.neuroimage.2023.120182).

Please also consider citing the these publications:

- Geuzaine, C., & Remacle, J.-F. (2009). Gmsh: a three-dimensional finite element mesh generator with built-in pre- and post-processing facilities. International Journal for Numerical Methods in Engineering, 79, 1309-1331.

- Andreux, M., Rodola, E., Aubry, M., & Cremers, D. (2014). Anisotropic Laplace-Beltrami operators for shape analysis. In European Conference on Computer Vision (pp. 299-312). Springer, Cham.

- Iglesias, J. E., Augustinack, J. C., Nguyen, K., Player, C. M., Player, A., Wright, M., ... & Fischl, B. (2015). A computational atlas of the hippocampal formation using ex vivo, ultra-high resolution MRI: application to adaptive segmentation of in vivo MRI. Neuroimage, 115, 117-137.

- Yushkevich, P. A., Pluta, J., Wang, H., Ding, S.L., Xie, L., Gertje, E., Mancuso, L., Kliot, D., Das, S. R., & Wolk, D.A. (2015). Automated Volumetry and Regional Thickness Analysis of Hippocampal Subfields and Medial Temporal Cortical Structures in Mild Cognitive Impairment. Human Brain Mapping, 36, 258-287.
