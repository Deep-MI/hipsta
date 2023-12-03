# Hippocampal shape and thickness analysis


## Purpose:

The hipsta package is a collection of scripts for hippocampal shape and
thickness analysis. This page provides usage information as well as a technical
description of processing steps and outputs.

## Usage:

### As a command-line script:

The hipsta scripts can be run from the command line as follows:

```
----------------------------------------
Hippocampal shape and thickness analysis
----------------------------------------

usage: run_hipsta [--filename <filename>] [--hemi <lh|rh>] [--lut <freesurfer|ashs|filename>] [--outputdir <directory>] [--no-cleanup] [--no-crop] [--upsample] [--upsample-size <float> <float> <float>]
                  [--no-merge-molecular-layer] [--automask-head] [--automask-tail] [--automask-head-margin <int>] [--automask-tail-margin <int>] [--no-gauss-filter]
                  [--gauss-filter-size <float> <float>] [--long-filter] [--long-filter-size <int>] [--no-close-mask] [--mca <mri_mc|skimage>] [--remesh] [--smooth <int>] [--cut-range <float> <float>]
                  [--aniso-alpha <float> [<float> ...]] [--aniso-smooth <int>] [--thickness-grid <float> <float> <float> <float> <float> <float> <float> <float> <float>] [--help] [--more-help]
                  [--version]

This program conducts a thickness analysis of the hippocampus, based on a FreeSurfer, ASHS, or custom hippocampal subfield segmentation.

Required arguments:
  --filename <filename>
                        Filename of a segmentation file.
  --hemi <lh|rh>        Hemisphere. Either 'lh' or 'rh'.
  --lut <freesurfer|ashs|filename>
                        Look-up table: a text file with numeric and verbal segmentation labels. 'freesurfer' and 'ashs' are keywords for built-in tables.
  --outputdir <directory>
                        Directory where the results will be written.

Optional arguments:
  --no-cleanup          Do not remove files that may be useful for diagnostic or debugging purposes, but are not necessary otherwise.
  --no-crop             Do not crop image.
  --upsample            Upsample to the smallest voxel edge length.
  --upsample-size <float> <float> <float>
                        Upsampling factors. Should be between 0 and 1. If all zeros, upsample to the smallest voxel edge length. Default: 0 0 0
  --no-merge-molecular-layer
                        Do not merge molecular layer (only applicable for FreeSurfer segmentations).
  --automask-head       Automated boundary detection for hippocampal head.
  --automask-tail       Automated boundary detection for hippocampal tail.
  --automask-head-margin <int>
                        Margin for automated boundary detection for hippocampal head. Default: 0
  --automask-tail-margin <int>
                        Margin for automated boundary detection for hippocampal tail. Default: 0
  --no-gauss-filter     Do not apply gaussian filter.
  --gauss-filter-size <float> <float>
                        Filter width and threshold for gaussian filtering. Default: 1 50.
  --long-filter         Apply filter along longitudinal axis, i.e. attempt to create smooth transitions between slices.
  --long-filter-size <int>
                        Size of longitudinal filter. Default: 5
  --no-close-mask       Do not apply closing operation to mask, i.e. do not attempt to close small holes.
  --mca <mri_mc|skimage>
                        Type of marching-cube algorithm. Either 'mri_mc' or 'skimage'. Default: 'skimage'
  --remesh              Apply remeshing operation to surface, i.e. create a regular, evenly spaced surface grid.
  --smooth <int>        Mesh smoothing iterations. Default: 5
  --cut-range <float> <float>
                        Range for tetrahedral boundary cutting. Default: -0.975, 0.975
  --aniso-alpha <float> [<float> ...]
                        Anisotropy parameter(s). Can be one or two or two numbers. Default: 40
  --aniso-smooth <int>  Anisotropy smoothing iterations. Default: 3
  --thickness-grid <float> <float> <float> <float> <float> <float> <float> <float> <float>
                        Extent and resolution of the grid used for thickness computation; three lists of three numbers: negative extent of x axis, positive extent of x axis, resolution on x axis.
                        Repeat for the y and z axes. Default: -0.9 0.9 41 -0.975 0.975 21 -0.9 0.9 11

Getting help:
  --help                Display this help message and exit
  --more-help           Display extensive help message and exit
  --version             Display version number and exit
  ```

### As a python package:

The hipsta scripts can also be run within a Python environment:

```python
import hipsta
hipsta.run_hipsta(filename="my_filename.mgz", hemi="lh", lut="freesurfer", outputdir="my_output_directory")
```

A full list of available arguments for the `run_hipsta` function can be obtained with:

```python
import hipsta
help(hipsta.run_hipsta)
```


## Description:

This script performs the following major processing steps:

1. process image
    1. convert format
    2. cropping (optional)
    3. upsampling (optional)
2. process labels
    1. automated identification of head and tail boundaries (optional)
    2. extract and merge subfield labels
    3. merge molecular layer (optional)
3. process mask
    1. binarize
    2. apply gaussian filter operations (optional)
    3. apply longitudinal filter operations (optional)
    3. fill holes using dilation and erosion (optional)
4. create surface
    1. extract initial surface using marching cube algorithm
    2. remesh surface (optional)
    3. smooth surface
    4. create QC plots
5. create tetrahedral mesh from triangular surface
6. create label files for tetra mesh
7. cut open tetrahedral mesh at anterior and posterior ends
    1. remove boundary mask
    2. mesh cutting
    3. check boundaries
8. create cube parametrization
    1. compute cube parametrization
    2. create QC plots
9. compute thickness and curvature
    1. compute thickness and curvature
    2. create QC plots
10. map subfield values (and other volume-based data, optional)
11. create supplementary files for visualization

## Outputs:

1. Thickness and curvature estimates

The primary output of the hippocampal shape and thickness package are surface
files and associated thickness values in the 'tickness' folder.

The thickness values will be stored in csv tables:

- <lh|rh>.grid-segments-<x|y|z>.csv

Here, x corresponds to the medial-->lateral dimension, y to the
posterior-->anterior dimension, and z to the exterior-->interior dimension.

The thickness values are therefore in the z files, whereas the x and y files
contain length estimates in the other two directions.

The thickness values will also be stored as mgh and psol overlay files that
can be overlaid onto the mid-surface vtk file.

- <lh|rh>.mid-surface.vtk
- <lh|rh>.mid-surface.thickness.<mgh|psol>

It is also possible to overlay the projected subfield boundary files onto
the midsurface:

- <lh|rh>.mid-surface_hsf.<mgh|psol>

Besides the thickness values, mean and gaussian curvature estimates are
provided for interior, mid, and exterior surfaces:

- <lh|rh>.<int|mid|ext>-surface.vtk
- <lh|rh>.<int|mid|ext>-surface.<mean|gauss>-curv.csv
- <lh|rh>.<int|mid|ext>-surface.<mean|gauss>-curv.mgh
- <lh|rh>.<int|mid|ext>-surface.<mean|gauss>-curv.psol

The other files within the 'thickness' folder files are intermediate files that
were created / used during the thickness computation:

- <lh|rh>.<int|mid|ext>-surface.csv
- <lh|rh>.mid-surface.csv
- <lh|rh>.grid-lines.csv
- <lh|rh>.grid-lines-<x|y|z>.vtk
- <lh|rh>.mid-surface_hsf.csv
- <lh|rh>.mid-surface_hsf-bnd.vtk

2. Logfile and intermediate volume and surface files

The logfile can be used to check if the processing finished without errors:

- logfile.txt

Also the intermediate volume and surface files can be useful for quality control or troubleshooting.

- Result of processing step 1 ('process image'): <lh|rh>.image.mgz
- Result of processing step 2 ('process labels'): <lh|rh>.labels.mgz
- Result of processing step 3 ('process mask'): <lh|rh>.mask.mgz
- Result of processing step 4 ('create surface'): <lh|rh>.surf.vtk
- Result of processing step 5 ('create tetrahedral mesh'): <lh|rh>.tetra.vtk
- Result of processing step 7 ('cut open tetrahedral mesh'): <lh|rh>.cut.vtk

3. A set of sub-folders with intermediate results:

    1. An 'image' folder with intermediate results, primarily from basic image
    processing.

    2. A 'labels' folder with intermediate results, primarily label files.

    3. A 'mask' folder with intermediate results, primarily files created
    during creation and post-processing of binary masks.

    4. A 'surface' folder with intermediate results, primarily files created
    during surface construction.

    5. A 'tetra-mesh' folder with intermediate results, primarily files
    created during tetrahedral mesh construction.

    6. A 'tetra-labels' folder with intermediate results, primarily files
    created during tetrahedral mesh construction.

    7. A 'tetra-cut' folder with intermediate results, primarily files
    created during mesh cutting.

    8. A 'tetra-cube' folder with intermediate results, primarily files
    created during cube parametrization.

    9. A 'thickness' folder with thickness overlays and mid-surfaces.

    10. A 'qc' folder with QC plots.


## Supported segmentations:

- The `freesurfer` segmentation should work as-is.
- If using the `ashs` segmentation, additional labels for the hippocampal head
  (label value: 20) and tail (label value: 5) labels are required. Use `--lut ashs`.
- If using a custom look-up table (`--lut <filename>`), the expected format for
  the file is: `<numeric ID> <name> <R> <G> <B> <A>`. `R`, `G`, `B`, `A` are
  numerical values for RGB colors and transparency (alpha) and will be ignored.
  For example, `236 Subiculum 255 0 0 0`. Each line may only contain a single
  anatomical structure. Do not use a header line. The following labels need to
  be present in the look-up table: `presubiculum`,  `subiculum`, `ca1`, `ca2`,
  `ca3`, `ca4`, `head`, and `tail`. Additional labels may be present, but will
  be ignored. Lines starting with `#` will be ignored (and can be used for
  comments).
- Multiple substructures can have the same numeric ID, e.g. `presubiculum` and
  `subiculum` can have the same numeric ID if these substructures are not
  distinguished in the segmentation. `head` and `tail` can have multiple labels
  if these are distinguished in the segmentation, but should be combined for
  the processing with this package.


## Requirements:

1. A FreeSurfer version (6.x or 7.x) must be sourced, i.e. FREESURFER_HOME must
   exist as an environment variable and point to a valid FreeSurfer installation.

2. A hippocampal subfield segmentation created by FreeSurfer 7.11 or later
   or the ASHS software. A custom segmentation is also permissible (some restrictions
   and settings apply; see [Supported Segmentations](#supported-segmentations)).

3. Python 3.8 or higher including the lapy, numpy, scipy, nibabel, pyvista, and
   pyacvd libraries. See `requirements.txt` for a full list.

4. The gmsh package (version 2.x; http://gmsh.info) must be installed. Can be
   downloaded from e.g. https://gmsh.info/bin/Linux/gmsh-2.16.0-Linux64.tgz
   or https://gmsh.info/bin/MacOSX/gmsh-2.16.0-MacOSX.dmg. The 'gmsh' binary
   must also be on the $PATH, i.e `export PATH=${PATH}:/path/to/my/gmsh`


## References:

Please cite the following publications if you use these scripts in your work:

- Diers, K., Baumeister, H., Jessen, F., Düzel, E., Berron, D., & Reuter, M. (2023).
  An automated, geometry-based method for hippocampal shape and thickness analysis.
  Neuroimage, 276:120182. doi: 10.1016/j.neuroimage.2023.120182.

Please also consider citing the these publications:

- Geuzaine, C., & Remacle, J.-F. (2009). Gmsh: a three-dimensional finite element mesh
  generator with built-in pre- and post-processing facilities. International Journal
  for Numerical Methods in Engineering, 79, 1309-1331.

- Andreux, M., Rodola, E., Aubry, M., & Cremers, D. (2014). Anisotropic Laplace-Beltrami
  operators for shape analysis. In European Conference on Computer Vision (pp. 299-312).
  Springer, Cham.

- Iglesias, J. E., Augustinack, J. C., Nguyen, K., Player, C. M., Player, A., Wright, M.,
  ... & Fischl, B. (2015). A computational atlas of the hippocampal formation using ex vivo,
  ultra-high resolution MRI: application to adaptive segmentation of in vivo MRI. Neuroimage,
  115, 117-137.

- Yushkevich, P. A., Pluta, J., Wang, H., Ding, S.L., Xie, L., Gertje, E., Mancuso, L.,
  Kliot, D., Das, S. R., & Wolk, D.A. (2015). Automated Volumetry and Regional Thickness
  Analysis of Hippocampal Subfields and Medial Temporal Cortical Structures in Mild Cognitive
  Impairment. Human Brain Mapping, 36, 258-287.
