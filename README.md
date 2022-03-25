# shapetools.py


## Purpose:

This is a script for the creation and analysis of hippocampal surfaces.


## Description:

This script performs the following major processing steps:

0. preprocess image (cropping, upsampling, automask)
1. create labels
   1. extract HSF labels from label image
   2. merge individual labels
2. attach molecular layer to neighboring regions (if present)
3. create mask
   1. fill holes using dilation and erosion
   2. create final mask
4. create surface
   1. create initial surface using marching cube algorithm
   2. topology fixing (optional)
   3. smooth (fixed) surface
   4. remesh surface (optional)
5. create tetrahedral mesh from triangular surface
6. create label files for tetra mesh
7. remove boundary mask (preprocessing for mesh cutting)
8. cut open tetrahedral mesh at anterior and posterior ends
9. create cube parametrization
10. compute thickness and curvature
11. map subfield values (and other volume-based data)
12. create supplementary files for visualization


## Usage:

`python3 shapetools.py --filename <filename> --hemi <hemi> --lut <lookup-table> [--outputdir <OUTDIR>] [further options]`

**Mandatory arguments:**   | **Description**
---------------------------|----------------------------------------------------
`--filename`               | A segmentation file.
`--hemi`                   | Hemisphere. Either 'lh' or 'rh'.
`--lut`                    | Look-up table. A valid filename or one of the following: 'fs711', 'ashs', 'ashs-ca2ca3', 'ashs-noca3', 'ashs-ctx', 'ashs-ctx-nohc', 'ashs-ent'.

**Optional arguments:**    | **Description**
---------------------------|----------------------------------------------------
`--outputdir`              | Directory where the results will be written. If not given, a subfolder within each subject's directory will be created.
`--no-cleanup`             | Do not delete files that may be useful for diagnostic or debugging purposes, but are not strictly necessary otherwise.

**Getting help:**    | **Description**
---------------------------|----------------------------------------------------
`--help`                   | Display this help and exit.
`--more-help`              | Display extensive help and exit.
`--version`                | Display version number and exit.

**Experimental arguments:**| **Description**
---------------------------|----------------------------------------------------
`--no-crop`                | Do not crop image.
`--upsample [ <float> <float> <float> ]` | A list of parameters to fine-tune image upsampling: specify three multiplicative factors to determine the final voxel size. If none are given, resample to the smallest voxel size.
`--automated-head`         | Automatically identify boundary towards hippocampal head.
`--automated-tail`         | Automatically identify boundary towards hippocampal tail.
`--margin-head <int>`      | Margin for automated head identification.
`--margin-tail <int>`      | Margin for automated tail identification.
`--dilation-erosion <int> <int>` | Settings for dilation-erosion procedure. Specify two numbers for width of dilation and erosion. Default is `0 0`.
`--no-filter`              | Do not filter image.
`--filter <int> <int>`     | A list of parameters to fine-tune image filtering. Specify two numbers for filter width and threshold. Default is `0.5 50`
`--long-filter`            | Filter image along longitudinal axis, i.e. attempt to create smooth transitions between slices.
`--close-mask [<string>]`  | Apply closing operation to mask, i.e. attempt to close small holes. Default is `regular`, alternative is `experimental`.
`--mca <string>`           | Marching-cube algorithm. Either `mri_mc` (default) or `mri_tessellate`.
`--mcc <integer>`          | Marching-cube connectivity. Only used for `mri_mc` algorithm. Default is `1`. See `mri_mc --help` for details.
`--topological-fixing`     | Use FreeSurfer's topology fixing program to refine and fix initial surfaces. Can only be used with FreeSurfer-processed data. Expects two files as input, brain.mgz and wm.mgz.
`--smooth <integer>`       | Mesh smoothing iterations. Default is `5`.
`--remesh [ <int> ]`       | Switch on remeshing, i.e. create a regular, evenly spaced surface grid. If a number is given, resample to that
`--no-check-surface`       | Do not check surface and proceed even if holes are present.
`--no-check-boundaries`    | Do not check boundaries and proceed if there are less / more than two continuous boundary loops.
`--no-qc`                  | Do not perform QC.
`--cut-params <float> <float>` | A list of parameters to fine-tune the cut operation. Default cutting range is `[-0.975, 0.975]`.
`--aniso-alpha <float> [ <float> ]` | Anisotropy parameters. Specify either one or two numbers.
`--aniso-smooth <int>`     | Anisotropy smoothing. Specify a number to control the number of smoothing iterations. Default is `10`.
`--thickness-params <...>` | A list of parameters to fine-tune the thickness computation. Specify three lists of three numbers: negative extent of x axis, positive extent of x axis, resolution on x axis. Repeat for the y and z axes. Default values are `-0.9 0.9 41 -0.975 0.975 21 -0.9 0.9 11`.
`--allow-ragged`           | Allow ragged mid-surfaces.
`--allow-ragged-triangles` | Allow triangles for ragged mid-surfaces.


## Examples:

`python3 shapetools.py --filename /my/segmentation/file --hemi lh --lut fs711 --topological-fixing /path/to/brain.mgz /path/to/wm.mgz --outputdir /my/output/directory`

`python3 shapetools.py --filename /my/segmentation/file --hemi lh --lut ashs --outputdir /my/output/directory`


## Outputs:

1. The primary output are thickness values that will be computed and stored as
   CSV or PSOL or MGH overlay files within the 'thickness' folder. The PSOL or
   MGH files can be overlaid onto the mid-surface vtk file that is also found
   within the 'thickness' folder.

2. Intermediate volume and surface files, each prefixed with the particular
   operation performed. The basic pattern for the filenames is:

   `<hemisphere>.<prefixes>.<hsf-labels>.<suffixes>`

Abbreviation   | Intpretation
-------------- | ----------------------------
*Hemisphere*   |
lh             | left hemisphere
rh             | right hemisphere
*Prefixes*     |
crop           | crop image
ups            | upsample image
am             | auto-mask for head and tail
ml             | merged molecular layer
de             | dilation and erosion
filt           | gaussian filtering and thresholding
cm             | close mask
lf             | longitudinal filter
mc             | marching cube algorithm
tf             | topological fixing
sm             | smoothing
rs             | remeshing
tet            | tetrahedral mesh
bnd            | boundary mesh
cut            | mesh cut at anterior/posterior ends
rm             | remove free vertices
open           | removed anterior/posterior ends to create an open mesh
uvw            | cube parametrization (xyz --> uvw)
*HSF-Labels*   |
234, 236, etc. | several numerics corresponding to the hippocampal subfield segmentation look-up table
*Suffixes*     |
assigned       | molecular layer split into subregions
merged         | molecular layer subregions merged with other regions

3. Several folders with intermediate results

Step | HSF label   | Prefix | Folder         | Contents
-----|-------------|--------|----------------|-----------------------------
1    | \<none\>    | \<none\> | labels       | files with volume labels
2    | HSFLABEL_01 | ml     | merge-ml       | files created during the merging of the molecular layer
3    | HSFLABEL_02 | de     | mask           | files created during creation and post-processing of binary masks
4a   | HSFLABEL_03 | mc     | \<none\>       | \<none\>
4b   | HSFLABEL_05 | tf     | fixed-surface  | files created during topological fixing (optional)
4c   | HSFLABEL_06 | rs     | \<none\>       | \<none\>
5    | HSFLABEL_07 | tet    | \<none\>       | \<none\>
6    | \<none\>    | \<none\> | tetra-label  | files created during tetrahedral mesh construction
7    | \<none\>    | \<none\> | tetra-cut    | files created during preprocessing for tetrahedral mesh cutting
8    | HSFLABEL_09 | cut    | tetra-cut      | files created during tetrahedral mesh cutting
9    | \<none\>    | bnd, open, uvw | tetra-cube | files created during cube parametrization
10+11| \<none\>    | \<none\> | thickness    | folder with thickness overlays and mid-surface files


## Custom segmentations:

- If using the `ashs` segmentation, additional labels for the hippocampal head
  (label value: 20) and tail (label value: 5) labels are required. Use `--lut ashs`.
- Topological fixing (`--topological-fixing <filename1> <filename2>`) should
  not be used unless working with FreeSurfer-processed data. `<filename1>` is
  the `brain.mgz` file and `<filename2>` is the `wm.mgz` file, both found within
  the `mri` subdirectory of an individual FreeSurfer output directory.
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
  the processing within the hippocampal shapetools.


## Installation:

Use the following code to download this package from its GitHub repository:

`git clone https://github.com/reuter-lab/hippocampal-shape-tools.git`

You can use the following code to install the required Python depencies into
your local python package directory:

`pip install --user -r requirements.txt`


## Requirements:

1. A FreeSurfer version (6.x or 7.x) must be sourced, i.e. FREESURFER_HOME must
exist as an environment variable and point to a valid FreeSurfer installation.

2. A hippocampal subfield segmentation created by FreeSurfer 7.11 (development
versions after 6.0 will also work) or a subfield segmentation obtained from the
ASHS IKND 7T Young Adults atlas. A custom segmentation is also permissible
(some restrictions and settings apply; see `Custom Segmentations`).

3. Python 3.5 or higher including the lapy, numpy, scipy, nibabel, pyvista, and
pyacvd libraries. See `requirements.txt` for a full list. The lapy package can
be obtained from https://github.com/Deep-MI/LaPy.

4. The gmsh package (verson 2.x; http://gmsh.info) must be installed. Can be
downloaded from e.g. http://gmsh.info/bin/Linux/older/gmsh-2.16.0-Linux64.tgz
or http://gmsh.info/bin/MacOSX/older/gmsh-2.16.0-MacOSX.dmg. The 'gmsh' binary
must also be on the $PATH, i.e `export PATH=${PATH}:/path/to/my/gmsh`

5. The PYTHONPATH environment variable should include the toolbox directory,
e.g. `export PYTHONPATH=${PYTHONPATH}:/path/to/hippocampal-shapetools`.


## References:

Please cite the following publications if you use these scripts in your work:

- Reuter, M., Wolter, F.-E., Shenton, M., & Niethammer, M. (2009). Laplace-Beltrami Eigenvalues and Topological Features of Eigenfunctions for Statistical Shape Analysis. Computer-Aided Design, 41, 739-755.

- Geuzaine, C., & Remacle, J.-F. (2009). Gmsh: a three-dimensional finite element mesh generator with built-in pre- and post-processing facilities. International Journal for Numerical Methods in Engineering, 79, 1309-1331.

- Andreux, M., Rodola, E., Aubry, M., & Cremers, D. (2014). Anisotropic Laplace-Beltrami operators for shape analysis. In European Conference on Computer Vision (pp. 299-312). Springer, Cham.

- Iglesias, J. E., Augustinack, J. C., Nguyen, K., Player, C. M., Player, A., Wright, M., ... & Fischl, B. (2015). A computational atlas of the hippocampal formation using ex vivo, ultra-high resolution MRI: application to adaptive segmentation of in vivo MRI. Neuroimage, 115, 117-137.

- Yushkevich, P. A., Pluta, J., Wang, H., Ding, S.L., Xie, L., Gertje, E., Mancuso, L., Kliot, D., Das, S. R., & Wolk, D.A. (2015). Automated Volumetry and Regional Thickness Analysis of Hippocampal Subfields and Medial Temporal Cortical Structures in Mild Cognitive Impairment. Human Brain Mapping, 36, 258-287.
