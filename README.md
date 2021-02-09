# shapetools.py


## Purpose:

This is a script for the creation and analysis of hippocampal surfaces.


## Description:

This script performs the following major processing steps:

1. create labels
   1. extract HSF labels from label image
   2. merge individual labels
2. attach molecular layer to neighboring regions
3. create mask
   1. fill holes using dilation and erosion
   2. create final mask
4. create surface
   1. create initial surface using marching cube algorithm
   2. topology fixing (optional)
   3. smooth (fixed) surface
5. create tetrahedral mesh from triangular surface
6. create label files for tetra mesh
7. remove boundary mask (preprocessing for mesh cutting)
8. cut open tetrahedral mesh at anterior and posterior ends
9. create cube parametrization
10. compute thickness and curvature
11. map subfield values (and other volume-based data)


## Usage:

`python3 shapetools.py --subjid <SUBJID> --hemi <hemi> --lut <lookup-table>
  [--subjdir <SUBJDIR>] [--outputdir <OUTDIR>] [further options]`

**Mandatory arguments:**   | **Description**
---------------------------|----------------------------------------------------
`--subjid`                 | Subject ID. Must match a subdirectory in Freesurfer's SUBJECTS_DIR or in a manually specified `--subjdir`
`--hemi`                   | Hemisphere. Either 'lh' or 'rh'.
`--lut`                    | Look-up table. Either 'fs711' or 'ashs' or a valid filename.

**Optional arguments:**    | **Description**
---------------------------|----------------------------------------------------
`--subjdir`                | Subjects directory. This directory must contain a folder that matches the subject ID and contains Freesurfer output including the hippocampal subfield segmentation. If not given,  Freesurfer's SUBJECTS_DIR environment variable will be used.
`--outputdir`              | Directory where the results will be written. If not given, a subfolder within each subject's directory will be created.
`--sfx`                    | Subfield segmentation suffix. These are different for the label files created with different versions and settings of the subfield segmentation (default: hippoAmygLabels-T1.v20).
`--no-topological-fixing`  | Disable topological fixing. This is a simple switch and no further parameter is required; if it is given, use Freesurfer's topology fixing program to refine and fix initial surfaces.
`--no-cleanup`             | Do not delete files that may be useful for diagnostic or debugging purposes, but are not strictly necessary otherwise.
`--help`                   | Display this help and exit.
`--version`                | Display version number and exit.


## Examples:

`python3 shapetools.py --subjid my_subject --hemi lh --lut fs711 --subjdir /my/subjects/directory --outputdir /my/output/directory --sfx hippoAmygLabels-T1.v20`

`python3 shapetools.py --subjid my_subject --hemi lh --lut ashs  --subjdir /my/subjects/directory --suffix my_suffix --no-topological-fixing`

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
ml             | merged molecular layer
de             | dilation and erosion
mc             | marching cube algorithm
tf             | topological fixing
rs             | remeshing and smoothing
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
1    | \<none\>    | \<none\> | labels         | files with volume labels
2    | HSFLABEL_01 | ml     | merge-ml       | files created during the merging of the molecular layer
3    | HSFLABEL_02 | de     | mask           | files created during creation and post-processing of binary masks
4a   | HSFLABEL_03 | mc     | \<none\>       | \<none\>
4b   | HSFLABEL_05 | tf     | fixed-surface  | files created during topological fixing (optional)
4c   | HSFLABEL_06 | rs     | \<none\>       | \<none\>
5    | HSFLABEL_07 | tet    | \<none\>       | \<none\>
6    | \<none\>    | \<none\> | tetra-label    | files created during tetrahedral mesh construction
7    | \<none\>    | \<none\> | tetra-cut      | files created during preprocessing for tetrahedral mesh cutting
8    | HSFLABEL_09 | cut    | tetra-cut      | files created during tetrahedral mesh cutting
9    | \<none\>    | bnd, open, uvw | tetra-cube     | files created during cube parametrization
10+11| \<none\>    | \<none\> | thickness      | folder with thickness overlays and mid-surface files


## Requirements:

1. A Freesurfer version (6.x or 7.x) must be sourced, i.e. FREESURFER_HOME must
exist as an environment variable and point to a valid Freesurfer installation.

2. A hippocampal subfield segmentation created by Freesurfer 7.11 (development
versions after 6.0 will also work) must be present in the 'mri' subdirectory of
the Freesurfer 6.0 or newer processed SUBJECTS_DIR/<subjid> directory.
Alternatively, a subfield segmentation obtained from the ASHS IKND 7T Young
Adults is also permissible (some restrictions on filenames and folders apply;
see CHANGES.md).

3. Python 3.5 or higher including the lapy, numpy, scipy, and nibabel
libraries.

4. The gmsh package (verson 2.x; http://gmsh.info) must be installed. Can be
downloaded from e.g. http://gmsh.info/bin/Linux/older/gmsh-2.16.0-Linux64.tgz
or http://gmsh.info/bin/MacOSX/older/gmsh-2.16.0-MacOSX.dmg. The 'gmsh' binary
must also be on the $PATH, i.e `export PATH=${PATH}:/path/to/my/gmsh`


## References:

Please cite the following publications if you use these scripts in your work:

- Reuter, M., Wolter, F.-E., Shenton, M., & Niethammer, M. (2009). Laplace-Beltrami Eigenvalues and Topological Features of Eigenfunctions for Statistical Shape Analysis. Computer-Aided Design, 41, 739-755.

- Geuzaine, C., & Remacle, J.-F. (2009). Gmsh: a three-dimensional finite element mesh generator with built-in pre- and post-processing facilities. International Journal for Numerical Methods in Engineering, 79, 1309-1331.

- Andreux, M., Rodola, E., Aubry, M., & Cremers, D. (2014). Anisotropic Laplace-Beltrami operators for shape analysis. In European Conference on Computer Vision (pp. 299-312). Springer, Cham.

- Iglesias, J. E., Augustinack, J. C., Nguyen, K., Player, C. M., Player, A., Wright, M., ... & Fischl, B. (2015). A computational atlas of the hippocampal formation using ex vivo, ultra-high resolution MRI: application to adaptive segmentation of in vivo MRI. Neuroimage, 115, 117-137.

- Yushkevich, P. A., Pluta, J., Wang, H., Ding, S.L., Xie, L., Gertje, E., Mancuso, L., Kliot, D., Das, S. R., & Wolk, D.A. (2015). Automated Volumetry and Regional Thickness Analysis of Hippocampal Subfields and Medial Temporal Cortical Structures in Mild Cognitive Impairment. Human Brain Mapping, 36, 258-287.
