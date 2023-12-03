## Changes

### v0.4.0-beta

- fundamental revision incl. changes in command-line interface

### v0.3.0-beta

- introduced logging
- full lapy integration
- name change: subfieldDNA --> hippocampal shapetools (briefly: shapetools)
- support for ASHS segmentations
- mandatory `lut` command line argument
- dropped `hsflist` command line argument
- changes in default names for output folders and files

### v0.2.0-beta

- re-introduced remeshing as a preprocessing step; requires additional python
  packages `pyvista` and `pyacvd`. See README.md.
- added sheet selection argument for `mapValues()` function, plus output of
  interior and exterior surfaces.
- starting with this version, the sign of the mean curvature has changed to
  comply with Freesurfer conventions; what was negative is now positive and
  vice versa.

### v0.1.1-beta

- standalone functionality for `mapValues()`  function
- curvature output for interior and exterior surfaces

### v0.1.0-beta

- direction of `x`, `y`, `z` dimensions for output files in the `thickness`
  folder has changed:
  - `x` is now from medial to lateral (i.e. Presubiculum to CA3)
  - `y` is now from posterior to anterior
  - `z` is now from inferior to superior
- curvature and subfield label overlays as additional default outputs in the
  `thickness` folder
- introduction of the `mapValues()` function to project volume-based data onto
  the hippocampal mid-surface (needs to be run separately)
- introduction of `lapy` dependency
- `hsflist` command line argument is no longer necessary
- changes and improvements in various algorithms and parameters
