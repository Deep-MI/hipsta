# shapetools.py

## Changes

### v0.2.0-beta

- name change: subfieldDNA --> hippocampal shapetools (briefly: shapetools)
- support for ASHS_IKND_7T segmentations
- mandatory `lut` command line argument
- dropped `hsflist` command line argument
- changes in default names for output folders and files

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

## Notes

- If using the `ashs` segmentation, additional labels for the hippocampal head
  (255) and tail (254) labels are required. The expected file / folder
  structure is: <SUBJECTS_DIR>/SUBJECT_ID>/mri/<lh|rh>.<SUFFIX>.mgz
  Use `--lut ashs`. Specify `--hemi` and `--suffix`. Do not use topological
  fixing (`--no-topological-fixing`).
- Multiple substructures can have same label, but same substructure cannot have
  multiple labels.
