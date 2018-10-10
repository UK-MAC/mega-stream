# Changelog
All notable changes to this project will be documented in this file.

## [Unreleased]

### Added
- Timing output of each sweep (for first and last rank) and logical communication operations
- Model bandwidths using just compute time (rather than total time)
- mega-sweep3d learned `--prof` flag to skip (serial) output checking when profiling

### Changed
- README updated
- Ensured `psi*` arrays are copied/zeroed in parallel in mega-sweep


## [v2.1] - 2018-02-07

### Added
- 3D version of mega-sweep.
- Performance model of 2D mega-sweep cache bandwidth.
- Footprint output of more arrays.
- Output of population of all groups in mega-sweep.
- Check for zero values in angular flux in mega-sweep as a validation for missing iterations in the sweep.

### Changed
- Default 2D mega-sweep problem size.

### Fixed
- Output formatting of mega-sweep footprint.


## [v2.0] - 2017-11-06

### Added
- The mega-sweep mini-app has been created. Currently a Fortran only implementation.

### Changed
- Restructured directory to contain baseline and optimised mega-stream alongside mega-sweep.


## [v1.0] - 2017-06-30
New kernel, using 5 nested loops, is the major change of this release.

### Added
- Added matplotlib graph generation scripts for exploring problem space.
- Output of total execution time.
- Intel optimisations of newest kernel.

### Changed
- Kernel completely updated to use 5 nested loops, in C and Fortran versions.
- Kernel pulled into own subroutine.
- Memory bandwidth model updated for new kernel

### Fixed
- Support memory alignment on Apple platforms.

### Removed
- Verification check, as no longer applies to new kernel. New strategy is to print out the final value as a sanity check. It is too difficult for the new kernel to calculate the what expected result is a priori.
- OpenMP 4 versions.


## [v0.3] - 2016-12-20

### Added
- OPTIONS variable to Fortran makefile.

### Changed
- Updated Fortran version to match C version with multi-dimensional arrays.


## [v0.2.2] - 2016-12-15

### Fixed
- Memory OpenMP 4 bandwidth calculation


## [v0.2.1] - 2016-12-15

### Fixed
- Memory bandwidth calculation


## [v0.2] - 2016-12-14

### Added
- SIMD reduction of result array over smallest dimension.
- Option to run benchmark kernel many times, default 100.
- OPTIONS variable in Makefile for additional compiler flags.
- Output of array sizes and memory footprint.
- README file added.
- Swap command line option.

### Changed
- Make arrays multi-dimensional.
- Align the arrays using C11 alignment allocators.
- Updated results checking, with a tolerance.

## [v0.1] - 2016-12-13
Initial release of mega-stream benchmark.

