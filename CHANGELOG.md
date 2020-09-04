## 0.2.2 (September 04, 2020)
  - Enabling markdown support in PyPi, see https://stackoverflow.com/a/26737258
  - Merge branch 'master' of https://github.com/jacopo-chevallard/PySpecLines
  - Update README.md

## 0.2.1 (February 14, 2019)
  - Add possibility to choose verbose level
  - De-redden the spectrum, only then de-redshift it; download the dust maps the first time the de-reddening option is used
  - Removed dust maps from the repo, as the files were too large to allow uploading package to PyPI. The maps are downlaoded the first time that the de-reddening option is used
  - Added more extended explanations to README

## 0.2.0 (February 14, 2019)
  - Added required dependencies to correct for Galactic extinction
  - Added '--deredden' option
  - Added de-reddenning of spectrum before computing fluxes/EWs
  - Added FITS file containing Galactic extinction maps
  - Adding MANIFEST to install additional files required by the package

## 0.1.6 (February 13, 2019)
  - You're now correctly plotting the Gaussian fit obtained with PyMC, when PyMC is used (using the posterior median of the parameters)

## 0.1.5 (February 13, 2019)
  - Added example file for UV lines
  - Added option to show plot each time a line (or group of lines) are fitted
  - Added another JSON file example, showing how to fit multiple kinematic components. Note that components with the same initial width estimate will have their widths 'tied' together during the fitting
  - Can now fit multiple kinematic components along with multiple lines
  - Added further options: can decide how many iterations fit continuum-fit lines to perform; can add an estimate of the relative accuracy with which the continuum can be estimated, hence placing a lower 'floor' on the final errors; can plot with log y-axis

## 0.1.4 (February 12, 2019)
  - Bug fix: removing spurious commands

## 0.1.3 (February 12, 2019)
  - Removed unused modules
  - Bug fix: added missing module astropy.log; removed unused modules
  - Update README.md

## 0.1.2 (February 12, 2019)
  - Removed python 2.6 support
  - Cleaned

## 0.1.1 (February 12, 2019)
  - Added CHANGELOG
  - Testing Shippable CI
