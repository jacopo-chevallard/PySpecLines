# PySpecLines

Some useful scripts to compute integrated line fluxes and equivalent widths using 3 different methods:

- numerical (trapezoidal) integration
- Gaussian fit using a "standard" Levenberg-Marquardt algorithm
- Gaussian fit using an MCMC algorithm


## Examples

The script **requires** one or more input FITS files containing the calibrated spectra and the following columns:
- ``wl``, wavelength array in ``ang``
- ``flux``, flux array in ``erg s^-1 cm^-2 ang^-1``
- ``err``, error array in  ``erg s^-1 cm^-2 ang^-1``

If the spectrum is provided in the observed frame, then you must provide a ``REDSHIFT`` keyword in the FITS header containing the object redshift.

- Compute the fluxes and EWs using numerical integration
  ``./get_integrated_fluxes.py --file my_spectrum.fits --json-file  emission_lines_EWs_config.json``

- Compute the fluxes and EWs using Gaussian fit (Levenberg-Marquardt)
  ``./get_integrated_fluxes.py --file my_spectrum.fits --json-file  emission_lines_EWs_config.json --gaussian-fit``

- Compute the fluxes and EWs using Gaussian fit (MCMC)
  ``./get_integrated_fluxes.py --file my_spectrum.fits --json-file  emission_lines_EWs_config.json --gaussian-fit --use-PyMC --MCMC-samples 5000``

## Requirements

- numpy
- astropy
- matplotlib
- pyspeckit
