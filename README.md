
# PySpecLines

Some useful scripts to compute integrated line fluxes and equivalent widths from galaxy spectra using 3 different methods:

- numerical (trapezoidal) integration
- Gaussian fit using a "standard" Levenberg-Marquardt algorithm
- Gaussian fit using an MCMC algorithm

## Installing the package

To install ``PySpecLines`` you can use ``pip``, which will take care of installing the required dependencies as well
```
pip install pyspeclines
```

To upgrade to the latest available version you can run
```
pip install pyspeclines --upgrade
```

## Data format

### FITS binary table

``PySpecLines`` works with FITS tables. Each spectrum should be contained in a separate FITS file, which **must** contain the following columns:

- ``wl``, wavelength array in ``ang``
- ``flux``, flux array in ``erg s^-1 cm^-2 ang^-1``
- ``err``, error array in  ``erg s^-1 cm^-2 ang^-1``

### FITS header

``PySpecLines`` also uses the following (**optional**) FITS header keywords:

- ``REDSHIFT``, used to de-redshift a spectrum from the observed-frame to the rest-frame
- ``RA``, right ascension of the object in ``deg``, in the ``ICRS`` frame, used to correct for Galactic absorption
- ``DEC``, declination of the object in ``deg``, in the ``ICRS`` frame, used to correct for Galactic absorption

When both keywords ``RA`` and ``DEC`` are present, ``PySpecLines`` can de-redden the spectrum (passing the option ``-deredden``) using the Galactic dust map of Schlegel, Finkbeiner & Davis (1998) recalibrated by Schlafly & Finkbeiner (2011), and the ``R_V=3.1`` extinction curve of Fitzpatrick (1999).  

### JSON configuration file

A JSON file allows you to select and configure the emission lines to be measured. Some example JSON files are provided in the [PySpecLines/files](https://github.com/jacopo-chevallard/PySpecLines/tree/master/PySpecLines/files) folder. The JSON file contains a dictionary of ``key : values`` where ``key`` labels the line (or group of lines) and ``values`` contains multiple entries.

Below we report some simple examples:

- single line: 
  ```json
      {
      "HeII4686" : {
        "wl_central": [4686.0], 
        "wl_range":[4680.0, 4681.0], 
        "continuum_left":[4672.0, 4679.0],
        "continuum_right":[4692.0, 4700.0]
      }
    }
  ```
  - ``wl_central`` is used as starting point for the line center in the Gaussian fitting
  - ``wl_range`` is the range of numerical integration of the line
  - ``continuum_left`` is the range of numerical integration of the left-continuum
  - ``continuum_right`` is the range of numerical integration of the right-continuum 
  - when using Gaussian fitting, the range over which the continuum is fitted is ``[continuum_left[0], continuum_right[1]]``

- line doublet: 
  ```json
    {
    "SII6716_SII6731" : {
        "wl_central": [6716.0, 6731.0], 
        "wl_range":[6710.0, 6723.0], 
        "exclude":[6710.0, 6721.0, 6726.0, 6740.0],
        "continuum_left":[6695.0, 6705.0],
        "continuum_right":[6740.0, 6750.0]
      }
    }
  ```
  - wrt to the example above, the ``key`` is composed of two labels separated by an underscore ``_``
  - ``exclude`` allows to define regions (``[exclude[0], exclude[1]], [exclude[2], exclude[3]]``) excluded from the continuum fitting 
  
- multiple kinematic components
  ```json
    {
    "OIII5007N_OIII5007B" : {
        "wl_central": [5007.0, 5007.0], 
        "width": [100.0, 400.0], 
        "wl_range":[5000.0, 5014.0], 
        "exclude":[4995.0, 5020.0],
        "continuum_left":[4990.0, 5000.0],
        "continuum_right":[5020.0, 5034.0]
      }
    }
  ```
  - ``width`` allows to define multiple kinematic components, in this case a "narrow" (labelled ``OIII5007N``) and a "broad" (labelled ``OIII5007B``) component, whose starting widths **must** be set to different values.
  
- multiple lines with multiple kinematic components
  ```json
    {
    "NII6548_HalphaN_HalphaB_NII6584" : {
        "wl_central": [6548.05, 6563.0, 6563.0, 6584.0], 
        "width": [100.0, 100.0, 400.0, 100.0], 
        "exclude":[6542.0, 6580.0],
        "wl_range":[6541.0, 6575.0], 
        "continuum_left":[6515.0,6542.0],
        "continuum_right":[6595.0,6610.0]
      }
    }
  ```
  - in this case we want to use the same width for different lines (``NII6548``, ``HalphaN`` and ``NII6584``) and a different width for ``HalphaB``. We thus use the same ``width`` value for ``NII6548``, ``HalphaN`` and ``NII6584``, as this will "tie" together their widths during the Gaussian fitting, while the ``width`` of ``HalphaB`` will be kept separate.
  


## Examples

You can see the different available options running
```
pyspeclines --help
```

If the spectrum is provided in the observed frame, then you must provide a ``REDSHIFT`` keyword in the FITS header containing the object redshift.

- Compute the fluxes and EWs using numerical integration
  ```
  pyspeclines --file my_spectrum.fits --json-file  emission_lines_EWs_config.json
  ```

- Compute the fluxes and EWs using Gaussian fit (Levenberg-Marquardt)
  ```
  pyspeclines --file my_spectrum.fits --json-file  emission_lines_EWs_config.json --gaussian-fit
  ```

- Compute the fluxes and EWs using Gaussian fit (MCMC)
  ```
  pyspeclines --file my_spectrum.fits --json-file  emission_lines_EWs_config.json --gaussian-fit --use-PyMC --MCMC-samples 5000
  ```
