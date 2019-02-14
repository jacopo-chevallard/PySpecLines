#! /usr/bin/env python

from collections import OrderedDict
import os
import json
from astropy.io import fits
from astropy import constants as const
from astropy import log
import dustmaps.sfd
import dustmaps.config 
from astropy.coordinates import SkyCoord
import astropy.units as units
import numpy as np
import pyspeckit
import matplotlib.pyplot as plt
import pkg_resources
import pyspeclines
import extinction

SEP = "_"
FLUX_PREFIX = "F" + SEP
EW_PREFIX = "EW" + SEP
ERR_SUFFIX = SEP + "err"

def sum_errors_in_quadrature(errors):

    _errors = np.array(errors)
    _err = np.sqrt(np.sum(_errors**2))
    return _err

def get_multiple_keys(key):

    _key_sp = key.split(SEP)
    _key_sp.append(key)

    return _key_sp

def velocity_to_pixel(velocity, wl_central, dwl):

    _dwl = velocity / const.c.to('km/s').value * wl_central

    pxl = _dwl / dwl

    return pxl

def pixel_to_velocity(pxl, wl_central, dwl):

    _dwl = pxl * dwl

    velocity = _dwl / wl_central * const.c.to('km/s').value

def compute_fluxes_EWs(file_name, json_file, args):

    # Load JSON structure containing list of lines (and corresponding wl,
    # continuum regions, ...) for which we will compute integrated fluxes and EWs
    with open(json_file) as f:
        lines = json.load(f, object_pairs_hook=OrderedDict)

    # Create empty dictionaries that will contain the relevant quantities
    EWs = OrderedDict()
    EWs_errors = OrderedDict()
    integrated_fluxes = OrderedDict()
    integrated_errors = OrderedDict()

    # Load the input spectrum
    hdulist = fits.open(file_name)
    wl = hdulist[1].data["wl"]
    spectrum = hdulist[1].data["flux"]
    error = hdulist[1].data["err"]

    # Apply correction for Galactic extinction
    if args.deredden:
        if 'RA' not in hdulist[1].header:
            raise ValueError("'RA' keyword (in degrees) must be present in the FITS header!")
        if 'DEC' not in hdulist[1].header:
            raise ValueError("'DEC' keyword (in degrees) must be present in the FITS header!")

        package_dir = pyspeclines.__path__[0]
        dust_map_path = os.path.join(os.path.dirname(package_dir), 'PySpecLines', 'files', 'dustmaps')
        dustmaps.config.config['data_dir'] = dust_map_path

        # If the dust maps are not present on the local machine, download them !
        dust_map_ngp = os.path.join(dust_map_path, 'SFD_dust_4096_ngp.fits')
        dust_map_sgp = os.path.join(dust_map_path, 'SFD_dust_4096_sgp.fits')
        if not os.path.isfile(dust_map_ngp) or not os.path.isfile(dust_map_sgp):
            dustmaps.sfd.fetch()
        
        # Get the color excerss E(B-V) at the location of the galaxy
        sfd = dustmaps.sfd.SFDQuery()
        ra, dec = hdulist[1].header['RA'], hdulist[1].header['DEC']
        coord = SkyCoord(ra, dec, frame='icrs', unit=(units.deg, units.deg))
        E_B_V = sfd(coord)
        log.info('Applying a Galactic extinction correction E(B-V) = ' + str(E_B_V))

        # Convert E(B-V) to A_V
        R_V = 3.1
        A_V = E_B_V * R_V

        # Get the Fitzpatrick (1999) extinction curve
        extinction_curve = extinction.fitzpatrick99(wl, -A_V, R_V) 

        # De-redden spectrum and error
        spectrum = extinction.apply(extinction_curve, spectrum)
        error = extinction.apply(extinction_curve, error)

    # De-redshift the spectrum
    redshift = None
    if 'redshift' in hdulist[1].header:
        redshift = hdulist[1].header["redshift"]
        log.info('De-redshifting the spectrum using z = ' + str(redshift))
        wl /= (1.+redshift)
        spectrum *= (1.+redshift)
        error *= (1.+redshift)

    # Cycle across all the lines 
    for key, value in lines.iteritems():

        print "\nkey: ", key

        if args.gaussian_fit:

            # Select the region over which you'll fit a (or multliple)
            # Gaussian + continuum. We exploit the continuum windows
            # provided by the user to define this region, considering N
            # pixels to the left of the left-region, and N pixels to the
            # right of the right region
            il0 = np.searchsorted(wl, value["continuum_left"][0]) ; il0 = max([0, il0-args.extend_region])
            ir1 = np.searchsorted(wl, value["continuum_right"][1]) ; ir1 = min([len(wl)+1, ir1+args.extend_region])

            _wl = wl[il0:ir1]
            _dwl = _wl[1] - _wl[0]  
            _spectrum = spectrum[il0:ir1]
            _error = error[il0:ir1]

            exclude = None
            if "exclude" in value:
                exclude = value["exclude"]

            # Main pyspeckit class that will enable you to do the analysis
            sp = pyspeckit.Spectrum(data=_spectrum, error=_error, xarr=_wl,
                        xarrkwargs={'unit':'AA'},
                        unit='$erg/s/cm^2/AA$')
            
            # Fit the baseline (continuum). In this first fit of the
            # continuum, pyspeckit will try to fit the lines as well, so
            # the continuum will be biased high. This is why we will refit
            # the continuum (and the lines) later on.
            sp.baseline(subtract=False, highlight_fitregion=False,
                    order=args.continuum_degree, excludefit=False, exclude=exclude)

            # Fit the lines using some input guesses for the amplitude,
            # width, and centroid of the line(s). As centroid(s) use the
            # user-provided wl of the transitions.  Note that here we tie
            # togther the width of the lines, i.e. multiple lines share the
            # same line-width.
            amplitude_guess = np.mean(spectrum)
            _width_guess = 1.
            guesses = list() ; tied = list()
            width_velocity = None
            for c, component in enumerate(value["wl_central"]):
                center_guess = component
                width_guess = _width_guess
                if "width" in value:
                    width_velocity = value["width"][c]
                    width_guess = velocity_to_pixel(width_velocity, center_guess, _dwl)

                guesses = guesses + [amplitude_guess, center_guess, width_guess]

                # The first line has line-width as a free parameter
                if c == 0:
                    tied = tied + ['', '', '']
                # For the successive lines the width can be fixed or variable
                else:
                    if width_velocity is None:
                        tied = tied + ['', '', 'p[2]']
                    else:
                        width_velocities = list(value["width"][0:c])
                        try:
                            p = width_velocities.index(width_velocity)
                            p = 2 + p * 3
                            tied = tied + ['', '', 'p[2]']
                        except:
                            tied = tied + ['', '', '']

            # Fit the line(s) with a Gaussian function
            sp.specfit(fittype='gaussian', guesses=guesses, tied=tied)

            for N in range(args.n_iter):
                # Fit again the continuum, this time *accounting* for the
                # presence of the emission lines we just fitted. The continuum
                # will not be biased high anymore.
                sp.baseline(subtract=False, highlight_fitregion=False,
                        order=args.continuum_degree, excludefit=True, exclude=exclude)

                # Fit again the lines (using the re-computed continuum)
                sp.specfit(fittype='gaussian', guesses=guesses, tied=tied)

            sp.plotter(errstyle='fill')
            if args.log_flux: 
                plt.yscale('log')

            # Retrieve the fitted continuum
            continuum = sp.baseline.basespec

            print "##########################"
            print "## Gaussian Fit results ##" 
            print "##########################"
            integrated_flux = list() ; integrated_error =  list()
            EW = list() ; EW_err = list()
            par_values = list()

            for c, component in enumerate(value["wl_central"]):
                c0, c1 = c*3, (c+1)*3

                # Retrieve the parameters of the Gaussian fit and corresponding errors
                amplitude, center, width = sp.specfit.parinfo.values[c0:c1]
                par_values.append(amplitude) ; par_values.append(center) ; par_values.append(width)
                amplitude_err, center_err, width_err = sp.specfit.parinfo.errors[c0:c1]

                # Compute the integrated flux (integral of a Gaussian)
                _integrated_flux = np.sqrt(2*np.pi)*width*amplitude
                integrated_flux.append(_integrated_flux)

                # Basic error propagation to obtain the error on the integrated flux
                _integrated_error = np.sqrt((amplitude_err/amplitude)**2 + (width_err/width)**2)
                integrated_error.append(_integrated_error*_integrated_flux)

                # Compute the median continuum around the fitted line center
                i0 = np.searchsorted(_wl, center-3.*width)
                i1 = np.searchsorted(_wl, center+3*width)
                _continuum = np.median(continuum[i0:i1+1])

                # Compute the EW and error
                _EW = _integrated_flux/_continuum
                EW.append(_EW)
                _EW_err = _integrated_error
                EW_err.append(_EW*_EW_err)

                if args.continuum_error is not None:
                    _continuum_error = args.continuum_error * _continuum
                    integrated_error[c] = sum_errors_in_quadrature([integrated_error[c], _continuum_error])
                    _continuum_error = args.continuum_error
                    EW_err[c] = sum_errors_in_quadrature([EW_err[c], _continuum_error])
                
                print "Flux, error, line center, width: ", integrated_flux[c], integrated_error[c], center, pixel_to_velocity(width, center, _dwl)
                print "EW, error: ", EW[c], EW_err[c]

            # Here we sum up the different components (useful to have also
            # the total flux of a line doublet)
            if len(integrated_flux) > 1:
                _integrated_flux = np.sum(integrated_flux)
                _integrated_error = np.sqrt(np.sum(np.array(integrated_error)**2))
                integrated_flux.append(_integrated_flux)
                integrated_error.append(_integrated_error)

                _EW = np.sum(EW)
                _EW_err = np.sqrt(np.sum(np.array(EW_err)**2))
                EW.append(_EW)
                EW_err.append(_EW_err)
                print "(Sum of components) Flux, error: ", integrated_flux[-1], integrated_error[-1]
                print "(Sum of components) EW, error: ", EW[-1], EW_err[-1]

            ########################################################################
            # Here we use MCMC to compute the line fluxes are corresponding errors #
            ########################################################################
            if args.use_pymc:

                # Adopt non-informative priors (I guess they mean uniform)
                MC_uninformed = sp.specfit.get_pymc()
                burn = int(args.frac_burn * args.n_samples)

                # Sample the parameter space with adaptive Metropolis-Hastings
                MC_uninformed.sample(args.n_samples, burn=burn, tune_interval=250)

                # Define quantiles that will be used to compute the
                # posterior central credible region (= errors)
                _perc = [0.5]
                for interval in args.credible_intervals:
                    _p_low, _p_up = 0.5*(1.-interval), 1.-0.5*(1.-interval)
                    _perc = _perc + [_p_low, _p_up]

                _perc = 100.*np.array(_perc)

                print "###############################"
                print "## MCMC Gaussian Fit results ##" 
                print "###############################"
                integrated_flux = list() ; integrated_error =  list()
                EW = list() ; EW_err = list()
                areas = list() ; areas_continuum = list()
                par_values = list()

                for c, component in enumerate(value["wl_central"]):

                    # Retrieve the parameters of the Gaussian fit and compute area of the Gaussian (= integrated flux)
                    _amplitude = MC_uninformed.trace('AMPLITUDE'+str(c))[:]
                    _center = MC_uninformed.trace('SHIFT'+str(c))[:]
                    _width =MC_uninformed.trace('WIDTH'+str(c))[:] 

                    amplitude, center, width = np.median(_amplitude), np.median(_center), np.median(_width) 
                    par_values.append(amplitude) ; par_values.append(center) ; par_values.append(width)

                    _area = np.sqrt(2*np.pi)*_width*_amplitude
                    areas.append(_area)

                    # Compute the quantiles
                    _integrated_flux, err_low, err_up = np.percentile(_area, _perc)
                    integrated_flux.append(_integrated_flux)
                    integrated_error.append(0.5 * (err_up-err_low))

                    # Compute the median continuum around the fitted line center
                    width = np.median(_width)
                    i0 = np.searchsorted(_wl, center-3.*width)
                    i1 = np.searchsorted(_wl, center+3*width)
                    _continuum = np.median(continuum[i0:i1+1])

                    # Compute the EW and error
                    __EW = _area/_continuum
                    areas_continuum.append(__EW)
                    _EW, err_low, err_up = np.percentile(__EW, _perc)
                    EW.append(_EW)
                    EW_err.append(0.5 * (err_up-err_low))

                    if args.continuum_error is not None:
                        _continuum_error = args.continuum_error * _continuum
                        integrated_error[c] = sum_errors_in_quadrature([integrated_error[c], _continuum_error])
                        _continuum_error = args.continuum_error
                        EW_err[c] = sum_errors_in_quadrature([EW_err[c], _continuum_error])

                    print "Flux, error: ", integrated_flux[c], integrated_error[c]
                    print "EW, error: ", EW[c], EW_err[c]

                # Here we sum up the different components (useful to have also
                # the total flux of a line doublet)
                if len(integrated_flux) > 1:
                    _area = np.zeros(len(areas[0]))
                    _area_continuum = np.zeros(len(areas[0]))
                    for _a, _a_c in zip(areas, areas_continuum):
                        _area = _area + _a
                        _area_continuum = _area_continuum + _a_c

                    _integrated_flux, err_low, err_up = np.percentile(_area, _perc)
                    integrated_flux.append(_integrated_flux)
                    integrated_error.append(0.5 * (err_up-err_low))

                    _EW, err_low, err_up = np.percentile(_area_continuum, _perc)
                    EW.append(_EW)
                    EW_err.append(0.5 * (err_up-err_low))

                    print "(Sum of components) Flux, error: ", integrated_flux[-1], integrated_error[-1]
                    print "(Sum of components) EW, error: ", EW[-1], EW_err[-1]


            sp.specfit.plot_fit(pars=par_values)

            integrated_fluxes[key] = integrated_flux
            integrated_errors[key] = integrated_error
            EWs[key] = EW
            EWs_errors[key] = EW_err

            if args.show_plot:
                plt.show()

            fig_name = key + ".pdf"
            sp.plotter.savefig(fig_name)

        else:

            integrated_flux = list() ; integrated_error =  list()
            EW = list() ; EW_err = list()

            # Compute the average left continuum
            il0 = np.searchsorted(wl, value["continuum_left"][0]) ; il0 -= 1
            il1 = np.searchsorted(wl, value["continuum_left"][1])
            if il0 == il1:
                il0 -= 1
            flux_left = np.trapz(spectrum[il0:il1+1], x=wl[il0:il1+1]) / (wl[il1]-wl[il0])
            _wlleft = 0.5*(wl[il0]+wl[il1])
            #print "wl_lwft, flux_left: ", _wlleft, flux_left

            # Compute the average right continuum
            ir0 = np.searchsorted(wl, value["continuum_right"][0]) ; ir0 -= 1
            ir1 = np.searchsorted(wl, value["continuum_right"][1])
            if ir0 == ir1:
                ir0 -= 1
            flux_right = np.trapz(spectrum[ir0:ir1+1], x=wl[ir0:ir1+1]) / (wl[ir1]-wl[ir0])
            _wlright = 0.5*(wl[ir0]+wl[ir1])
            #print "wl_lwft, flux_left: ", _wlright, flux_right

            # Approximate the continuum with a straght line
            grad = (flux_right-flux_left)/(_wlright-_wlleft)
            intercept = flux_right - grad*_wlright
        
            #### Compute EW
            i0 = np.searchsorted(wl, value["wl_range"][0]) ; i0 -= 1
            i1 = np.searchsorted(wl, value["wl_range"][1])
            n_wl = i1-i0+1

            # Interpolate the _spectrum at the edges to have to right integration limits
            _wl = np.copy(wl[i0:i1+1])
            _spectrum = np.copy(spectrum[i0:i1+1])
            #for i in range(len(_wl)):
            #    print "wl, flux: ", _wl[i]*(1.+redshift), _spectrum[i]/(1.+redshift)
            _error = np.copy(error[i0:i1+1])
            _spectrum[0] = _spectrum[0] + (_spectrum[1]-_spectrum[0])/(_wl[1]-_wl[0]) * (value["wl_range"][0]-_wl[0])
            _spectrum[-1] = _spectrum[-2] + (_spectrum[-2]-_spectrum[-1])/(_wl[-2]-_wl[-1]) * (value["wl_range"][1]-_wl[-2])

            # Build the actual wl array over which you will perform the integration
            _wl = np.zeros(n_wl)
            _wl[0] = value["wl_range"][0]
            _wl[-1] = value["wl_range"][1]
            _wl[1:-1] = wl[i0+1:i1]

            relative_error = _error/_spectrum
            _integrated_relative_error = np.sqrt(np.sum(relative_error**2))

            integrand = 1.0 - _spectrum/(grad*_wl+intercept)
            _EW = -np.trapz(integrand, x=_wl)
            EW.append(_EW)
            _EW_err = abs(_EW) * _integrated_relative_error
            EW_err.append(_EW_err)

            EWs[key] = EW
            EWs_errors[key] = EW

            #print "_wl:  ", _wl
            #print "_spectrum:  ", _spectrum
            #print "continuum:  ", grad*_wl+intercept
            integrand = _spectrum-(grad*_wl+intercept)
            _integrated_flux = np.trapz(integrand, x=_wl)
            integrated_flux.append(_integrated_flux)
            _integrated_error = abs(_integrated_flux) * _integrated_relative_error
            integrated_error.append(_integrated_error)

            integrated_fluxes[key] = integrated_flux
            integrated_errors[key] = integrated_error

            print  "EW, flux, rel_err: ", EW, integrated_flux, integrated_error


    new_hdulist = fits.HDUList([fits.PrimaryHDU()])
    ##################################################################
    # Create a list of columns from the dictionary containing the EWs
    ##################################################################
    cols = list()
    for key, value in EWs.iteritems():

        if len(value) > 1:
            _key_sp = get_multiple_keys(key)
            for k in _key_sp:
                _key = EW_PREFIX + k
                col = fits.Column(name=_key, format='E')
                cols.append(col)

                _key = _key + ERR_SUFFIX
                col = fits.Column(name=_key, format='E')
                cols.append(col)
        else:
            _key = EW_PREFIX + key
            col = fits.Column(name=_key, format='E')
            cols.append(col)

            _key = _key + ERR_SUFFIX
            col = fits.Column(name=_key, format='E')
            cols.append(col)

    # Create a new binary table HDU 
    columns = fits.ColDefs(cols)
    new_hdu = fits.BinTableHDU.from_columns(columns, nrows=1)
    new_hdu.name = 'EQUIVALENT WIDTHS'

    for key, value in EWs.iteritems():
        if len(value) > 1:
            _key_sp = get_multiple_keys(key)
            for c, k in enumerate(_key_sp):
                _key = EW_PREFIX + k
                new_hdu.data[_key] = EWs[key][c]

                _key = _key + ERR_SUFFIX
                new_hdu.data[_key] = EWs_errors[key][c]
        else:
            _key = EW_PREFIX + key
            new_hdu.data[_key] = EWs[key]

            _key = _key + ERR_SUFFIX
            new_hdu.data[_key] = EWs_errors[key]

    # Add the new HDU to the Beagle file
    if new_hdu.name in new_hdulist:
        new_hdulist[new_hdu.name] = new_hdu
    else:
        new_hdulist.append(new_hdu)

    ##################################################################
    # Create a list of columns from the dictionary containing the integrated fluxes
    ##################################################################
    cols = list()
    for key, value in integrated_fluxes.iteritems():

        if len(value) > 1:
            _key_sp = get_multiple_keys(key)
            for k in _key_sp:
                _key = FLUX_PREFIX + k
                col = fits.Column(name=_key, format='E')
                cols.append(col)

                _key = _key + ERR_SUFFIX
                col = fits.Column(name=_key, format='E')
                cols.append(col)
        else:
            _key = FLUX_PREFIX + key
            col = fits.Column(name=_key, format='E')
            cols.append(col)

            _key = _key + ERR_SUFFIX
            col = fits.Column(name=_key, format='E')
            cols.append(col)

    # Create a new binary table HDU 
    columns = fits.ColDefs(cols)
    new_hdu = fits.BinTableHDU.from_columns(columns, nrows=1)
    new_hdu.name = 'INTEGRATED FLUXES'

    for key, value in integrated_fluxes.iteritems():
        if len(value) > 1:
            _key_sp = get_multiple_keys(key)
            for c, k in enumerate(_key_sp):
                _key = FLUX_PREFIX + k
                new_hdu.data[_key] = integrated_fluxes[key][c]

                _key = _key + ERR_SUFFIX
                new_hdu.data[_key] = integrated_errors[key][c]
        else:
            _key = FLUX_PREFIX + key
            new_hdu.data[_key] = integrated_fluxes[key]

            _key = _key + ERR_SUFFIX
            new_hdu.data[_key] = integrated_errors[key]

    # Add the new HDU to the Beagle file
    if new_hdu.name in new_hdulist:
        new_hdulist[new_hdu.name] = new_hdu
    else:
        new_hdulist.append(new_hdu)

    if args.gaussian_fit:
        _file_name = file_name.split(".fits")[0] + "_gaussian_fluxes.fits"
    else:
        _file_name = file_name.split(".fits")[0] + "_numerical_fluxes.fits"

    new_hdulist.writeto(_file_name, overwrite=True)




