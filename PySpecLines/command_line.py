#!/usr/bin/env python
import argparse
import numpy as np
from matplotlib import rc
from pathos.multiprocessing import ProcessingPool 
from astropy import log

from _version import __version__
from pyspeclines import *


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--file',
        help="Name of the FITS file containing the observed spectrum for which EWs and integrated fluxes must be computed.",
        action="store", 
        type=str, 
        nargs="+",
        dest="file", 
        required=True
    )

    parser.add_argument(
        '--json-file',
        help="JSON file containing the list of emission lines for which EWs will be computed.",
        action="store", 
        type=str, 
        dest="json_file", 
        required=True
    )

    parser.add_argument(
        '--gaussian-fit',
        help="Use Gaussian line fitting instead of numerical (trapezoidal) integration to "\
              "compute fluxes and EWs.",
        action="store_true", 
        dest="gaussian_fit"
    )

    parser.add_argument(
        '--verbose',
        help="Choose the verbose level",
        action="store",
        type=str,
        dest="verbose",
        default="WARNING",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"]
    )

    parser.add_argument(
        '--deredden',
        help="Deredden the spectrum before computing the fluxes, using the dust maps of "\
             "Schlafly & Finkbeiner (2011) and the R_V=3.1 extinction curve of Fitzpatrick (1999)",
        action="store_true", 
        dest="deredden"
    )

    parser.add_argument(
        '--iterate-fit',
        help="Number of iterations fit continuum - fit lines.",
        action="store", 
        type=int, 
        default=3,
        dest="n_iter"
    )

    parser.add_argument(
        '--extend-continuum-region',
        help="Number of pixels that will be added to the left and to the right of the continuum "\
                "region to define the region of the continuum+Gaussian fit",
        action="store", 
        type=int, 
        default=0,
        dest="extend_region"
    )

    parser.add_argument(
        '--continuum-degree',
        help="Degree of the polynomial used to fit the continuum",
        action="store", 
        type=int, 
        default=1,
        dest="continuum_degree"
    )

    parser.add_argument(
        '--continuum-error',
        help="Relative error of the continuum determination",
        action="store", 
        type=float, 
        dest="continuum_error"
    )

    parser.add_argument(
        '--use-PyMC',
        help="Use PyMC to perform an MCMC exploration of the parameter space",
        action="store_true", 
        dest="use_pymc"
    )

    parser.add_argument(
        '--log-flux',
        help="Plot logarithmic y-axis scale",
        action="store_true", 
        dest="log_flux"
    )

    parser.add_argument(
        '--show-plot',
        help="Show plot of each fitted line(s)",
        action="store_true", 
        dest="show_plot"
    )

    parser.add_argument(
        '--credible-intervals',
        help="Credible intervals to be used to compute the error on the fluxes",
        action="store", 
        nargs="+",
        type=float, 
        default=[0.68],
        dest="credible_intervals"
    )

    parser.add_argument(
        '--MCMC-samples',
        help="Number of MCMC samples to be drawn",
        action="store", 
        type=int, 
        default=10000,
        dest="n_samples"
    )

    parser.add_argument(
        '--MCMC-burn_in',
        help="Fraction of samples to be discarded for burn-in",
        action="store", 
        type=float, 
        default=0.2,
        dest="frac_burn"
    )

    parser.add_argument(
        '-np',
        help="Number of processors to use",
        action="store", 
        type=int, 
        dest="n_proc",
        default=1
    )

    # Add package version
    parser.add_argument('-v', '--version', 
            action='version', 
            version='%(prog)s ' + __version__ + ' - Author: Jacopo Chevallard'
            )

    # Get parsed arguments
    args = parser.parse_args()    

    log.setLevel(args.verbose)

    # Create "pool" of processes
    if args.n_proc > 1:
        pool = ProcessingPool(nodes=args.n_proc)

    file_names = args.file
    json_file = args.json_file

    if args.n_proc > 1:
        pool.map(compute_fluxes_EWs, 
                file_names, 
                [json_file,]*len(file_names),
                [args,]*len(file_names))

    else:
        for file_name in file_names:
            compute_fluxes_EWs(file_name, json_file, args)
