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
        '--extend-continuum-region',
        help="Number of pixels that will be added to the left and to the right of the continuum "\
                "region to define the region of the continuum+Gaussian fit",
        action="store", 
        type=int, 
        default=10,
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
        '--use-PyMC',
        help="Use PyMC to perform an MCMC exploration of the parameter space",
        action="store_true", 
        dest="use_pymc"
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


    log.setLevel('WARNING')

    # Add package version
    parser.add_argument('-v', '--version', 
            action='version', 
            version='%(prog)s ' + __version__ + ' - Author: Jacopo Chevallard'
            )

    # Get parsed arguments
    args = parser.parse_args()    

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
