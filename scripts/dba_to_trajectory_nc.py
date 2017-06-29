#!/usr/bin/env python

import os
import sys
from gnc.TrajectoryNetCDFWriter import *
import logging
import argparse

DELAYED_MODE_EXTENSIONS = ['dbd',
    'ebd']
REALTIME_MODE_EXTENSIONS = ['sbd',
    'mbd',
    'nbd',
    'tbd']
    
def main(args):
    """Write one or more Slocum glider ascii dba files to a CF-compliant Trajectory
    NetCDF file
    """
    
    status = 0
    
    # Set up the erddapfoo.lib.m2m.M2mClient logger
    log_level = getattr(logging, args.loglevel.upper())
    log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
    logging.basicConfig(format=log_format, level=log_level)
    
    if not os.path.isdir(args.config_path):
        logging.error('Invalid configuration directory: {:s}'.format(args.config_dir))
        return 1
        
    if not args.output_path:
        args.output_path = os.path.realpath(os.curdir)
        logging.info('No NetCDF destination specified. Using cwd: {:s}'.format(args.output_path))
        
    if not os.path.isdir(args.output_path):
        logging.error('Invalid output_dir: {:s}'.format(args.output_dir))
        return 1
        
    # Create the Trajectory NetCDF writer
    ncw = TrajectoryNetCDFWriter(args.config_path, clobber=True)
    ## Set the NetCDF destination
    #ncw.output_path = args.output_path
    # Write the NetCDF files
    destination_nc_files = ncw.dbas_to_nc(args.dba_files, args.output_path, clobber=args.clobber)

    if not destination_nc_files:
        destination_nc_files = []
        
    # Print the list of files created
    for dest_nc_file in destination_nc_files:
        sys.stdout.write('{:s}\n'.format(dest_nc_file))
        
    return status
    
if __name__ == '__main__':
    
    arg_parser = argparse.ArgumentParser(description=main.__doc__)
    
    arg_parser.add_argument('config_path',
        help='Location of deployment configuration files')
        
    arg_parser.add_argument('dba_files',
        help='Source ASCII dba files to process',
        nargs='+')
        
    arg_parser.add_argument('-o', '--output_path',
        help='NetCDF destination directory.  Must exist.  <Default=cwd()>')
        
    arg_parser.add_argument('-c', '--clobber',
        help='Clobber existing NetCDF files if they exist',
        action='store_true')
        
    arg_parser.add_argument('-m', '--mode',
        type=str,
        help='String typically used to denote realtime or delayed mode status. If not specified, is inferred from source file type(s)')
        
    arg_parser.add_argument('-l', '--loglevel',
        help='Verbosity level <Default=info>',
        type=str,
        choices=['debug', 'info', 'warning', 'error', 'critical'],
        default='info')
    
    parsed_args = arg_parser.parse_args()

    sys.exit(main(parsed_args))