#!/usr/bin/env python

import os
import sys
from gnc.ProfileNetCDFWriter import ProfileNetCDFWriter, GliderNetCDFWriterException
from gnc.readers.dba import create_llat_dba_reader
from gnc.yo import find_profiles
import logging
import argparse
import tempfile
import shutil
import numpy as np
from dateutil import parser
import datetime

DELAYED_MODE_EXTENSIONS = ['dbd',
    'ebd']
REALTIME_MODE_EXTENSIONS = ['sbd',
    'mbd',
    'nbd',
    'tbd']
    
def main(args):
    """Write one or more Slocum glider ascii dba files to a CF-compliant Profile
    NetCDF file.
    """
    
    status = 0
    
    # Set up the erddapfoo.lib.m2m.M2mClient logger
    log_level = getattr(logging, args.loglevel.upper())
    log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
    logging.basicConfig(format=log_format, level=log_level)
    
    if not os.path.isdir(args.config_dir):
        logging.error('Invalid configuration directory: {:s}'.format(args.config_dir))
        return 1
        
    if not args.output_dir:
        args.output_dir = os.path.realpath(os.curdir)
        
    if not os.path.isdir(args.output_dir):
        logging.error('Invalid output_dir: {:s}'.format(args.output_dir))
        return 1
        
    # Temp directory
    tmpdir = tempfile.mkdtemp()
    logging.debug('Temporary NetCDF directory: {:s}'.format(tmpdir))
    
    move_pairs = []

    # Initialize the NetCDF writer   
    ncw = ProfileNetCDFWriter(args.config_dir, profile_id=args.profilestart)
    ncw.clobber = True
    
    # Create the deployment/trajectory name
    try:
        trajectory_dt = parser.parse(ncw.attributes['deployment']['trajectory_date'])
    except ValueError as e:
        logging.error('Error parsing deployment trajectory_date: {:s} ({:s})'.format(ncw.attributes['deployment']['trajectory_date'], e))
        return 1
        
    #profile_count = args.profile_start
    
    for dba_file in args.dba_files:
        
        if not os.path.isfile(dba_file):
            logging.error('Invalid dba file specified: {:s}'.format(dba_file))
            continue
            
        logging.debug('Processing dba file: {:s}'.format(dba_file))
        
        # Parse the dba file
        dba = create_llat_dba_reader(dba_file)
        if len(dba['data']) == 0:
            logging.warning('Skipping empty dba file: {:s}'.format(dba_file))
            continue
            
        # Create the yo for profile indexing find the profile minima/maxima
        try:
            profile_times = find_profiles(dba['data'])
        except ValueError as e:
            logging.error('{:s}: {:s}'.format(dba_file, e))
            
        if len(profile_times) == 0:
            logging.debug('No profiles indexed: {:s}'.format(dba_file))
            continue
            
        # All timestamps from stream
        ts = [r['drv_timestamp'] for r in dba['data'] if 'drv_timestamp' in r]
        
        for profile_interval in profile_times:
            
            # Profile start time
            p0 = profile_interval[0]
            # Profile end time
            p1 = profile_interval[-1]
            # Find all rows in ts that are between p0 & p1
            p_inds = np.flatnonzero(np.logical_and(ts >= p0, ts <= p1))
            profile_stream = dba['data'][p_inds[0]:p_inds[-1]]
            
            if args.mode:
                file_mode = args.mode
            else:
                try:
                    if dba['dbd_meta']['filename_extension'] in DELAYED_MODE_EXTENSIONS:
                        file_mode = 'delayed'
                    elif dba['dbd_meta']['filename_extension'] in REALTIME_MODE_EXTENSIONS:
                        file_mode = 'rt'
                    else:
                        logging.warning('Skipping {:s}: Unknown mode filetype: {:s}'.format(dba_file, dba['dbd_meta']['filename_extension']))
                        continue
                except KeyError as e:
                    logging.error(e)
                    status = 1
                    continue
                    
            # Calculate and convert profile mean time to a datetime
            pro_mean_dt = datetime.datetime.utcfromtimestamp(np.mean(profile_interval))
            
            # Create the output NetCDF path
            pro_mean_ts = pro_mean_dt.strftime('%Y%m%dT%H%M%SZ')
            profile_filename = '{:s}-{:s}-{:s}'.format(ncw.attributes['deployment']['glider'], pro_mean_ts, file_mode)
            # Path to temporarily hold file while we create it
            _, tmp_nc = tempfile.mkstemp(dir=tmpdir, suffix='.nc', prefix=os.path.basename(profile_filename))
        
            out_nc_file = os.path.join(args.output_dir, '{:s}.nc'.format(profile_filename))
        
            try:
                ncw.init_nc(tmp_nc)
            except (GliderNetCDFWriterException, IOError) as e:
                logging.error(e)
                status = 1
                continue
            
            try:
                ncw.open_nc()
                # Add command line call used to create the file
                ncw.update_history('{:s} {:s}'.format(sys.argv[0], dba_file))
            except (GliderNetCDFWriterException, IOError) as e:
                logging.error(e)
                status = 1
                shutil.remove(out_nc_file)
                continue
            
            # Create the trajectory string and set the trajectory variable
            trajectory_string = '{:s}-{:s}-rt'.format(ncw.attributes['deployment']['glider'],
                trajectory_dt.strftime('%Y%m%dT%H%M'),
                file_mode)
            ncw.set_trajectory_id(trajectory_string)
            # Update the global title attribute
            ncw._nc.title = 'Slocum Glider Profile: {:s}'.format('{:s}-{:s}-{:s}'.format(ncw.attributes['deployment']['glider'], pro_mean_ts, file_mode))
            
            # Create the source file scalar variable
            ncw.set_source_file_var(dba['dbd_meta']['filename_label'], dba['dbd_meta'])
            
            # Add the derived sensor definitions
            dba_sensors = [s['sensor'] for s in dba['sensors']]
            if 'drv_timestamp' in dba_sensors and 'drv_timestamp' in ncw.nc_sensor_defs:
                ncw.update_sensor_def('drv_timestamp', dba['sensors'][dba_sensors.index('drv_timestamp')])
            if 'drv_m_gps_lat' in dba_sensors and 'drv_m_gps_lat' in ncw.nc_sensor_defs:
                ncw.update_sensor_def('drv_m_gps_lat', dba['sensors'][dba_sensors.index('drv_m_gps_lat')])
            if 'drv_m_gps_lon' in dba_sensors and 'drv_m_gps_lon' in ncw.nc_sensor_defs:
                ncw.update_sensor_def('drv_m_gps_lon', dba['sensors'][dba_sensors.index('drv_m_gps_lon')])
            if 'drv_pressure' in dba_sensors and 'drv_pressure' in ncw.nc_sensor_defs:
                ncw.update_sensor_def('drv_pressure', dba['sensors'][dba_sensors.index('drv_pressure')])
            if 'drv_depth' in dba_sensors and 'drv_depth' in ncw.nc_sensor_defs:
                ncw.update_sensor_def('drv_depth', dba['sensors'][dba_sensors.index('drv_depth')])
            if 'drv_interp_m_gps_lat' in dba_sensors and 'drv_interp_m_gps_lat' in ncw.nc_sensor_defs:
                ncw.update_sensor_def('drv_interp_m_gps_lat', dba['sensors'][dba_sensors.index('drv_interp_m_gps_lat')])
            if 'drv_interp_m_gps_lon' in dba_sensors and 'drv_interp_m_gps_lon' in ncw.nc_sensor_defs:
                ncw.update_sensor_def('drv_interp_m_gps_lon', dba['sensors'][dba_sensors.index('drv_interp_m_gps_lon')])
            
            # Write the data to the NetCDF file
            for r in profile_stream:
                ncw.stream_dict_insert(r)
            
            # Write scalar profile variable and permanently close the NetCDF file 
            nc_file = ncw.finish_nc()
            
            # Add the output NetCDF file name to the list of those to be moved to args.output_dir
            if nc_file:
                move_pairs.append([tmp_nc, out_nc_file])
            
    # Move all created NetCDF files to args.output_dir
    destination_nc_files = []
    for tmp_nc,out_nc in move_pairs:
        if os.path.isfile(out_nc):
            if not args.clobber:
                logging.info('Skipping existing NetCDF file: {:s}'.format(out_nc))
                continue
            else:
                logging.info('Clobbering existing NetCDF file: {:s}'.format(out_nc))
                try:
                    os.remove(out_nc)
                except OSError as e:
                    logging.error(e)
                    continue
                
        # Move the tmp_nc to out_nc
        try:
            shutil.move(tmp_nc, out_nc)
            destination_nc_files.append(out_nc)
        except:
            logging.error('Error moving {:s}: {:s}'.format(tmp_nc, e))
            status = 1
            
    # Delete the temporary directory once files have been moved
    try:
        logging.debug('Removing temporary directory: {:s}'.format(tmpdir))
        shutil.rmtree(tmpdir)
    except OSError as e:
        logging.error(e)
            
    # Print the list of files created
    for dest_nc_file in destination_nc_files:
        sys.stdout.write('{:s}\n'.format(dest_nc_file))
        
    return status
    
if __name__ == '__main__':
    
    arg_parser = argparse.ArgumentParser(description=main.__doc__)
    
    arg_parser.add_argument('config_dir',
        help='Location of deployment configuration files')
        
    arg_parser.add_argument('dba_files',
        help='Source ASCII dba files to process',
        nargs='+')
        
    arg_parser.add_argument('-p', '--profilestart',
        help='Number specifying the starting profile id <Default=1>',
        type=int,
        default=1
    )
        
    arg_parser.add_argument('-o', '--output_dir',
        help='NetCDF destination directory.  Must exist.  <Default=pwd>')
        
    arg_parser.add_argument('-c', '--clobber',
        help='Clobber existing NetCDF files if they exist',
        action='store_true')
        
    arg_parser.add_argument('-m', '--mode',
        type=str,
        help='String typically used to denote realtime or delayed mode status. If not specified, is inferred from source file type(s)')
        
    arg_parser.add_argument('-l', '--loglevel',
        help='Verbosity level <Default=warning>',
        type=str,
        choices=['debug', 'info', 'warning', 'error', 'critical'],
        default='warning')
    
    parsed_args = arg_parser.parse_args()

    sys.exit(main(parsed_args))
