import sys
import os
import logging
import json
import datetime
import shutil
import glob
from copy import deepcopy
import numpy as np
import tempfile
from netCDF4 import Dataset, stringtoarr
from netCDF4 import default_fillvals as NC_FILL_VALUES
from shapely.geometry import Polygon
from dateutil import parser
from gnc.readers.dba import create_llat_dba_reader

NETCDF_FORMATS = ['NETCDF3_CLASSIC',
    'NETCDF4_CLASSIC',
    'NETCDF4']
DELAYED_MODE_EXTENSIONS = ['dbd',
    'ebd']
REALTIME_MODE_EXTENSIONS = ['sbd',
    'mbd',
    'nbd',
    'tbd']
    
_SECONDS_PER_YEAR = 60*60*24*365
_SECONDS_PER_DAY = 60*60*24
    
class TrajectoryNetCDFWriter(object):
    
    def __init__(self, config_path, format='NETCDF4_CLASSIC', comp_level=1, clobber=False):
        
        # Initilize private properties
        self._config_path = None
        self._out_nc = None
        self._format = None
        self._comp_level = None
        self._clobber = None
        self._logger = logging.getLogger(os.path.basename(__file__))
        self._glider_base_dir = os.path.dirname(__file__)
        self._default_sensor_defs_path = os.path.realpath('{:s}/../resources/trajectory-sensor_defs.json'.format(self._glider_base_dir))
        self._default_global_atts_path = os.path.realpath('{:s}/../resources/global_attributes.json'.format(self._glider_base_dir))
        
        self._deployment_path = None
        self._global_attributes_path = None
        
        # Sensor definitions
        self._default_sensor_defs = None
        self._config_sensor_defs = None
        self._nc_sensor_defs = None
        
        #self._tmp_dir = None
        #self._has_tmp_nc = None
        #self._move_pairs = []
        #self._output_nc_files = []
        
        self._attributes = {'deployment' : {},
            'global_attributes' : {},
            'instruments' : {}}
        self._nc = None
        self._record_dimension = None
        self._stream_index = 0
        
        # Set properties from _init_
        self.config_path = config_path
        self.format = format
        self.comp_level = comp_level
        self.clobber = clobber
        
        
    @property
    def config_path(self):
        
        return self._config_path
        
    @config_path.setter
    def config_path(self, config_path):
        
        # Do now allow changing of the configuration path if self._nc is not None
        if self._nc:
            raise GliderNetCDFWriterException('Cannot reconfigure with existing netCDF4.Dataset: {:s}'.format(self._nc))
            
        if not os.path.isdir(config_path):
            raise OSError('Invalid configuration path: {:s}'.format(config_path))
            
        # Make sure config path has deployment.json and global_attributes.json 
        self._deployment_path = os.path.join(config_path, 'deployment.json')
        if not os.path.isfile(self._deployment_path):
            raise GliderNetCDFWriterException('Missing deployment.json file: {:s}'.format(self._deployment_pat))
        self._global_attributes_path = os.path.join(config_path, 'global_attributes.json')
        if not os.path.isfile(self._global_attributes_path):
            raise GliderNetCDFWriterException('Missing global_attributes.json file: {:s}'.format(self._global_attributes_path))
        self._instruments_path = os.path.join(config_path, 'instruments.json')
        if not os.path.isfile(self._global_attributes_path):
            raise GliderNetCDFWriterException('Missing global_attributes.json file: {:s}'.format(self._instruments_path))
            
        # Check for config path sensor defs used to override the defaults
        self._config_sensor_defs_path = os.path.join(config_path, 'sensor_defs.json')
        
        # Set the configuration path    
        self._config_path = config_path
        
        # Load the default sensor definitions
        self._load_default_sensor_defs()
        
        # Load the configuration path sensor definitions
        self._load_config_sensor_defs()
        
        # Update self._default_sensor_defs with self._config_sensor_defs
        self._update_sensor_defs()
        
        # Read in NetCDF attribute files
        self._read_attributes()
        
    #@property
    #def output_path(self):
    #    
    #    return self._output_path
    #    
    #@output_path.setter
    #def output_path(self, output_path):
    #    
    #    if not os.path.isdir(output_path):
    #        self._logger.warning('Invalid output_path specified: {:s}'.format(output_path))
    #        return
    #        
    #    self._output_path = output_path
        
    @property
    def format(self):
        
        return self._format
        
        
    @format.setter
    def format(self, format):
        
        if format not in NETCDF_FORMATS:
            raise ValueError('Invalid NetCDF format: {:s}'.format(format))
            
        self._format = format
        
    @property
    def comp_level(self):
        
        return self._comp_level
        
    
    @comp_level.setter
    def comp_level(self, comp_level):
        
        if comp_level not in range(11):
            raise ValueError('Compression level must be a value from 0 - 10')
            
        self._comp_level = comp_level
        
    @property
    def clobber(self):
        
        return self._clobber
        
        
    @clobber.setter
    def clobber(self, clobber):
        
        if type(clobber) != bool:
            raise ValueError('Clobber value must be boolean')
            
        self._clobber = clobber
            
    @property
    def default_sensor_defs_path(self):
        
        return self._default_sensor_defs_path
            
    
    @property
    def default_sensor_defs(self):
        
        return self._default_sensor_defs
        
        
    @property
    def config_sensor_defs(self):
        
        return self._config_sensor_defs
        
    
    @property
    def nc_sensor_defs(self):
        
        return self._nc_sensor_defs
        
        
    @property
    def attributes(self):
        
        return self._attributes
        
       
    def dbas_to_nc(self, dba_files, output_path, clobber=False, mode=None):
        
        #if not self._output_path:
        #    self._logger.warning('No NetCDF output_path specified')
        #    return
        if not os.path.isdir(output_path):
            self._logger.error('Invalid output_path specified: {:s}'.format(output_path))
            return
        
        # Create the deployment/trajectory name
        try:
            trajectory_dt = parser.parse(self._attributes['deployment']['trajectory_date'])
        except ValueError as e:
            logging.error('Error parsing deployment trajectory_date: {:s} ({:s})'.format(self._attributes['deployment']['trajectory_date'], e))
            return
            
        # Create a temporary directory for creating/writing NetCDF prior to 
        # moving them to output_path
        tmp_dir = tempfile.mkdtemp()
        self._logger.debug('Temporary NetCDF directory: {:s}'.format(tmp_dir))
            
        # Write one NetCDF file for each input file
        output_nc_files = []
        for dba_file in dba_files:
        
            if not os.path.isfile(dba_file):
                logging.error('Invalid dba file specified: {:s}'.format(dba_file))
                continue
                
            self._logger.info('Processing dba file: {:s}'.format(dba_file))
            
            # Parse the dba file
            dba = create_llat_dba_reader(dba_file)
            if len(dba['data']) == 0:
                logging.warning('Skipping empty dba file: {:s}'.format(dba_file))
                continue
                
            # Split the filename and extension
            dba_filename, dba_ext = os.path.splitext(os.path.basename(dba_file))
            
            # Guess at the realtime/delayed status, based on file type, if not specified
            # by the user
            if mode:
                file_mode = mode
            else:
                try:
                    if dba['dbd_meta']['filename_extension'] in DELAYED_MODE_EXTENSIONS:
                        if mode:
                            file_mode = mode
                        else:
                            file_mode = 'delayed'
                    elif dba['dbd_meta']['filename_extension'] in REALTIME_MODE_EXTENSIONS:
                        if mode:
                            file_mode = mode
                        else:
                            file_mode = 'rt'
                    else:
                        logging.error('No realtime/delayed mode specified and unable to guess: {:s}'.format(dba_file))
                        continue
                except KeyError as e:
                    logging.error(e)
                    continue
            
            # Create the output NetCDF path

            out_nc_file = os.path.join(output_path, '{:s}-{:s}.nc'.format(dba_filename, file_mode))
                
            # Clobber existing files as long as self._clobber == True.  If not, skip
            # this file
            if os.path.isfile(out_nc_file):
                if self._clobber:
                    self._logger.info('Clobbering existing file: {:s}'.format(out_nc_file))
                else:
                    self._logger.warning('Skipping existing NetCDF file: {:s}'.format(out_nc_file))
                    continue
                    
            # Path to hold file while we create it
            _, tmp_nc = tempfile.mkstemp(dir=tmp_dir, suffix='.nc', prefix=os.path.basename(__file__))
            
            try:
                self.init_nc(tmp_nc)
            except (GliderNetCDFWriterException, IOError) as e:
                logging.error('Error initializing {:s}: {:s}'.format(tmp_nc, e))
                continue
            
            try:
                self.open_nc()
                # Add command line call used to create the file
                self.update_history('{:s} {:s}'.format(sys.argv[0], dba_file))
            except (GliderNetCDFWriterException, IOError) as e:
                logging.error('Error opening {:s}: {:s}'.format(tmp_nc, e))
                shutil.remove(tmp_nc)
                continue
            
            # Create and set the trajectory
            trajectory_string = '{:s}-{:s}-rt'.format(self.attributes['deployment']['glider'],
                trajectory_dt.strftime('%Y%m%dT%H%M'),
                file_mode)
            self.set_trajectory_id(trajectory_string)
            # Update the global title attribute
            self._nc.title = 'Slocum Glider dba file: {:s}'.format(trajectory_string)
            
            # Create the source file scalar variable
            self.set_source_file_var(dba['dbd_meta']['filename_label'], dba['dbd_meta'])
            
            # Add the derived sensor definitions
            dba_sensors = [s['sensor'] for s in dba['sensors']]
            if 'drv_timestamp' in dba_sensors and 'drv_timestamp' in self.nc_sensor_defs:
                self.update_sensor_def('drv_timestamp', dba['sensors'][dba_sensors.index('drv_timestamp')])
            if 'drv_m_gps_lat' in dba_sensors and 'drv_m_gps_lat' in self.nc_sensor_defs:
                self.update_sensor_def('drv_m_gps_lat', dba['sensors'][dba_sensors.index('drv_m_gps_lat')])
            if 'drv_m_gps_lon' in dba_sensors and 'drv_m_gps_lon' in self.nc_sensor_defs:
                self.update_sensor_def('drv_m_gps_lon', dba['sensors'][dba_sensors.index('drv_m_gps_lon')])
            if 'drv_pressure' in dba_sensors and 'drv_pressure' in self.nc_sensor_defs:
                self.update_sensor_def('drv_pressure', dba['sensors'][dba_sensors.index('drv_pressure')])
            if 'drv_depth' in dba_sensors and 'drv_depth' in self.nc_sensor_defs:
                self.update_sensor_def('drv_depth', dba['sensors'][dba_sensors.index('drv_depth')])
            if 'drv_interp_m_gps_lat' in dba_sensors and 'drv_interp_m_gps_lat' in self.nc_sensor_defs:
                self.update_sensor_def('drv_interp_m_gps_lat', dba['sensors'][dba_sensors.index('drv_interp_m_gps_lat')])
            if 'drv_interp_m_gps_lon' in dba_sensors and 'drv_interp_m_gps_lon' in self.nc_sensor_defs:
                self.update_sensor_def('drv_interp_m_gps_lon', dba['sensors'][dba_sensors.index('drv_interp_m_gps_lon')])
            
            # Write the data to the NetCDF file
            for r in dba['data']:
                self.stream_dict_insert(r)
            
            # Permanently close the NetCDF file after writing it  
            nc_file = self.finish_nc()
            
            # Add the output NetCDF file name to the list of those to be moved to args.output_dir
            if nc_file:
                shutil.move(tmp_nc, out_nc_file)
                
            output_nc_files.append(out_nc_file)
            
                
        #        self._move_pairs.append([tmp_nc, out_nc_file])
        #        self._has_tmp_nc = True
        #        
        ## Check for tmp NetCDF files. If none, delete the temporary directory
        #if not self._has_tmp_nc:
        #    shutil.rmtree(self._tmp_dir)
        #    self._tmp_dir = None
        #    return
        #        
        ## Move all created NetCDF files to args.output_dir
        #self._output_nc_files = []
        #for tmp_nc,out_nc in self._move_pairs:
        #    if os.path.isfile(out_nc):
        #        if not self.clobber:
        #            self._logger.info('Skipping existing NetCDF file: {:s}'.format(out_nc))
        #            continue
        #        else:
        #            self._logger.info('Clobbering existing NetCDF file: {:s}'.format(out_nc))
        #            try:
        #                os.remove(out_nc)
        #            except OSError as e:
        #                self._logger.error(e)
        #                continue
        #            
        #    # Move the tmp_nc to out_nc
        #    try:
        #        shutil.move(tmp_nc, out_nc)
        #        self._output_nc_files.append(out_nc)
        #    except:
        #        self._logger.error('Error moving {:s}: {:s}'.format(tmp_nc, e))
        #        continue
                
        # Delete the temporary directory once files have been moved
        try:
            self._logger.debug('Removing temporary directory: {:s}'.format(tmp_dir))
            shutil.rmtree(tmp_dir)
        except OSError as e:
            logging.error(e)
        
        return output_nc_files
        

    def init_nc(self, out_nc):
        
        if not self._record_dimension:
            raise GliderNetCDFWriterException('No record dimension found in sensor definitions')
            
        if self._nc:
            raise GliderNetCDFWriterException('Existing netCDF4.Dataset: {:s}'.format(self._nc))
            
        self._nc = Dataset(out_nc, mode='w', clobber=True, format=self._format)
        
        # Create the record dimension
        self._nc.createDimension(self._record_dimension['name'], size=self._record_dimension['dimension_length'])
        
        # Store the NetCDF destination name
        self._out_nc = out_nc
        
        # Write global attributes
        # Add date_created, date_modified, date_issued globals
        nc_create_ts = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')
        self._attributes['global']['date_created'] = nc_create_ts
        self._attributes['global']['date_issued'] = nc_create_ts
        self._attributes['global']['date_modified'] = nc_create_ts
        # Add history attribute if not present in self._attributes['global']
        if 'history' not in self._attributes['global']:
            self._attributes['global']['history'] = ' '
        if 'id' not in self._attributes['global']:
            self._attributes['global']['id'] = ' '
        
        # Add the global cdm_data_type attribute
        # MUST be 'Trajectory'
        self._attributes['global']['cdm_data_type'] = 'Trajectory'
        # Add the global featureType attribute
        # MUST be 'trajectory'
        self._attributes['global']['featureType'] = 'trajectory'
        
        # Write the NetCDF global attributes
        self.set_global_attributes()
        
        # Update global history attribute
        self.update_history(self)
        
        # Create platform container variable
        self.set_platform()
        
        # Create instrument container variables
        self.set_instruments()
        
        self._nc.close()

        
    def open_nc(self):
        """Open the current NetCDF file in append mode
        """
        
        if not self._out_nc:
            self._logger.error('The NetCDF file has not been initialized')
            return
            
        if self._nc and self._nc.isopen():
            raise GliderNetCDFWriterException('netCDF4.Dataset is already open: {:s}'.format(self._nc))
            
        # Open the NetCDF in append mode
        self._nc = Dataset(self._out_nc, mode='a')
        
        # Starting index of record dimension
        self._stream_index = self._get_record_dim_len()
    
    def finish_nc(self):
        """Close the NetCDF file permanently and delete instance properties preventing
        new NetCDF files from being created
        """
        
        if not self._out_nc or not os.path.isfile(self._out_nc):
            self._logger.error('No output NetCDF file specified')
            return
            
        if not self._nc:
            self._logger.error('The NetCDF file has not been initialized')
            return
            
        if not self._nc.isopen():
            self._logger.warning('The NetCDF file is already closed: {:s}'.format(self._output_path))
            return
            
        # Update global geospatial attributes
        self._update_geospatial_global_attributes()
        # Update global time_coverage attributes
        self._update_time_coverage_global_attributes()
        
        self._nc.close()
        
        #output_nc = self._output_path
        
        self._nc = None
        #self._output_path = None
        
        return self._out_nc
            
        
    def update_history(self, message):
        """ Updates the global history attribute with the message appended to
        and ISO8601:2004 timestamp
        """

        # Get timestamp for this access
        now_time_ts = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
        history_string = '{:s}: {:s}\n'.format(now_time_ts, message)
        if 'history' not in self._nc.ncattrs():
            self._nc.setncattr('history', history_string)
            return
        
        previous_history = self._nc.history.strip()
        if not previous_history:
            self._nc.history = history_string
        else:
            self._nc.history += history_string
        
        
    def set_platform(self):
        """ Creates a variable that describes the glider
        """

        self.set_scalar('platform')
        for key, value in sorted(self._attributes['deployment']['platform'].items()):
            self._nc.variables['platform'].setncattr(key, value)
        
    def set_global_attributes(self):
        """ Sets a dictionary of values as global attributes
        """

        for key, value in sorted(self._attributes['global'].items()):
            self._nc.setncattr(key, value)
            
    
    def _update_time_coverage_global_attributes(self):
        """Update all global time_coverage attributes.  The following global
        attributes are created/updated:
            time_coverage_start
            time_coverage_end
            time_coverage_duration
        """
        
        time_var_name = self.sensor_def_exists('drv_timestamp')
        if not time_var_name:
            self._logger.warning('Skipping set global time_coverage_start/end attributes')
        else:
            min_timestamp = self._nc.variables[time_var_name][:].min()
            max_timestamp = self._nc.variables[time_var_name][:].max()
            try:
                dt0 = datetime.datetime.utcfromtimestamp(min_timestamp)
            except ValueError as e:
                self._logger.error('Error parsing min {:s}: {:s} ({:s})'.format(time_var_name, min_timestamp, e))
            try:
                dt1 = datetime.datetime.utcfromtimestamp(max_timestamp)
            except ValueError as e:
                self._logger.error('Error parsing max {:s}: {:s} ({:s})'.format(time_var_name, max_timestamp, e))
                
            self._nc.setncattr('time_coverage_start', dt0.strftime('%Y-%m-%dT%H:%M:%SZ'))
            self._nc.setncattr('time_coverage_end', dt1.strftime('%Y-%m-%dT%H:%M:%SZ'))
            self._nc.setncattr('time_coverage_duration', self.delta_to_iso_duration(dt1 - dt0))
                    
    def _update_geospatial_global_attributes(self):
        """Update all global geospatial_ min/max attributes.  The following global
        attributes are created/updated:
            geospatial_lat_min
            geospatial_lat_max
            geospatial_lon_min
            geospatial_lon_max
            geospatial_bounds
            geospatial_vertical_min
            geospatial_vertical_max
        """
        
        lat_var_name = self.sensor_def_exists('drv_interp_m_gps_lat')
        lon_var_name = self.sensor_def_exists('drv_interp_m_gps_lon')
        if not lat_var_name or not lon_var_name:
            self._logger.warning('Skipping set global geospatial_lat/lon attributes')
            return
        else:
            
            if lat_var_name in self._nc.variables and lon_var_name in self._nc.variables:
                min_lat = self._nc.variables[lat_var_name][:].min()
                max_lat = self._nc.variables[lat_var_name][:].max()
                min_lon = self._nc.variables[lon_var_name][:].min()
                max_lon = self._nc.variables[lon_var_name][:].max()
                
                # Create polygon WKT and set geospatial_bounds
                coords = ((max_lat, min_lon),
                    (max_lat, max_lon),
                    (min_lat, max_lon),
                    (min_lat, min_lon),
                    (max_lat, min_lon))
                polygon = Polygon(coords)
                polygon_wkt = polygon.wkt
            else:
                min_lat = np.nan
                max_lat = np.nan
                min_lon = np.nan
                max_lon = np.nan
                polygon_wkt = u'POLYGON EMPTY'
            
            # Set the global attributes
            self._nc.setncattr('geospatial_lat_min', min_lat)
            self._nc.setncattr('geospatial_lat_max', max_lat)
            self._nc.setncattr('geospatial_lon_min', min_lon)
            self._nc.setncattr('geospatial_lon_max', max_lon)
            self._nc.setncattr('geospatial_bounds', polygon_wkt)
            
        depth_var_name = self.sensor_def_exists('drv_depth')
        if not depth_var_name:
            self._logger.warning('Skipping set global geospatial_vertical attributes')
        else:
            if depth_var_name in self._nc.variables:
                min_depth = self._nc.variables[depth_var_name][:].min()
                max_depth = self._nc.variables[depth_var_name][:].max()
            else:
                min_depth = np.nan
                max_depth = np.nan
                
            self._nc.setncattr('geospatial_vertical_min', min_depth)
            self._nc.setncattr('geospatial_vertical_max', max_depth)
            

    def sensor_def_exists(self, sensor):
        """Return the NetCDF variable name from the sensor definition, if it
        exists.  Returns None if not found"""
        
        if sensor not in self._nc_sensor_defs:
            self._logger.warning('No {:s} sensor definition found'.format(sensor))
            return None
            
        return self._nc_sensor_defs[sensor]['name']
        
            
    def update_sensor_def(self, sensor, new_def):
        """Updates the sensor definition with the sensor key with the key,value
        pairs in new_def"""
        
        if sensor not in self._nc_sensor_defs:
            self._nc_sensor_defs[sensor] = new_def
            return
            
        for k,v in new_def.items():
            if type(v) == dict:
                if k not in self._nc_sensor_defs[sensor]:
                    self._nc_sensor_defs[sensor][k] = {}
                for ak,av in v.items():
                    self._nc_sensor_defs[sensor][k][ak] = av
            else:
                self._nc_sensor_defs[sensor][k] = v
            
            
    def check_datatype_exists(self, key):
        if key not in self._nc_sensor_defs:
            raise KeyError('Unknown datatype {:s} cannot '
                           'be inserted to NetCDF: No sensor definition found'.format(key))

        datatype = self._nc_sensor_defs[key]
        if datatype['name'] not in self._nc.variables:
            self.set_datatype(key, datatype)

        return datatype
        
        
    def set_datatype(self, key, desc):
        """ Sets up a datatype description for the dataset
        """
        
        #if 'is_dimension' in desc and desc['is_dimension']:
        #    if key not in self._nc.dimensions.keys():
        #        try:
        #            self._nc.createDimension(desc['name'], desc['dimension_length'])
        #        except RuntimeError as e:
        #            raise GliderNetCDFWriterException('{:s}: {:s}->dim={:s}'.format(e, key, desc['dimension']))

        if len(desc) == 0:
            return  # Skip empty configurations

        if desc['name'] in self._nc.variables:
            return  # This variable already exists

        if desc['dimension'] is None:
            dimension = ()
        else:
            dimension = (desc['dimension'],)

        datatype = self._nc.createVariable(
            desc['name'],
            desc['type'],
            dimensions=dimension,
            zlib=True,
            complevel=self._comp_level,
            fill_value=NC_FILL_VALUES[desc['type']]
        )

        # Add an attribute to note the variable name used in the source data file
        desc['attrs']['source_variable'] = key
        desc['attrs']['coordinates'] = 'drv_m_gps_lon drv_m_gps_lat drv_depth drv_timestamp'
        if 'long_name' not in desc['attrs'] or not desc['attrs']['long_name'].strip():
            desc['attrs']['long_name'] = key
        for k, v in sorted(desc['attrs'].items()):
            datatype.setncattr(k, v)

        #if 'status_flag' in desc:
        #    status_flag = desc['status_flag']
        #    status_flag_name = self.get_status_flag_name(desc['name'])
        #    datatype.setncattr('ancillary_variables', status_flag_name)
        #    status_flag_var = self.nc.createVariable(
        #        status_flag_name,
        #        'i1',
        #        dimension,
        #        zlib=True,
        #        complevel=self.COMP_LEVEL,
        #        fill_value=NC_FILL_VALUES['i1']
        #    )
        #    # Append defaults
        #    sf_standard_name = desc['attrs']['standard_name'] + ' status_flag'
        #    status_flag['attrs'].update({
        #        'standard_name': sf_standard_name,
        #        'flag_meanings': self.QC_FLAG_MEANINGS,
        #        'valid_min': self.QC_FLAGS[0],
        #        'valid_max': self.QC_FLAGS[-1],
        #        'flag_values': self.QC_FLAGS
        #    })
        #    for key, value in sorted(status_flag['attrs'].items()):
        #        status_flag_var.setncattr(key, value)
                
                
    def set_scalar(self, key, value=None):
        datatype = self.check_datatype_exists(key)

        # Set None or NaN values to _FillValue
        if value is None or np.isnan(value):
            value = NC_FILL_VALUES[datatype['type']]

        self._nc.variables[datatype['name']].assignValue(value)

        #if "status_flag" in datatype:
        #    status_flag_name = self.get_status_flag_name(datatype['name'])
        #    flag = self.perform_qaqc(key, value)
        #    self.nc.variables[status_flag_name].assignValue(flag)
            
    
    def set_instruments(self):
        """ Adds a list of instrument descriptions to the dataset
        """

        for description in self._attributes['instruments']:
            self._set_instrument(
                description['name'],
                description['type'],
                description['attrs']
            )
            
            
    def set_trajectory_id(self, trajectory_string):
        """ Sets the trajectory dimension and variable for the dataset and the
        global id attribute

        Input:
            - glider: Name of the glider deployed.
            - deployment_date: String or DateTime of when glider was
                first deployed.
        """
            
        if 'trajectory' not in self._nc.variables:
            # Setup Trajectory Dimension
            self._nc.createDimension('traj_strlen', len(trajectory_string))

            # Setup Trajectory Variable
            trajectory_var = self._nc.createVariable(
                u'trajectory',
                'S1',
                ('traj_strlen',),
                zlib=True,
                complevel=self._comp_level
            )

            attrs = {
                'cf_role': 'trajectory_id',
                'long_name': 'Trajectory/Deployment Name',  # NOQA
                'comment': 'A trajectory is a single deployment of a glider and may span multiple data files.'  # NOQA
            }
            for key, value in sorted(attrs.items()):
                trajectory_var.setncattr(key, value)
        else:
            trajectory_var = self._nc.variables['trajectory']

        # Set the trajectory variable data
        trajectory_var[:] = stringtoarr(trajectory_string, len(trajectory_string))
        
        if not self._nc.getncattr('id').strip():
            self._nc.id = trajectory_string  # Global id variable
            
    def set_source_file_var(self, source_file_string, attrs=None):
        """ Sets the trajectory dimension and variable for the dataset and the
        global id attribute

        Input:
            - glider: Name of the glider deployed.
            - deployment_date: String or DateTime of when glider was
                first deployed.
        """
            
        if 'source_file' not in self._nc.variables:
            # Setup Trajectory Dimension
            self._nc.createDimension('source_file_strlen', len(source_file_string))

            # Setup Trajectory Variable
            source_file_var = self._nc.createVariable(
                u'source_file',
                'S1',
                ('source_file_strlen',),
                zlib=True,
                complevel=self._comp_level
            )

            if attrs:
                attrs['long_name'] = 'Source data file'
                attrs['comment'] = 'Name of the source data file and associated file metadata'
                for key, value in sorted(attrs.items()):
                    source_file_var.setncattr(key, value)
        else:
            source_file_var = self._nc.variables['source_file']

        # Set the trajectory variable data
        source_file_var[:] = stringtoarr(source_file_string, len(source_file_string))
        
        if not self._nc.getncattr('source').strip():
            self._nc.source = 'Observational Slocum glider data from source dba file {:s}'.format(source_file_string)  # Global source variable
        
        
    #def stream_dict_insert(self, line, qaqc_methods={}):
    def stream_dict_insert(self, line):
        """ Adds a data point glider_binary_data_reader library to NetCDF

        Input:
        - line: A dictionary of values where the key is a given
                <value name>-<units> pair that matches a description
                in the datatypes.json file.
        """
        
        if not self._nc_sensor_defs:
            raise GliderNetCDFWriterException('No sensor definitions defined')
            
        for name, value in line.items():
            #if name == 'timestamp':
            #    continue  # Skip timestamp, inserted above

            try:
                datatype = self.check_datatype_exists(name)
            except KeyError:
                self._logger.debug("Datatype {:s} does not exist".format(name))
                continue

            datatype = self._nc_sensor_defs[name]
            if datatype['dimension'] == self._record_dimension['name']:
                self.set_array_value(name, self._stream_index, value)
            else:
                self.set_scalar(name, value)
                #if name == "m_water_vx-m/s":
                #    self.fill_uv_vars(line)

        self._stream_index += 1
        
        
    def set_array_value(self, key, index, value=None):
        datatype = self.check_datatype_exists(key)

        # Set None or NaN values to _FillValue
        if value is None or np.isnan(value):
            value = NC_FILL_VALUES[datatype['type']]

        self._nc.variables[datatype['name']][index] = value
                
        
    def contains(self, datatype_key):
        if datatype_key in self.datatypes:
            field_name = self.datatypes[datatype_key]['name']
            return field_name in self.nc.variables
        else:
            return False

    def get_scalar(self, datatype_key):
        if self.contains(datatype_key):
            field_name = self.datatypes[datatype_key]['name']
            return self.nc.variables[field_name].getValue()

    def copy_field(self, src_glider_nc, datatype_key):
        datatype = self.check_datatype_exists(datatype_key)
        field_name = datatype['name']

        if src_glider_nc.contains(field_name):
            src_variable = src_glider_nc.nc.variables[field_name]

            if 'time' in src_variable.dimensions:
                self.set_array(datatype_key, src_variable[:])
            else:
                self.set_scalar(datatype_key, src_variable.getValue())

        else:
            raise KeyError(
                'Field not found in source glider NetCDF: %s' % field_name
            )
            

    def copy_glider_datatypes(self, src_glider_nc, datatype_keys):
        for datatype_key in datatype_keys:
            try:
                self.copy_field(src_glider_nc, datatype_key)
            except KeyError:
                self._logger.error("Copy field failed")
                
    
    def delta_to_iso_duration(self, time_delta):
        """Parse a datetime.timedelta object and return a ISO 8601:2004 duration formatted
        string
        """
        
        seconds = time_delta.total_seconds()
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        days, hours = divmod(hours, 24)
        days, hours, minutes = map(int, (days, hours, minutes))
        seconds = round(seconds, 6)
    
        ## build date
        date = ''
        if days:
            date = '%sD' % days
    
        ## build time
        time = u'T'
        # hours
        bigger_exists = date or hours
        if bigger_exists:
            time += '{:02}H'.format(hours)
        # minutes
        bigger_exists = bigger_exists or minutes
        if bigger_exists:
            time += '{:02}M'.format(minutes)
        # seconds
        if seconds.is_integer():
            seconds = '{:02}'.format(int(seconds))
        else:
            # 9 chars long w/leading 0, 6 digits after decimal
            seconds = '%09.6f' % seconds
        # remove trailing zeros
        seconds = seconds.rstrip('0')
        time += '{}S'.format(seconds)
        return u'P' + date + time
                

    def _get_record_dim_len(self):
        
        if not self._record_dimension:
            self._logger.error('No record dimension defined')
            return
            
        if self._record_dimension['name'] in self._nc.variables:
            return len(self._nc.variables[self._record_dimension['name']])
        else:
            return 0
            
            
    def _netcdf_to_np_op(self, variable_data, operation):
        array = np.array(variable_data)
        array[array == NC_FILL_VALUES['f8']] = float('nan')

        result = operation(array)
        if result == np.nan:
            result = NC_FILL_VALUES['f8']
        return result
        
        
    def _set_instrument(self, name, var_type, attrs):
        """ Adds a description for a single instrument
        """

        if name not in self._nc.variables:
            self._nc.createVariable(
                name,
                var_type,
                fill_value=NC_FILL_VALUES[var_type]
            )

        for key, value in sorted(attrs.items()):
            self._nc.variables[name].setncattr(key, value)
            
            
    def _load_default_sensor_defs(self):
        """Load the default sensor definitions from the package"""
        
        if not os.path.isfile(self._default_sensor_defs_path):
            self._logger.warning('No default sensor definitions found: {:s}'.format(self._default_sensor_defs_path))
            return
            
        self._logger.debug('Loading default sensor definitions: {:s}'.format(self._default_sensor_defs_path))
        
        try:
            with open(self._default_sensor_defs_path, 'r') as fid:
                self._default_sensor_defs = json.load(fid)
        except ValueError as e:
            self._logger.error('Error parsing default sensor definitions: {:s} ({:s})'.format(self._default_sensor_defs_path, e))
            self._nc_sensor_defs = None
            self._default_sensor_defs = None
            
    def _load_config_sensor_defs(self):
        """Load the configuration path sensor definitions if they exist"""
        
        if not os.path.isfile(self._config_sensor_defs_path):
            return
            
        self._logger.debug('Loading deployment sensor definitions: {:s}'.format(self._config_sensor_defs_path))
        
        try:
            with open(self._config_sensor_defs_path, 'r') as fid:
                self._config_sensor_defs = json.load(fid)
        except ValueError as e:
            self._logger.error('Error parsing user-defined sensor definitions: {:s} ({:s})'.format(self._config_sensor_defs_path, e))
            self._nc_sensor_defs = None
            self._default_sensor_defs = None
            
        
    def _update_sensor_defs(self):
        """Updates the sensor definition dicts in self._default_sensor_defs
        with key,value pairs in self._config_sensor_defs if the key is missing or
        if the value mapped to the key is changed"""
        
        if not self._default_sensor_defs:
            self._logger.debug('Instance contains no default sensor definitions') 
            self._nc_sensor_defs = None
            return
            
        self._nc_sensor_defs = deepcopy(self._default_sensor_defs)
        
        if not self._config_sensor_defs:
            self._logger.debug('Instance contains no configured sensor definitions') 
            self._configure_record_dimension()
            return
        
        for sensor,config_sensor_def in self._config_sensor_defs.items():
            self.update_sensor_def(sensor, config_sensor_def)
        
        # Set record dimension
        self._configure_record_dimension()
        
        
    def _configure_record_dimension(self):
        """Check self._nc_sensor_defs for the record dimension
        """
        
        # Check for the unlimited record dimension after all sensor defs have been
        # updated
        dims = [self._nc_sensor_defs[s] for s in self._nc_sensor_defs if 'is_dimension' in self._nc_sensor_defs[s] and self._nc_sensor_defs[s]['is_dimension']]
        if not dims:
            self._logger.warning('No record dimension specified in sensor definitions')
            self._logger.warning('Cannot write NetCDF data until a record dimension is defined')
            return
            
        if len(dims) != 1:
            self._logger.warning('Multiple record dimensions specified:')
            for dim in dims:
                self._logger.warning('Record dimension: {:s}'.format(dim['name']))
            self._logger.warning('Only one record dimension is allowed')
            
        self._record_dimension = dims[0]
        
    def _read_attributes(self):
        """Load in configurations from self._config_path"""
        
        # Read in default global attributes if the file is there
        if os.path.isfile(self._default_global_atts_path):
            with open(self._default_global_atts_path, 'r') as fid:
                self._attributes['global'] = json.load(fid)
        else:
            self._logger.warning('Default global attributes files does not exist: {:s}'.format(self._default_global_atts_path))
            
        with open(self._deployment_path, 'r') as fid:
            self._attributes['deployment'] = json.load(fid)
        with open(self._global_attributes_path, 'r') as fid:
            self._attributes['global'].update(json.load(fid))
        with open(self._instruments_path, 'r') as fid:
            self._attributes['instruments'] = json.load(fid)
            
        # Add the deployment['global_attributes'] to the globals
        self._attributes['global'].update(self._attributes['deployment']['global_attributes'])
        
        
    def _cleanup(self):
        
        if self._tmp_dir and os.path.isdir(self._tmp_dir):
            self._logger.debug('Cleaning up...')
            self._logger.debug('Removing temp NetCDF directory: {:s}'.format(self._tmp_dir))
            try:
                shutil.rmtree(self._tmp_dir)
            except OSError as e:
                self._logger.error(e)
            
            
    def __repr__(self):
        
        return '<NetCDFWriter(config_path={:s}, format={:s})>'.format(self._config_path, self._format)
        
    
class GliderNetCDFWriterException(Exception):
    
    def __init__(self, message, errors={}):
        self.message = message
        self.errors = errors
        
    
    def __str__(self):
        return self.message
