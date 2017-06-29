
import numpy as np
import os
from gnc.TrajectoryNetCDFWriter import TrajectoryNetCDFWriter, GliderNetCDFWriterException


class ProfileNetCDFWriter(TrajectoryNetCDFWriter):
    
    # We want to inherit from TrajectoryNetCDFWriter but also add a few more properties,
    # specific to writing profile NetCDFs, so we have to call 
    # TrajectoryNetCDFWriter.__init__(self, ...)
    def __init__(self, config_path, format='NETCDF4_CLASSIC', comp_level=1, clobber=False, profile_id=1):
        
        TrajectoryNetCDFWriter.__init__(self, config_path, format=format, comp_level=1, clobber=clobber)
        self._profile_id = profile_id
        self._default_sensor_defs_path = os.path.realpath('{:s}/../resources/profile-sensor_defs.json'.format(self._glider_base_dir))
        
        # Clear self._nc_sensor_defs to remove the trajectory sensor defs
        self._nc_sensor_defs = None
        
        # Need to reset the config_path to trigger reading of profile-based sensor
        # definitions from self._default_sensor_defs_path
        self.config_path = config_path
        
    @property
    def profile_id(self):
        return self._profile_id
        
    @profile_id.setter
    def profile_id(self, profile_id):
        self._profile_id = profile_id
    
    
    def _set_profile_var(self):
        """ Sets Profile ID in NetCDF File

        """

        self.set_scalar('profile_id', self._profile_id)
        
        self._profile_id += 1
            
    
    def finish_nc(self):
        """Close the NetCDF file permanently and delete instance properties preventing
        new NetCDF files from being created
        """
        
        if not self._nc:
            self._logger.error('The NetCDF file has not been initialized')
            return
            
        if not self._nc.isopen():
            self._logger.warning('The NetCDF file is already closed: {:s}'.format(self._output_path))
            return
            
        # Set profile variables
        self._update_profile_vars()
        # Update global geospatial attributes
        self._update_geospatial_global_attributes()
        # Update global time_coverage attributes
        self._update_time_coverage_global_attributes()
        
        self._nc.close()
        
        output_nc = self._output_path
        
        self._nc = None
        self._output_path = None
        
        return output_nc
        
        
    def _update_profile_vars(self):
        """ Internal function that updates all profile variables
        before closing a file
        """

        self._set_profile_var()
        
        time_var_name = self.sensor_def_exists('drv_timestamp')
        if not time_var_name:
            self._logger.warning('Skipping creation of profile_time variable')
        else:
            if time_var_name in self._nc.variables:
                profile_time = self._netcdf_to_np_op(
                    self._nc.variables[time_var_name][:],
                    np.nanmin
                )
                self.set_scalar('profile_time', profile_time)
            else:
                self._logger.warning('Cannot set profile_time (missing {:s} variable)'.format(time_var_name))

        lat_var_name = self.sensor_def_exists('drv_interp_m_gps_lat')
        lon_var_name = self.sensor_def_exists('drv_interp_m_gps_lon')
        if lon_var_name in self._nc.variables:
            profile_lon = self._netcdf_to_np_op(
                self._nc.variables[lon_var_name][:],
                np.average
            )
            self.set_scalar('profile_lon', profile_lon)
        else:
            self._logger.warning('Cannot set profile_lon (missing {:s} variable)'.format(lon_var_name))

        if lat_var_name in self._nc.variables:
            profile_lat = self._netcdf_to_np_op(
                self._nc.variables[lat_var_name][:],
                np.average
            )
            self.set_scalar('profile_lat', profile_lat)
        else:
            self._logger.warning('Cannot set profile_lat (missing {:s} variable)'.format(lat_var_name))

    def __repr__(self):
        
        return '<ProfileNetCDFWriter(config_path={:s}, format={:s})>'.format(self._config_path, self._format)