import numpy as np
import os
import sys
from dateutil import parser
import re
from gsw import z_from_p
import logging
from gnc.gps import interpolate_gps, get_decimal_degrees

months = ['Jan',
    'Feb',
    'Mar',
    'Apr',
    'May',
    'Jun',
    'Jul',
    'Aug',
    'Sep',
    'Oct',
    'Nov',
    'Dec']
    
VALID_DBA_TIMESTAMP_SENSORS = ['m_present_time',
    'sci_m_present_time']
VALID_DBA_PRESSURE_SENSORS = ['sci_water_pressure',
    'm_water_pressure',
    'm_pressure']
    
logger = logging.getLogger(os.path.basename(__file__))
    
    
def create_llat_dba_reader(dba_file, timesensor=None, pressuresensor=None):
    
    if not os.path.isfile(dba_file):
        logger.error('dba file does not exist: {:s}'.format(dba_file))
        return
    
    dba = create_dba_reader(dba_file)
    
    if len(dba['data']) == 0:
        logger.warning('No data records parsed: {:s}'.format(dba_file))
        return
       
    timestamp_sensor = None
    pressure_sensor = None
     
    dba_sensors = [s['sensor'] for s in dba['sensors']]
    if timesensor and timesensor not in dba_sensors:
        logger.warning('Specified timesensor {:s} not found in dba: {:s}'.format(timesensor, dba_file))
        return dba
    else:
        for t in VALID_DBA_TIMESTAMP_SENSORS:
            if t in dba_sensors:
                timestamp_sensor = t
                sensor_def = {'sensor' : 'drv_timestamp',
                    'attrs' : {}}
                sensor_def['attrs']['dba_sensor'] = timestamp_sensor
                sensor_def['attrs']['comment'] = 'Alias for {:s}'.format(timestamp_sensor)
                dba['sensors'].append(sensor_def)
                break
                
    if pressuresensor and pressuresensor not in dba_sensors:
        logger.warning('Specified pressuresensor {:s} not found in dba: {:s}'.format(pressuresensor, dba_file))
    else:
        for p in VALID_DBA_PRESSURE_SENSORS:
            if p in dba_sensors:
                pressure_sensor = p
                sensor_def = {'sensor' : 'drv_pressure',
                    'attrs' : {}}
                sensor_def['attrs']['dba_sensor'] = pressure_sensor
                sensor_def['attrs']['comment'] = 'Alias for {:s}'.format(pressure_sensor)
                dba['sensors'].append(sensor_def)
                break
        
    # Create the aliases
    counter = 0
    raw_gps = np.empty((len(dba['data']),4)) * np.nan
    for r in dba['data']:
        row_sensors = r.keys()
        if timestamp_sensor in row_sensors:
            r['drv_timestamp'] = r[timestamp_sensor]
            
            # Fill in the timestamps
            raw_gps[counter,0] = r['drv_timestamp']
            
            # Fill in raw_gps
            if 'm_gps_lat' in row_sensors and 'm_gps_lon' in row_sensors:
                
                # Ignore if the m_gps_lat and/or m_gps_lon value is the default
                # masterdata value
                if r['m_gps_lat'] <= 9000 or r['m_gps_lon'] < 18000:
                    r['drv_m_gps_lat'] = get_decimal_degrees(r['m_gps_lat'])
                    r['drv_m_gps_lon'] = get_decimal_degrees(r['m_gps_lon'])
                    raw_gps[counter,1] = r['drv_m_gps_lat']
                    raw_gps[counter,2] = r['drv_m_gps_lon']
                else:
                    logger.warning('Ignoring m_gps_lat/m_gps_lon default masterdata values: {:s}'.format(dba_file))
            
        if pressure_sensor and pressure_sensor in row_sensors:
            r['drv_pressure'] = r[pressure_sensor]
            raw_gps[counter,3] = r[pressure_sensor]
            
        counter += 1
    
    # Interpolate m_gps_lat/lon and add records to dataset['data'] 
    #try:
    interp_lat,interp_lon = interpolate_gps(raw_gps[:,0], raw_gps[:,1], raw_gps[:,2])
    #except IndexError as e:
    #    logger.error('{:s}: {:s}'.format(dba_file, e))
        #return
        
    # Calculate depth from pressure (multiplied by 10 to get to decibars) and latitude
    # Negate the results so that increasing values note increasing depths
    depths = -z_from_p(raw_gps[:,3]*10, interp_lat)
    depth_def = {'sensor' : 'drv_depth', 'attrs' : {}}
    depth_def['attrs']['comment'] = 'Depth calculated from pressure*10 (decibars) and raw latitude. Negative values denote increased depths'
    depth_def['attrs']['ancillary_variables'] = '{:s},drv_m_gps_lat'.format(pressure_sensor)
    dba['sensors'].append(depth_def)
        
    counter = 0
    for r in dba['data']:
        if not np.isnan(interp_lat[counter]):
            r['drv_interp_m_gps_lat'] =  interp_lat[counter]
        if not np.isnan(interp_lon[counter]):
            r['drv_interp_m_gps_lon'] = interp_lon[counter]
        if not np.isnan(depths[counter]):
            r['drv_depth'] = depths[counter]
        
        counter += 1
            
    # Create and add the converted m_gps_lat/lon sensor defs
    lat_dd_def = {'sensor' : 'drv_m_gps_lat', 'attrs' : {}}
    lat_dd_def['attrs']['comment'] = 'm_gps_lat converted to decimal degrees'
    lat_dd_def['attrs']['dba_sensor'] = 'm_gps_lat'
    dba['sensors'].append(lat_dd_def)
        
    lon_dd_def = {'sensor' : 'drv_m_gps_lon', 'attrs' : {}}
    lon_dd_def['attrs']['comment'] = 'm_gps_lon converted to decimal degrees'
    lon_dd_def['attrs']['dba_sensor'] = 'm_gps_lon'
    dba['sensors'].append(lon_dd_def)
        
    
    # Create and add the converted & interpolated m_gps_lat/lon sensor defs
    ilat_dd_def = {'sensor' : 'drv_interp_m_gps_lat', 'attrs' : {}}
    ilat_dd_def['attrs']['comment'] = 'm_gps_lat converted to decimal degrees and interpolated'
    ilat_dd_def['attrs']['dba_sensor'] = 'm_gps_lat'
    dba['sensors'].append(ilat_dd_def)
    ilon_dd_def = {'sensor' : 'drv_interp_m_gps_lon', 'attrs' : {}}
    ilon_dd_def['attrs']['comment'] =  'm_gps_lon converted to decimal degrees and interpolated'
    ilon_dd_def['attrs']['dba_sensor'] = 'm_gps_lon'
    dba['sensors'].append(ilon_dd_def)
        
    return dba
        
    
def create_dba_reader(dba_file):
    
    if not os.path.isfile(dba_file):
        logger.error('dba file does not exist: {:s}'.format(dba_file))
        return
        
    dba = parse_dba(dba_file)
    
    # Create the array of dicts mapping native sensor name to sensor value
    dba['data'] = [{s['sensor']:s['value'] for s in r} for r in dba['data']]
        
    return dba
    
def parse_dba(dba_file):
    '''Parse a Slocum Dinkum Binary ASCII (dba) file
    
    Args:
        dbd_file: the Slocum Dinkum Binary Data file to parse
    Returns:
        A dictionary representation of the data file'''
    
    if not os.path.isfile(dba_file):
        logger.error('dba file does not exist: {:s}'.format(dba_file))
        return
        
    dbd_meta = parse_dbd_header(dba_file)
    # Add the full path to the dba_file
    dbd_meta['source_file'] = os.path.realpath(dba_file)
    # Add the dba file size
    dba_stat = os.stat(dba_file)
    dbd_meta['file_size_bytes'] = dba_stat.st_size
    sensors_meta = parse_dbd_sensor_meta(dba_file)
    data = parse_dbd_data(dba_file)
    
    dbd_json = {'dbd_meta' : dbd_meta,
        'sensors' : sensors_meta,
        'data' : data}
        
    return dbd_json
    

def fileopen_time_to_datetime(fileopen_time):
    
    # Convert dbd_obj['dbd_meta']['fileopen_time'] to datetime
    t = fileopen_time.replace('__', '_').split('_')
    
    try:
        dt = parser.parse(
            '{:s}-{:d}-{:s} {:s}Z'.format(t[4], months.index(t[1])+1, t[2], t[3]))
    except ValueError as e:
        sys.stderr.write('{:s}: {:s}\n'.format(fileopen_time, e))
        return
        
    return dt
        
def dbd2docs(dbd_obj):
    
    doc = {'rows' : []}
    
    # Convert dbd_obj['dbd_meta']['fileopen_time'] to datetime
    t = dbd_obj['dbd_meta']['fileopen_time'].replace('__', '_').split('_')

    try:
        fileopen_time = parser.parse(
            '{:s}-{:d}-{:s} {:s}'.format(t[4], months.index(t[1])+1, t[2], t[3]))
    except ValueError as e:
        sys.stderr.write('{:s}: {:s}\n'.format(dbd_obj['dbd_meta']['filename_label'], e))
        sys.stderr.flush()
        return doc
        
    # Parse the glider name from dbd_obj['dbd_meta']['filename']
    dbd_re = re.compile('^(.*)\-\d{4}\-\d{1,3}\-\d{1,3}\-\d{1,3}$')
    m = dbd_re.search(dbd_obj['dbd_meta']['filename'])
    if not m:
        sys.stderr.write('Failed to parse glider name: {:s}\n'.format(dbd_obj['dbd_meta']['filename']))
        return doc
        
    glider = m.groups()[0]
    
    dbd_obj['dbd_meta']['fileopen_time'] = fileopen_time
    #doc['dbd'] = dbd_obj['dbd_meta']
    
    doc['rows'] = [{'dbd' : dbd_obj['dbd_meta'], 'glider' : glider, 'sensors' : r} for r in dbd_obj['data']]
    
    return doc
    
def dbd2doc_json(dbd_obj):
    
    doc = {'dbd' : None,
        'rows' : []}
    
    # Convert dbd_obj['dbd_meta']['fileopen_time'] to datetime
    #t = dbd_obj['dbd_meta']['fileopen_time'].split('_')
    #fileopen_time = parser.parse(
    #    '{:s}-{:d}-{:s} {:s}'.format(t[4], months.index(t[1])+1, t[2], t[3]))
    # Parse the glider name from dbd_obj['dbd_meta']['filename']
    dbd_re = re.compile('^(.*)\-\d{4}\-\d{1,3}\-\d{1,3}\-\d{1,3}$')
    m = dbd_re.search(dbd_obj['dbd_meta']['filename'])
    if not m:
        sys.stderr.write('Failed to parse glider name: {:s}\n'.format(dbd_obj['dbd_meta']['filename']))
        return doc
        
    glider = m.groups()[0]
    
    doc['rows'] = [{'dbd_id' : None, 'glider' : glider, 'sensors' : r} for r in dbd_obj['data']]
    
    #dbd_obj['dbd_meta']['fileopen_time'] = fileopen_time
    doc['dbd'] = dbd_obj['dbd_meta']
    
    return doc

    
def parse_dbd_header(dbd_file):
    '''Parse an ascii Slocum Dinkum Binary Data file for file metadata
    Args:
        dbd_file: the Slocum Dinkum Binary Data file to parse
    Returns:
        A dictionary containing the file's metadata'''
    
    if not os.path.exists(dbd_file):
        sys.stderr.write('File does not exist: {:s}\n'.format(dbd_file))
        return
        
    num_regex = re.compile('^\d{1,7}$')
    try:
        with open(dbd_file, 'r') as fid:
            dbd_header = {}
            for f in fid:
                
                tokens = f.strip().split(': ')
                if len(tokens) != 2:
                    return dbd_header
                    
                v = tokens[1]
                match = num_regex.match(v)
                if match:
                    v = float(v)
                
                dbd_header[tokens[0]] = v
    except IOError as e:
        sys.stderr.write('{:s}\n'.format(e))
        return
    
    return dbd_header
            
def parse_dbd_sensor_meta(dbd_file):
    '''Parse a Slocum Dinkum Binary Data file for sensor metadata
    Args:
        dbd_file: the Slocum Dinkum Binary Data file to parse
    Returns:
        An array of dictionaries containing the file's sensor metadata'''
    
    if not os.path.exists(dbd_file):
        sys.stderr.write('File does not exist: {:s}\n'.format(dbd_file))
        return {}
        
    # Parse the file header lines
    dbd_meta = parse_dbd_header(dbd_file)
    num_header_lines = int(dbd_meta['num_ascii_tags'])
    
    fid = open(dbd_file, 'r')
    
    sensors = []
    line_count = 0
    while line_count < num_header_lines:
        
        fid.readline()
        line_count += 1
    
    # Get the sensor names line
    sensors_line = fid.readline().strip()
    # Get the sensor units line
    units_line = fid.readline().strip()
    # Get the datatype byte storage information
    bytes_line = fid.readline().strip()
    
    fid.close()
    
    sensors = sensors_line.split()
    units = units_line.split()
    datatype_bytes = bytes_line.split()
    
    return [{'sensor': sensors[i], 'units' : units[i], 'bytes' : int(datatype_bytes[i])} for i in range(len(sensors))]
            
def parse_dbd_data(dbd_file):
    '''Parse a Slocum Dinkum Binary Data file for sensor data values
    Args:
        dbd_file: the Slocum Dinkum Binary Data file to parse
    Returns:
        An array containing one or more arrays of dictionaries containing sensor
        data for each row in the data file.  Sensor data containing NaN values is 
        discarded'''
    
    if not os.path.exists(dbd_file):
        sys.stderr.write('File does not exist: {:s}\n'.format(dbd_file))
        return
        
    # Parse the file header lines
    dbd_meta = parse_dbd_header(dbd_file)
    num_header_lines = int(dbd_meta['num_ascii_tags'])
    num_label_lines = int(dbd_meta['num_label_lines'])
    # Total number of header lines before the data matrix starts
    total_header_lines = num_header_lines + num_label_lines
    
    # Parse the sensor header lines
    sensors = parse_dbd_sensor_meta(dbd_file)
    
    
    data_table = np.loadtxt(dbd_file, skiprows=total_header_lines)
    
    data = []
    for r in range(data_table.shape[0]):
        row = [{'sensor' : sensors[i]['sensor'], 'units' : sensors[i]['units'], 'value' : data_table[r,i]} for i in range(data_table[r,:].size) if not np.isnan(data_table[r,i])]
        data.append(row)
        
    return data
