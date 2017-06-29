/***********************************************************************************/
/*  Copyright 2014 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <meteoio/plugins/libncpp.h>
#include <meteoio/MathOptim.h>
#include <meteoio/ResamplingAlgorithms2D.h>
#include <meteoio/dataClasses/Coords.h>
#include <meteoio/dataClasses/CoordsAlgorithms.h>
#include <meteoio/IOUtils.h>
#include <meteoio/IOExceptions.h>

#include <algorithm>

using namespace std;
using namespace mio;  // for the IOExceptions and IOUtils

namespace ncpp {

//
// NetCDF C Library wrappers
//
void open_file(const std::string& filename, const int& omode, int& ncid)
{
	const int status = nc_open(filename.c_str(), omode, &ncid);
	if (status != NC_NOERR)
		throw IOException("Could not open netcdf file '" + filename + "': " + nc_strerror(status), AT);
}

void create_file(const std::string& filename, const int& cmode, int& ncid)
{
	const int status = nc_create(filename.c_str(), cmode, &ncid);
	if (status != NC_NOERR)
		throw IOException("Could not create netcdf file '" + filename + "': " + nc_strerror(status), AT);
}

void get_variable(const int& ncid, const std::string& varname, int& varid)
{
	const int status = nc_inq_varid(ncid, varname.c_str(), &varid);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve varid for variable '" + varname + "': " + nc_strerror(status), AT);
}

void get_unlimited_dimname(const int& ncid, std::string& dimname)
{
	int unlim_id;
	const int status = nc_inq_unlimdim(ncid, &unlim_id);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve unlimited dimension: " + std::string(nc_strerror(status)), AT);
	if (unlim_id==-1) { //no unlimited dimension has been defined
		dimname="";
		return;
	}

	char name[NC_MAX_NAME+1];
	const int status2 = nc_inq_dimname (ncid, unlim_id, name);
	if (status2 != NC_NOERR)
		throw IOException("Could not retrieve the name of the unlimited dimension: " + std::string(nc_strerror(status)), AT);

	dimname = std::string( name );
}

void get_dimension(const int& ncid, const std::string& dimname, int& dimid)
{
	const int status = nc_inq_dimid(ncid, dimname.c_str(), &dimid);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve dimid for dimension '" + dimname + "': " + nc_strerror(status), AT);
}

void get_dimension(const int& ncid, const std::string& dimname, int& dimid, size_t& dimlen)
{
	get_dimension(ncid, dimname, dimid);

	const int status = nc_inq_dimlen(ncid, dimid, &dimlen);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve length for dimension '" + dimname + "': " + nc_strerror(status), AT);
}

std::string get_DimAttribute(const int& ncid, const std::string& dimname, const std::string& attr_name)
{
	int dimid;
	get_dimension(ncid, dimname, dimid);
	
	size_t attr_len;
	int status = nc_inq_attlen (ncid, dimid, attr_name.c_str(), &attr_len);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve attribute '" + attr_name + "' for var '" + dimname + "': " + nc_strerror(status), AT);

	char* value = new char[attr_len + 1]; // +1 for trailing null
	status = nc_get_att_text(ncid, dimid, attr_name.c_str(), value);
	if (status != NC_NOERR)
		throw IOException("Could not read attribute '" + attr_name + "' for var '" + dimname + "': " + nc_strerror(status), AT);

	value[attr_len] = '\0';
	const std::string attr_value(value);
	delete[] value;
	return attr_value;
}

void get_VarAttribute(const int& ncid, const std::string& varname, const std::string& attr_name, std::string& attr_value)
{
	int varid;
	get_variable(ncid, varname, varid);
	
	size_t attr_len;
	int status = nc_inq_attlen (ncid, varid, attr_name.c_str(), &attr_len);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve attribute '" + attr_name + "' for var '" + varname + "': " + nc_strerror(status), AT);

	char* value = new char[attr_len + 1]; // +1 for trailing null
	status = nc_get_att_text(ncid, varid, attr_name.c_str(), value);
	if (status != NC_NOERR)
		throw IOException("Could not read attribute '" + attr_name + "' for var '" + varname + "': " + nc_strerror(status), AT);

	value[attr_len] = '\0';
	attr_value = string(value);

	delete[] value;
}

void get_attribute(const int& ncid, const std::string& varname, const int& varid, const std::string& attr_name, std::string& attr_value)
{
	size_t attr_len;

	int status = nc_inq_attlen (ncid, varid, attr_name.c_str(), &attr_len);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve attribute '" + attr_name + "' for var '" + varname + "': " + nc_strerror(status), AT);

	char* value = new char[attr_len + 1]; // +1 for trailing null
	status = nc_get_att_text(ncid, varid, attr_name.c_str(), value);
	if (status != NC_NOERR)
		throw IOException("Could not read attribute '" + attr_name + "' for var '" + varname + "': " + nc_strerror(status), AT);

	value[attr_len] = '\0';
	attr_value = std::string(value);

	delete[] value;
}

void get_attribute(const int& ncid, const std::string& varname, const int& varid, const std::string& attr_name, double& attr_value)
{
	size_t attr_len;

	int status = nc_inq_attlen (ncid, varid, attr_name.c_str(), &attr_len);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve attribute '" + attr_name + "' for var '" + varname + "': " + nc_strerror(status), AT);

	status = nc_get_att_double(ncid, varid, attr_name.c_str(), &attr_value);
	if (status != NC_NOERR)
		throw IOException("Could not read attribute '" + attr_name + "' for var '" + varname + "': " + nc_strerror(status), AT);
}

bool check_attribute(const int& ncid, const int& varid, const std::string& attr_name)
{
	size_t attr_len;
	const int status = nc_inq_attlen (ncid, varid, attr_name.c_str(), &attr_len);

	if (status != NC_NOERR) return false;

	return true;
}

bool check_variable(const int& ncid, const std::string& varname)
{
	int varid;
	const int status = nc_inq_varid(ncid, varname.c_str(), &varid);

	if (status != NC_NOERR) return false;

	return true;
}

bool check_dim_var(const int& ncid, const std::string& dimname)
{
	int dimid;
	const int status = nc_inq_dimid(ncid, dimname.c_str(), &dimid);
	if (status != NC_NOERR) return false;

	return check_variable(ncid, dimname);
}

// Retrieve all variables with a certain set of dimensions
void get_variables(const int& ncid, const std::vector<std::string>& dimensions, std::vector<std::string>& variables)
{
	int nr_of_variables = -1;
	int status = nc_inq_nvars(ncid, &nr_of_variables);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve variables for dataset: " + string(nc_strerror(status)), AT);

	// Variable IDs in a NetCDF file are consecutive integers starting with 0
	for (int ii=0; ii<nr_of_variables; ++ii) {
		char name[NC_MAX_NAME+1];
		const int stat = nc_inq_varname(ncid, ii, name);
		if (stat != NC_NOERR) throw IOException(nc_strerror(stat), AT);

		const std::string varname(name);
		const bool check = check_dimensions(ncid, varname, ii, dimensions);

		if (check) variables.push_back(varname);
	}
}

// NOTE: the dimension names in the vector 'names' have to have the same order
//       as the dimensions ids retrieved for the variable
bool check_dimensions(const int& ncid, const std::string& varname, const int& varid, const std::vector<std::string>& names)
{
	int dimids[NC_MAX_VAR_DIMS], ndimsp;

	const int status = nc_inq_var(ncid, varid, NULL, NULL, &ndimsp, dimids, NULL);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve dimensions for variable '" + varname + "': " + nc_strerror(status), AT);

	if ((int)names.size() != ndimsp) return false;

	for (int ii=0; ii<ndimsp; ++ii) {
		char name[NC_MAX_NAME+1];

		const int stat = nc_inq_dimname(ncid, dimids[ii], name);
		if (stat != NC_NOERR) throw IOException(nc_strerror(stat), AT);

		const std::string dimname(name);
		const bool exists = (dimname == names[ii]); //(find(names.begin(), names.end(), dimname) != names.end());

		if (!exists) return false;
	}

	return true; // dimension check successfull
}

void get_dimension(const int& ncid, const std::string& varname, const int& varid,
                   std::vector<int>& dimid, std::vector<int>& dim_varid, std::vector<std::string>& dimname, std::vector<size_t>& dimlen)
{
	dimid.clear(); dim_varid.clear(); dimname.clear(); dimlen.clear();

	int dimids[NC_MAX_VAR_DIMS], ndimsp;

	int status = nc_inq_var(ncid, varid, NULL, NULL, &ndimsp, dimids, NULL);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve dimensions for variable '" + varname + "': " + nc_strerror(status), AT);

	for (int ii=0; ii<ndimsp; ++ii) {
		char name_cstr[NC_MAX_NAME+1];
		status = nc_inq_dimname(ncid, dimids[ii], name_cstr);
		if (status != NC_NOERR) throw IOException(nc_strerror(status), AT);
		const std::string name( name_cstr );

		size_t length=0;
		status = nc_inq_dimlen(ncid, dimids[ii], &length);
		if (status != NC_NOERR) throw IOException("Could not read dimension length for '" + name  + "':" + nc_strerror(status), AT);

		int dimvarid;
		status = nc_inq_varid(ncid, name_cstr, &dimvarid);
		if (status != NC_NOERR)
			throw IOException("Could not retrieve varid for variable '" + name + "': " + nc_strerror(status), AT);

		dimid.push_back(dimids[ii]);
		dim_varid.push_back(dimvarid);
		dimname.push_back(name);
		dimlen.push_back(length);
	}
}

void get_wrf_dimension(const int& ncid, const std::string& varname, const int& varid,
                   std::vector<int>& dimid, std::vector<int>& dim_varid, std::vector<std::string>& dimname, std::vector<size_t>& dimlen)
{
	dimid.clear(); dim_varid.clear(); dimname.clear(); dimlen.clear();

	int dimids[NC_MAX_VAR_DIMS], ndimsp;
	int status = nc_inq_var(ncid, varid, NULL, NULL, &ndimsp, dimids, NULL);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve dimensions for variable '" + varname + "': " + nc_strerror(status), AT);

	for (int ii=0; ii<ndimsp; ++ii) {
		char name_cstr[NC_MAX_NAME+1];
		status = nc_inq_dimname(ncid, dimids[ii], name_cstr);
		if (status != NC_NOERR) throw IOException(nc_strerror(status), AT);
		const std::string name( name_cstr );

		size_t length=0;
		status = nc_inq_dimlen(ncid, dimids[ii], &length);
		if (status != NC_NOERR) throw IOException("Could not read dimension length for '" + name  + "':" + nc_strerror(status), AT);

		int dimvarid;
		if (name=="Time") {
			status = nc_inq_varid(ncid, "Times", &dimvarid);
		} else if (name=="south_north") {
			status = nc_inq_varid(ncid, "XLAT", &dimvarid);
		} else if (name=="west_east") {
			status = nc_inq_varid(ncid, "XLONG", &dimvarid);
		} else
			status = nc_inq_varid(ncid, name_cstr, &dimvarid);

		if (status != NC_NOERR)
			throw IOException("Could not retrieve varid for variable '" + name + "': " + nc_strerror(status), AT);

		dimid.push_back(dimids[ii]);
		dim_varid.push_back(dimvarid);
		dimname.push_back(name);
		dimlen.push_back(length);
	}
}

void read_data_2D(const int& ncid, const std::string& varname, const int& varid,
                  const size_t& record, const size_t& nr_of_records, const size_t& length, double*& data)
{
	const size_t start[] = {record, 0};
	const size_t count[] = {nr_of_records, length};

	const int status = nc_get_vara_double(ncid, varid, start, count, data);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve data for variable '" + varname + "': " + nc_strerror(status), AT);
}

void read_value(const int& ncid, const std::string& varname, const int& varid, double& data)
{
	read_value(ncid, varname, varid, 0, data);
}

void read_value(const int& ncid, const std::string& varname, const int& varid, const size_t& pos, double& data)
{
	const size_t index[] = {pos};
	const int status = nc_get_var1_double(ncid, varid, index, &data);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve data for variable '" + varname + "': " + nc_strerror(status), AT);
}

void read_data(const int& ncid, const std::string& varname, const int& varid,
               const size_t& pos, const size_t& latlen, const size_t& lonlen, double*& data)
{
	const size_t start[] = {pos, 0, 0};
	const size_t count[] = {1, latlen, lonlen};

	const int status = nc_get_vara_double(ncid, varid, start, count, data);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve data for variable '" + varname + "': " + nc_strerror(status), AT);
}

void read_data(const int& ncid, const std::string& varname, const int& varid, double*& data)
{
	const int status = nc_get_var_double(ncid, varid, data);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve data for variable '" + varname + "': " + nc_strerror(status), AT);
}

void read_wrf_latlon(const int& ncid, const int& latid, const int& lonid, const size_t& latlen, const size_t& lonlen, double*& lat, double*& lon)
{
	int ndims;
	int status = nc_inq_var (ncid, latid, 0, 0, &ndims, 0, 0);
	if ( status != NC_NOERR)
		throw IOException("Could not retrieve number of dimensions of latitudes: " + std::string(nc_strerror(status)), AT);
	if (ndims!=3)
		throw IOException("The WRF schema has been selected, but it seems that the NetCDF file does not follow this schema", AT);

	const size_t start[] = {0, 0, 0};
	const size_t count_lat[] = {1, latlen, 1};
	status = nc_get_vara_double(ncid, latid, start, count_lat, lat);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve latitudes: " + std::string(nc_strerror(status)), AT);

	const size_t count_lon[] = {1, 1, lonlen};
	status = nc_get_vara_double(ncid, lonid, start, count_lon, lon);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve latitudes: " + std::string(nc_strerror(status)), AT);
}

void write_data(const int& ncid, const std::string& varname, const int& varid, const double * const data)
{
	const int status = nc_put_var_double(ncid, varid, data);
	if (status != NC_NOERR)
		throw IOException("Could not write data for variable '" + varname + "': " + nc_strerror(status), AT);
}

void write_data(const int& ncid, const std::string& varname, const int& varid, const size_t& nrows, const size_t& ncols,
                const size_t& pos_start, const double * const data)
{
	const size_t start[] = {pos_start, 0, 0};
	const size_t count[] = {1, nrows, ncols};

	const int status = nc_put_vara_double(ncid, varid, start, count, data);
	if (status != NC_NOERR) {
		throw IOException("Could not write variable '" + varname + "': " + string(nc_strerror(status)), AT);
	}
}

void write_data(const int& ncid, const std::string& varname, const int& varid, const int * const data)
{
	const int status = nc_put_var_int(ncid, varid, data);
	if (status != NC_NOERR)
		throw IOException("Could not write data for variable '" + varname + "': " + nc_strerror(status), AT);
}

void write_data(const int& ncid, const std::string& varname, const int& varid, const size_t& nrows, const size_t& ncols,
                const size_t& pos_start, const int * const data)
{
	const size_t start[] = {pos_start, 0, 0};
	const size_t count[] = {1, nrows, ncols};

	const int status = nc_put_vara_int(ncid, varid, start, count, data);
	if (status != NC_NOERR) {
		throw IOException("Could not write variable '" + varname + "': " + string(nc_strerror(status)), AT);
	}
}

// Adding a record value (e.g. timestamp), in case it doesn't already exist and
// that the value is greater than the last record variable value. For example,
// timestamps have to be strictly monotonically increasing or already existent.
size_t add_record(const int& ncid, const std::string& varname, const int& varid, const double& data)
{
	int dimid;
	size_t dimlen;

	get_dimension(ncid, varname, dimid, dimlen);

	//check if record already exists
	if (dimlen > 0) {
		double last_value = IOUtils::nodata;
		read_value(ncid, varname, varid, dimlen-1, last_value);

		if (last_value == data) return (dimlen - 1); //The timestamp already exists

		if (last_value > data) {
			const size_t pos = find_record(ncid, varname, dimid, data); // Search for a possible match

			if (pos != IOUtils::npos) {
				return pos;
			} else {
				throw IOException("The variable '" + varname + "' has to be linearly increasing", AT);
			}
		}
	}

	write_record(ncid, varname, varid, dimlen, 1, &data);
	return dimlen;
}

std::vector<Date> read_wrf_Time(const int& ncid, const int& dimid, const size_t& dimlen)
{
	static const size_t DateStrLen = 19; //HACK DateStrLen = 19, defined in Dimensions

	char *record_value = (char*)calloc(dimlen, sizeof(char)*DateStrLen);
	const int status = nc_get_var_text(ncid, dimid, record_value);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve data for Time variable: " + std::string(nc_strerror(status)), AT);

	std::vector<Date> dates(dimlen);
	for(size_t ii=0; ii<dimlen; ii++) {
		std::string tmp(DateStrLen, '\0');
		for(size_t jj=0; jj<DateStrLen; jj++) {
			const char c = record_value[ii*DateStrLen+jj];
			tmp[jj] = (c!='_')? c : 'T';
		}
		IOUtils::convertString(dates[ii], tmp, 0.); //assuming GMT
	}
	free( record_value );
	return dates;
}

bool get_dimensionMinMax(const int& ncid, const std::string& varname, double &min, double &max)
{
	int dimid;
	size_t dimlen;
	get_dimension(ncid, varname, dimid, dimlen);
	if (dimlen<=0) return false; //data not found
	
	double *record_value = new double[dimlen];
	read_data(ncid, varname, dimid, record_value);
	min = record_value[0];
	max =  record_value[dimlen-1];

	delete[] record_value;
	return true;
}

bool get_wrf_dimensionMinMax(const int& ncid, const std::string& varname, double &min, double &max)
{
	int dimid;
	size_t dimlen;
	get_dimension(ncid, varname, dimid, dimlen);
	if (dimlen<=0) return false;

	if (varname=="Time") {
		const std::vector<Date> dates( read_wrf_Time(ncid, dimid, dimlen) );
		min = dates.front().getJulian();
		max = dates.back().getJulian();
	} else {
		double *record_value = new double[dimlen];
		read_data(ncid, varname, dimid, record_value);
		min = record_value[0];
		max =  record_value[dimlen-1];

		delete[] record_value;
	}

	return true;
}

bool get_recordMinMax(const int& ncid, const std::string& varname, const int& varid, double &min, double &max)
{
	int dimid;
	size_t dimlen;

	get_dimension(ncid, varname, dimid, dimlen);
	if (dimlen<=0) return false;

	//check if record already exists
	double *record_value = new double[dimlen];
	read_data(ncid, varname, varid, record_value);

	min = record_value[0];
	max =  record_value[dimlen-1];

	delete[] record_value;

	return true;
}

// Finding a certain record variable value (e.g. timestamp) by retrieving all
// record values and then performing a linear search
size_t find_record(const int& ncid, const std::string& varname, const double& data)
{
	int dimid;
	size_t dimlen;
	get_dimension(ncid, varname, dimid, dimlen);
	if (dimlen<=0) return IOUtils::npos; // data not found

	//check if record already exists
	double *record_value = new double[dimlen];
	read_data(ncid, varname, dimid, record_value);

	for (size_t ii=0; ii<dimlen; ii++) {
		if (record_value[ii] == data) {
			delete[] record_value;
			return ii;
		}
	}

	delete[] record_value;
	return IOUtils::npos; // data not found
}

size_t find_wrf_record(const int& ncid, const std::string& varname, const double& data)
{
	int dimid;
	size_t dimlen;
	get_dimension(ncid, varname, dimid, dimlen);
	if (dimlen<=0) return IOUtils::npos; // data not found

	if (varname=="Time") {
		const std::vector<Date> dates( read_wrf_Time(ncid, dimid, dimlen) );

		for (size_t ii=0; ii<dimlen; ii++) {
			if (dates[ii].getJulian() == data) return ii;
		}
	} else {
		double *record_value = new double[dimlen];
		read_data(ncid, varname, dimid, record_value);

		for (size_t ii=0; ii<dimlen; ii++) {
			if (record_value[ii] == data) {
				delete[] record_value;
				return ii;
			}
		}

		delete[] record_value;
	}

	return IOUtils::npos; // data not found
}

size_t find_record(const int& ncid, const std::string& varname, const int& varid, const double& data)
{
	int dimid;
	size_t dimlen;
	get_dimension(ncid, varname, dimid, dimlen);
	if (dimlen<=0) return IOUtils::npos; // data not found

	//check if record already exists
	double *record_value = new double[dimlen];
	read_data(ncid, varname, varid, record_value);

	for (size_t ii=0; ii<dimlen; ii++) {
		if (record_value[ii] == data) {
			delete[] record_value;
			return ii;
		}
	}

	delete[] record_value;
	return IOUtils::npos; // data not found
}

// In case the dimension length of the record variable is less than start_pos
// values will be added (containing the _FillValue) until a length of start_pos-1
// has been reached. Finally the length amount elements from start_pos and on
// will be added.
void write_record(const int& ncid, const std::string& varname, const int& varid, const size_t& start_pos, const size_t& length, const double * const data)
{
	const size_t start[] = {start_pos};
	const size_t count[] = {length};

	const int status = nc_put_vara_double(ncid, varid, start, count, data);
	if (status != NC_NOERR)
		throw IOException("Could not write data for record variable '" + varname + "': " + nc_strerror(status), AT);
}

void write_record(const int& ncid, const std::string& varname, const int& varid, const size_t& start_pos, const size_t& length, const int * const data)
{
	const size_t start[] = {start_pos};
	const size_t count[] = {length};

	const int status = nc_put_vara_int(ncid, varid, start, count, data);
	if (status != NC_NOERR)
		throw IOException("Could not write data for record variable '" + varname + "': " + nc_strerror(status), AT);
}

void add_dimension(const int& ncid, const std::string& dimname, const size_t& length, int& dimid)
{
	const int status = nc_def_dim(ncid, dimname.c_str(), length, &dimid);
	if (status != NC_NOERR)
		throw IOException("Could not define dimension '" + dimname + "': " + nc_strerror(status), AT);
}

void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const double& attr_value)
{
	const int status = nc_put_att_double(ncid, varid, attr_name.c_str(), NC_DOUBLE, 1, &attr_value);
	if (status != NC_NOERR)
		throw IOException("Could not add attribute '" + attr_name + "': " + nc_strerror(status), AT);
}

void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const std::string& attr_value)
{
	const int status = nc_put_att_text(ncid, varid, attr_name.c_str(), attr_value.size(), attr_value.c_str());
	if (status != NC_NOERR)
		throw IOException("Could not add attribute '" + attr_name + "': " + nc_strerror(status), AT);
}

void add_0D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, int& varid)
{
	int dimid;
	const int status = nc_def_var(ncid, varname.c_str(), xtype, 0, &dimid, &varid);
	if (status != NC_NOERR)
		throw IOException("Could not define variable '" + varname + "': " + nc_strerror(status), AT);
}

void add_1D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, const int& dimid, int& varid)
{
	const int status = nc_def_var(ncid, varname.c_str(), xtype, 1, &dimid, &varid);
	if (status != NC_NOERR)
		throw IOException("Could not define variable '" + varname + "': " + nc_strerror(status), AT);
}

void add_2D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, const int& dimid1, const int& dimid2, int& varid)
{
	vector<int> dimids;
	dimids.push_back(dimid1);
	dimids.push_back(dimid2);

	const int status = nc_def_var(ncid, varname.c_str(), xtype, 2, &dimids[0], &varid);
	if (status != NC_NOERR)
		throw IOException("Could not define variable '" + varname + "': " + nc_strerror(status), AT);
}

void add_3D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, const int& dimid_record, const int& dimid1, const int& dimid2, int& varid)
{
	vector<int> dimids;
	dimids.push_back(dimid_record); // has to be the first one, the slowest changing index
	dimids.push_back(dimid1);
	dimids.push_back(dimid2);


	const int status = nc_def_var(ncid, varname.c_str(), xtype, 3, &dimids[0], &varid);
	if (status != NC_NOERR)
		throw IOException("Could not define variable '" + varname + "': " + nc_strerror(status), AT);
}

void start_definitions(const std::string& filename, const int& ncid)
{
	const int status = nc_redef(ncid);
	if (status != NC_NOERR)
		throw IOException("Could not open define mode for file '" + filename + "': " + nc_strerror(status), AT);

}

void end_definitions(const std::string& filename, const int& ncid)
{
	const int status = nc_enddef(ncid);
	if (status != NC_NOERR)
		throw IOException("Could not close define mode for file '" + filename + "': " + nc_strerror(status), AT);

}

void close_file(const std::string& filename, const int& ncid)
{
	const int status = nc_close(ncid);
	if (status != NC_NOERR)
		throw IOException("Could not close netcdf file  '" + filename + "': " + nc_strerror(status), AT);

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Wrappers to MeteoIO's data classes
void copy_grid(const std::string& coordin, const std::string& coordinparam, const size_t& latlen, const size_t& lonlen, const double * const lat, const double * const lon,
                         const double * const grid, const double& nodata, mio::Grid2DObject& grid_out)
{
	mio::Coords llcorner(coordin, coordinparam);
	llcorner.setLatLon(lat[0], lon[0], IOUtils::nodata);
	double resampling_factor_x = IOUtils::nodata, resampling_factor_y=IOUtils::nodata;
	const double cellsize = calculate_cellsize(latlen, lonlen, lat, lon, resampling_factor_x, resampling_factor_y);
	grid_out.set(lonlen, latlen, cellsize, llcorner);
	
	//Handle the case of llcorner/urcorner swapped
	if (lat[0]<=lat[latlen-1]) {
		for (size_t kk=0; kk < latlen; kk++) {
			const size_t row = kk*lonlen;
			if (lon[0]<=lon[lonlen-1]) {
				for (size_t ll=0; ll < lonlen; ll++)
					grid_out(ll, kk) = mio::IOUtils::standardizeNodata(grid[row + ll], nodata);
			} else {
				for (size_t ll=0; ll < lonlen; ll++)
					grid_out(ll, kk) = mio::IOUtils::standardizeNodata(grid[row + (lonlen -1) - ll], nodata);
			}
		}
	} else {
		for (size_t kk=0; kk < latlen; kk++) {
			const size_t row = ((latlen-1) - kk)*lonlen;
			if (lon[0]<=lon[lonlen-1]) {
				for (size_t ll=0; ll < lonlen; ll++)
					grid_out(ll, kk) = mio::IOUtils::standardizeNodata(grid[row + ll], nodata);
			} else {
				for (size_t ll=0; ll < lonlen; ll++)
					grid_out(ll, kk) = mio::IOUtils::standardizeNodata(grid[row + (lonlen -1) - ll], nodata);
			}
		}
	}
	
	if (resampling_factor_x != mio::IOUtils::nodata || resampling_factor_y != mio::IOUtils::nodata) {
		grid_out.grid2D = mio::ResamplingAlgorithms2D::BilinearResampling(grid_out.grid2D, resampling_factor_x, resampling_factor_y);
	}
}

/* The Grid2DObject holds data and meta data for quadratic cells. However the NetCDF file
 * stores the grid as discrete latitude and longitude values. It is necessary to calculate
 * the distance between the edges of the grid and determine the cellsize. This cellsize may
 * be different for X and Y directions. We then choose one cellsize for our grid and
 * determine a factor that will be used for resampling the grid to likewise consist of
 * quadratic cells.
 */
double calculate_cellsize(const size_t& latlen, const size_t& lonlen, const double * const lat_array, const double * const lon_array,
                                    double& factor_x, double& factor_y)
{
	double alpha;
	const double cntr_lat = .5*(lat_array[0]+lat_array[latlen-1]);
	const double cntr_lon = .5*(lon_array[0]+lon_array[lonlen-1]);
	const double lat_length = CoordsAlgorithms::VincentyDistance(cntr_lat-.5, cntr_lon, cntr_lat+.5, cntr_lon, alpha);
	const double lon_length = CoordsAlgorithms::VincentyDistance(cntr_lat, cntr_lon-.5, cntr_lat, cntr_lon+.5, alpha);

	const double distanceX = (lon_array[lonlen-1] - lon_array[0]) * lon_length;
	const double distanceY = (lat_array[latlen-1] - lat_array[0]) * lat_length;

	//round to 1cm precision for numerical stability
	const double cellsize_x = static_cast<double>(Optim::round( distanceX / static_cast<double>(lonlen)*100. )) / 100.;
	const double cellsize_y = static_cast<double>(Optim::round( distanceY / static_cast<double>(latlen)*100. )) / 100.;

	if (cellsize_x == cellsize_y) {
		return cellsize_x;
	} else {
		const double cellsize = std::min(cellsize_x, cellsize_y);
		factor_x =  cellsize_x / cellsize;
		factor_y =  cellsize_y / cellsize;
		return cellsize;
	}
}

/* Fill the arrays of lat/lon with the lat/lon intervals
 * as calculated from the cellsize of the grid object.
 */
void calculate_dimensions(const mio::Grid2DObject& grid, double*& lat_array, double*& lon_array)
{
	//There is a trick here: walking along a line of constant northing does NOT lead to a constant latitude. Both grids
	//are shifted (even if a little), which means that the center of lat/lon is != center of east./north..
	//So, in order to find the center of the domain, we do a few iteration to converge toward a reasonnable approximation
	double alpha;
	double lat_length, lon_length, cntr_lat=grid.llcorner.getLat(), cntr_lon=grid.llcorner.getLon();
	for(size_t ii=0; ii<5; ii++) {
		lat_length = CoordsAlgorithms::VincentyDistance(cntr_lat-.5, cntr_lon, cntr_lat+.5, cntr_lon, alpha);
		lon_length = CoordsAlgorithms::VincentyDistance(cntr_lat, cntr_lon-.5, cntr_lat, cntr_lon+.5, alpha);
		cntr_lat = (.5*static_cast<double>(grid.getNy())*grid.cellsize) / lat_length + grid.llcorner.getLat();
		cntr_lon = (.5*static_cast<double>(grid.getNx())*grid.cellsize) / lon_length + grid.llcorner.getLon();
	}

	const double min_lat =  cntr_lat - (0.5*static_cast<double>(grid.getNy())*grid.cellsize) / lat_length;
	const double min_lon = cntr_lon - (0.5*static_cast<double>(grid.getNx())*grid.cellsize) / lon_length;
	const double max_lat = cntr_lat + (0.5*static_cast<double>(grid.getNy())*grid.cellsize) / lat_length;
	const double max_lon = cntr_lon + (0.5*static_cast<double>(grid.getNx())*grid.cellsize) / lon_length;
	const double lat_interval = abs(max_lat - min_lat);
	const double lon_interval = abs(max_lon - min_lon);

	for (size_t ii=0; ii<grid.getNy(); ++ii) {
		lat_array[ii] = min_lat + (lat_interval * static_cast<double>(ii)) / (static_cast<double>(grid.getNy())-1);
	}
	for (size_t ii=0; ii<grid.getNx(); ++ii) {
		lon_array[ii] = min_lon + (lon_interval * static_cast<double>(ii)) / (static_cast<double>(grid.getNx()-1));
	}
}

// Fill a NetCDF 2D array with the data from a Grid2DObject
void fill_grid_data(const mio::Grid2DObject& grid, double*& data)
{
	const size_t nrows = grid.getNy(), ncols = grid.getNx();
	for (size_t kk=0; kk<nrows; ++kk) {
		for (size_t ll=0; ll<ncols; ++ll) {
			data[kk*ncols + ll] = grid.grid2D(ll,kk);
		}
	}
}

void fill_grid_data(const mio::Grid2DObject& grid, const double& new_nodata, int*& data)
{
	const size_t nrows = grid.getNy(), ncols = grid.getNx();
	for (size_t kk=0; kk<nrows; ++kk) {
		for (size_t ll=0; ll<ncols; ++ll) {
			const double val = grid.grid2D(ll,kk);
			if (val!=IOUtils::nodata) 
				data[kk*ncols + ll] = static_cast<int>( Optim::round(val) );
			else 
				data[kk*ncols + ll] = static_cast<int>( new_nodata );
		}
	}
}


} //end namespace

