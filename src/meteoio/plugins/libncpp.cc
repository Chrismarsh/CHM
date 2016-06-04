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
	attr_value = string(value);

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

		const string varname(name);
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

		const string dimname(name);
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
		int dimvarid;
		size_t length=0;
		char name[NC_MAX_NAME+1];

		status = nc_inq_dimname(ncid, dimids[ii], name);
		if (status != NC_NOERR) throw IOException(nc_strerror(status), AT);

		status = nc_inq_dimlen(ncid, dimids[ii], &length);
		if (status != NC_NOERR) throw IOException("Could not read dimension length for '" + string(name)  + "':" + nc_strerror(status), AT);

		status = nc_inq_varid(ncid, name, &dimvarid);
		if (status != NC_NOERR)
			throw IOException("Could not retrieve varid for variable '" + string(name) + "': " + nc_strerror(status), AT);

		dimid.push_back(dimids[ii]);
		dim_varid.push_back(dimvarid);
		dimname.push_back(string(name));
		dimlen.push_back(length);
	}
}

void read_data_2D(const int& ncid, const std::string& varname, const int& varid,
                  const size_t& record, const size_t& nr_of_records, const size_t& length, double*& data)
{
	size_t start[] = {record, 0};
	size_t count[] = {nr_of_records, length};

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
	size_t index[] = {pos};

	const int status = nc_get_var1_double(ncid, varid, index, &data);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve data for variable '" + varname + "': " + nc_strerror(status), AT);

}

void read_data(const int& ncid, const std::string& varname, const int& varid,
               const size_t& pos, const size_t& latlen, const size_t& lonlen, double*& data)
{
	size_t start[] = {pos, 0, 0};
	size_t count[] = {1, latlen, lonlen};

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

void write_data(const int& ncid, const std::string& varname, const int& varid, const double * const data)
{
	const int status = nc_put_var_double(ncid, varid, data);
	if (status != NC_NOERR)
		throw IOException("Could not write data for variable '" + varname + "': " + nc_strerror(status), AT);
}

void write_data(const int& ncid, const std::string& varname, const int& varid, const size_t& nrows, const size_t& ncols,
                const size_t& pos_start, const double * const data)
{
	size_t start[] = {pos_start, 0, 0};
	size_t count[] = {1, nrows, ncols};

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
	size_t start[] = {pos_start, 0, 0};
	size_t count[] = {1, nrows, ncols};

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

bool get_recordMinMax(const int& ncid, const std::string& varname, const int& varid, double &min, double &max)
{
	int dimid;
	size_t dimlen;

	get_dimension(ncid, varname, dimid, dimlen);
	
	//check if record already exists
	if (dimlen > 0) {
		double *record_value = new double[dimlen];
		read_data(ncid, varname, varid, record_value);

		min = record_value[0];
		max =  record_value[dimlen-1];

		delete[] record_value;
	} else 
		return false; // data not found
		
	return true;
}

// Finding a certain record variable value (e.g. timestamp) by retrieving all
// record values and then performing a linear search
size_t find_record(const int& ncid, const std::string& varname, const int& varid, const double& data)
{
	int dimid;
	size_t dimlen;

	get_dimension(ncid, varname, dimid, dimlen);

	//check if record already exists
	if (dimlen > 0) {
		double *record_value = new double[dimlen];
		read_data(ncid, varname, varid, record_value);

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

// In case the dimension length of the record variable is less than start_pos
// values will be added (containing the _FillValue) until a length of start_pos-1
// has been reached. Finally the length amount elements from start_pos and on
// will be added.
void write_record(const int& ncid, const std::string& varname, const int& varid, const size_t& start_pos, const size_t& length, const double * const data)
{
	size_t start[] = {start_pos};
	size_t count[] = {length};

	const int status = nc_put_vara_double(ncid, varid, start, count, data);
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
	double resampling_factor_x = IOUtils::nodata, resampling_factor_y=IOUtils::nodata;
	const double cellsize = calculate_cellsize(latlen, lonlen, lat, lon, resampling_factor_x, resampling_factor_y);
	const double cntr_lat = .5*(lat[0]+lat[latlen-1]);
	const double cntr_lon = .5*(lon[0]+lon[lonlen-1]);

	mio::Coords cntr(coordin, coordinparam);
	cntr.setLatLon(cntr_lat, cntr_lon, IOUtils::nodata); //it will be moved to llcorner later, after correcting the aspect ratio
	grid_out.set(lonlen, latlen, cellsize, cntr);
	
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
	
	if (resampling_factor_x != mio::IOUtils::nodata) {
		grid_out.grid2D = mio::ResamplingAlgorithms2D::BilinearResampling(grid_out.grid2D, resampling_factor_x, resampling_factor_y);
	}
	
	//computing lower left corner by using the center point as reference, AFTER we corrected for the aspect ratio
	grid_out.llcorner.moveByXY(-.5*(double)grid_out.getNx()*cellsize, -.5*(double)grid_out.getNy()*cellsize);
}

/* The Grid2DObject holds data and meta data for quadratic cells. However the NetCDF file
 * stores the grid as discrete latitude and longitude values. It is necessary to calculate
 * the distance between the edges of the grid and determine the cellsize. This cellsize may
 * be different for X and Y directions. We then choose one cellsize for our grid and
 * determine a factor that will be used for resampling the grid to likewise consist of
 * quadratic cells.
 */
double calculate_cellsize(const size_t& latlen, const size_t& lonlen, const double * const lat, const double * const lon,
                                    double& factor_x, double& factor_y)
{
	const double cntr_lat = .5*(lat[0]+lat[latlen-1]);
	const double cntr_lon = .5*(lon[0]+lon[lonlen-1]);
	double alpha;

	const double distanceX = mio::Coords::VincentyDistance(cntr_lat, lon[0], cntr_lat, lon[lonlen-1], alpha);
	const double distanceY = mio::Coords::VincentyDistance(lat[0], cntr_lon, lat[latlen-1], cntr_lon, alpha);

	// lonlen, latlen are decremented by 1; n linearly connected points have (n-1) connections
	const double cellsize_x = distanceX / static_cast<double>(lonlen-1);
	const double cellsize_y = distanceY / static_cast<double>(latlen-1);

	// round to 1cm precision for numerical stability
	const double cellsize = static_cast<double>(Optim::round( std::min(cellsize_x, cellsize_y)*100. )) / 100.;

	if (cellsize_x == cellsize_y) {
		return cellsize_x;
	} else {
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
	lat_array[0] = grid.llcorner.getLat();
	lon_array[0] = grid.llcorner.getLon();

	// The idea is to use the difference in coordinates of the upper right and the lower left
	// corner to calculate the lat/lon intervals between cells
	Coords urcorner(grid.llcorner);
	urcorner.setGridIndex(static_cast<int>(grid.getNx() - 1), static_cast<int>(grid.getNy() - 1), IOUtils::nodata, true);
	grid.gridify(urcorner);

	const double lat_interval = (urcorner.getLat() - lat_array[0]) / static_cast<double>(grid.getNy()-1);
	const double lon_interval = (urcorner.getLon() - lon_array[0]) / static_cast<double>(grid.getNx()-1);

	// The method to use interval*ii is consistent with the corresponding
	// calculation of the Grid2DObject::gridify method -> numerical stability
	for (size_t ii=1; ii<grid.getNy(); ++ii) {
		lat_array[ii] = lat_array[0] + lat_interval*static_cast<double>(ii);
	}

	for (size_t ii=1; ii<grid.getNx(); ++ii) {
		lon_array[ii] = lon_array[0] + lon_interval*static_cast<double>(ii);
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

