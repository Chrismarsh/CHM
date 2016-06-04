/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/meteoLaws/Meteoconst.h> //for PI
#include <meteoio/MathOptim.h>

#include <string.h>
#include <algorithm>

#include <meteoio/plugins/ARPSIO.h>

using namespace std;

namespace mio {
/**
 * @page arps ARPSIO
 * @section arps_format Format
 * This is for reading grid data in the ARPS grid format (it transparently supports both true ARPS ascii grids and grids modified by the ARPSGRID utility). DEM reading works well while reading meteo parameters might be a rough ride (since ARPS files do not always contain a consistent set of meteo fields).
 *
 * @section arps_units Units
 * All units are assumed to be MKSA.
 *
 * @section arps_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
 * - DEMFILE: path and file containing the DEM; [Input] section
 * - ARPS_XCOORD: x coordinate of the lower left corner of the grids; [Input] section
 * - ARPS_YCOORD: y coordinate of the lower left corner of the grids; [Input] section
 * - GRID2DPATH: path to the input directory where to find the arps files to be read as grids; [Input] section
 * - GRID2DEXT: arps file extension, or <i>none</i> for no file extension (default: .asc)
 */

const double ARPSIO::plugin_nodata = -999.; //plugin specific nodata value
const char* ARPSIO::default_ext=".asc"; //filename extension

ARPSIO::ARPSIO(const std::string& configfile)
        : cfg(configfile),
          fin(NULL), filename(),  coordin(), coordinparam(), coordout(), coordoutparam(),
          grid2dpath_in(), ext(default_ext), dimx(0), dimy(0), dimz(0), cellsize(0.),
          xcoord(IOUtils::nodata), ycoord(IOUtils::nodata), zcoord(), is_true_arps(true)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	setOptions();
}

ARPSIO::ARPSIO(const Config& cfgreader)
        : cfg(cfgreader),
          fin(NULL), filename(),  coordin(), coordinparam(), coordout(), coordoutparam(),
          grid2dpath_in(), ext(default_ext), dimx(0), dimy(0), dimz(0), cellsize(0.),
          xcoord(IOUtils::nodata), ycoord(IOUtils::nodata), zcoord(), is_true_arps(true)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	setOptions();
}

ARPSIO& ARPSIO::operator=(const ARPSIO& source) {
	if(this != &source) {
		fin = NULL;
		filename = source.filename;
		coordin = source.coordin;
		coordinparam = source.coordinparam;
		coordout = source.coordout;
		coordoutparam = source.coordoutparam;
		grid2dpath_in = source. grid2dpath_in;
		ext = source.ext;
		dimx = source.dimx;
		dimy = source.dimy;
		dimz = source.dimz;
		cellsize = source.cellsize;
		xcoord = source.xcoord;
		ycoord = source.ycoord;
		zcoord = source.zcoord;
		is_true_arps = source.is_true_arps;
	}
	return *this;
}

void ARPSIO::setOptions()
{
	string tmp;
	cfg.getValue("GRID2D", "Input", tmp, IOUtils::nothrow);
	if (tmp == "ARPS") { //keep it synchronized with IOHandler.cc for plugin mapping!!
		cfg.getValue("GRID2DPATH", "Input", grid2dpath_in);
	}

	cfg.getValue("ARPS_XCOORD", "Input", xcoord, IOUtils::dothrow);
	cfg.getValue("ARPS_YCOORD", "Input", ycoord, IOUtils::dothrow);

	//default value has been set in constructor
	cfg.getValue("GRID2DEXT", "Input", ext, IOUtils::nothrow);
	if(ext=="none") ext.clear();
}

ARPSIO::~ARPSIO() throw()
{
	cleanup();
}

void ARPSIO::read2DGrid(Grid2DObject& grid_out, const std::string& i_name)
{
	const std::string _filename = grid2dpath_in +"/" + i_name;

	openGridFile(_filename);

	const unsigned int layer=2;
	if(is_true_arps)
		readGridLayer("zp coordinat", layer, grid_out);
	else
		readGridLayer("zp_coordinat", layer, grid_out);

	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARPSIO::read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	std::string date_str = date.toString(Date::ISO);
	std::replace( date_str.begin(), date_str.end(), ':', '.');
	const std::string name = grid2dpath_in + "/" + date_str + ext;
	openGridFile(name);

	//Radiation parameters
	if(parameter==MeteoGrids::ISWR) readGridLayer("radsw", 2, grid_out);
	if(parameter==MeteoGrids::RSWR) {
		Grid2DObject net;
		readGridLayer("radsw", 2, grid_out);
		readGridLayer("radswnet", 2, net);
		grid_out.grid2D -= net.grid2D;
	}
	if(parameter==MeteoGrids::ILWR) readGridLayer("radlwin", 2, grid_out);
	if(parameter==MeteoGrids::ALB) {
		Grid2DObject rswr, iswr;
		readGridLayer("radsw", 2, iswr);
		readGridLayer("radswnet", 2, grid_out); //net radiation
		rswr.grid2D = iswr.grid2D - grid_out.grid2D;
		grid_out.grid2D = iswr.grid2D/rswr.grid2D;
	}

	//Wind grids
	if(parameter==MeteoGrids::U) readGridLayer("u", 2, grid_out);
	if(parameter==MeteoGrids::V) readGridLayer("v", 2, grid_out);
	if(parameter==MeteoGrids::W) readGridLayer("w", 2, grid_out);
	if(parameter==MeteoGrids::VW) {
		Grid2DObject V;
		readGridLayer("u", 2, grid_out); //U
		readGridLayer("v", 2, V);
		for(size_t jj=0; jj<grid_out.getNy(); jj++) {
			for(size_t ii=0; ii<grid_out.getNx(); ii++) {
				grid_out(ii,jj) = sqrt( Optim::pow2(grid_out(ii,jj)) + Optim::pow2(V(ii,jj)) );
			}
		}
	}
	if(parameter==MeteoGrids::DW) {
		Grid2DObject V;
		readGridLayer("u", 2, grid_out); //U
		readGridLayer("v", 2, V);
		for(size_t jj=0; jj<grid_out.getNy(); jj++) {
			for(size_t ii=0; ii<grid_out.getNx(); ii++) {
				grid_out(ii,jj) = fmod( atan2( grid_out(ii,jj), V(ii,jj) ) * Cst::to_deg + 360., 360.); // turn into degrees [0;360)
			}
		}
	}

	//Basic meteo parameters
	if(parameter==MeteoGrids::P) readGridLayer("p", 2, grid_out);
	if(parameter==MeteoGrids::TSG) readGridLayer("tsoil", 2, grid_out); //or is it tss for us?
	/*if(parameter==MeteoGrids::RH) {
		//const double epsilon = Cst::gaz_constant_dry_air / Cst::gaz_constant_water_vapor;
		readGridLayer("qv", 2, grid_out); //water vapor mixing ratio
		//Atmosphere::waterSaturationPressure(T);
		//HACK: compute relative humidity out of it!
		//through potential temperature -> local temperature?
	}*/

	//Hydrological parameters
	if(parameter==MeteoGrids::HS) readGridLayer("snowdpth", 2, grid_out);
	if(parameter==MeteoGrids::PSUM) {
		readGridLayer("prcrate1", 2, grid_out); //in kg/m^2/s
		grid_out.grid2D *= 3600.; //we need kg/m^2/h
	}

	//DEM
	std::string dem_marker="zp coordinat";
	if(!is_true_arps) dem_marker="zp_coordinat";
	if(parameter==MeteoGrids::DEM) readGridLayer(dem_marker, 2, grid_out);
	if(parameter==MeteoGrids::SLOPE) {
		DEMObject dem;
		dem.setUpdatePpt(DEMObject::SLOPE);
		readGridLayer(dem_marker, 2, dem);
		dem.update();
		grid_out.set(dem.cellsize, dem.llcorner, dem.slope);
	}
	if(parameter==MeteoGrids::AZI) {
		DEMObject dem;
		dem.setUpdatePpt(DEMObject::SLOPE);
		readGridLayer(dem_marker, 2, dem);
		dem.update();
		grid_out.set(dem.cellsize, dem.llcorner, dem.azi);
	}

	if(grid_out.empty()) {
		ostringstream ss;
		ss << "No suitable data found for parameter " << MeteoGrids::getParameterName(parameter) << " ";
		ss << "at time step " << date.toString(Date::ISO) << " in file \"" << name << "\"";
		throw NoAvailableDataException(ss.str(), AT);
	}

	rewind(fin);
}

void ARPSIO::read3DGrid(Grid3DObject& grid_out, const std::string& i_name)
{
	const std::string _filename = grid2dpath_in + "/" + i_name;
	openGridFile(_filename);

	//resize the grid just in case
	grid_out.grid3D.resize(dimx, dimy, dimz);

	// Read until the parameter is found
	std::string parameter; //HACK
	moveToMarker(parameter);

	//read the data we are interested in
	for (size_t ix = 0; ix < dimx; ix++) {
		for (size_t iy = 0; iy < dimy; iy++) {
			for (size_t iz = 0; iz < dimz; iz++) {
				double tmp;
				if(fscanf(fin," %16lf%*[\n]",&tmp)==1) {
					grid_out.grid3D(ix,iy,iz) = tmp;
				} else {
					cleanup();
					throw InvalidFormatException("Failure in reading 3D grid in file "+_filename, AT);
				}
			}
		}
	}

	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARPSIO::readDEM(DEMObject& dem_out)
{
	std::string _filename;
	cfg.getValue("DEMFILE", "Input", _filename);
	openGridFile(_filename);
	if(is_true_arps) {
		readGridLayer(std::string("zp coordinat"), 2 ,dem_out);
	} else {
		readGridLayer(std::string("zp_coordinat"), 2 ,dem_out);
	}
}

void ARPSIO::readLanduse(Grid2DObject& /*landuse_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARPSIO::readAssimilationData(const Date& /*date_in*/, Grid2DObject& /*da_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARPSIO::readStationData(const Date&, std::vector<StationData>& /*vecStation*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARPSIO::readMeteoData(const Date& /*dateStart*/, const Date& /*dateEnd*/,
                           std::vector< std::vector<MeteoData> >& /*vecMeteo*/,
                           const size_t&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARPSIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& /*vecMeteo*/,
                            const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARPSIO::readPOI(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARPSIO::write2DGrid(const Grid2DObject& /*grid_in*/, const std::string& /*name*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARPSIO::write2DGrid(const Grid2DObject&, const MeteoGrids::Parameters&, const Date&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARPSIO::initializeGRIDARPS()
{
	double v1, v2;

	//go to read the sizes
	moveToMarker("nnx");
	//finish reading the line and move to the next one
	if(fscanf(fin,"%*[^\n]")!=0) {
		cleanup();
		throw InvalidFormatException("Error in file format of file "+filename, AT);
	}
	if (fscanf(fin," %u %u %u \n",&dimx,&dimy,&dimz)!=3) {
		cleanup();
		throw InvalidFormatException("Can not read dimx, dimy, dimz from file "+filename, AT);
	}
	if (dimx==0 || dimy==0 || dimz==0) {
		cleanup();
		throw IndexOutOfBoundsException("Invalid dimx, dimy, dimz from file "+filename, AT);
	}

	//initializing cell size
	moveToMarker("x_coordinate");
	if (fscanf(fin,"%lg %lg",&v1,&v2)!=2) {
		cleanup();
		throw InvalidFormatException("Can not read first two x coordinates from file "+filename, AT);
	}
	const double cellsize_x = v2 - v1;
	moveToMarker("y_coordinate");
	if (fscanf(fin,"%lg %lg",&v1,&v2)!=2) {
		cleanup();
		throw InvalidFormatException("Can not read first two y coordinates from file "+filename, AT);
	}
	const double cellsize_y = v2 - v1;
	if(cellsize_x!=cellsize_y) {
		cleanup();
		throw InvalidFormatException("Only square cells currently supported! Non compliance in file "+filename, AT);
	}
	cellsize = cellsize_y;

	//HACK: zcoords must be read from zp. But they are NOT constant...

}

void ARPSIO::initializeTrueARPS(const char curr_line[ARPS_MAX_LINE_LENGTH])
{
	double v1, v2;

	//go to read the sizes
	if (sscanf(curr_line," nx = %u, ny = %u, nz = %u ",&dimx,&dimy,&dimz)!=3) {
		cleanup();
		throw InvalidFormatException("Can not read dimx, dimy, dimz from file "+filename, AT);
	}
	if (dimx==0 || dimy==0 || dimz==0) {
		cleanup();
		throw IndexOutOfBoundsException("Invalid dimx, dimy, dimz from file "+filename, AT);
	}

	//initializing cell size
	moveToMarker("x coordinate");
	if (fscanf(fin,"%lg %lg",&v1,&v2)!=2) {
		cleanup();
		throw InvalidFormatException("Can not read first two x coordinates from file "+filename, AT);
	}
	const double cellsize_x = v2 - v1;
	moveToMarker("y coordinate");
	if (fscanf(fin,"%lg %lg",&v1,&v2)!=2) {
		cleanup();
		throw InvalidFormatException("Can not read first two y coordinates from file "+filename, AT);
	}
	const double cellsize_y = v2 - v1;
	if(cellsize_x!=cellsize_y) {
		cleanup();
		throw InvalidFormatException("Only square cells currently supported! Non compliance in file "+filename, AT);
	}
	cellsize = cellsize_y;

	moveToMarker("z coordinate");
	while (fscanf(fin,"%lg",&v1)==1) {
		zcoord.push_back( v1 );
	}
	if(zcoord.size()!=dimz) {
		ostringstream ss;
		ss << "Expected " << dimz << " z coordinates in file \""+filename+"\", found " << zcoord.size();
		cleanup();
		throw InvalidFormatException(ss.str(), AT);
	}
}

void ARPSIO::openGridFile(const std::string& in_filename)
{
	unsigned int v1;
	filename = in_filename;

	if (!IOUtils::fileExists(filename)) throw FileAccessException(filename, AT); //prevent invalid filenames
	if((fin=fopen(filename.c_str(),"r")) == NULL) {
		cleanup();
		throw FileAccessException("Can not open file "+filename, AT);
	}

	//identify if the file is an original arps file or a file modified by ARPSGRID
	char dummy[ARPS_MAX_LINE_LENGTH];
	for (int j=0; j<5; j++) {
		//the first easy difference in the structure happens at line 5
		if(fgets(dummy,ARPS_MAX_STRING_LENGTH,fin)==NULL) {
			cleanup();
			throw InvalidFormatException("Fail to read header lines of file "+filename, AT);
		}
	}
	if (sscanf(dummy," nx = %u, ny = ", &v1)<1) {
		//this is an ASCII file modified by ARPSGRID
		is_true_arps=false;
		initializeGRIDARPS();
	} else {
		//this is a true ARPS file
		initializeTrueARPS(dummy);
	}

	//come back to the begining of the file
	rewind(fin);
}

void ARPSIO::cleanup() throw()
{
	if (fin!=NULL) {//close fin if open
		fclose(fin);
		fin=NULL;
	}
	/*if (fin.is_open()) {//close fin if open
		fin.close();
	}*/

	zcoord.clear();
}

/** @brief Read a specific layer for a given parameter from the ARPS file
 * @param parameter The parameter to extract. This could be any of the following:
 *        - x_coordinate for getting the X coordinates of the mesh
 *        - y_coordinate for getting the Y coordinates of the mesh
 *        - zp_coordinat for getting the Z coordinates of the mesh
 *        - u for getting the u component of the wind field
 *        - v for getting the v component of the wind field
 *        - w for getting the w component of the wind field
 * @param layer     Index of the layer to extract (1 to dimz)
 * @param grid      [out] grid containing the values. The grid will be resized if necessary.
*/
void ARPSIO::readGridLayer(const std::string& parameter, const unsigned int& layer, Grid2DObject& grid)
{
	if(layer<1 || layer>dimz) {
		cleanup();
		ostringstream tmp;
		tmp << "Layer " << layer << " does not exist in ARPS file " << filename << " (nr layers=" << dimz << ")";
		throw IndexOutOfBoundsException(tmp.str(), AT);
	}

	//resize the grid just in case
	Coords llcorner(coordin, coordinparam);
	llcorner.setXY(xcoord, ycoord, IOUtils::nodata);
	grid.set(dimx, dimy, cellsize, llcorner);

	// Read until the parameter is found
	moveToMarker(parameter);

	// move to the begining of the layer of interest
	if(layer>1) {
		double tmp;
		const size_t jmax=dimx*dimy*(layer-1);
		for (size_t j = 0; j < jmax; j++)
			if(fscanf(fin," %16lf%*[\n]",&tmp)==EOF) {
				cleanup();
				throw InvalidFormatException("Fail to skip data layers in file "+filename, AT);
			}
	}

	//read the data we are interested in
	for (size_t iy = 0; iy < dimy; iy++) {
		for (size_t ix = 0; ix < dimx; ix++) {
			double tmp;
			if(fscanf(fin," %16lf%*[\n]",&tmp)==1) {
				grid(ix,iy) = tmp;
			} else {
				cleanup();
				throw InvalidFormatException("Fail to read data layer in file "+filename, AT);
			}
		}
	}
}

void ARPSIO::moveToMarker(const std::string& marker)
{
	char dummy[ARPS_MAX_LINE_LENGTH];
	int nb_elems=0;
	do {
		nb_elems=fscanf(fin," %[^\t\n] ",dummy); //HACK: possible buffer overflow
	} while (!feof(fin) && strcmp(dummy,marker.c_str()) != 0 && nb_elems!=0);
	if(feof(fin)) {
		cleanup();
		const std::string message = "End of file "+filename+" should NOT have been reached when looking for "+marker;
		throw InvalidFormatException(message, AT);
	}
	if(nb_elems==0) {
		cleanup();
		const std::string message = "Matching failure in file "+filename+" when looking for "+marker;
		throw InvalidFormatException(message, AT);
	}
}

} //namespace
