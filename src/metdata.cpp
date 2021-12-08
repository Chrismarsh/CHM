//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//

#include "metdata.hpp"

metdata::metdata(std::string mesh_proj4)
{
    _nc = nullptr;
    _use_netcdf = false;
    _n_timesteps = 0;
    _mesh_proj4 = mesh_proj4;
    is_first_timestep = true;

    OGRSpatialReference srs;
    srs.importFromProj4(_mesh_proj4.c_str());
    _is_geographic = srs.IsGeographic();
}

metdata::~metdata()
{

}

void metdata::load_from_netcdf(const std::string& path,std::map<std::string, boost::shared_ptr<filter_base> > filters)
{
    if(_mesh_proj4 == "")
        BOOST_THROW_EXCEPTION(forcing_error() << errstr_info( "Met loader not initialized with proj4 string" ));

    LOG_DEBUG << "Found NetCDF file";

    _use_netcdf = true;
    _nc = std::make_unique<netcdf>();

    // make a copy of the filters and track what they provide
    for(auto& itr : filters)
    {
        _netcdf_filters[itr.first] = itr.second;
        for(auto& p : itr.second->provides())
        {
            _provides_from_nc_filters.insert(p);
        }
    }

    // spatial reference conversions to ensure the virtual station coordinates are the same as the meshes'
    OGRSpatialReference insrs, outsrs;
    insrs.SetWellKnownGeogCS("EPSG:4326");
    bool err = outsrs.importFromProj4(_mesh_proj4.c_str());

    if(err)
    {
        BOOST_THROW_EXCEPTION(forcing_error() << errstr_info( "Failure importing mesh proj4 string" ));
    }
    OGRCoordinateTransformation* coordTrans =  nullptr;
    if(!_is_geographic)
    {
        coordTrans = OGRCreateCoordinateTransformation(&insrs, &outsrs);
        if(!coordTrans)
        {
            CHM_THROW_EXCEPTION(forcing_error,"Error creating CRS transform in Met loader");
        }
    }


    try
    {
        _nc->open_GEM(path);

        _variables = _nc->get_variable_names();

        _variables.insert(_provides_from_nc_filters.begin(),_provides_from_nc_filters.end());

        _start_time = _nc->get_start();
        _end_time = _nc->get_end();
        _n_timesteps = _nc->get_ntimesteps();

        _nstations = _nc->get_xsize() * _nc->get_ysize();
        _stations.resize(_nstations);

        LOG_DEBUG << "Grid is (y)" << _nc->get_ysize() << " by (x)" << _nc->get_xsize();

        LOG_DEBUG << "Loading lat/long grid...";

        auto lat = _nc->get_lat();
        auto lon = _nc->get_lon();
        auto e = _nc->get_z();

        LOG_DEBUG << "Initializing datastructure";

        // #pragma omp parallel for
        // hangs, unclear why, critical sections around the json and gdal calls
        // don't seem to help. Probably _dD_tree is not thread safe
        for (size_t y = 0; y < _nc->get_ysize(); y++)
        {
            for (size_t x = 0; x < _nc->get_xsize(); x++)
            {
                double latitude = 0;
                double longitude = 0 ;
                double z = 0;

                latitude = lat[y][x];;
                longitude = lon[y][x];
                z = e[y][x];

                if( std::isnan(z))
                {
                    CHM_THROW_EXCEPTION(forcing_error,
                                        "Elevation for x,y=" + std::to_string(x) + ","+ std::to_string(y)+
                                            " is NaN. Regardless of the timestep the model is started from, it looks for timestep = 0 to define the elevations. Ensure it is defined then.");
                }

                size_t index = x + y * _nc->get_xsize();
                std::string station_name = std::to_string(index); // these don't really have names

                //need to convert the input lat/long into the coordinate system our mesh is in
                if (!_is_geographic)
                {
                    //CRS created with the “EPSG:4326” or “WGS84” strings use the latitude first, longitude second axis order.
                    if (!coordTrans->Transform(1, &latitude, &longitude))
                    {
                        BOOST_THROW_EXCEPTION(forcing_error() << errstr_info(
                            "Station=" + station_name + ": unable to convert coordinates to mesh format."));
                    }
                }

                double elevation = z;

                std::shared_ptr<station> s = std::make_shared<station>(station_name,
                    longitude, latitude, elevation, _variables);

                s->_nc_x = x;
                s->_nc_y = y;

                //index this linear array as if it were 2D to make the lazy load in the main run() loop easier.
                //it will allow us to pull out the station for a specific x,y more easily.
                _stations.at(index) = s;

                _dD_tree.insert( boost::make_tuple(Kernel::Point_2(s->x(),s->y()),s) );
            }
        }

    } catch(netCDF::exceptions::NcException& e)
    {
        OGRCoordinateTransformation::DestroyCT(coordTrans);
        BOOST_THROW_EXCEPTION(forcing_error() << errstr_info(e.what()));
    }

    OGRCoordinateTransformation::DestroyCT(coordTrans);
    _dt = _nc->get_dt();

    _current_ts = _start_time;
}

void metdata::load_from_ascii(std::vector<ascii_metdata> stations, int utc_offset)
{
    if(_mesh_proj4 == "")
    {
        CHM_THROW_EXCEPTION(config_error,"Metdata has not been initialized with the mesh and has no CRS.");
    }

    // a set of the ids we've loaded, ensure there are no duplicated IDs as there is some assumption we are not loading the same thing twic
    std::set<std::string> loaded_ids;

    for(auto& itr: stations)
    {
        if( (itr.latitude > 90 || itr.latitude < -90) ||
            (itr.longitude > 180 || itr.longitude < -180) )
        {
            BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("Station " + itr.id + " coordinate is invalid."));
        }


        if (!_is_geographic)
        {
            // spatial reference conversions to ensure the virtual station coordinates are the same as the meshes'
            OGRSpatialReference insrs, outsrs;
            insrs.SetWellKnownGeogCS("EPSG:4326");
            bool err = outsrs.importFromProj4(_mesh_proj4.c_str());
            if(err)
            {
                BOOST_THROW_EXCEPTION(forcing_error() << errstr_info( "Failure importing mesh proj4 string" ));
            }

            OGRCoordinateTransformation* coordTrans =  OGRCreateCoordinateTransformation(&insrs, &outsrs);
            //CRS created with the “EPSG:4326” use the latitude first, longitude second axis order.
            if (!coordTrans->Transform(1, &itr.latitude, &itr.longitude))
            {
                BOOST_THROW_EXCEPTION(forcing_error() << errstr_info(
                    "Station=" + itr.id + ": unable to convert coordinates to mesh format."));
            }
            // ensure we clean this up correctly
            OGRCoordinateTransformation::DestroyCT(coordTrans);
        }

        std::shared_ptr<station> s = std::make_shared<station>();
        s->ID(itr.id);
        s->x(itr.longitude);
        s->y(itr.latitude);
        s->z(itr.elevation);

        if(loaded_ids.find(s->ID()) == loaded_ids.end())
            loaded_ids.insert(s->ID());
        else
            CHM_THROW_EXCEPTION(forcing_error, "Stations with duplicated ID (" + s->ID() + ") inserted.");

        // load the ascii data into the timeseries object
        _ascii_stations.insert( std::make_pair(s->ID(), std::make_unique<ascii_data>()));
        _ascii_stations[s->ID()]->_obs.open(itr.path);

        // computes dt
        if(_ascii_stations[s->ID()]->_obs.get_date_timeseries().size() == 1)
        {
            CHM_THROW_EXCEPTION(model_init_error,"Unable to determine model timestep from only 1 input timestep.");
        }

        _ascii_stations[s->ID()]->_itr =  _ascii_stations[s->ID()]->_obs.begin();
        _ascii_stations[s->ID()]->id = s->ID();

        _variables = _ascii_stations[s->ID()]->_obs.list_variables();

        std::set<std::string> provides;

        for (auto filt : itr.filters)
        {
            filt->init();
            _ascii_stations[s->ID()]->filters.push_back(filt);         //copy the ptr for the filters, we own these now

            auto p = filt->provides();
            provides.insert(p.begin(),p.end());
        }

        _variables.insert(provides.begin(),provides.end());
        // init the datastore with timeseries variables + anything from the filters
        s->init(_variables);
        _stations.push_back(s);

        _dD_tree.insert( boost::make_tuple(Kernel::Point_2(s->x(),s->y()),s) );
    }

    // compute the dt for all stations and ensure they match
    std::vector<boost::posix_time::time_duration> dts;
    for(auto& itr: _ascii_stations)
    {
        auto t0 = itr.second->_obs.get_date_timeseries().at(0);
        auto t1 = itr.second->_obs.get_date_timeseries().at(1);
        auto dt = (t1 - t0);
        dts.push_back(dt);
    }

    if(dts.empty() || !std::all_of(dts.begin(),dts.end(),
                     [&](const auto& r) {return r==dts.front();}) )
    {
        BOOST_THROW_EXCEPTION(model_init_error() << errstr_info("Input forcing file timesteps are not consistent"));
    }

    _dt = dts.front();

    std::tie(_start_time,_end_time) = find_unified_start_end();
    subset(_start_time,_end_time); //subset assumes we have a valid dt

    _current_ts = _start_time;
    _n_timesteps = _ascii_stations.begin()->second->_obs.get_date_timeseries().size(); //grab the first timeseries, they are all the same period now
    _nstations = _ascii_stations.size();

}

boost::posix_time::ptime metdata::start_time()
{
    return _start_time;
}
boost::posix_time::ptime metdata::end_time()
{
    return _end_time;
}

boost::posix_time::time_duration metdata::dt()
{
    return _dt;
}
size_t metdata::dt_seconds()
{
    return dt().total_seconds();
}

size_t metdata::n_timestep()
{
    return _n_timesteps;
}
boost::posix_time::ptime metdata::current_time()
{
    return _current_ts;
}

std::string metdata::current_time_str()
{
    return boost::posix_time::to_iso_string(current_time());
}
std::string metdata::start_time_str()
{
    return boost::posix_time::to_iso_string(start_time());
}
std::string metdata::end_time_str()
{
    return boost::posix_time::to_iso_string(end_time());
}

std::set<std::string> metdata::list_variables()
{
    return _variables;
}

void metdata::check_ts_consistency()
{
    //ensure all the stations have the same start and end times
    // per-timestep agreeent happens during runtime.

    for (auto& itr : _ascii_stations)
    {
        if (itr.second->_obs.get_date_timeseries().at(0) != _start_time ||
            itr.second->_obs.get_date_timeseries().back() != _end_time)
        {
            BOOST_THROW_EXCEPTION(forcing_timestep_mismatch()
                                      <<
                                      errstr_info("Timestep mismatch at station: " + itr.second->id));
        }
    }
}
void metdata::subset(boost::posix_time::ptime start, boost::posix_time::ptime end)
{
    if( _dt.total_seconds() == 0)
    {
        CHM_THROW_EXCEPTION(forcing_error,"dt = 0");
    }
    // the netcdf files are simple and don't need this subsetting
    if(!_use_netcdf)
    {
        for(auto& itr : _ascii_stations)
        {
            itr.second->_obs.subset(start, end);
            itr.second->_itr = itr.second->_obs.begin();
        }
    }

    _start_time = start;
    _end_time = end;
    _current_ts = _start_time;
    _n_timesteps = ( (_end_time+_dt) - _start_time).total_seconds() / _dt.total_seconds(); // need to add +dt so that we are inclusive of the last timestep
}
std::pair<boost::posix_time::ptime,boost::posix_time::ptime> metdata::start_end_time()
{
    return std::make_pair(_start_time,_end_time);
}
std::pair<boost::posix_time::ptime,boost::posix_time::ptime> metdata::find_unified_start_end()
{
    if(!_use_netcdf)
    {

        _start_time = _ascii_stations.begin()->second->_obs.get_date_timeseries().at(0);
        _end_time = _ascii_stations.begin()->second->_obs.get_date_timeseries().back();

        // find the latest start time and the earliest end time
        for (auto& itr: _ascii_stations)
        {
            auto tmp = itr.second->_obs.get_date_timeseries().front();

            if (tmp > _start_time)
                _start_time = tmp;

            tmp = itr.second->_obs.get_date_timeseries().back();
            if (tmp < _end_time)
                _end_time = tmp;
        }
    }
    else
    {
        _start_time = _nc->get_start();
        _end_time = _nc->get_end();
    }
    return std::make_pair(_start_time,_end_time);
}

size_t metdata::nstations()
{
    return _nstations;
}

void metdata::write_stations_to_ptv(const std::string& path)
{
    assert(_nstations > 0);
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkStringArray> labels = vtkSmartPointer<vtkStringArray>::New();
    labels->SetName("Station name");
    labels->SetNumberOfValues(_nstations);
    for(size_t i = 0; i < _nstations;i++)
    {
        auto s = _stations.at(i);
        points->InsertNextPoint(s->x(), s->y(), s->z());
        labels->SetValue(i, s->ID() );
    }
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    polydata->GetPointData()->AddArray(labels);
    // Write the file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(path.c_str());

#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(polydata);
#else
    writer->SetInputData(polydata);
#endif

    writer->Write();
}

std::shared_ptr<station> metdata::at(size_t idx)
{
    return _stations.at(idx);
}

bool metdata::next()
{
    bool has_next = false;

    // allows for doing first timestep loading without incrementing the timestep
    if(!is_first_timestep)
        _current_ts = _current_ts + _dt;

    if(_use_netcdf)
    {
        has_next = next_nc();
    }
    else
    {
        has_next = next_ascii();
    }

    is_first_timestep = false;
    return has_next;
}

bool metdata::next_ascii()
{

    for(size_t i = 0; i < nstations();i++)
    {
        auto s = _stations.at(i);
        auto& proxy = _ascii_stations[s->ID()];

        //the very first timestep needs to handle loading the data without incrementing the internal iterators
        if(!is_first_timestep)
        {
            ++proxy->_itr;
            if (proxy->_itr == proxy->_obs.end())
                return false;
        }


        if(proxy->_itr->get_posix() != _current_ts)
        {
            CHM_THROW_EXCEPTION(forcing_error,
                "Mismatch between model timestep and ascii file timestep. Current model = " +
                boost::posix_time::to_simple_string(_current_ts) + ", ascii was:"+
                boost::posix_time::to_simple_string(proxy->_itr->get_posix()) +" @station id="+s->ID());
        }

        // don't use the stations variable map as it'll contain anything inserted by a filter which won't exist in the ascii file
        for (auto &v: proxy->_obs.list_variables() )
        {
            (*s)[v] = proxy->_itr->get(v);

        }

        //get the list of filters to run for this station
        auto filters = proxy->filters;

        for (auto& filt : filters)
        {
            filt->process(s);
        }

        s->set_posix(_current_ts);

    }

    return true;
}
bool metdata::next_nc()
{
    if(_current_ts > _end_time) //_current_ts is already ++ from the next() call
    {
        return false; // we've run out of data, we done
    }


    //The call to netCDF isn't thread safe. It is protected by a critical section but it's costly, and not running this
    // in parallel is about 2x faster  #pragma omp parallel for
    for(size_t i = 0; i < nstations();i++)
    {
        auto s = _stations.at(i);
        s->set_posix(_current_ts);

        // don't use the stations variable map as it'll contain anything inserted by a filter which won't exist in the nc file
        for (auto &v: _nc->get_variable_names() )
        {
            double d = _nc->get_var(v, _current_ts, s->_nc_x, s->_nc_y);
            (*s)[v] = d;

        }

        // run all the filters for this station
        for (auto& f : _netcdf_filters)
        {
            f.second->process(s);
        }
    }

    return true;

}

std::vector< std::shared_ptr<station> > metdata::get_stations_in_radius(double x, double y, double radius )
{
    // define exact circular range query  (fuzziness=0)
    Kernel::Point_2 center(x, y);
    Fuzzy_circle exact_range(center, radius);

    std::vector<boost::tuple<Kernel::Point_2, std::shared_ptr<station> > > result;
    _dD_tree.search(std::back_inserter(result), exact_range);

    std::vector< std::shared_ptr<station> > stations;
    stations.reserve(result.size());
    for (auto& itr : result)
    {
        stations.push_back( boost::get<1>(itr));
    }
    return stations;

}

std::vector< std::shared_ptr<station> > metdata::nearest_station(double x, double y,unsigned int N)
{
    Kernel::Point_2 query(x,y);
    Neighbor_search search(_dD_tree, query, N);

    std::vector< std::shared_ptr<station> > stations;
    for (auto itr : search)
    {
        stations.push_back( boost::get<1>(itr.first));
    }
    return stations;

}

void metdata::prune_stations(std::unordered_set<std::string>& station_ids)
{
    _stations.erase(
        std::remove_if(std::begin(_stations), std::end(_stations),
        [&](auto const& it)
        {
          return (station_ids.find(it->ID()) != std::end(station_ids));
        }),
        std::end(_stations));

    _nstations = _stations.size();
}

std::vector< std::shared_ptr<station>>& metdata::stations()
{
    return _stations;
}