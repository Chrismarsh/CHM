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

metdata::metdata()
{
    _nc = nullptr;
    _use_netcdf = false;
}

metdata::metdata(boost::shared_ptr< triangulation > mesh) :
 metdata()
{
    _mesh = mesh;
}

void metdata::load_from_netcdf(const std::string& path,std::map<std::string, boost::shared_ptr<filter_base> > filters)
{
    if(!_mesh)
        BOOST_THROW_EXCEPTION(forcing_error() << errstr_info( "metdata class not initialized with a mesh ptr" ));

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
    insrs.SetWellKnownGeogCS("WGS84");
    bool err = outsrs.importFromProj4(_mesh->proj4().c_str());
    if(err)
    {
        BOOST_THROW_EXCEPTION(forcing_error() << errstr_info( "Failure importing mesh proj4 string" ));
    }
    OGRCoordinateTransformation* coordTrans =  nullptr;
    if(!_mesh->is_geographic())
        coordTrans = OGRCreateCoordinateTransformation(&insrs, &outsrs);

    try
    {
        _nc->open_GEM(path);

        auto variables = _nc->get_variable_names();

        variables.insert(_provides_from_nc_filters.begin(),_provides_from_nc_filters.end());

        auto date_vec = _nc->get_datevec();

        _start_time = date_vec.at(0);
        _end_time = date_vec.back();

        _nstations = _nc->get_xsize() * _nc->get_ysize();
        LOG_DEBUG << "Grid is (y)" << _nc->get_ysize() << " by (x)" << _nc->get_xsize();

        LOG_DEBUG << "Loading lat/long grid...";

        auto lat = _nc->get_lat();
        auto lon = _nc->get_lon();
        auto e = _nc->get_z();

        LOG_DEBUG << "Initializing datastructure";

// #pragma omp parallel for
// hangs, unclear why, critical sections around the json and gdal calls don't seem to help
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

                size_t index = x + y * _nc->get_xsize();
                std::string station_name = std::to_string(index); // these don't really have names

                //need to convert the input lat/long into the coordinate system our mesh is in
                if (!_mesh->is_geographic())
                {
                    if (!coordTrans->Transform(1, &longitude, &latitude))
                    {
                        BOOST_THROW_EXCEPTION(forcing_error() << errstr_info(
                            "Station=" + station_name + ": unable to convert coordinates to mesh format."));
                    }
                }

                double elevation = z;

                std::shared_ptr<station> s = std::make_shared<station>(station_name,
                    longitude, latitude, elevation, variables);

                //index this linear array as if it were 2D to make the lazy load in the main run() loop easier.
                //it will allow us to pull out the station for a specific x,y more easily.
                _stations.at(index) = s;
            }
        }

    } catch(netCDF::exceptions::NcException& e)
    {
        delete coordTrans;
        BOOST_THROW_EXCEPTION(forcing_error() << errstr_info(e.what()));
    }

    delete coordTrans;
    _dt = _nc->get_dt();

    _current_ts = _start_time;

}

void metdata::load_from_ascii(std::vector<ascii_metdata> stations, int utc_offset)
{
    for(auto& itr: stations)
    {
        if( (itr.latitude > 90 || itr.latitude < -90) ||
            (itr.longitude > 180 || itr.longitude < -180) )
        {
            BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("Station " + itr.id + " coordinate is invalid."));
        }

        // spatial reference conversions to ensure the virtual station coordinates are the same as the meshes'
        OGRSpatialReference insrs, outsrs;
        insrs.SetWellKnownGeogCS("WGS84");
        bool err = outsrs.importFromProj4(_mesh->proj4().c_str());
        if(err)
        {
            BOOST_THROW_EXCEPTION(forcing_error() << errstr_info( "Failure importing mesh proj4 string" ));
        }
        OGRCoordinateTransformation* coordTrans =  nullptr;
        if(!_mesh->is_geographic())
            coordTrans = OGRCreateCoordinateTransformation(&insrs, &outsrs);

        if (!_mesh->is_geographic())
        {
            if (!coordTrans->Transform(1, &itr.longitude, &itr.latitude))
            {
                BOOST_THROW_EXCEPTION(forcing_error() << errstr_info(
                    "Station=" + itr.id + ": unable to convert coordinates to mesh format."));
            }
        }

        std::shared_ptr<station> s = std::make_shared<station>();
        s->ID(itr.id);
        s->x(itr.longitude);
        s->y(itr.latitude);
        s->z(itr.elevation);

        // load the ascii data into the timeseries object
        _ascii_stations.insert( std::make_pair(s->ID(), std::make_unique<ascii_data>()));
        _ascii_stations[s->ID()]->_obs.open(itr.path);

        // computes dt
        if(_ascii_stations[s->ID()]->_obs.get_date_timeseries().size() == 1)
        {
            CHM_THROW_EXCEPTION(model_init_error,"Unable to determine model timestep from only 1 input timestep.");
        }

        _ascii_stations[s->ID()]->_itr =  _ascii_stations[s->ID()]->_obs.begin();

        std::set<std::string> provides;
        auto obs_variables = _ascii_stations[s->ID()]->_obs.list_variables();
        provides.insert(obs_variables.begin(), obs_variables.end());
        //copy the ptr for the filters, we own these now
        for (auto filt : itr.filters)
        {
            filt->init();
            _ascii_stations[s->ID()]->filters.push_back(filt);

            auto p = filt->provides();
            provides.insert(p.begin(),p.end());
        }

        // init the datastore with timeseries variables + anything from the filters
        s->init(provides);
    }


    std::tie(_start_time,_end_time) = start_end_times();

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
    _current_ts = _start_time;
}

boost::posix_time::ptime metdata::start_time()
{
    return _start_time;
}
boost::posix_time::ptime metdata::end_time()
{
    return _end_time;
}

size_t metdata::dt()
{

    return _dt.total_seconds();
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
                                      errstr_info("Timestep mismatch at station: " + itr.second.id));
        }
    }
}
void metdata::subset(boost::posix_time::ptime start, boost::posix_time::ptime end)
{

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
}
std::pair<boost::posix_time::ptime,boost::posix_time::ptime> metdata::start_end_times()
{
    if(!_use_netcdf)
    {
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
    _current_ts = _current_ts + _dt;

    if(_use_netcdf)
    {
        return next_nc();
    }
    else
    {
        return next_ascii();
    }

    return false;
}

bool metdata::next_ascii()
{

    for(size_t i = 0; i < nstations();i++)
    {
        auto s = _stations.at(i);
        auto& proxy = _ascii_stations[s->ID()];

        ++proxy->_itr;
        if (proxy->_itr == proxy->_obs.end())
            return false;

        if(proxy->_itr->get_posix() != _current_ts)
        {
            CHM_THROW_EXCEPTION(forcing_error,
                std::to_string("Mismatch between model timestep and ascii file timestep. Current model = ")+
                _current_ts + ", ascii was:"+proxy->_itr->get_posix()+" @station id="+s->ID());
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

    }

    return true;
}
bool metdata::next_nc()
{
    if(_current_ts > _end_time) //_current_ts is already ++ from the next() call
    {
        return false; // we've run out of data, we done
    }
    // don't use the stations variable map as it'll contain anything inserted by a filter which won't exist in the nc file
    for (auto &v: _nc->get_variable_names() )
    {
        // Need to decide how best to do this as it may be too memory expensive to load the entire thing into ram
        // auto data = _nc->get_var(itr, t);

#pragma omp parallel for
        for (size_t y = 0; y < _nc->get_ysize(); y++)
        {
            for (size_t x = 0; x < _nc->get_xsize(); x++)
            {
                size_t index = x + y * _nc->get_xsize();
                auto s = _stations.at(index);

                //sanity check that we are getting the right station for this xy pair
                if (s->ID() != std::to_string(index))
                {
                    CHM_THROW_EXCEPTION(forcing_error, "Station=" + s->ID() + ": wrong ID");
                }

                double d =  _nc->get_var(v, _current_ts, x, y);
                (*s)[v] = d;

                for (auto &f : _netcdf_filters)
                {
                    f.second->process(s);
                }
            }
        }
    }

    return true;


}