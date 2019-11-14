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

                //this ctor will init an empty timeseries for us
                std::shared_ptr<station> s = std::make_shared<station>(station_name, longitude, latitude,
                                                                       elevation);
                s->raw_timeseries()->init(variables, date_vec);
                s->reset_itrs();

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

void metdata::load_from_ascii(std::vector<ascii_metadata> stations, int utc_offset)
{
    for(auto& itr: stations)
    {
        if( (itr.latitude > 90 || itr.latitude < -90) ||
            (itr.longitude > 180 || itr.longitude < -180) )
        {
            BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("Station " + itr.name + " coordinate is invalid."));
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

        std::shared_ptr<station> s = std::make_shared<station>();

        if (!_mesh->is_geographic())
        {
            if (!coordTrans->Transform(1, &itr.longitude, &itr.latitude))
            {
                BOOST_THROW_EXCEPTION(forcing_error() << errstr_info(
                    "Station=" + itr.name + ": unable to convert coordinates to mesh format."));
            }
        }

        s->x(itr.longitude);
        s->y(itr.latitude);
        s->z(itr.elevation);

        s->open(itr.path);

    }

    std::tie(_start_time,_end_time) = start_end_times();

    // computes dt
    if(_stations.at(0)->date_timeseries().size() == 1)
    {
        BOOST_THROW_EXCEPTION(model_init_error() << errstr_info("Unable to determine model timestep from only 1 input timestep."));
    }

    auto t0 = _stations.at(0)->date_timeseries().at(0);
    auto t1 = _stations.at(0)->date_timeseries().at(1);
    _dt = (t1 - t0);

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

    for (size_t i = 0; i < _stations.size(); i++)
    {
        if (_stations.at(i)->date_timeseries().at(0) != _start_time ||
            _stations.at(i)->date_timeseries().back() != _end_time)
        {
            BOOST_THROW_EXCEPTION(forcing_timestep_mismatch()
                                      <<
                                      errstr_info("Timestep mismatch at station: " + _stations.at(i)->ID()));
        }
    }
}
void metdata::subset(boost::posix_time::ptime start, boost::posix_time::ptime end)
{

    // the netcdf files are simple and don't need this subsetting
    if(!_use_netcdf)
    {
#pragma omp for
        for(size_t i = 0; i < _stations.size(); ++i)
        {
            _stations.at(i)->raw_timeseries()->subset(start, end);
            _stations.at(i)->reset_itrs();
        }
    }

    _start_time = start;
    _end_time = end;
    _current_ts = _start_time;
}
std::pair<boost::posix_time::ptime,boost::posix_time::ptime> metdata::start_end_times()
{
// find the latest start time and the earliest end time
    for (size_t i = 0; i < _stations.size(); ++i)
    {
        auto tmp = _stations.at(i)->date_timeseries().at(0);

        if (tmp > _start_time)
            _start_time = tmp;

        tmp = _stations.at(i)->date_timeseries().back();
        if (tmp < _end_time)
            _end_time = tmp;
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

void metdata::next()
{
    _current_ts = _current_ts + _dt;

    if(_use_netcdf)
    {
        next_nc();
    }
}

void metdata::next_nc()
{
    // don't use the stations variable map as it'll contain anything inserted by a filter which won't exist in the nc file
    for (auto &itr: _nc->get_variable_names() )
    {
//        auto data = _nc->get_var(itr, t);

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

//                double d = data[y][x];
                double d =  _nc->get_var(itr, _current_ts, x, y);
                (*s)[itr] = d;
//                s->now().set(itr, d);

                for (auto &f : _netcdf_filters)
                {
                    f.second->process(s);
                }
            }
        }
    }



}