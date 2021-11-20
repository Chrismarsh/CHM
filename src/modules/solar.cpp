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

#include "solar.hpp"
REGISTER_MODULE_CPP(solar);

solar::solar(config_file cfg)
        : module_base("solar", parallel::data, cfg)
{
    provides("solar_el");
    provides("solar_az");

    provides_parameter("svf");
}
solar::~solar()
{

}
void solar::run(mesh_elem &face)
{
    double Lon = 0;
    double Lat = 0;
    if(global_param->is_geographic())
    {
        Lon = face->center().x();
        Lat = face->center().y();
    }
    else{
        auto& data = face->get_module_data<solar::data>(ID);
        Lon = data.lng;
        Lat = data.lat;
    }


    //Following the RA DEC to Az Alt conversion sequence explained here:
    //http://www.stargazing.net/kepler/altaz.html

    //UTC offset. Don't know how to use datetime's UTC converter yet....
    boost::posix_time::time_duration UTC_offset = boost::posix_time::hours(global_param->_utc_offset);
    std::tm tm = boost::posix_time::to_tm(global_param->posix_time()+UTC_offset);
    double year =  tm.tm_year + 1900.; //convert from epoch
    double month =  tm.tm_mon + 1.;//conert jan == 0
    double day =   tm.tm_mday; //starts at 1, ok
    double hour = tm.tm_hour; // 0 = midnight, ok
    double min = tm.tm_min; // 0, ok
    double sec = tm.tm_sec; // [0,60] in c++11, ok http://en.cppreference.com/w/cpp/chrono/c/tm
    double Alt = face->center().z();//0.; //TODO: fix this?

    if (month <= 2.0)
    {
        year = year -1.0;
        month = month +12.0;
    }

    double jd = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - \
        floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5 + \
        (hour + min/60. + sec/3600.)/24.;

    double d = jd-2451543.5;
    // Keplerian Elements for the Sun (geocentric)
    double w = 282.9404+4.70935*pow(10,-5)*d; //    (longitude of perihelion degrees)
    double a = 1.000000;          //    (mean distance, a.u.)
    double e = 0.016709- 1.151*pow(10.,-9.)*d;  //    (eccentricity)
    double M = fmod(356.0470+0.9856002585*d,360.0); //  (mean anomaly degrees)
    double L = w + M;                     //(Sun's mean longitude degrees)
    double oblecl = 23.4393-3.563e-7*d;  //(Sun's obliquity of the ecliptic)

    //auxiliary angle
    double E = M+(180./M_PI)*e*sin(M*(M_PI/180.))*(1+e*cos(M*(M_PI/180.)));

    //rectangular coordinates in the plane of the ecliptic (x axis toward
    //perhilion)
    double x = cos(E*(M_PI/180.))-e;
    double y = sin(E*(M_PI/180.))*sqrt(1.-e*e);

    //find the distance and true anomaly
    double r = sqrt(x*x + y*y);
    double v = atan2(y,x)*(180./M_PI);


    //find the longitude of the sun
    double lon = v + w;

    //compute the ecliptic rectangular coordinates
    double xeclip = r*cos(lon*(M_PI/180.));
    double yeclip = r*sin(lon*(M_PI/180.));
    double zeclip = 0.0;

    //rotate these coordinates to equitorial rectangular coordinates
    double xequat = xeclip;
    double yequat = yeclip*cos(oblecl*(M_PI/180.))+zeclip*sin(oblecl*(M_PI/180.));
    double zequat = yeclip*sin(23.4406*(M_PI/180.))+zeclip*cos(oblecl*(M_PI/180.));

    //convert equatorial rectangular coordinates to RA and Decl:
    r = sqrt(xequat*xequat + yequat*yequat + zequat*zequat)-(Alt/149598000.0); //roll up the altitude correction
    double RA = atan2(yequat,xequat)*(180./M_PI);
    double delta = asin(zequat/r)*(180./M_PI);

    //Find the J2000 value
    double J2000 = jd - 2451545.0;

    double UTH = hour+min/60.0+sec/3600.0;   //Calculate local siderial time
    double GMST0=fmod(L+180.,360.)/15.;
    double SIDTIME = GMST0 + UTH + Lon/15.;

    //Replace RA with hour angle HA
    double HA = (SIDTIME*15. - RA);

    //convert to rectangular coordinate system
    x = cos(HA*(M_PI/180.))*cos(delta*(M_PI/180.));
    y = sin(HA*(M_PI/180.))*cos(delta*(M_PI/180.));
    double z = sin(delta*(M_PI/180.));

    //rotate this along an axis going east-west.
    double xhor = x*cos((90.-Lat)*(M_PI/180.))-z*sin((90.-Lat)*(M_PI/180.));
    double yhor = y;
    double zhor = x*sin((90.-Lat)*(M_PI/180.))+z*cos((90.-Lat)*(M_PI/180.));

    //Find the h and AZ
    double Az = atan2(yhor,xhor)*(180./M_PI) + 180.;
    double El = asin(zhor)*(180./M_PI);

    (*face)["solar_az"_s]=Az;
    (*face)["solar_el"_s]=El;

}
void solar::init(mesh& domain)
{

    //number of steps along the search vector to check for a higher point
    int steps = cfg.get("svf.steps",10);
    //max distance to search
    double max_distance = cfg.get("svf.max_distance",1000.0);

    //size of the step to take
    double size_of_step = max_distance / steps;

    //number of azimuthal sections
    int N = cfg.get("svf.nsectors", 12);


    double azimuthal_width = 360./(double)N; // in degrees


    bool svf_compute = cfg.get("svf.compute",true);

    OGRSpatialReference monUtm;
    OGRSpatialReference monGeo;
    OGRCoordinateTransformation* coordTrans = nullptr;



    // we are UTM and need to convert internally to lat long to calc the solar position
    if(!domain->is_geographic())
    {
        monUtm.importFromProj4(domain->proj4().c_str());
        monUtm.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

        monGeo.SetWellKnownGeogCS("WGS84");
        monGeo.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

        coordTrans = OGRCreateCoordinateTransformation(&monUtm, &monGeo);
    }

    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {

	       auto face = domain->face(i);

	       // we are UTM and need to convert internally to lat long to calc the solar position
	       if(!domain->is_geographic())
	       {
		   double x = face->center().x();
		   double y = face->center().y();
		   int reprojected = coordTrans->Transform(1, &x, &y);

		   auto& d = face->make_module_data<solar::data>(ID);
		   d.lat = y;
		   d.lng = x;
	       }

	       double svf = 0.0;

	       if(svf_compute)
	       {
               Point_3 me = face->center();
               auto cosSlope = cos(face->slope());
               auto sinSlope = sin(face->slope());

               //for each search azimuthal sector
               for (int k = 0; k < N; k++)
               {
                   double phi = 0.;
                   // search along each azimuth in j step increments to find horizon angle
                   for (int j = 1; j <= steps; ++j)
                   {
                   double distance = j * size_of_step;

                   auto f = domain->find_closest_face(math::gis::point_from_bearing(me, k * azimuthal_width, distance));

                   double z_diff = (f->center().z() - me.z());
                   if (z_diff > 0)
                   {
                       double dist = math::gis::distance(f->center(), me);
                       phi = std::max(atan(z_diff / dist), phi);
                   }
                   }

                   auto cosPhi = cos(phi);
                   auto sinPhi = sin(phi);
                   auto azi_in_rad = (k * azimuthal_width * M_PI / 180.);

                   svf += cosSlope * cosPhi * cosPhi +
                 sinSlope * cos(azi_in_rad - face->aspect()) * (M_PI_2 - phi - sinPhi * cosPhi);

               }

               svf /= (double)N;
	       } else{
		        svf = 1.;
	       }
	       face->parameter("svf"_s) = std::max(0.0, svf);

    }


    delete coordTrans;

}
