#include "solar.hpp"

solar::solar(config_file cfg)
        :module_base(parallel::data)
{
    provides("solar_el");
    provides("solar_az");
}
solar::~solar()
{

}
void solar::run(mesh_elem &face, boost::shared_ptr <global> global_param)
{
    double Lon = 0;
    double Lat = 0;
    if(global_param->is_geographic())
    {
        Lon = face->center().x();
        Lat = face->center().y();

    }
    else{
        auto data = face->get_module_data<solar::data>(ID);
        Lon = data->lng;
        Lat = data->lat;
    }

    //approximate UTC offset
    int _utc_offset=6;

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
    double Alt = 0.; //TODO: fix this?

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

    face->set_face_data("solar_az",Az);
    face->set_face_data("solar_el",El);

}
void solar::init(mesh domain, boost::shared_ptr<global> global_param)
{
    if(!domain->is_geographic())
    {

        OGRSpatialReference monUtm;
        monUtm.SetWellKnownGeogCS("WGS84");

        int UTM = domain->UTM_zone();
        monUtm.SetUTM( abs(UTM) , UTM > 0); //UTM > 0 => northern hemisphere. SetUTM takes 2nd arg = true for Northern hemisphere

        OGRSpatialReference monGeo;
        monGeo.SetWellKnownGeogCS("WGS84");

        // we are UTM and need to convert internally to lat long to calc the solar position
#pragma omp parallel for
        for (size_t i = 0; i < domain->size_faces(); i++)
        {
            auto face = domain->face(i);

            OGRCoordinateTransformation* coordTrans = OGRCreateCoordinateTransformation(&monUtm, &monGeo);

            double x = face->center().x();
            double y = face->center().y();
            int reprojected = coordTrans->Transform(1, &x, &y);
            delete coordTrans;

            auto d = face->make_module_data<solar::data>(ID);
            d->lat = y;
            d->lng = x;
        }
    }



}
