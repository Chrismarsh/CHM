#include "global.hpp"

void global::update()
{
    solar_el_az();
}
int global::year()
{
    return _current_date.date().year();
}
int global::day()
{
    return _current_date.date().day();
}
int global::month()
{
    return _current_date.date().month();
}
int global::hour()
{
    return boost::posix_time::to_tm(_current_date).tm_hour;
}
int global::min()
{
    return boost::posix_time::to_tm(_current_date).tm_min;
}
int global::sec()
{
    return boost::posix_time::to_tm(_current_date).tm_sec;
}
boost::posix_time::ptime global::posix_time()
{
    return _current_date;
}

double global::solar_el()
{
    return _solar_el;
}

double global::solar_az()
{
    return _solar_az;  
}

void global::solar_el_az()
{
    //UTC offset. Don't know how to use datetime's UTC converter yet....
    boost::posix_time::time_duration UTC_offset = boost::posix_time::hours(_utc_offset);
                
    std::tm tm = boost::posix_time::to_tm(_current_date+UTC_offset);
    double year =  tm.tm_year + 1900; //convert from epoch
    double month =  tm.tm_mon + 1;//conert jan == 0
    double day =   tm.tm_mday; //starts at 1, ok
    double hour = tm.tm_hour; // 0 = midnight, ok
    double min = tm.tm_min; // 0, ok
    double sec = tm.tm_sec; // [0,60] in c++11, ok http://en.cppreference.com/w/cpp/chrono/c/tm
    double Alt = 0; //TODO: fix this?
    
    if (month <= 2.0)
    {
        year = year -1.0;
        month = month +12.0;
    }
    
    double jd = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - \
        floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5 + \
        (hour + min/60 + sec/3600)/24;

    double d = jd-2451543.5;
    // Keplerian Elements for the Sun (geocentric)
    double w = 282.9404+4.70935*10*pow(10,-5)*d; //    (longitude of perihelion degrees)
    double a = 1.000000;          //    (mean distance, a.u.)
    double e = 0.016709- 1.151*pow(10,-9)*d;  //    (eccentricity)
    double M = fmod(356.0470+0.9856002585*d,360.0); //  (mean anomaly degrees)
    double L = w + M;                     //(Sun's mean longitude degrees)
    double oblecl = 23.4393-3.563e-7*d;  //(Sun's obliquity of the ecliptic)

    double pi = 3.1415926;
    //auxiliary angle
    double E = M+(180/pi)*e*sin(M*(pi/180))*(1+e*cos(M*(pi/180))); 
    
    //rectangular coordinates in the plane of the ecliptic (x axis toward
    //perhilion)
    double x = cos(E*(pi/180))-e;
    double y = sin(E*(pi/180))*sqrt(1-e*e);
    
    //find the distance and true anomaly
    double r = sqrt(x*x + y*y);
    double v = atan2(y,x)*(180/pi);
    
    double Lon = _lon;
    double Lat = _lat;
    //find the longitude of the sun
    double lon = v + w;

        //compute the ecliptic rectangular coordinates
    double xeclip = r*cos(lon*(pi/180));
    double yeclip = r*sin(lon*(pi/180));
    double zeclip = 0.0;

    //rotate these coordinates to equitorial rectangular coordinates
    double xequat = xeclip;
    double yequat = yeclip*cos(oblecl*(pi/180))+zeclip*sin(oblecl*(pi/180));
    double zequat = yeclip*sin(23.4406*(pi/180))+zeclip*cos(oblecl*(pi/180));

        //convert equatorial rectangular coordinates to RA and Decl:
     r = sqrt(xequat*xequat + yequat*yequat + zequat*zequat)-(Alt/149598000); //roll up the altitude correction
    double RA = atan2(yequat,xequat)*(180/pi);
    double delta = asin(zequat/r)*(180/pi);


    //Following the RA DEC to Az Alt conversion sequence explained here:
    //http://www.stargazing.net/kepler/altaz.html

    //Find the J2000 value
    double J2000 = jd - 2451545.0;

     double   UTH = hour;   //Calculate local siderial time
    double GMST0=fmod(L+180,360)/15;
    double SIDTIME = GMST0 + UTH + Lon/15;

        //Replace RA with hour angle HA
    double HA = (SIDTIME*15 - RA);

    //convert to rectangular coordinate system
    x = cos(HA*(pi/180))*cos(delta*(pi/180));
    y = sin(HA*(pi/180))*cos(delta*(pi/180));
    double z = sin(delta*(pi/180));

    //rotate this along an axis going east-west.
    double xhor = x*cos((90-Lat)*(pi/180))-z*sin((90-Lat)*(pi/180));
    double yhor = y;
    double zhor = x*sin((90-Lat)*(pi/180))+z*cos((90-Lat)*(pi/180));

    //Find the h and AZ 
    double Az = atan2(yhor,xhor)*(180/pi) + 180;
    double El = asin(zhor)*(180/pi);
    
    _solar_az = Az;
    _solar_el = El;
}


std::string global::get_variable(std::string variable)
{
    return _variables(variable);
    
}