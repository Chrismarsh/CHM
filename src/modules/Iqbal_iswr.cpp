//
// Created by chris on 06/11/15.
//

#include "Iqbal_iswr.h"

Iqbal_iswr::Iqbal_iswr(config_file cfg)
        :module_base(parallel::data)
{

    depends("t");
    depends("rh");
    depends("cloud_frac");
    depends("solar_el");

    provides("iswr");
    provides("iswr_direct");
    provides("iswr_diffuse");
    provides("atm_trans");
}

Iqbal_iswr::~Iqbal_iswr()
{

}

void Iqbal_iswr::run(mesh_elem &face)
{
    double pressure = mio::Atmosphere::stdAirPressure(face->get_z());//101325.0;
    double altitude = face->get_z();
    double sun_elevation = face->face_data("solar_el");
    sun_elevation = sun_elevation < 0? 0. : sun_elevation;

    if (sun_elevation < 3)
    {
        face->set_face_data("iswr",0);
        face->set_face_data("iswr_direct",0);
        face->set_face_data("iswr_diffuse",0);
        return;
    }

    double ta = face->face_data("t")+273.15;
    double rh = face->face_data("rh")/100.0;
    double R_toa = 1375;
    double R_direct=0;
    double R_diffuse=0;
    double ground_albedo = 0.1;


    //these pow cost us a lot here, but replacing them by fastPow() has a large impact on accuracy (because of the exp())
    const double olt = 0.32;   //ozone layer thickness (cm) U.S.standard = 0.34 cm
    const double w0 = 0.9;     //fraction of energy scattered to total attenuation by aerosols (Bird and Hulstrom(1981))
    const double fc = 0.84;    //fraction of forward scattering to total scattering (Bird and Hulstrom(1981))
    const double alpha = 1.3;  //wavelength exponent (Iqbal(1983) p.118). Good average value: 1.3+/-0.5. Related to the size distribution of the particules
    const double beta = 0.03;  //amount of particules index (Iqbal(1983) p.118). Between 0 & .5 and above.
    const double zenith = 90. - sun_elevation; //this is the TRUE zenith because the elevation is the TRUE elevation
    const double cos_zenith = cos(zenith*mio::Cst::to_rad); //this uses true zenith angle

    // relative optical air mass Young (1994), see http://en.wikipedia.org/wiki/Airmass
    //const double mr = 1. / (cos_zenith + 0.50572 * pow( 96.07995-zenith , -1.6364 )); //pbl: this should use apparent zenith angle, and we only get true zenith angle here...
    // relative optical air mass, Young, A. T. 1994. Air mass and refraction. Applied Optics. 33:1108–1110.
    const double mr = ( 1.002432*cos_zenith*cos_zenith + 0.148386*cos_zenith + 0.0096467) /
                      ( cos_zenith*cos_zenith*cos_zenith + 0.149864*cos_zenith*cos_zenith
                        + 0.0102963*cos_zenith +0.000303978);

    // actual air mass: because mr is applicable for standard pressure
    // it is modified for other pressures (in Iqbal (1983), p.100)
    // pressure in Pa
    const double ma = mr * (pressure/101325.);

    // the equations for all the transmittances of the individual atmospheric constituents
    // are from Bird and Hulstrom (1980, 1981) and can be found summarized in Iqbal (1983)
    // on the quoted pages

    // broadband transmittance by Rayleigh scattering (Iqbal (1983), p.189)
    const double taur = exp( -0.0903 * pow(ma,0.84) * (1. + ma - pow(ma,1.01)) );

    // broadband transmittance by ozone (Iqbal (1983), p.189)
    const double u3 = olt * mr; // ozone relative optical path length
    const double alpha_oz = 0.1611 * u3 * pow(1. + 139.48 * u3, -0.3035) -
                            0.002715 * u3 / ( 1. + 0.044  * u3 + 0.0003 * u3 * u3); //ozone absorbance
    const double tauoz = 1. - alpha_oz;

    // broadband transmittance by uniformly mixed gases (Iqbal (1983), p.189)
    const double taug = exp( -0.0127 * pow(ma, 0.26) );

    // saturation vapor pressure in Pa
    //const double Ps = exp( 26.23 - 5416./ta ); //as used for the parametrization
    const double Ps = mio::Atmosphere::waterSaturationPressure(ta);

    // Leckner (1978) (in Iqbal (1983), p.94), reduced precipitable water
    const double w = 0.493 * rh * Ps / ta;

    // pressure corrected relative optical path length of precipitable water (Iqbal (1983), p.176)
    // pressure and temperature correction not necessary since it is included in its numerical constant
    const double u1 = w * mr;

    // broadband transmittance by water vapor (in Iqbal (1983), p.189)
    const double tauw = 1. - 2.4959 * u1  / (pow(1.0 + 79.034 * u1, 0.6828) + 6.385 * u1);

    // broadband total transmittance by aerosols (in Iqbal (1983), pp.189-190)
    // using Angstroem's turbidity formula Angstroem (1929, 1930) for the aerosol thickness
    // in Iqbal (1983), pp.117-119
    // aerosol optical depth at wavelengths 0.38 and 0.5 micrometer
    const double ka1 = beta * pow(0.38, -alpha);
    const double ka2 = beta * pow(0.5, -alpha);

    // broadband aerosol optical depth:
    const double ka  = 0.2758 * ka1 + 0.35 * ka2;

    // total aerosol transmittance function for the two wavelengths 0.38 and 0.5 micrometer:
    const double taua = exp( -pow(ka, 0.873) * (1. + ka - pow(ka, 0.7088)) * pow(ma, 0.9108) );

    // broadband transmittance by aerosols due to absorption only (Iqbal (1983) p. 190)
    const double tauaa = 1. - (1. - w0) * (1. - ma + pow(ma, 1.06)) * (1. - taua);

    // broadband transmittance function due to aerosols scattering only
    // Iqbal (1983) p. 146 (Bird and Hulstrom (1981))
    const double tauas = taua / tauaa;

    // direct normal solar irradiance in range 0.3 to 3.0 micrometer (Iqbal (1983) ,p.189)
    // 0.9751 is for this wavelength range.
    // Bintanja (1996) (see Corripio (2002)) introduced a correction beta_z for increased
    // transmittance with altitude that is linear up to 3000 m and than fairly constant up to 5000 - 6000 m
    const double beta_z = (altitude<3000.)? 2.2*1.e-5*altitude : 2.2*1.e-5*3000.;

    //Now calculating the radiation
    //Top of atmosphere radiation (it will always be positive, because we check for sun elevation before)
    const double tau_commons = tauoz * taug * tauw * taua;

    // Diffuse radiation from the sky
    const double factor = 0.79 * R_toa * tau_commons / (1. - ma + pow( ma,1.02 ));  //avoid recomputing pow() twice
    // Rayleigh-scattered diffuse radiation after the first pass through atmosphere (Iqbal (1983), p.190)
    const double Idr = factor * 0.5 * (1. - taur );

    // aerosol scattered diffuse radiation after the first pass through atmosphere (Iqbal (1983), p.190)
    const double Ida = factor * fc  * (1. - tauas);

    // cloudless sky albedo Bird and Hulstrom (1980, 1981) (in Iqbal (1983) p. 190)
    //in Iqbal, it is recomputed with ma=1.66*pressure/101325.; and alb_sky=0.0685+0.17*(1.-taua_p)*w0;
    const double alb_sky = 0.0685 + (1. - fc) * (1. - tauas);


    //Now, we compute the direct and diffuse radiation components
    //Direct radiation. All transmitances, including Rayleigh scattering (Iqbal (1983), p.189)
    R_direct = 0.9751*( taur * tau_commons + beta_z ) * R_toa ;

    // multiple reflected diffuse radiation between surface and sky (Iqbal (1983), p.154)
    const double Idm = (Idr + Ida + R_direct) * ground_albedo * alb_sky / (1. - ground_albedo * alb_sky);
    R_diffuse = (Idr + Ida + Idm)*cos_zenith; //Iqbal always "project" diffuse radiation on the horizontal


    double elevation_threshold = 2.0 * mio::Cst::to_rad;

    if( sun_elevation < elevation_threshold ) {
        //if the Sun is too low on the horizon, we put all the radiation as diffuse
        //the splitting calculation that might take place later on will reflect this
        //instead point radiation, it becomes the radiation of a horizontal sky above the domain
        R_diffuse += R_direct*cos_zenith; //HACK
       // R_diffuse =0.0;
        R_direct = 0.;
    }


    double cf = face->face_data("cloud_frac");
    double dir = R_direct  * (0.6 + 0.2*cos_zenith) * (1.0-cf);


    face->set_face_data("iswr_direct",dir);
    face->set_face_data("iswr_diffuse",R_diffuse);
    face->set_face_data("iswr",dir+R_diffuse);
    face->set_face_data("atm_trans",(dir+R_diffuse)/1375.);

}
