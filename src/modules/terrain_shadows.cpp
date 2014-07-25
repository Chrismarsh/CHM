
#include "solar.hpp"

terrain_shadow::terrain_shadow(std::string ID)
{
    _provides->push_back("is_shadow");
    this->ID = ID;
    _parallel_type = parallel::domain;
    LOG_DEBUG << "Successfully instantiated module " << this->ID;

}

void terrain_shadow::run(mesh domain, boost::shared_ptr<global> global_param)
{


    double A = global_param->solar_az() * M_PI / 180.0;
    double E = global_param->solar_el() * M_PI / 180.0;

    if (E < 5)
    {
        for (triangulation::Finite_faces_iterator fit = domain->finite_faces_begin(); fit != domain->finite_faces_end(); ++fit)
        {
            //interpolate the station data to the current element
            triangulation::Face_handle face = fit;
            face->set_face_data("terrain_shadow", 0); //unshadowed
        }
        return;
    }
    //euler rotation matrix K
    arma::mat K;
    // eqns(6) & (7) in Montero
    double z0 = M_PI - A;
    double q0 = M_PI / 2.0 - E;

    K << cos(z0) << sin(z0) << 0 << arma::endr
      << -cos(q0) * sin(z0) << cos(q0) * cos(z0) << sin(q0) << arma::endr
      << sin(q0) * sin(z0) << -cos(z0) * sin(q0) << cos(q0) << arma::endr;

    Delaunay prj(*domain);
    for (triangulation::Finite_faces_iterator fit = prj->finite_faces_begin(); fit != prj->finite_faces_end(); ++fit)
    {
        //interpolate the station data to the current element
        triangulation::Face_handle face = fit;
        
        
	arma::vec coord(3);
        coord(0) = face->vertex(0)->point().x();
        coord(1) = face->vertex(0)->point().y();
        coord(2) = face->vertex(0)->point().z();
        
        coord = K*coord;
        
        face->
    }
    


}

Solar::~Solar()
{


}



