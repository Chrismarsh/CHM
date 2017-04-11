
#include "fast_shadow.hpp"

fast_shadow::fast_shadow(config_file cfg)
        : module_base(parallel::data)
{
    depends("solar_az");
    depends("solar_el");
    provides("shadow");

    //number of steps along the search vector to check for a higher point
    steps = cfg.get("steps",10);
    //max distance to search
    max_distance = cfg.get("max_distance",1000.0);

    //size of the step to take
    size_of_step = max_distance / steps;
}

fast_shadow::~fast_shadow()
{


}

void fast_shadow::run(mesh_elem& face)
{
        //currentl degrees
        double solar_az = face->face_data("solar_az") ;
        double solar_el = face->face_data("solar_el") ;

        face->set_face_data("shadow", 0);
        //bail early
        if (solar_el < 5)
            return;

        Point_3 me = face->center();
        auto cosSlope = cos(face->slope());
        auto sinSlope = sin(face->slope());

        double phi = 0.;
        // search along each azimuth in j step increments to find horizon angle
        for (int j = 1; j <= steps; ++j)
        {
            double distance = j * size_of_step;

            auto f = face->find_closest_face(solar_az, distance);

            double z_diff = f->center().z() - me.z() ;
            if (z_diff > 0)
            {
                double dist = math::gis::distance(f->center(), me);
                phi = std::max(atan(z_diff / dist), phi);
            }
        }

        if (phi > (solar_el *M_PI / 180.) )
        {
            face->set_face_data("shadow", 1);
        }

}

