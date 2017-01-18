#include "snow_slide.hpp"

snow_slide::snow_slide(config_file cfg)
        : module_base(parallel::domain)
{
    depends("snowdepthavg");
    depends("swe");

    provides("swe");
    provides("snowdepthavg");
    provides("minDepth");

}

snow_slide::~snow_slide()
{

}

void snow_slide::run(mesh domain)
{

    // Sort faces by elevation + snowdepth
    


    // Loop through each face, from highest to lowest triangle
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
            auto face = domain->face(i);

	    auto data = face->get_module_data<snow_slide::data>(ID);

	    // Loop through each face, from highest to lowest triangle
	    double minDepth = data->minDepth;
	    double snowdepthavg = face->face_data("snowdepthavg");
	    double swe = face->face_data("swe");

	    // testing
            if (snowdepthavg > minDepth) {
		LOG_DEBUG << "avalanche! " << snowdepthavg << " " << minDepth;
		snowdepthavg = 0;
		swe = 0;	
	    } 

	    // Save model states
	    face->set_face_data("snowdepthavg",snowdepthavg); 
	    face->set_face_data("swe",swe);
	    // testing states
	    face->set_face_data("minDepth",minDepth);

    }
}

void snow_slide::init(mesh domain)
{
    // Snow Slide paramters
    //double S_m = 25; // Minimum slope angle  (Degrees)
   
    for(size_t i=0;i<domain->size_faces();i++)
    {
        auto face = domain->face(i);

        auto d = face->make_module_data<snow_slide::data>(ID);

	// Parameterize the Minimum snow holding depth
	double slopeDeg = face->slope()*180/M_PI; // radians to degres
	d->minDepth = 3178.4 * pow(slopeDeg,-1.998); // (m??) Estimate min depth that avanlanch occurs 
    }
}
