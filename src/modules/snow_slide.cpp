#include "snow_slide.hpp"

snow_slide::snow_slide(config_file cfg)
        : module_base(parallel::domain)
{
    depends("snowdepthavg");
    depends("swe");

    provides("swe");
    provides("snowdepthavg");
    provides("maxDepth");

}

snow_slide::~snow_slide()
{

}

void snow_slide::run(mesh domain)
{

    // Make a vector of pairs (elevation, pointer to face)
    std::vector< std::pair<double, mesh_elem> > sorted_z;
    for(int i = 0; i  <domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        sorted_z.push_back( std::make_pair( face->center().z(), face) );

    }

    // Sort faces by elevation (TODO: add snowpack height)
    std::sort(sorted_z.begin(), sorted_z.end(), [](std::pair<double, mesh_elem> &a, std::pair<double, mesh_elem> &b) {
        return b.first < a.first;
    });


    // Loop through each face, from highest to lowest triangle
    for (size_t i = 0; i < sorted_z.size(); i++)
    {
    	auto face = sorted_z[i].second;

	auto data = face->get_module_data<snow_slide::data>(ID);

	// Loop through each face, from highest to lowest triangle
	double maxDepth = data->maxDepth;
	double snowdepthavg = face->face_data("snowdepthavg");
	double swe = face->face_data("swe");

	// Check if face snowdepth have exceeded maxDepth
        if (snowdepthavg > maxDepth) {
		LOG_DEBUG << "avalanche! " << snowdepthavg << " " << maxDepth;
                // Get elevations of neighbors
		double min_z = 99999.0; // A large number
                int min_i = 10; // something larger than 2
		for(int i = 0; i < 3; ++i) {
		    auto n = face->neighbor(i);
    
		    // If not null (edge case) check is less than lowFace
    		    if(n != nullptr) {
			if(n->center().z() < min_z) { 
			    min_z = n->center().z();
			    min_i = i;
			}
		    }
	        }
		auto n0 = face->neighbor(min_i);
		//LOG_DEBUG << "neighbor elev " << n0->center().z() << ". min_z " << min_z;
		
		// Route snow to neighbor
		n0->set_face_data("snowdepthavg",n0->face_data("snowdepthavg")+snowdepthavg);
                n0->set_face_data("swe",n0->face_data("swe")+swe);

		// Remove snow from initial face
		face->set_face_data("snowdepthavg",0.0);
                face->set_face_data("swe",0.0);

        } 

	// Do nothing if threshold maxDepth not exceded 
        // testing states
        face->set_face_data("maxDepth",maxDepth);

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

	// Parametrize the Minimum snow holding depth
    double slopeDeg = std::max(10.0,face->slope()*180/M_PI);  // radians to degres, limit to >10 degrees to avoid inf
	
d->maxDepth = 3178.4 * pow(slopeDeg,-1.998); // (m??) Estimate min depth that avanlanch occurs

    }
}
