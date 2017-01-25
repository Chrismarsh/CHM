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

    // Make a vector of pairs (elevation + snowdepth, pointer to face)
    std::vector< std::pair<double, mesh_elem> > sorted_z;
    for(size_t i = 0; i  <domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        sorted_z.push_back( std::make_pair( face->center().z() + face->face_data("snowdepthavg"), face) );

    }

    // Sort faces by elevation + snowdepth
    std::sort(sorted_z.begin(), sorted_z.end(), [](std::pair<double, mesh_elem> &a, std::pair<double, mesh_elem> &b) {
        return b.first < a.first;
    });

    // Loop through each face, from highest to lowest triangle surface
    for (size_t i = 0; i < sorted_z.size(); i++)
    {
    	auto face = sorted_z[i].second; // Get pointer to face

	auto data = face->get_module_data<snow_slide::data>(ID); // Get stored data for face

	// Loop through each face, from highest to lowest triangle
	double maxDepth = data->maxDepth;
	double snowdepthavg = face->face_data("snowdepthavg");
	double swe = face->face_data("swe");

	// Check if face snowdepth have exceeded maxDepth
        if (snowdepthavg >  maxDepth) {
		LOG_DEBUG << "avalanche! " << snowdepthavg << " " << maxDepth;
                // Calc weights for routing snow
		
		double z_s = face->center().z() + snowdepthavg; // Current face elevation + snowdepth
                //LOG_DEBUG << "current height " << z_s;
		std::vector<double> w = {0,0,0}; // Weights for each face neighbor to route snow to
                double w_dem = 0; // Denomenator for weights (sum of all elev diffs)
                
                for(int i = 0; i < 3; ++i) {
		    auto n = face->neighbor(i);
    
		    // If not null (edge case) check is less than lowFace
    		    if(n != nullptr) {
                        // Calc weighting based on height diff 
                        // (std::max insures that if one neighbor is higher, its weight will be zero)
			w[i] = std::max(0.0, z_s -  (n->center().z()+n->face_data("snowdepthavg")) );
                        //LOG_DEBUG << "weight of " << i << " is " << w[i];
                        w_dem += w[i]; // Store weight denominator 
		    }
		}
                // Divide by sum height differences to create weights that sum to unityi
                //LOG_DEBUG << "weight is " << w_dem;
                if(w_dem != 0) { // prevent divide by zero
                std::transform(w.begin(), w.end(), w.begin(),
                   [w_dem](double cw) { return cw/w_dem; });
                }
                //LOG_DEBUG << "weights are " << w[0] << " " << w[1] << " "<< w[2] << " ";	
		
		// Route snow to each neighbor based on weights
                // TODO: Include vegetation heigh in roughting calcs
		for(int i = 0; i < 3; ++i) {
                    auto n = face->neighbor(i);
                    if(n != nullptr) {
                        n->set_face_data("snowdepthavg",n->face_data("snowdepthavg")+snowdepthavg*w[i]);
                        n->set_face_data("swe",n->face_data("swe")+swe*w[i]);
                    }
		}
                // Remove snow from initial face (TODO: should we only remove above maxDepth??)
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
