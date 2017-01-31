#include "snow_slide.hpp"

snow_slide::snow_slide(config_file cfg)
        : module_base(parallel::domain)
{
    depends("snowdepthavg");
    depends("swe");

    provides("delta_swe");
    provides("delta_snowdepthavg");
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
        auto face = domain->face(i); // Get face
        // Make copy of snowdepthavg and swe to modify within snow_slide (not saved)
        auto data = face->get_module_data<snow_slide::data>(ID); // Get data
        data->snowdepthavg_copy = face->face_data("snowdepthavg"); // Store copy of snowdepth for snow_slide use
        data->swe_copy = face->face_data("swe");
        // Initalize snow transport to zero
        data->delta_snowdepthavg = 0.0;
        data->delta_swe = 0.0;
        sorted_z.push_back( std::make_pair( face->center().z() + face->face_data("snowdepthavg"), face) );

    }

    // Sort faces by elevation + snowdepth
    std::sort(sorted_z.begin(), sorted_z.end(), [](std::pair<double, mesh_elem> &a, std::pair<double, mesh_elem> &b) {
        return b.first < a.first;
    });


    // Loop through each face, from highest to lowest triangle surface
    for (size_t i = 0; i < sorted_z.size(); i++) {
    	auto face = sorted_z[i].second; // Get pointer to face
        double cen_area = face->get_area(); // Area of center triangle
        auto data = face->get_module_data<snow_slide::data>(ID); // Get stored data for face

        // Loop through each face, from highest to lowest triangle
        double maxDepth = data->maxDepth;
        double snowdepthavg = data->snowdepthavg_copy;
        double swe = data->swe_copy;

        // Check if face snowdepth have exceeded maxDepth
        if (snowdepthavg > 0.3 ) { //maxDepth) {
            LOG_DEBUG << "avalanche! " << snowdepthavg << " " << maxDepth;

            double z_s = face->center().z() + snowdepthavg; // Current face elevation + snowdepth
            std::vector<double> w = {0,0,0}; // Weights for each face neighbor to route snow to
            double w_dem = 0; // Denomenator for weights (sum of all elev diffs)

            // Calc weights for routing snow
            for(int i = 0; i < 3; ++i) {
                auto n = face->neighbor(i); // Pointer to neighbor face

                // If not null (edge case) check is less than lowFace
                if(n != nullptr) {
                    auto n_data = n->get_module_data<snow_slide::data>(ID); // pointer to face's data
                    // Calc weighting based on height diff
                    // (std::max insures that if one neighbor is higher, its weight will be zero)
                    w[i] = std::max(0.0, z_s -  (n->center().z()+n_data->snowdepthavg_copy) );
                    w_dem += w[i]; // Store weight denominator
                }
            }
            // Divide by sum height differences to create weights that sum to unity
            if(w_dem != 0) { // prevent divide by zero
            std::transform(w.begin(), w.end(), w.begin(),
               [w_dem](double cw) { return cw/w_dem; });
            }

            // Route snow to each neighbor based on weights
            // TODO: Include vegetation height in routing calcs
            for(int i = 0; i < 3; ++i) {
                auto n = face->neighbor(i);
                if(n != nullptr) {
                    auto n_data = n->get_module_data<snow_slide::data>(ID); // pointer to face's data

                    // Update copy of snowdepth and swe of neighbor
                    n_data->snowdepthavg_copy += snowdepthavg*w[i]; // (m)
                    n_data->swe_copy += swe*w[i]; // (m)

                    // Update mass transport of neighbor
                    n_data->delta_snowdepthavg += snowdepthavg*w[i]* cen_area; // Fraction of snowdepth (m) * center triangle area (m^2) = volune of snow depth (m^3)
                    n_data->delta_swe += swe*w[i] * cen_area; // (m) * (m^2) = (m^3) of swe
                }
            }
            // Remove snow from initial face (TODO: should we only remove above maxDepth??)
            data->snowdepthavg_copy = 0; // data refers to current/center cell
            data->swe_copy = 0;
            // Update mass transport
            data->delta_snowdepthavg -= snowdepthavg * cen_area;
            data->delta_swe -= swe * cen_area;
        }

        // Save state variables at end of time step
        face->set_face_data("maxDepth",maxDepth);
        face->set_face_data("delta_snowdepthavg",data->delta_snowdepthavg);
        face->set_face_data("delta_swe",data->delta_swe);

    } // End of each face

    // Mass balance check (Note: because of the way delta values are stored, they should NOT sum to zero over the domain)
    //double mass_bal = 0;
    //for (size_t i = 0; i < domain->size_faces(); i++) {
    //    auto face = domain->face(i);
    //    mass_bal += face->face_data("delta_swe");
    //}
    //LOG_DEBUG << "Mass balance = " << mass_bal; 
}

void snow_slide::init(mesh domain)
{
    // Initilaze for each triangle
    for(size_t i=0;i<domain->size_faces();i++)
    {
        auto face = domain->face(i);

        auto d = face->make_module_data<snow_slide::data>(ID);

	// Parametrize the Minimum snow holding depth
        double slopeDeg = std::max(10.0,face->slope()*180/M_PI);  // radians to degres, limit to >10 degrees to avoid inf
        d->maxDepth = 3178.4 * pow(slopeDeg,-1.998); // (m??) Estimate min depth that avanlanch occurs
    }
}
