#include "snow_slide.hpp"

snow_slide::snow_slide(config_file cfg)
        : module_base(parallel::domain)
{
    depends("snowdepthavg");
    depends("swe");

    provides("delta_swe");
    provides("delta_snowdepthavg");
    provides("delta_snowdepthavg_depthTESTING");
    provides("maxDepth");

}

snow_slide::~snow_slide()
{

}

void snow_slide::run(mesh domain)
{
    // TODO: swe is stored in mm through CHM, we convert it to m for consisent use in snowslide
    // TODO: Need to make it m throughout CHM.

    // Make a vector of pairs (elevation + snowdepth, pointer to face)
    std::vector< std::pair<double, mesh_elem> > sorted_z;
    for(size_t i = 0; i  <domain->size_faces(); i++)
    {
        auto face = domain->face(i); // Get face
        // Make copy of snowdepthavg and swe to modify within snow_slide (not saved)
        auto data = face->get_module_data<snow_slide::data>(ID); // Get data
        data->snowdepthavg_copy = face->face_data("snowdepthavg"); // Store copy of snowdepth for snow_slide use
        data->swe_copy = face->face_data("swe")/1000; // mm to m
        // Initalize snow transport to zero
        data->delta_snowdepthavg = 0.0;
        data->delta_swe = 0.0; // m
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

        // Get current triangle snow info at beginning of time step
        double maxDepth = data->maxDepth;
        double snowdepthavg = data->snowdepthavg_copy; // m
        double swe = data->swe_copy; // m

        // Check if face snowdepth have exceeded maxDepth
        if (snowdepthavg > maxDepth) {
            //LOG_DEBUG << "avalanche! " << snowdepthavg << " " << maxDepth;
            double del_depth = snowdepthavg - maxDepth; // Amount to be removed (positive) [m]
            double del_swe   = swe * (1 - maxDepth / snowdepthavg); // Amount of swe to be removed (positive) [m]
            double orig_mass = del_swe * cen_area;

            double z_s = face->center().z() + snowdepthavg; // Current face elevation + snowdepth
            std::vector<double> w = {0, 0, 0}; // Weights for each face neighbor to route snow to
            double w_dem = 0; // Denomenator for weights (sum of all elev diffs)

            // Calc weights for routing snow
            for (int i = 0; i < 3; ++i) {
                auto n = face->neighbor(i); // Pointer to neighbor face

                // If not null (edge case) check is less than lowFace
                if (n != nullptr) {
                    auto n_data = n->get_module_data<snow_slide::data>(ID); // pointer to face's data
                    // Calc weighting based on height diff
                    // (std::max insures that if one neighbor is higher, its weight will be zero)
                    w[i] = std::max(0.0, z_s - (n->center().z() + n_data->snowdepthavg_copy));
                    w_dem += w[i]; // Store weight denominator
                }
            }
            // Check for case where all weights are zero (edge cell with higher neighbors)
            if(w_dem==0) {
                LOG_DEBUG << "above maxDepth but no cells to dump too. Cell " << i;
                continue; // Don't do anything to this cell or neighbors
            }
            // Divide by sum height differences to create weights that sum to unity
            if (w_dem != 0) { // prevent divide by zero
                std::transform(w.begin(), w.end(), w.begin(),
                               [w_dem](double cw) { return cw / w_dem; });
            }

            double out_mass = 0; // DEBUG check
            // Route snow to each neighbor based on weights
            for (int j = 0; j < 3; ++j) {
                auto n = face->neighbor(j);
                if (n != nullptr) {
                    double n_area = n->get_area(); // Area of neighbor triangle
                    auto   n_data = n->get_module_data<snow_slide::data>(ID); // pointer to face's data
                    // // Update neighbor snowdepth and swe (copies only for internal snowSlide use)
                    // Update copy of snowdepth and swe of neighbor
                    n_data->snowdepthavg_copy += del_depth * (cen_area/n_area) * w[j]; // (m)
                    // Here we must make an assumption of the pack density (because we do not have access to
                    // layer information (if exists (i.e. snowpack is running), or it doesn't (i.e. snobal is running))
                    // Assume uniform density.
                    //
                    // Remember (swe, snowdepthavg, and maxDepth refer to the tri being avalanched, while
                    // n_data->swe_copy refer to the neighbor tri)
                    n_data->swe_copy += del_swe * (cen_area/n_area) * w[j]; // (m)

                    // Update mass transport to neighbor
                    n_data->delta_snowdepthavg += del_depth * cen_area * w[j]; // Fraction of snowdepth (m) *
                    // center triangle area (m^2) = volune of snow depth (m^3)
                    n_data->delta_swe +=  del_swe * cen_area * w[j]; // (m) * (m^2) = (m^3) of swe
                    out_mass += del_swe * cen_area * w[j];
                }
            }
            // Remove snow from initial face
            data->snowdepthavg_copy = maxDepth; // data refers to current/center cell
            data->swe_copy = swe * maxDepth / snowdepthavg;

            // Update mass transport (m^3)
            data->delta_snowdepthavg -= del_depth * cen_area;
            data->delta_swe -= del_swe * cen_area;

            // Check mass transport balances for current avalanche cell
            if (std::abs(orig_mass-out_mass)>0.0001) {
                for (int j = 0; j < 3; ++j) {
                    LOG_DEBUG << "weight is " << w[j];
                }
                LOG_DEBUG << "Moved mass total is " << out_mass;
                LOG_DEBUG << "diff = " << orig_mass-out_mass;
                LOG_DEBUG << " something bad happend";
            }

        }

        // Save state variables at end of time step
        face->set_face_data("maxDepth", maxDepth);
        face->set_face_data("delta_snowdepthavg", data->delta_snowdepthavg);
        face->set_face_data("delta_snowdepthavg_depthTESTING", data->delta_snowdepthavg / cen_area); // m
        face->set_face_data("delta_swe", data->delta_swe);

    } // End of each face


}

void snow_slide::init(mesh domain)
{
    // Initilaze for each triangle
    for(size_t i=0;i<domain->size_faces();i++)
    {
        auto face = domain->face(i);

        auto d = face->make_module_data<snow_slide::data>(ID);

        double Z_CanTop;
        if(face->has_parameter("landcover") )
        {
            int LC   = face->get_parameter("landcover");
            Z_CanTop = global_param->parameters.get<double>("landcover." + std::to_string(LC) + ".CanopyHeight");
        } else {
            Z_CanTop = 0.0;
        }

	    // Parametrize the Minimum snow holding depth
        double slopeDeg = std::max(10.0,face->slope()*180/M_PI);  // radians to degres, limit to >10 degrees to avoid inf
        d->maxDepth = std::max(3178.4 * pow(slopeDeg,-1.998),Z_CanTop); // (m) Estimate min depth that avanlanch occurs
        // Max of either veg height or derived max holding snow depth.
    }
}
