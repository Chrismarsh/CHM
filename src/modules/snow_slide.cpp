#include "snow_slide.hpp"

snow_slide::snow_slide(config_file cfg)
        : module_base(parallel::domain)
{
    depends("snowdepthavg");
    depends("swe");

    provides("delta_avalanche_mass");
    provides("delta_avalanche_snowdepth");
    provides("maxDepth");

}

snow_slide::~snow_slide()
{

}

void snow_slide::checkpoint(mesh domain,  netcdf& chkpt)
{

    chkpt.create_variable1D("snow_slide:delta_avalanche_snowdepth", domain->size_faces());
    chkpt.create_variable1D("snow_slide:delta_avalanche_mass", domain->size_faces());

    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        chkpt.put_var1D("snow_slide:delta_avalanche_snowdepth",i,face->get_module_data<data>(ID)->delta_avalanche_snowdepth);
        chkpt.put_var1D("snow_slide:delta_avalanche_mass",i,face->get_module_data<data>(ID)->delta_avalanche_mass);
    }
}

void snow_slide::load_checkpoint(mesh domain,  netcdf& chkpt)
{
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        face->get_module_data<data>(ID)->delta_avalanche_snowdepth = chkpt.get_var1D("snow_slide:delta_avalanche_snowdepth",i);
        face->get_module_data<data>(ID)->delta_avalanche_mass = chkpt.get_var1D("snow_slide:delta_avalanche_mass",i);
    }
}

void snow_slide::run(mesh domain)
{

    // Make a vector of pairs (elevation + snowdepth, pointer to face)
    tbb::concurrent_vector< std::pair<double, mesh_elem> > sorted_z(domain->size_faces());

#pragma omp parallel for
    for(size_t i = 0; i  < domain->size_faces(); i++)
    {
        auto face = domain->face(i); // Get face
        // Make copy of snowdepthavg and swe to modify within snow_slide (not saved)
        auto data = face->get_module_data<snow_slide::data>(ID); // Get data
        data->snowdepthavg_copy = face->face_data("snowdepthavg"); // Store copy of snowdepth for snow_slide use
        data->swe_copy = face->face_data("swe")/1000; // mm to m
        // Initalize snow transport to zero
        data->delta_avalanche_snowdepth = 0.0;
        data->delta_avalanche_mass = 0.0; // m
        sorted_z.at(i) = std::make_pair( face->center().z() + face->face_data("snowdepthavg"), face) ;

    }

    // Sort faces by elevation + snowdepth
//    std::sort(sorted_z.begin(), sorted_z.end(), [](const std::pair<double, mesh_elem> &a,const std::pair<double, mesh_elem> &b) {
//        return b.first < a.first;
//    });
    tbb::parallel_sort(sorted_z.begin(), sorted_z.end(), [](const std::pair<double, mesh_elem> &a,const std::pair<double, mesh_elem> &b) {
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
            bool edge_flag = false; // Flag for determiing if current cell is an edge (handel routing differently)

            // Calc weights for routing snow
            // Possible Cases:
            //      1) edge cell, then edge_flag is true, and snow is dumpped off mesh
            //      2) non-edge cell, w_dem is greater than 0 -> there is atleast one lower neighbor, route so to it/them
            //      3) non-edge cell, w_dem = 0, "sink" case. Don't route any snow.
            for (int i = 0; i < 3; ++i) {
                auto n = face->neighbor(i); // Pointer to neighbor face

                // Check if not-null (null indicates edge cell)
                if (n != nullptr) {
                    auto n_data = n->get_module_data<snow_slide::data>(ID); // pointer to face's data
                    // Calc weighting based on height diff
                    // (std::max insures that if one neighbor is higher, its weight will be zero)
                    w[i] = std::max(0.0, z_s - (n->center().z() + n_data->snowdepthavg_copy));
                    w_dem += w[i]; // Store weight denominator
                } else { // It is an edge cell, set flag
                    edge_flag = true;
                }
            }

            // Case 1) Edge cell
            if(edge_flag) {
                // Special case: dump snow out of domain (loosing mass) by just removing from current edge cell.
                // Don't route and exit loop.

                // Remove snow from initial face
                data->snowdepthavg_copy = maxDepth;
                data->swe_copy = swe * maxDepth / snowdepthavg;
                // Update mass transport (m^3)
                data->delta_avalanche_snowdepth -= del_depth * cen_area;
                data->delta_avalanche_mass -= del_swe * cen_area;

                // Save state variables
                face->set_face_data("delta_avalanche_snowdepth", data->delta_avalanche_snowdepth);
                face->set_face_data("delta_avalanche_mass", data->delta_avalanche_mass);
            }

            // Case 2) Non-Edge cell, but w_dem=0, "sink" cell. Don't route snow.
            if(w_dem==0) {
                continue; // Restart to next loop.
            }

            // Must be Case 3), Divide by sum height differences to create weights that sum to unity
            if (w_dem != 0) {
                std::transform(w.begin(), w.end(), w.begin(),
                               [w_dem](double cw) { return cw / w_dem; });
            }

            // Case 3), Non-Edge cell, w_dem>0, route snow to down slope neighbor(s).
            double out_mass = 0; // Mass balance check
            // Route snow to each neighbor based on weights
            for (int j = 0; j < 3; ++j) {
                auto n = face->neighbor(j);
                if (n != nullptr) {
                    double n_area = n->get_area(); // Area of neighbor triangle
                    auto   n_data = n->get_module_data<snow_slide::data>(ID); // pointer to face's data

                    // // Update neighbor snowdepth and swe (copies only for internal snowSlide use)
                    // Here we must make an assumption of the pack density (because we do not have access to
                    // layer information (if exists (i.e. snowpack is running), or it doesn't (i.e. snobal is running))
                    // Therefore, we assume uniform density.
                    // The (cen_area/n_area) term converts depth change from orig cell to volume, then back to a depth term
                    // using the neighbor's area.
                    n_data->snowdepthavg_copy += del_depth * (cen_area/n_area) * w[j]; // (m)
                    n_data->swe_copy += del_swe * (cen_area/n_area) * w[j]; // (m)

                    // Update mass transport to neighbor
                    n_data->delta_avalanche_snowdepth += del_depth * cen_area * w[j]; // Fraction of snowdepth (m) *
                    // center triangle area (m^2) = volune of snow depth (m^3)
                    n_data->delta_avalanche_mass +=  del_swe * cen_area * w[j]; // (m) * (m^2) = (m^3) of swe
                    out_mass += del_swe * cen_area * w[j];
                }
            }
            // Remove snow from initial face
            data->snowdepthavg_copy = maxDepth; // data refers to current/center cell
            data->swe_copy = swe * maxDepth / snowdepthavg; // Uses ratio of depth change to calc new swe
            // This relys on the assumption of uniform density.

            // Update mass transport (m^3)
            data->delta_avalanche_snowdepth -= del_depth * cen_area;
            data->delta_avalanche_mass -= del_swe * cen_area;

            // Check mass transport balances for current avalanche cell
            if (std::abs(orig_mass-out_mass)>0.0001) {
                LOG_DEBUG << "Moved mass total is " << out_mass;
                LOG_DEBUG << "diff = " << orig_mass-out_mass;
                LOG_DEBUG << "Mass balance of avalanche times step was not conserved.";
            }

        }

        // Save state variables at end of time step
        face->set_face_data("delta_avalanche_snowdepth", data->delta_avalanche_snowdepth);
        face->set_face_data("delta_avalanche_mass", data->delta_avalanche_mass);

    } // End of each face


}

void snow_slide::init(mesh domain)
{
    // Get Parameters that control maxDepth function
    double avalache_mult = cfg.get("avalache_mult",3178.4);
    double avalache_pow  = cfg.get("avalache_pow",-1.998);

    // Initialize for each triangle
    for(size_t i=0;i<domain->size_faces();i++)
    {
        auto face = domain->face(i);

        auto d = face->make_module_data<snow_slide::data>(ID);

        double Z_CanTop;
        if(face->has_vegetation())
        {

            Z_CanTop = face->veg_attribute("CanopyHeight");
        } else {
            Z_CanTop = 0.0;
        }

	    // Parametrize the Minimum snow holding depth
        double slopeDeg = std::max(10.0,face->slope()*180/M_PI);  // radians to degres, limit to >10 degrees to avoid inf
        d->maxDepth = std::max(avalache_mult * pow(slopeDeg,avalache_pow),Z_CanTop); // (m) Estimate min depth that avalanche occurs
        // Max of either veg height or derived max holding snow depth.
        face->set_face_data("maxDepth", d->maxDepth);
    }
}
