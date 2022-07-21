//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//

#include "snow_slide.hpp"
REGISTER_MODULE_CPP(snow_slide);

snow_slide::snow_slide(config_file cfg)
        : module_base("snow_slide", parallel::domain, cfg)
{
    depends("snowdepthavg",SpatialType::neighbor);
    depends("swe",SpatialType::neighbor);

    use_vertical_snow = cfg.get("use_vertical_snow",true);

    provides("delta_avalanche_mass");
    provides("delta_avalanche_snowdepth");

    provides("delta_avalanche_mass_sum");
    provides("delta_avalanche_snowdepth_sum");

    provides("maxDepth");

    provides("ghost_ss_snowdepthavg_vert_copy");
    provides("ghost_ss_snowdepthavg_to_xfer");
    provides("ghost_ss_swe_to_xfer");

    provides("ghost_ss_delta_avalanche_snowdepth");
    provides("ghost_ss_delta_avalanche_swe");

    provides("ghost_ss_sum_swe_xfer");

    provides("test");

}

snow_slide::~snow_slide()
{

}

void snow_slide::checkpoint(mesh& domain,  netcdf& chkpt)
{

    chkpt.create_variable1D("snow_slide:delta_avalanche_snowdepth", domain->size_faces());
    chkpt.create_variable1D("snow_slide:delta_avalanche_mass", domain->size_faces());

    chkpt.create_variable1D("snow_slide:delta_avalanche_snowdepth_sum", domain->size_faces());
    chkpt.create_variable1D("snow_slide:delta_avalanche_mass_sum", domain->size_faces());

    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        chkpt.put_var1D("snow_slide:delta_avalanche_snowdepth",i,face->get_module_data<data>(ID).delta_avalanche_snowdepth);
        chkpt.put_var1D("snow_slide:delta_avalanche_mass",i,face->get_module_data<data>(ID).delta_avalanche_mass);

        chkpt.put_var1D("snow_slide:delta_avalanche_snowdepth_sum",i,  (*face)["delta_avalanche_snowdepth_sum"_s]);
        chkpt.put_var1D("snow_slide:delta_avalanche_mass_sum",i, (*face)["delta_avalanche_mass_sum"_s]);

    }
}

void snow_slide::load_checkpoint(mesh& domain,  netcdf& chkpt)
{
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        face->get_module_data<data>(ID).delta_avalanche_snowdepth = chkpt.get_var1D("snow_slide:delta_avalanche_snowdepth",i);
        face->get_module_data<data>(ID).delta_avalanche_mass = chkpt.get_var1D("snow_slide:delta_avalanche_mass",i);

        (*face)["delta_avalanche_snowdepth_sum"_s] = chkpt.get_var1D("snow_slide:delta_avalanche_snowdepth_sum",i);
        (*face)["delta_avalanche_mass_sum"_s] = chkpt.get_var1D("snow_slide:delta_avalanche_mass_sum",i);

    }
}

void snow_slide::run(mesh& domain)
{

    int done = 1; // int because we run a global reduce on it to determine a min global state

    int iterations = 0; // number of iterations we've run
    do
    {

        int this_iter_moved_snow = false;
        // Make a vector of pairs (elevation + snowdepth, pointer to face)
        tbb::concurrent_vector<std::pair<double, mesh_elem>> sorted_z(domain->size_faces());

        // only init these on the first iteration
        if(iterations ==0)
        {
#pragma omp parallel for
            for (size_t i = 0; i < domain->size_faces(); i++)
            {
                auto face = domain->face(i); // Get face
                auto& data = face->get_module_data<snow_slide::data>(ID); // Get data

                    // Make copy of snowdepthavg and swe to modify within snow_slide (not saved)
                    data.snowdepthavg_copy = (*face)["snowdepthavg"_s]; // Store copy of snowdepth for snow_slide use
                    data.snowdepthavg_vert_copy = (*face)["snowdepthavg_vert"_s]; // Vertical snow depth
                    data.swe_copy = (*face)["swe"_s] / 1000.0;                    // mm to m
                    data.slope = face->slope();                                   // slope in rad

                    // Initalize snow transport to zero
                    data.delta_avalanche_snowdepth = 0.0;
                    data.delta_avalanche_mass = 0.0; // m

                    for (int j = 0; j < 3; ++j)
                    {
                        auto n = face->neighbor(j);
                        if (n!=nullptr && n->is_ghost)
                        {
                             #pragma omp critical
                            {
                                (*n)["ghost_ss_sum_swe_xfer"] = 0;
                            }
                        }
                    }

            }
        }


        //#ifdef USE_MPI
//        // update our ghosts from other ranks values
////        // first iteration this will be 0
////        domain->ghost_neighbors_communicate_variable("ghost_ss_snowdepthavg_to_xfer"_s);
//        domain->ghost_neighbors_communicate_variable("ghost_ss_snowdepthavg_vert_copy"_s);
//#endif

#pragma omp parallel for
        for (size_t i = 0; i < domain->size_faces(); i++)
        {
            auto face = domain->face(i); // Get face
            auto& data = face->get_module_data<snow_slide::data>(ID); // Get data

            // Use these placeholders to store the copy as we can then access it on the neighbours
            // vert missing as we can reconstruct it easily
            (*face)["ghost_ss_snowdepthavg_to_xfer"] = 0; // snowdepth to transfer this timestep
            (*face)["ghost_ss_swe_to_xfer"] = 0;

            (*face)["ghost_ss_snowdepthavg_vert_copy"] = data.snowdepthavg_vert_copy;
            (*face)["ghost_ss_delta_avalanche_snowdepth"] = 0;
            (*face)["ghost_ss_delta_avalanche_swe"] = 0;

            for (int j = 0; j < 3; ++j)
            {
                auto n = face->neighbor(j);
                if (n!=nullptr && n->is_ghost)
                {
                    #pragma omp critical
                    {
                        (*n)["ghost_ss_snowdepthavg_to_xfer"] = 0;
                        (*n)["ghost_ss_swe_to_xfer"] = 0;
                        (*n)["ghost_ss_delta_avalanche_snowdepth"] = 0;
                        (*n)["ghost_ss_delta_avalanche_swe"] = 0;
                    }
                }
            }

            sorted_z.at(i) = std::make_pair(
                face->center().z() + data.snowdepthavg_vert_copy, face);

        }

#ifdef USE_MPI
        // update everyone's ghosts with our values
        domain->ghost_neighbors_communicate_variable("ghost_ss_snowdepthavg_vert_copy"_s);
#endif

        // Sort faces by elevation + snowdepth
        tbb::parallel_sort(sorted_z.begin(), sorted_z.end(),
                           [](const std::pair<double, mesh_elem>& a, const std::pair<double, mesh_elem>& b)
                           { return b.first < a.first; });

        // Loop through each face, from highest to lowest triangle surface
        for (size_t i = 0; i < sorted_z.size(); i++)
        {
            auto face = sorted_z[i].second;                           // Get pointer to face
            double cen_area = face->get_area();                       // Area of center triangle
            auto& data = face->get_module_data<snow_slide::data>(ID); // Get stored data for face

            // Get current triangle snow info at beginning of time step
            double maxDepth = data.maxDepth;

            double snowdepthavg = data.snowdepthavg_copy;           // m - Snow depth perpendicular to the surface
            double snowdepthavg_vert = data.snowdepthavg_vert_copy; // m - Vertical snow depth
            double swe = data.swe_copy;                             // m

            // Check if face normal snowdepth have exceeded normal maxDepth
            if (snowdepthavg > maxDepth)
            {
                done = 0;
                this_iter_moved_snow = true;

                double del_depth = snowdepthavg - maxDepth;           // Amount to be removed (positive) [m]
                double del_swe = swe * (1 - maxDepth / snowdepthavg); // Amount of swe to be removed (positive) [m]
                double orig_mass = del_swe * cen_area;

                double z_s = face->center().z() + snowdepthavg_vert; // Current face elevation + vertical snowdepth
                std::vector<double> w = {0, 0, 0};                   // Weights for each face neighbor to route snow to
                double w_dem = 0;                                    // Denomenator for weights (sum of all elev diffs)


                // Calc weights for routing snow
                // Possible Cases:
                //      1) edge cell, then edge_flag is true, and snow is dumped off mesh
                //      2) non-edge cell, w_dem is greater than 0 -> there is at least one lower neighbor, route so to it/them 3) non-edge cell, w_dem = 0, "sink" case. Don't route any snow.
                // Calc weighting based on height diff
                // std::max insures that if one neighbor is higher, its weight will be zero
                for (int i = 0; i < 3; ++i)
                {
                    auto n = face->neighbor(i); // Pointer to neighbor face

                    // this is a domain edge
                    if (n == nullptr)
                    {
                        // pretend our missing face has the same elevation as us, but has no snow so it take can some transport
                        w[i] = std::max(0.0, z_s - face->center().z());

                    }
                    else if (n->is_ghost)
                    {
                        w[i] = std::max(0.0, z_s - (n->center().z() + (*n)["ghost_ss_snowdepthavg_vert_copy"_s]));
                    }
                    // Only non-ghost will have these
                    else
                    {
                        auto& n_data = n->get_module_data<snow_slide::data>(ID); // pointer to face's data
                        w[i] = std::max(0.0, z_s - (n->center().z() + n_data.snowdepthavg_vert_copy));
                    }
                    w_dem += w[i]; // Store weight denominator
                }

                // Case 2) Non-Edge cell, but w_dem=0, "sink" cell. Don't route snow.
                if (w_dem == 0)
                {
                    continue; // Restart to next iteration.
                }

                // Must be Case 3), Divide by sum height differences to create weights that sum to unity
                if (w_dem != 0)
                {
                    std::transform(w.begin(), w.end(), w.begin(), [w_dem](double cw) { return cw / w_dem; });
                }

                // Case 3), Non-Edge cell, w_dem>0, route snow to down slope neighbor(s).
                double out_mass = 0; // Mass balance check
                // Route snow to each neighbor based on weights
                for (int j = 0; j < 3; ++j)
                {
                    auto n = face->neighbor(j);
                    if (n == nullptr)
                    {
                        // Special case: dump snow out of domain (loosing mass) by just removing from current edge cell.
                        out_mass += del_swe * cen_area * w[j];
                    }
                    else //move the mass to a non domain edge neighbour triangle
                    {
                        double n_area = n->get_area(); // Area of neighbor triangle

                        double delta_sd_avg = del_depth * (cen_area / n_area) * w[j]; // (m)
                        double delta_swe = del_swe * (cen_area / n_area) * w[j]; // (m)

                        // Fraction of snowdepth (m) *center triangle area (m^2) = volume of snow depth (m^3)
                        double delta_sd_avg_m3 = del_depth * cen_area * w[j]; // (m3)
                        double delta_swe_m3 = del_swe * cen_area * w[j]; // (m3)

                        if (n->is_ghost)
                        {
                            // amounts to move to ghosts. SD Vert is calculated for normal sd
                            (*n)["ghost_ss_snowdepthavg_to_xfer"_s] += delta_sd_avg;
                            (*n)["ghost_ss_swe_to_xfer"_s] += delta_swe;

                            (*n)["ghost_ss_delta_avalanche_snowdepth"] += delta_sd_avg_m3;
                            (*n)["ghost_ss_delta_avalanche_swe"] += delta_swe_m3;
                        }
                        else
                        {
                            auto& n_data = n->get_module_data<snow_slide::data>(ID); // pointer to face's data

                            // Update neighbor snowdepth and swe (copies only for internal snowSlide use)
                            // Here we must make an assumption of the pack density (because we do not have access to
                            // layer information (if exists (i.e. snowpack is running), or it doesn't
                            // (i.e. snobal is running)) Therefore, we assume uniform density. The (cen_area/n_area)
                            // term converts depth change from orig cell to volume, then back to a depth term using
                            // the neighbor's area.
                            n_data.snowdepthavg_copy += delta_sd_avg; // (m)
                            n_data.swe_copy += delta_swe;            // (m)
                            // Update vertical snow depth
                            n_data.snowdepthavg_vert_copy = n_data.snowdepthavg_copy / std::max(0.001, cos(face->slope()));

                            // Update mass transport to neighbor
                            // Fraction of snowdepth (m) *center triangle area (m^2) = volune of snow depth (m^3)
                            n_data.delta_avalanche_snowdepth += delta_sd_avg_m3;
                            n_data.delta_avalanche_mass += delta_swe_m3; // (m) * (m^2) = (m^3) of swe
                          }


                        out_mass += del_swe * cen_area * w[j];
                    }

                }
                // Remove snow from initial face
                // here we are using the snowmodel normal depth and then covert it to a vert equivalent
                data.snowdepthavg_copy = maxDepth; // data refers to current/center cell
                data.snowdepthavg_vert_copy = data.snowdepthavg_copy / std::max(0.001, cos(face->slope()));
                data.swe_copy = swe * maxDepth / snowdepthavg; // Uses ratio of depth change to calc new swe
                // This relies on the assumption of uniform density.

                // Update mass transport (m^3)
                data.delta_avalanche_snowdepth -= del_depth * cen_area;
                data.delta_avalanche_mass -= del_swe * cen_area;

                // Check mass transport balances for current avalanche cell
                if (std::abs(orig_mass - out_mass) > 0.0001)
                {
                    LOG_DEBUG << "Moved mass total is " << out_mass;
                    LOG_DEBUG << "diff = " << orig_mass - out_mass;
                    LOG_DEBUG << "Mass balance of avalanche times step was not conserved.";
                }
            } // end if snowdepth > maxdepth
        } // End of each face

#ifdef USE_MPI
        // At this point we've set values on the our ghost faces. These correspond with actual faces on other ranks
        // So we need to send these data back
        domain->ghost_to_neighbors_communicate_variable("ghost_ss_snowdepthavg_to_xfer"_s);
        domain->ghost_to_neighbors_communicate_variable("ghost_ss_swe_to_xfer"_s);
        domain->ghost_to_neighbors_communicate_variable("ghost_ss_delta_avalanche_snowdepth"_s);
        domain->ghost_to_neighbors_communicate_variable("ghost_ss_delta_avalanche_swe"_s);
#endif

//#pragma omp parallel for
//        for (size_t i = 0; i < domain->size_faces(); i++)
//        {
//            auto face = domain->face(i); // Get face
//            auto val = (*face)["ghost_ss_snowdepthavg_to_xfer"_s];
//            (*face)["ghost_ss_sum_swe_xfer"] += (*face)["ghost_ss_snowdepthavg_to_xfer"];
//            if (val > 0)
//            {
//                LOG_DEBUG << "Face got from ghost = " << val;
//            }
//        }


        size_t ghost_transport = 0;
        #pragma omp parallel for
        for (size_t i = 0; i < domain->size_faces(); i++)
        {
            auto face = domain->face(i); // Get face                        // Get pointer to face
            auto& data = face->get_module_data<snow_slide::data>(ID); // Get stored data for face

            // these will be != on the triangles owned as ghosts by another rank
            data.snowdepthavg_copy += (*face)["ghost_ss_snowdepthavg_to_xfer"];
            data.snowdepthavg_vert_copy += (*face)["ghost_ss_snowdepthavg_to_xfer"] / std::max(0.001, cos(face->slope()));
            data.swe_copy += (*face)["ghost_ss_swe_to_xfer"_s];

            data.delta_avalanche_snowdepth += (*face)["ghost_ss_delta_avalanche_snowdepth"_s];
            data.delta_avalanche_mass += (*face)["ghost_ss_delta_avalanche_swe"_s]; // (m) * (m^2) = (m^3) of swe

            if((*face)["ghost_ss_delta_avalanche_snowdepth"_s] > 0)
            {
#pragma omp atomic
                ++ghost_transport;
//                LOG_DEBUG << (*face)["ghost_ss_delta_avalanche_snowdepth"_s];
            }

            (*face)["ghost_ss_snowdepthavg_vert_copy"] = data.snowdepthavg_vert_copy;

            // Save state variables at end of time step
            (*face)["delta_avalanche_snowdepth"_s] = data.delta_avalanche_snowdepth;
            (*face)["delta_avalanche_mass"_s] = data.delta_avalanche_mass;

            (*face)["delta_avalanche_snowdepth_sum"_s] += data.delta_avalanche_snowdepth;
            (*face)["delta_avalanche_mass_sum"_s] += data.delta_avalanche_mass;

        }

        // only do another iteration if we have incoming mass transport from the ghosts or if we have moved mass this itr
        // this algorithm tends to need a couple passes to make sure there are no straglers
        if(ghost_transport > 0 || this_iter_moved_snow)
            done = 0; //TODO: put this back to 0
        else
            done = 1;

        ++iterations;

        // a global all reduce to determine the minimum value across all ranks. Min = 0 implies we are not done and need to iterate again
        int global_done=1;
#ifdef USE_MPI
        boost::mpi::all_reduce(domain->_comm_world, done, global_done, boost::mpi::minimum<int>());
#endif

        done = global_done;

        if(!done && iterations > 500)
        {
            //bail
            done = 1;
            LOG_ERROR << "SnowSlide did not converge after 500 iterations";
        }

#ifdef USE_MPI
        if(!done)
            // because we don't have access to ndata.snowdepthavg_copy, pass it through here
            // only used for obtaining transport weights, but only do this comms if we are expecting another iter
            domain->ghost_neighbors_communicate_variable("ghost_ss_snowdepthavg_vert_copy"_s);
#endif

    }while(!done);

    LOG_DEBUG << "[SnowSlide] needed " << iterations << " iterations";

}

void snow_slide::init(mesh& domain)
{
    // Get Parameters that control maxDepth function
    double avalache_mult = cfg.get("avalache_mult",3178.4);
    double avalache_pow  = cfg.get("avalache_pow",-1.998);

    // Initialize for each triangle
    for(size_t i=0;i<domain->size_faces();i++)
    {
        auto face = domain->face(i);

        auto& d = face->make_module_data<snow_slide::data>(ID);

        double Z_CanTop = 0.0;
        if(face->has_vegetation())
        {
            Z_CanTop = face->veg_attribute("CanopyHeight");
        }

        // Parametrize the Minimum snow holding depth (taken vertically)
        double slopeDeg = std::max(10.0,face->slope()*180/M_PI);  // radians to degres, limit to >10 degrees to avoid inf

        // (m) Estimate min depth that avalanche occurs
        // The max depth parameterization is defined for the vertical case, so we need to correct it to be valid for comparing
        // against the "normal" snow depth
        d.maxDepth = std::max(avalache_mult * pow(slopeDeg,avalache_pow),Z_CanTop) * std::max(0.001,cos(face->slope()));

        // Max of either veg height or derived max holding snow depth.
        (*face)["maxDepth"_s]= d.maxDepth;
        (*face)["delta_avalanche_snowdepth_sum"_s] = 0;
        (*face)["delta_avalanche_mass_sum"_s] = 0;
    }
}
