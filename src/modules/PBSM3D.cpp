; //
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a
// novel modular unstructured mesh based approach for hydrological modelling
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
// Canadian Hydrological Model is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//

#include "PBSM3D.hpp"

REGISTER_MODULE_CPP(PBSM3D);

// starting at row_start until row_end find offset for col
inline unsigned int offset(const unsigned int& row_start, const unsigned int& row_end, const unsigned int* col_buffer,
                           const unsigned int& col)
{
    for (unsigned int i = row_start; i < row_end; ++i)
    {
        if (col_buffer[i] == col)
            return i;
    }
    return -1; // wrap it and index garbage
}

struct my_fill_topo_params
{
    double a1;
    double b1;
    double a2;
    double b2;
    double m_tpi;
    double s_tpi;
};

double my_fill_topo(double x, void* p)
{
    struct my_fill_topo_params* params = (struct my_fill_topo_params*)p;
    double a1 = (params->a1);
    double b1 = (params->b1);
    double a2 = (params->a2);
    double b2 = (params->b2);
    double m_tpi = (params->m_tpi);
    double s_tpi = (params->s_tpi);

    double x2 = x + 0.25;
    if (x2 <= 0)
        return (1. - a1 * tanh(b1 * x2)) * gsl_ran_gaussian_pdf(x - m_tpi, s_tpi);
    else
        return (1. - a2 * tanh(b2 * x2)) * gsl_ran_gaussian_pdf(x - m_tpi, s_tpi);
}

struct my_fill_topo2_params
{
    double a1;
    double b1;
    double a2;
    double b2;
    double m_tpi;
    double s_tpi;
    double hs;
    double ii;
};

double my_fill_topo2(double x, void* p)
{
    struct my_fill_topo2_params* params = (struct my_fill_topo2_params*)p;
    double a1 = (params->a1);
    double b1 = (params->b1);
    double a2 = (params->a2);
    double b2 = (params->b2);
    double m_tpi = (params->m_tpi);
    double s_tpi = (params->s_tpi);
    double hs = (params->hs);
    double ii = (params->ii);

    double x2 = x + 0.25;
    if (x2 <= 0)
        return (1. - a1 * tanh(b1 * x2)) * hs / ii * gsl_ran_gaussian_pdf(x - m_tpi, s_tpi);
    else
        return (1. - a2 * tanh(b2 * x2)) * hs / ii * gsl_ran_gaussian_pdf(x - m_tpi, s_tpi);
}

struct my_fill_topo3_params
{
    double m_tpi;
    double s_tpi;
    double fac_fill;
};

double my_fill_topo3(double x, void* p)
{
    struct my_fill_topo3_params* params = (struct my_fill_topo3_params*)p;
    double m_tpi = (params->m_tpi);
    double s_tpi = (params->s_tpi);
    double fac_fill = (params->fac_fill);

    return -fac_fill * x * gsl_ran_gaussian_pdf(x - m_tpi, s_tpi);
}

PBSM3D::PBSM3D(config_file cfg) : module_base("PBSM3D", parallel::domain, cfg)
{
    depends("U_2m_above_srf");
    depends("vw_dir");
    depends("swe");
    depends("t");
    depends("rh");
    depends("U_R");

    depends("p_snow");
    depends("p");

    // Determine if we want fetch, and if so, which one. Default is to use Pomeroy
    // Tanh
    // if we use Liston fetch, we don't need hours since snowfall. If we use
    // Essery 1999, then we need hours since snowfall we have to do this here so
    // we can correctly setup the depends for module linking.
    use_exp_fetch = cfg.get("use_exp_fetch", false);
    use_tanh_fetch = cfg.get("use_tanh_fetch", true);
    use_PomLi_probability = cfg.get("use_PomLi_probability", false);
    z0_ustar_coupling = cfg.get("z0_ustar_coupling", false);

    // Determine if we account for sub-grid topography impact on snow redistribution
    use_subgrid_topo = cfg.get("use_subgrid_topo", false);
    use_subgrid_topo_V2 = cfg.get("use_subgrid_topo_V2", false);

    if (use_exp_fetch && use_tanh_fetch)
        BOOST_THROW_EXCEPTION(module_error() << errstr_info("PBSM3d: Cannot specify both exp_fetch and tanh_fetch"));

    if (use_exp_fetch || use_tanh_fetch)
        depends("fetch");
    else
        depends("p_snow_hours");
    use_R94_lambda = cfg.get("use_R94_lambda", true);

    if (!use_R94_lambda)
    {
        N = cfg.get("N", 1);
        dv = cfg.get("dv", 0.8);
    }

    provides("blowingsnow_probability");
    debug_output = cfg.get("debug_output", false);

    if (debug_output)
    {
        nLayer = cfg.get("nLayer", 5);
        for (int i = 0; i < nLayer; ++i)
        {
            provides("K" + std::to_string(i));
            provides("c" + std::to_string(i));
            provides("cz" + std::to_string(i));
            provides("rm" + std::to_string(i));

            provides("csubl" + std::to_string(i));
            provides("settling_velocity" + std::to_string(i));
            provides("u_z" + std::to_string(i));
        }

        provides("is_drifting");
        provides("Km_coeff");
        provides("Qsusp_pbsm");
        provides("height_diff");
        provides("suspension_mass");
        provides("saltation_mass");
        //        provides("Ti");

        provides("w");
        provides("hs");
        provides("ustar");
        provides("l");
        provides("z0");
        provides("lambda");

        provides("U_10m");

        provides("csalt");
        provides("csalt_orig");
        provides("csalt_reset");
        provides("mass_qsalt");
        provides("c_salt_fetch_big");

        provides("u*_th");
        provides("u*_n");
        provides("tau_n_ratio");

        provides("dm/dt");
        provides("mm");

        provides("Qsubl");
        provides("Qsubl_mass");
        provides("sum_subl");
    }

    provides("drift_mass"); // kg/m^2
    provides("Qsusp");
    provides("Qsalt");

    provides("sum_drift");

    if (use_subgrid_topo)
    {
        provides("frac_contrib");
        provides("frac_contrib_nosnw");
        provides("hold_topo");
    }

    if (use_subgrid_topo_V2)
    {
        provides("frac_contrib");
        provides("hold_topo");
        provides("test_int");
        provides("test_err");
        provides("tpi_lim");
    }
}

void PBSM3D::init(mesh& domain)
{
    nLayer = cfg.get("nLayer", 5);

    susp_depth = 5;                      // 5m as per pomeroy
    v_edge_height = susp_depth / nLayer; // height of each vertical prism
    l__max = 40;                         // mixing length for diffusivity calculations

    do_fixed_settling = cfg.get("do_fixed_settling", false);

    // settling_velocity is used if the user chooses a fixed settling velocity
    // (do_fixed_settling=true)
    settling_velocity = cfg.get("settling_velocity",
                                0.5); // m/s, Lehning, M., H. Löwe, M. Ryser, and N. Raderschall
                                      // (2008), Inhomogeneous precipitation distribution and snow
                                      // transport in steep terrain, Water Resour. Res., 44(7),
                                      // 1–19, doi:10.1029/2007WR006545.

    if (settling_velocity < 0)
        BOOST_THROW_EXCEPTION(module_error() << errstr_info("PBSM3D settling velocity must be positive"));

    do_sublimation = cfg.get("do_sublimation", true);
    do_lateral_diff = cfg.get("do_lateral_diff", true);
    eps = cfg.get("smooth_coeff", 820);
    limit_mass = cfg.get("limit_mass", false);
    min_mass_for_trans = cfg.get("min_mass_for_trans", 5);
    min_sd_trans = cfg.get("min_sd_trans", 0.1); // m Snow holding capacity for flat terrain

    cutoff = cfg.get("cutoff",
                     0.3); // cutoff veg-snow diff (m) that we inhibit saltation entirely

    snow_diffusion_const = cfg.get("snow_diffusion_const",
                                   0.9); // Beta * K, this is beta and scales the eddy diffusivity
    rouault_diffusion_coeff = cfg.get("rouault_diffusion_coef", false);

    enable_veg = cfg.get("enable_veg", true);

    iterative_subl = cfg.get("iterative_subl", false);

    if (rouault_diffusion_coeff)
    {
        LOG_WARNING << "rouault_diffusion_coef overrides const "
                       "snow_diffusion_const values.";
    }

    n_non_edge_tri = 0;

    // use this to build the sparsity pattern for suspension matrix
    size_t ntri = domain->size_faces();
    std::vector<std::map<unsigned int, vcl_scalar_type>> C(ntri * nLayer);

    // sparsity pattern for drift
    std::vector<std::map<unsigned int, vcl_scalar_type>> A(ntri);

    LOG_DEBUG << "#face=" << ntri;

#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {

        auto face = domain->face(i);
        auto d = face->make_module_data<data>(ID);

        if (!face->has_vegetation() && enable_veg)
        {
            LOG_ERROR << "Vegetation is enabled, but no vegetation parameter was found.";
        }
        if (face->has_vegetation() && enable_veg)
        {
            d->CanopyHeight = face->veg_attribute("CanopyHeight");

            // only grab LAI if we are using the R90 lambda formulation
            if (use_R94_lambda)
                d->LAI = face->veg_attribute("LAI");
            else
                d->LAI = 0;
        }
        else
        {
            d->CanopyHeight = 0;
            d->LAI = 0;
            enable_veg = false;
        }

        d->F_fill.function = &my_fill_topo;
        d->F_fill2.function = &my_fill_topo2;
        d->F_fill3.function = &my_fill_topo3;

        // pre alloc for the windpseeds
        d->u_z_susp.resize(nLayer);

        auto& m = d->m;
        // edge unit normals
        m[0].set_size(3);
        m[0](0) = face->edge_unit_normal(0).x();
        m[0](1) = face->edge_unit_normal(0).y();
        m[0](2) = 0;

        m[1].set_size(3);
        m[1](0) = face->edge_unit_normal(1).x();
        m[1](1) = face->edge_unit_normal(1).y();
        m[1](2) = 0;

        m[2].set_size(3);
        m[2](0) = face->edge_unit_normal(2).x();
        m[2](1) = face->edge_unit_normal(2).y();
        m[2](2) = 0;

        // top
        m[3].set_size(3);
        m[3].fill(0);
        m[3](2) = 1;

        // bottom
        m[4].set_size(3);
        m[4].fill(0);
        m[4](2) = -1;

        // face areas
        for (int j = 0; j < 3; ++j)
            d->A[j] = face->edge_length(j) * v_edge_height;

        // top, bottom
        d->A[3] = d->A[4] = face->get_area();

        d->is_edge = false;
        // which faces have neighbours? Ie, are we an edge?
        for (int a = 0; a < 3; ++a)
        {
            A[i][i] = -9999;

            auto neigh = face->neighbor(a);
            if (neigh == nullptr || neigh->_is_ghost)
            {
                d->face_neigh[a] = false;
                d->is_edge = true;
            }
            else
            {
                d->face_neigh[a] = true;

                A[i][neigh->cell_local_id] = -9999;
            }
        }
        if (!d->is_edge)
        {
            d->cell_local_id = n_non_edge_tri;
            ++n_non_edge_tri;
        }

        d->sum_drift = 0;
        d->sum_subl = 0;

        // iterate over the vertical layers
        for (int z = 0; z < nLayer; ++z)
        {
            size_t idx = ntri * z + face->cell_local_id;
            for (int f = 0; f < 3; f++)
            {
                if (d->face_neigh[f])
                {
                    size_t nidx = ntri * z + face->neighbor(f)->cell_local_id;
                    C[idx][idx] = -9999;
                    C[idx][nidx] = -9999;
                }
                else
                {
                    C[idx][idx] = -9999;
                }
            }

            if (z == 0)
            {
                C[idx][idx] = -9999;
                C[idx][ntri * (z + 1) + face->cell_local_id] = -9999;
            }
            else if (z == nLayer - 1)
            {
                C[idx][idx] = -9999;
                C[idx][ntri * (z - 1) + face->cell_local_id] = -9999;
            }
            else // middle layers
            {
                C[idx][idx] = -9999;
                C[idx][ntri * (z + 1) + face->cell_local_id] = -9999;
                C[idx][ntri * (z - 1) + face->cell_local_id] = -9999;
            }
        }
    }

    viennacl::copy(C, vl_C); // copy C -> vl_C, sets up the sparsity pattern
    viennacl::copy(A, vl_A); // copy A -> vl_A, sets up the sparsity pattern

    b.resize(ntri * nLayer);
    bb.resize(ntri);
    nnz = vl_C.nnz();
    nnz_drift = vl_A.nnz();
}

void PBSM3D::run(mesh& domain)
{

    // needed for linear system offsets
    size_t ntri = domain->size_faces();

    // vcl_scalar_type is defined in the main CMakeLists.txt file.
    // Some GPUs do not have double precision so the run will fail if the wrong
    // precision is used

#ifdef VIENNACL_WITH_OPENCL
    viennacl::context host_ctx(viennacl::MAIN_MEMORY);
    vl_C.switch_memory_context(host_ctx);
    vl_A.switch_memory_context(host_ctx);

    b.switch_memory_context(host_ctx);
    bb.switch_memory_context(host_ctx);
#endif

    // zero CSR vector in vl_C
    viennacl::vector_base<vcl_scalar_type> init_temporary(
        vl_C.handle(), viennacl::compressed_matrix<vcl_scalar_type>::size_type(nnz + 1), 0, 1);
    // write:
    init_temporary = viennacl::zero_vector<vcl_scalar_type>(
        viennacl::compressed_matrix<vcl_scalar_type>::size_type(nnz + 1), viennacl::traits::context(vl_C));

    // zero-fill RHS
    b.clear();

    // get row buffer
    unsigned int const* row_buffer =
        viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(vl_C.handle1());
    unsigned int const* col_buffer =
        viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(vl_C.handle2());
    vcl_scalar_type* elements =
        viennacl::linalg::host_based::detail::extract_raw_pointer<vcl_scalar_type>(vl_C.handle());

#pragma omp parallel
    {
        // Helpers for the u* iterative solver
        // - needs to be here in thread pool, otherwise there are thread consistency
        // issues with the solver
        boost::uintmax_t max_iter = 500;
        auto tol = [](double a, double b) -> bool { return fabs(a - b) < 1e-8; };

        // ice density
        double rho_p = PhysConst::rho_ice;
#pragma omp for
        for (size_t i = 0; i < domain->size_faces(); i++)
        {

            auto face = domain->face(i);

            auto id = face->cell_local_id;

            auto d = face->get_module_data<data>(ID);
            auto& m = d->m;

            double fetch = 1000;
            if (use_exp_fetch || use_tanh_fetch)
                fetch = (*face)["fetch"_s];

            double frac_contrib = 1.;       // Default value for the fraction of the grid contributing to snow transport
            double frac_contrib_nosnw = 1.; // Default value for the fraction of the grid contributing to snow transport
            double min_sd_trans_avg = min_sd_trans; // Grid-averaged value for the topographic subgrid holding capacity

            // get wind from the face
            double uref = (*face)["U_R"_s];
            double snow_depth = (*face)["snowdepthavg"_s];
            snow_depth = is_nan(snow_depth) ? 0 : snow_depth;

            double u2 = (*face)["U_2m_above_srf"_s];
            double z10; // 10-m height above the snow surface
            z10 = 10. + snow_depth;

            double u10;
            if (z10 < Atmosphere::Z_U_R)
            {
                u10 = Atmosphere::log_scale_wind(uref, Atmosphere::Z_U_R, z10,
                                                 snow_depth); // used by the pom probability forumuation, so don't
                                                              // hide behide debug output
            }
            else
            {
                u10 = uref; // Extreme case (avalanche gone crazy case)
            }

            if (debug_output)
                (*face)["U_10m"_s] = u10;

            double swe = (*face)["swe"_s]; // mm   -->    kg/m^2
            swe = is_nan(swe) ? 0 : swe;   // handle the first timestep where swe won't have been
                                           // updated if we override the module order

            // height difference between snowcover and veg
            double height_diff = std::max(0.0, d->CanopyHeight - snow_depth);
            if (!enable_veg)
                height_diff = 0;
            if (debug_output)
                (*face)["height_diff"_s] = height_diff;

            // Topographic holding capacity associated with subgrid topographic features
            // This method combines the distribution of TPI with a filling function to obtain
            // an estimation of the subgrid snow depth distribution per triangle
            // A filling criteria is then used to detertmine which fraction of the triangle
            // contributes to snow transport (area of positive TPI + filled gullies)
            // and which fraction does not (gullies which are not filled).

            if (use_subgrid_topo_V2)
            {
                // Areas with negative TPI are assumed to be filled when SD = fac_fill * TPI.
                double fac_fill = 0.8;

                // Default values for the TPI threshold above which gullies are filled.
                double tpi_lim = -min_sd_trans;

                if (!is_nan(face->parameter("TPI_std"_s)))
                { // Std value of TPI is defined

                    double moy_tpi = std::max(-5.0, std::min(5.0, face->parameter("TPI_mean"_s)));
                    double std_tpi = std::min(5.0, std::max(0.1, face->parameter("TPI_std"_s)));

                    if (snow_depth > min_sd_trans)
                    {

                        // Coefficient for the filling function that give normalized snow depth as a function of TPI

                        double a1 = 1.5;
                        double b1 = 0.3;
                        double a2 = 0.6;
                        double b2 = 0.55;
                        if (snow_depth > 0.75 and snow_depth < 1.25)
                        {
                            double a1 = 1.15;
                        }
                        else if (snow_depth < 1.75)
                        {
                            double a1 = 1.1;
                            double b2 = 0.4;
                        }
                        else if (snow_depth < 2.25)
                        {
                            double a1 = 0.9;
                            double b2 = 0.4;
                        }
                        else if (snow_depth < 2.75)
                        {
                            double a1 = 0.85;
                            double b2 = 0.4;
                        }
                        else if (snow_depth < 3.25)
                        {
                            double a1 = 0.75;
                            double b2 = 0.4;
                        }
                        else
                        {
                            double a1 = 0.6;
                            double b2 = 0.35;
                        }

                        // Compute normalization factor
                        struct my_fill_topo_params params = {a1, b1, a2, b2, moy_tpi, std_tpi};
                        d->F_fill.params = &params;

                        gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
                        double result, error;
                        int code = gsl_integration_qags(&d->F_fill, -50, 50, 0, 1e-7, 1000, w, &result, &error);
                        gsl_integration_workspace_free(w);

                        (*face)["test_int"_s] = result;

                        // Determine TPI threshold above which gullies are considered as filled.
                        auto frootFn = [&](double xx) -> double {
                            return (1 - a1 * tanh(b1 * (xx + 0.25))) * snow_depth / result + fac_fill * xx;
                        };
                        try
                        {
                            auto r =
                                boost::math::tools::bracket_and_solve_root(frootFn, -1.0, 1.0, true, tol, max_iter);
                            tpi_lim = r.first + (r.second - r.first) / 2.0;
                        }
                        catch (...)
                        {
                            // Didn't converge
                        }

                        // Determine area-averaged snow depth which is stored in the non-filled gullies
                        struct my_fill_topo2_params params2 = {a1, b1, a2, b2, moy_tpi, std_tpi, snow_depth, result};
                        d->F_fill2.params = &params2;

                        gsl_integration_workspace* w2 = gsl_integration_workspace_alloc(1000);
                        double h1;
                        int code2 = gsl_integration_qags(&d->F_fill2, -50, tpi_lim, 0, 1e-7, 1000, w2, &h1, &error);
                        gsl_integration_workspace_free(w2);

                        // Determine area-averaged snow depth which is stored in the filled gullies
                        struct my_fill_topo3_params params3 = {moy_tpi, std_tpi, fac_fill};
                        d->F_fill3.params = &params3;

                        gsl_integration_workspace* w3 = gsl_integration_workspace_alloc(1000);
                        double h2;
                        int code3 = gsl_integration_qags(&d->F_fill3, tpi_lim, -min_sd_trans / fac_fill, 0, 1e-7, 1000,
                                                         w3, &h2, &error);
                        gsl_integration_workspace_free(w3);

                        // Determine area-averaged snow depth hold in the area of positive TPI
                        double h3 = min_sd_trans * gsl_cdf_gaussian_Q(-min_sd_trans / fac_fill - moy_tpi, std_tpi);

                        // Compute total holding capacity
                        min_sd_trans_avg = std::min(h1 + h2 + h3, snow_depth);
                    }

                    // Determine fraction of the triangle that contributes to snow transport
                    frac_contrib = gsl_cdf_gaussian_Q(tpi_lim - moy_tpi, std_tpi);
                }
                (*face)["tpi_lim"_s] = tpi_lim;
                (*face)["frac_contrib"_s] = frac_contrib;
                (*face)["hold_topo"_s] = min_sd_trans_avg;
            }

            // Topographic holding capacity associated with subgrid topographic features
            // This method uses the distribution of negative TPI and the mean snow depth to determine the
            // fraction of the triangle which contributes to snow transport and the associated topographic
            // holding capacity

            if (use_subgrid_topo)
            {

                if (!is_nan(face->parameter("TPI_neg_frac"_s))) // Fraction of negative TPI is defined
                {
                    double frac_neg = face->parameter("TPI_neg_frac"_s); // fraction of the grid covered by negative TPI
                    if (frac_neg > 0.)
                    {
                        double moy_tpi_neg = std::max(-10.0, std::min(-0.05, face->parameter("TPI_neg_mean"_s)));
                        double std_tpi_neg = std::min(5.0, std::max(0.1, face->parameter("TPI_neg_std"_s)));

                        //              // Derive parameters of the gamma distribution representing the distrib of TPI
                        double shape_gam = pow(moy_tpi_neg, 2.0) / pow(std_tpi_neg, 2.0);
                        double scale_gam = -pow(std_tpi_neg, 2.0) / moy_tpi_neg;
                        //
                        frac_contrib = (1. - frac_neg) + // f_TPI>0.
                                       frac_neg * gsl_cdf_gamma_P(snow_depth, shape_gam, scale_gam);

                        frac_contrib_nosnw = (1. - frac_neg);

                        if (snow_depth > min_sd_trans)
                        {

                            //   double fac = shape_gam/(gsl_sf_gamma(shape_gam)*pow(scale_gam,shape_gam));
                            double hint =
                                shape_gam * scale_gam -
                                scale_gam * (snow_depth * gsl_ran_gamma_pdf(snow_depth, shape_gam, scale_gam) -
                                             min_sd_trans * gsl_ran_gamma_pdf(min_sd_trans, shape_gam, scale_gam) +
                                             shape_gam * gsl_cdf_gamma_P(min_sd_trans, shape_gam, scale_gam) +
                                             shape_gam * gsl_cdf_gamma_Q(snow_depth, shape_gam, scale_gam));

                            min_sd_trans_avg =
                                (1. - frac_neg) * min_sd_trans +
                                frac_neg * (min_sd_trans * gsl_cdf_gamma_P(min_sd_trans, shape_gam, scale_gam) + hint +
                                            snow_depth * gsl_cdf_gamma_Q(snow_depth, shape_gam, scale_gam));
                        }
                    }
                }
                (*face)["frac_contrib"_s] = frac_contrib;
                (*face)["frac_contrib_nosnw"_s] = frac_contrib_nosnw;
                (*face)["hold_topo"_s] = min_sd_trans_avg;
            }

            double ustar = 1.3; // placeholder

            // The strategy here is as follows:
            // 0) If the exposed vegetation is above $cutoff, inhibit saltation and
            // use the classical vegheight*0.12=z0 and use that for calculating u*
            // 1) Wait for vegetation to fill up until it is within $cutoff of the
            // top of the veg
            //     then use Pomeroy & Li 2000 eqn 4 to calculate an iterative
            //     solution to u* under blowing snow conditions This effectively
            //     allows wind to blow snow out of the vegetation
            // 2) Calculate a u* that uses the blowing snow z0. Then test this
            // against the blowing snow u* threshold 3) If blowing snow isn't
            // happening, recalculate u* using the normal z0.

            // This is the lamdba from Li and Pomeroy eqn 4 that is used to include
            // exposed vegetation w/ the z0 estimate
            double lambda = 0;

            d->saltation = false; // default case

            // threshold friction velocity. Compute here as it's used below as well
            // Pomeroy and Li, 2000
            // Eqn 7
            double T = (*face)["t"_s];
            double u_star_saltation_threshold =
                0.35 + (1.0 / 150.0) * T + (1.0 / 8200.0) * T * T; // saltation threshold m/s
            if (debug_output)
                (*face)["u*_th"_s] = u_star_saltation_threshold;

            // we don't have too high of veg. Check for blowing snow
            //   if (height_diff <= cutoff && swe >= min_mass_for_trans &&
            //       !is_water(face))
            if (height_diff <= cutoff && snow_depth >= min_sd_trans_avg && !is_water(face))
            {

                // lambda -> 0 when height_diff ->, such as full or no veg
                if (use_R94_lambda)
                    // LAI/2.0 suggestion from Raupach 1994 (DOI:10.1007/BF00709229)
                    // Section 3(a)
                    lambda = 0.5 * d->LAI * height_diff;
                else
                    lambda = N * dv * height_diff; // Pomeroy formulation

                if (debug_output)
                    (*face)["lambda"_s] = lambda;

                if (z0_ustar_coupling)
                {
                    // Calculate the new value of z0 to take into account partially filled
                    // vegetation and the momentum sink
                    auto ustarFn = [&](double ustar) -> double {
                        // Li and Pomeroy 2000, eqn 5.
                        // This formulation has the following coeffs built in
                        // c_2 = 1.6;
                        // c_3 = 0.07519;
                        // c_4 = 0.5;
                        // g   = 9.81;

                        return u2 * PhysConst::kappa / log(2.0 / (0.6131702345e-2 * ustar * ustar + .5 * lambda)) -
                               ustar;
                    };
                    try
                    {
                        auto r = boost::math::tools::bracket_and_solve_root(ustarFn, 1.0, 1.0, false, tol, max_iter);
                        ustar = r.first + (r.second - r.first) / 2.0;
                    }
                    catch (...)
                    {
                        // Didn't converge
                        d->saltation = false;
                    }
                }
                else
                {
                    // follow PBSM (Pom & Li 2000; Alpine3D) and don't calculate the feedback of z0 on u*
                    ustar = u2 * PhysConst::kappa / log(2.0 / Snow::Z0_SNOW);
                }

                if (ustar >= u_star_saltation_threshold)
                {
                    d->saltation = true;

                    if (z0_ustar_coupling)
                    {
                        // Update z0 for blowing snow conditions
                        // Li and Pomeroy 2000, eqn 5.
                        // This formulation has the following coeffs built in
                        // c_2 = 1.6;
                        // c_3 = 0.07519;
                        // c_4 = 0.5;
                        // g   = 9.81;
                        d->z0 = 0.6131702345e-2 * ustar * ustar + .5 * lambda; // pom and li 2000, eqn 4
                    }
                    else
                    {
                        d->z0 = Snow::Z0_SNOW;
                    }
                }
            }

            if (!d->saltation)
            {
                // we still need a u* for spatial K estimation later
                d->z0 = Snow::Z0_SNOW;
                ustar = std::max(0.01, PhysConst::kappa * uref / log(Atmosphere::Z_U_R / d->z0));
            }

            // sanity checks
            d->z0 = std::max(Snow::Z0_SNOW, d->z0);
            ustar = std::max(0.01, ustar);
            if (debug_output)
                (*face)["ustar"_s] = ustar;
            if (debug_output)
                (*face)["z0"_s] = d->z0;

            // depth of saltation layer
            double hs = 0;
            if (d->saltation)
                hs = 0.08436 * pow(ustar, 1.27); // pomeroy

            d->hs = hs;
            if (debug_output)
                (*face)["hs"_s] = hs;
            if (debug_output)
                (*face)["is_drifting"_s] = 0;
            if (debug_output)
                (*face)["Qsusp_pbsm"_s] = 0; // for santiy checks against pbsm

            double Qsalt = 0;
            double c_salt = 0;
            double t = (*face)["t"_s] + 273.15;

            // Check if we can blow snow in this triagnle
            // Are we above saltation threshold?
            // Do we have enough mass in this triangle?
            // Has saltation been disabled because there is too much veg?
            if (d->saltation)
            {

                double rho_f =
                    mio::Atmosphere::stdDryAirDensity(face->get_z(),
                                                      t); // air density kg/m^3, comment in mio is wrong.1.225;

                if (debug_output)
                    (*face)["blowingsnow_probability"_s] = 0; // default to 0%

                if (debug_output)
                {
                    double pbsm_qsusp = pow(u10, 4.13) / 674100.0;
                    (*face)["Qsusp_pbsm"_s] = pbsm_qsusp;
                }

                if (debug_output)
                    (*face)["is_drifting"_s] = 1;

                // Pomeroy and Li 2000, eqn 8
                double Beta = 202.0; // 170.0;
                double m = 0.16;

                // tau_n_ratio = ustar_n^2 / ustar^2 from MacDonald 2009 eq 3;
                // we need (ustar_n / ustar)^2 which the original derivation gives
                // so this is is correctly squared
                double tau_n_ratio = (m * Beta * lambda) / (1.0 + m * Beta * lambda);

                if (debug_output)
                    (*face)["tau_n_ratio"_s] = tau_n_ratio;

                // Pomeroy 1992, eqn 12, see note above for ustar_n calc, but ustar_n
                // is correctly squared already
                c_salt =
                    rho_f / (3.29 * ustar) *
                    (1.0 - tau_n_ratio - (u_star_saltation_threshold * u_star_saltation_threshold) / (ustar * ustar));

                // occasionally happens to happen at low wind speeds where the
                // parameterization breaks.
                if (c_salt < 0 || std::isnan(c_salt))
                {
                    c_salt = 0;
                    d->saltation = false;
                }

                if (debug_output)
                    (*face)["c_salt_fetch_big"_s] = c_salt;

                // exp decay of Liston, eq 10
                // 95% of max saltation occurs at fetch = 500m
                // Liston, G., & Sturm, M. (1998). A snow-transport model for complex
                // terrain. Journal of Glaciology.
                if (use_exp_fetch && fetch < 500)
                {
                    double fetch_ref = 500;
                    double mu = 3.0;
                    c_salt *= 1.0 - exp(-mu * fetch / fetch_ref);
                }
                else if (use_tanh_fetch && fetch <= 300.) // use Pomeroy & Male 1986 tanh fetch
                {
                    double fetch_ref = 300;
                    double Lc = 0.5 * tanh(0.1333333333e-1 * fetch_ref - 2.0) + 0.5;

                    c_salt *= Lc;
                }

                // consider the temporal non-steady effects
                if (use_PomLi_probability) // Pomeroy and Li 2000 upscaled
                                           // probability
                {
                    // Essery, Li, and Pomeroy 1999
                    // Probability of blowing snow
                    double A = (*face)["p_snow_hours"_s];                              // hours since last snowfall
                    double u_mean = 11.2 + 0.365 * T + 0.00706 * T * T + 0.9 * log(A); // eqn 10  T -> air temp, degC
                    double delta = 0.145 * T + 0.00196 * T * T + 4.3;                  // eqn 11
                    double Pu10 = 1.0 / (1.0 + exp((sqrt(M_PI) * (u_mean - u10)) / delta)); // eqn 12
                    (*face)["blowingsnow_probability"_s] = Pu10;

                    // decrease the saltation by the probability amount
                    c_salt *= Pu10;
                }

                // consider subgrid topographic effect
                if (use_subgrid_topo or use_subgrid_topo_V2)
                {
                    c_salt *= frac_contrib;
                }

                // wind speed in the saltation layer Pomeroy and Gray 1990
                double uhs = 2.8 * u_star_saltation_threshold; // eqn 7

                // kg/(m*s)
                Qsalt = c_salt * uhs * hs; // integrate over the depth of the saltation layer, kg/(m*s)

                double mass = 0;
                double phi = (*face)["vw_dir"_s];
                Vector_2 v = -math::gis::bearing_to_cartesian(phi);

                // setup wind vector
                arma::vec uvw(3);
                uvw(0) = v.x(); // U_x
                uvw(1) = v.y(); // U_y
                uvw(2) = 0;
                double V = face->get_area();
                double udotm[3] = {0, 0, 0};
                double E[3] = {0, 0, 0};
                for (int j = 0; j < 3; ++j)
                {
                    udotm[j] = arma::dot(uvw, d->m[j]);
                    E[j] = face->edge_length(j);
                    mass += -E[j] * Qsalt * udotm[j];
                }

                mass /= V * global_param->dt();

                if (debug_output)
                {
                    (*face)["csalt_orig"_s] = c_salt;
                    (*face)["mass_qsalt"_s] = mass;
                }

                //          if( mass > swe)
                //          {
                //              c_salt = -swe*V/(hs*uhs*(E[0]*udotm[0]+E[1]*udotm[1]+E[2]*udotm[2])*global_param->dt());
                //              // kg/(m*s)
                //              Qsalt =
                //                      c_salt * uhs *
                //                      hs; // integrate over the depth of the saltation layer, kg/(m*s)
                //
                //              if (debug_output)
                //                  (*face)["csalt_reset"_s] = c_salt;
                //          }
            }

            if (debug_output)
                (*face)["csalt"_s] = c_salt;

            (*face)["Qsalt"_s] = Qsalt;

            double rh = (*face)["rh"_s] / 100.;
            double es = mio::Atmosphere::saturatedVapourPressure(t);
            double ea = rh * es / 1000.; // ea needs to be in kpa

            double v = 1.88e-5; // kinematic viscosity of air, below eqn 13 in Pomeroy 1993

            // iterate over the vertical layers
            for (int z = 0; z < nLayer; ++z)
            {
                // height in the suspension layer, floats above the snow surface
                double cz = z * v_edge_height + hs + v_edge_height / 2.; // cell center height

                // compute new U_z at this height in the suspension layer
                double u_z = 0;

                // Height above the ground (snow+free) of the suspension layer
                double hz = cz + snow_depth;

                // the suspension layer discretization 'floats' on top of the snow
                // surface so height_diff = d->CanopyHeight - snowdepth which is
                // looking to see if cz is within this part of the canopy
                if (d->saltation && cz < height_diff)
                {
                    // saltating so used the z0 with veg, but we are in the canopy so
                    // use the saltation vel

                    // wind speed in the saltation layer Pomeroy and Gray 1990
                    u_z = 2.8 * u_star_saltation_threshold; // eqn 7
                }
                else if (cz < height_diff)
                {
                    // we/re in a canopy, but not saltating, just do nothing
                    u_z = 0.01; // essentially do nothing when we are in sub canopy

                    //                // LAI used as attenuation coefficient introduced
                    //                by Inoue (1963) and increases with canopy density
                    //                double LAI =
                    //                std::max(0.01,face->veg_attribute("LAI"));
                    //                //bring wind down to canopy top
                    //                double u_cantop = std::max(0.01,
                    //                Atmosphere::log_scale_wind(uref,
                    //                Atmosphere::Z_U_R, d->CanopyHeight, 0 , d->z0));
                    //
                    //                u_z = Atmosphere::exp_scale_wind(u_cantop,
                    //                d->CanopyHeight, cz, LAI);
                }
                else
                {
                    if (hz < Atmosphere::Z_U_R)
                    {
                        u_z =
                            std::max(0.01, Atmosphere::log_scale_wind(uref, Atmosphere::Z_U_R, hz, snow_depth, d->z0));
                    }
                    else
                    {
                        u_z = std::max(0.01, uref);
                    }
                }

                d->u_z_susp.at(z) = u_z;

                // calculate dm/dt from
                // equation 13 from Pomeroy and Li 2000
                // To do so, use equations 12 - 16 in Pomeroy et al 2000
                // Pomeroy, J. W., and L. Li (2000), Prairie and arctic areal snow
                // cover mass balance using a blowing snow model, J. Geophys. Res.,
                // 105(D21), 26619–26634, doi:10.1029/2000JD900149. [online] Available
                // from: http://www.agu.org/pubs/crossref/2000/2000JD900149.shtml

                // these are from
                // Pomeroy, J. W., D. M. Gray, and P. G. Landine (1993), The prairie
                // blowing snow model: characteristics, validation, operation, J.
                // Hydrol., 144(1–4), 165–192.

                // eqn 18, mean particle radius
                // This is 'r_r' in Pomeroy and Gray 1995, eqn 53
                double rm = 4.6e-5 * pow(cz, -0.258);
                if (debug_output)
                {
                    (*face)["rm" + std::to_string(z)] = rm;
                    (*face)["cz" + std::to_string(z)] = cz;
                }

                // calculate mean mass, eqn 23, 24 in Pomeroy 1993 (PBSM)
                // 52, 53 P&G 1995
                double mm_alpha = 4.08 + 12.6 * cz; // 24
                double mm = 4. / 3. * M_PI * rho_p * rm * rm * rm *
                            (1.0 + 3.0 / mm_alpha + 2. / (mm_alpha * mm_alpha)); // mean mass, eqn 23

                // mean radius of mean mass particle
                double r_z = pow((3.0 * mm) / (4 * M_PI * rho_p), 0.3333333); // 50 in p&g 1995
                if (debug_output)
                    (*face)["mm"_s] = mm;

                double xrz = 0.005 * pow(u_z, 1.36); // eqn 16

                double omega = settling_velocity;
                ; // Settling velocity
                if (!do_fixed_settling)
                {
                    omega = 1.1e7 * pow(r_z, 1.8); // eqn 15 settling velocity
                }

                if (debug_output)
                    (*face)["settling_velocity" + std::to_string(z)] = omega;
                double Vr = omega + 3.0 * xrz * cos(M_PI / 4.0); // eqn 14

                double v = 1.88e-5;             // kinematic viscosity of air, below eqn 13 in
                                                // Pomeroy 1993
                double Re = 2.0 * r_z * Vr / v; // eqn  55 in p&g 1995

                double Nu, Sh;
                Nu = Sh = 1.79 + 0.606 * pow(Re, 0.5); // eqn 12

                // define above, T is in C, t is in K

                // (A.6)
                double D = 2.06e-5 * pow(t / 273.15,
                                         1.75); // diffusivity of water vapour in air, t in K,
                                                // eqn A-7 in Liston 1998 or Harder 2013 A.6

                // (A.9)
                double lambda_t =
                    0.000063 * t + 0.00673; //  thermal conductivity, user Harder 2013 A.9, Pomeroy's
                                            //  is off by an order of magnitude, this matches this
                                            //  https://www.engineeringtoolbox.com/air-properties-d_156.html

                // Standard constant value, e.g.,
                // https://link.springer.com/referenceworkentry/10.1007%2F978-90-481-2642-2_329
                //          double L = 2.38e6; // Latent heat of sublimation, J/kg
                double L = 2.838e6; // Latent heat of sublimation, J/kg, Corrected value

                double dmdtz = 0;

                // use Pomeroy and Li 2000 iterative sol'n for Schmidt's equation
                if (iterative_subl)
                {
                    /*
                     * The *1000 and /1000 are important unit conversions. Doesn't quite
                     * match the harder paper, but Phil assures me it is correct.
                     */
                    double mw = 0.01801528 * 1000.0; //[kg/mol]  ---> g/mol
                    double R = 8.31441 / 1000.0;     // [J mol-1 K-1]

                    double rho = (mw * ea) / (R * t);

                    // use Harder 2013 (A.5) Formulation, but Pa formulation for e
                    auto fx = [=](double Ti) {
                        return boost::math::make_tuple(
                            T +
                                D * L *
                                    (rho / (1000.0) -
                                     .611 * mw * exp(17.3 * Ti / (237.3 + Ti)) / (R * (Ti + 273.15) * (1000.0))) /
                                    lambda_t -
                                Ti,
                            D * L *
                                    (-0.6110000000e-3 * mw * (17.3 / (237.3 + Ti) - 17.3 * Ti / pow(237.3 + Ti, 2)) *
                                         exp(17.3 * Ti / (237.3 + Ti)) / (R * (Ti + 273.15)) +
                                     0.6110000000e-3 * mw * exp(17.3 * Ti / (237.3 + Ti)) / (R * pow(Ti + 273.15, 2))) /
                                    lambda_t -
                                1);
                    };

                    double guess = T;
                    double min = -50;
                    double max = 0;
                    int digits = 6;

                    double Ti = boost::math::tools::newton_raphson_iterate(fx, guess, min, max, digits);
                    double Ts = Ti + 273.15; // dmdtz expects in K

                    // now use equation 13 with our solved Ts to compute dm/dt(z)
                    dmdtz = 2.0 * M_PI * rm * lambda_t / L * Nu * (Ts - (t + 273.15)); // eqn 13 in Pomeroy and Li 2000
                }
                else // Use PBSM. Eqn 11 Pomeroy, Gray, Ladine, 1993  "﻿The
                     // prairie blowing snow model: characteristics, validation,
                     // operation"
                {
                    double M = 18.01; // molecular weight of water kg kmol-1
                    double R = 8313;  // universal fas constant J mol-1 K-1

                    double sigma =
                        (rh - 1.0) * (1.019 + 0.027 * log(cz)); // undersaturation, Pomeroy and Li 2000, eqn 14

                    double rho = (M * es) / (R * t); // saturation vapour density at t
                    // radiative energy absorebed by the particle -- take from CRHM's
                    // PBSM implimentation
                    double Qr = 0.9 * M_PI * rm * rm * 120.0; // 120.0 = PBSM_constants::Qstar (Solar Radiation
                                                              // Input), 0.9 comes from Schmidt (1972) assuming a
                                                              // snow particle albedo of 0.5 and a snow surface
                                                              // albedo of 0.8
                                                              // rm ise used here as in Liston and Sturm (1998)

                    // eqn 11 in PGL 1993, r_z is used here as in PG95 and Liston ans Sturm (1998)
                    dmdtz = Sh * rho * D *
                            (6.283185308 * Nu * R * r_z * sigma * t * t * lambda_t - L * M * Qr + Qr * R * t) /
                            (D * L * Sh * (L * M - R * t) * rho + lambda_t * t * t * Nu * R);
                }
                if (debug_output)
                    (*face)["dm/dt"_s] = dmdtz;

                if (debug_output)
                    (*face)["mm"_s] = mm;
                double csubl = dmdtz / mm; // EQN 21 POMEROY 1993 (PBSM)

                //      if (debug_output)
                //        (*face)["csubl" + std::to_string(z)] = dmdtz;

                // eddy diffusivity (m^2/s)
                // 0,1,2 will all be K = 0, as no horizontal diffusion process
                double K[5] = {0, 0, 0, 0, 0};

                // holds A_ * K_ / h_
                // _0 -> _2 are the horizontal sides
                // _3 -> is the top of the prism
                // _4 -> is the bottom the prism
                double alpha[5] = {0, 0, 0, 0, 0};

                // compute alpha and K for edges
                if (do_lateral_diff)
                {
                    for (int a = 0; a < 3; ++a)
                    {
                        auto neigh = face->neighbor(a);
                        alpha[a] = d->A[a];

                        // do just very low horz diffusion for numerics
                        K[a] = 0.00001;
                        alpha[a] *= K[a];
                    }
                }
                // Li and Pomeroy 2000
                double l = PhysConst::kappa * (cz + d->z0) * l__max / (PhysConst::kappa * (cz + d->z0) + l__max);
                if (debug_output)
                    (*face)["l"_s] = l;

                double w = omega; // settling_velocity;
                if (debug_output)
                    (*face)["w"_s] = w;

                double diffusion_coeff = snow_diffusion_const; // snow_diffusion_const is a shared param so
                                                               // need a seperate copy here we can
                                                               // overwrite
                if (rouault_diffusion_coeff)
                {
                    double c2 = 1.0;
                    double dc = 1.0 / (1.0 + (c2 * w * w) / (1.56 * ustar * ustar));
                    diffusion_coeff = dc; // nope, snow_diffusion_const is shared, use a new
                }
                if (debug_output)
                    (*face)["Km_coeff"_s] = diffusion_coeff;

                // snow_diffusion_const is pretty much a calibration constant. At 1 it
                // seems to over predict transports.
                // with pomeroy fall velocity, 0.3 gives good agreement w/ published
                // Qsusp values. Low value compensates for low fall velocity
                K[3] = K[4] = diffusion_coeff * ustar * l;

                if (debug_output)
                    (*face)["K" + std::to_string(z)] = K[3];
                // top
                alpha[3] = d->A[3] * K[3] / v_edge_height;
                // bottom
                alpha[4] = d->A[4] * K[4] / v_edge_height;

                double phi = (*face)["vw_dir"_s]; // wind direction
                Vector_2 vwind = -math::gis::bearing_to_cartesian(phi);

                // setup wind vector
                arma::vec uvw(3);
                uvw(0) = vwind.x(); // U_x
                uvw(1) = vwind.y(); // U_y
                uvw(2) = 0;

                // above we have just the direction so it's unit vector. scale it to
                // have the same magnitude as u_z
                uvw *= u_z / arma::norm(uvw, 2);

                // now we can add in the settling_velocity
                uvw(2) = -w;

                if (debug_output)
                    (*face)["u_z" + std::to_string(z)] = u_z;

                // negate as direction it's blowing instead of where it is from!!
                Vector_3 v3(-uvw(0), -uvw(1), uvw(2));
                if (debug_output)
                    face->set_face_vector("uvw" + std::to_string(z), v3);

                // holds wind velocity dot face normal
                double udotm[5];
                for (int j = 0; j < 5; ++j)
                {
                    udotm[j] = arma::dot(uvw, m[j]);
                }
                // lateral
                size_t idx = ntri * z + face->cell_local_id;

                //[idx][idx]
                size_t idx_idx_off = offset(row_buffer[idx], row_buffer[idx + 1], col_buffer, idx);
                b[idx] = 0;
                double V = face->get_area() * v_edge_height;

                // the sink term is added on for each edge check, which isn't right
                // and ends up double counting it so / by 5 for csubl and V so it's
                // not 5x counted.
                //   csubl /= 5.0;
                V /= 5.0;
                if (!do_sublimation)
                {
                    csubl = 0.0;
                }

                if (debug_output)
                    (*face)["csubl" + std::to_string(z)] = csubl;

                for (int f = 0; f < 3; f++)
                {
                    if (udotm[f] > 0)
                    {

                        if (d->face_neigh[f])
                        {
                            size_t nidx = ntri * z + face->neighbor(f)->cell_local_id;

                            elements[idx_idx_off] += V * csubl - d->A[f] * udotm[f] - alpha[f];

                            size_t idx_nidx_off = offset(row_buffer[idx], row_buffer[idx + 1], col_buffer, nidx);

                            elements[idx_nidx_off] += alpha[f];
                        }
                        else // missing neighbour case
                        {
                            // no mass in
//                            elements[ idx_idx_off ] += V*csubl-d->A[f]*udotm[f]-alpha[f];

                            // allow mass into the domain from ghost cell
                            elements[idx_idx_off] += -0.1e-1 * alpha[f] - 1. * d->A[f] * udotm[f] + csubl * V;
                        }
                    }
                    else
                    {
                        if (d->face_neigh[f])
                        {
                            size_t nidx = ntri * z + face->neighbor(f)->cell_local_id;

                            size_t idx_nidx_off = offset(row_buffer[idx], row_buffer[idx + 1], col_buffer, nidx);

                            elements[idx_idx_off] += V * csubl - alpha[f];
                            elements[idx_nidx_off] += -d->A[f] * udotm[f] + alpha[f];
                        }
                        else
                        {
                            // No mass in
//                            elements[ idx_idx_off ] +=  V*csubl-alpha[f];

                            // allow mass in
                            elements[idx_idx_off] += -0.1e-1 * alpha[f] - .99 * d->A[f] * udotm[f] + csubl * V;
                        }
                    }
                }

                // vertical layers
                if (z == 0)
                {

                    double alpha4 = d->A[4] * K[4] / (hs / 2.0 + v_edge_height / 2.0);

                    // bottom face, only turbulent diffusion
                    //              elements[idx_idx_off] += V * csubl - alpha4;

                    //            // includes advection term
                    elements[idx_idx_off] += V * csubl - d->A[4] * udotm[4] - alpha4;

                    b[idx] += -alpha4 * c_salt;

                    // ntri * (z + 1) + face->cell_local_id
                    size_t idx_nidx_off =
                        offset(row_buffer[idx], row_buffer[idx + 1], col_buffer, ntri * (z + 1) + face->cell_local_id);

                    if (udotm[3] > 0)
                    {
                        elements[idx_idx_off] += V * csubl - d->A[3] * udotm[3] - alpha[3];
                        elements[idx_nidx_off] += alpha[3];
                    }
                    else
                    {
                        elements[idx_idx_off] += V * csubl - alpha[3];
                        elements[idx_nidx_off] += -d->A[3] * udotm[3] + alpha[3];
                    }
                }
                else if (z == nLayer - 1) // top z layer
                {
                    //(kg/m^2/s)/(m/s)  ---->  kg/m^3
                    double cprecip = 0; //(*face)["p_snow"_s]/global_param->dt()/w;

                    // (*face)["p_snow"_s]=0;
                    // (*face)["p"_s]=0;

                    if (udotm[3] > 0)
                    {
                        elements[idx_idx_off] += V * csubl - d->A[3] * udotm[3] - alpha[3];
                        b[idx] += -alpha[3] * cprecip;
                    }
                    else
                    {
                        elements[idx_idx_off] += V * csubl - alpha[3];
                        b[idx] += d->A[3] * cprecip * udotm[3] - alpha[3] * cprecip;
                    }

                    // ntri * (z - 1) + face->cell_local_id
                    size_t idx_nidx_off =
                        offset(row_buffer[idx], row_buffer[idx + 1], col_buffer, ntri * (z - 1) + face->cell_local_id);
                    if (udotm[4] > 0)
                    {
                        elements[idx_idx_off] += V * csubl - d->A[4] * udotm[4] - alpha[4];
                        elements[idx_nidx_off] += alpha[4];
                    }
                    else
                    {
                        elements[idx_idx_off] += V * csubl - alpha[4];
                        elements[idx_nidx_off] += -d->A[4] * udotm[4] + alpha[4];
                    }
                }
                else // middle layers
                {
                    // ntri * (z + 1) + face->cell_local_id
                    size_t idx_nidx_off =
                        offset(row_buffer[idx], row_buffer[idx + 1], col_buffer, ntri * (z + 1) + face->cell_local_id);

                    if (udotm[3] > 0)
                    {
                        elements[idx_idx_off] += V * csubl - d->A[3] * udotm[3] - alpha[3];
                        elements[idx_nidx_off] += alpha[3];
                    }
                    else
                    {
                        elements[idx_idx_off] += V * csubl - alpha[3];
                        elements[idx_nidx_off] += -d->A[3] * udotm[3] + alpha[3];
                    }

                    idx_nidx_off =
                        offset(row_buffer[idx], row_buffer[idx + 1], col_buffer, ntri * (z - 1) + face->cell_local_id);
                    if (udotm[4] > 0)
                    {
                        elements[idx_idx_off] += V * csubl - d->A[4] * udotm[4] - alpha[4];
                        elements[idx_nidx_off] += alpha[4];
                    }
                    else
                    {
                        elements[idx_idx_off] += V * csubl - alpha[4];
                        elements[idx_nidx_off] += -d->A[4] * udotm[4] + alpha[4];
                    }
                }
            } // end z iter

        } // end face iter

    } // end pragma omp parallel thread pool

    // setup the compressed matrix on the compute device, if available
#ifdef VIENNACL_WITH_OPENCL
    viennacl::context gpu_ctx(viennacl::OPENCL_MEMORY);
    vl_C.switch_memory_context(gpu_ctx);
    b.switch_memory_context(gpu_ctx);
#endif

    // This solves the steady-state suspension layer concentration

    // configuration of preconditioner:
    viennacl::linalg::chow_patel_tag chow_patel_ilu_config;
    chow_patel_ilu_config.sweeps(3);       //  nonlinear sweeps
    chow_patel_ilu_config.jacobi_iters(2); //  Jacobi iterations per triangular 'solve' Rx=r
    viennacl::linalg::chow_patel_ilu_precond<viennacl::compressed_matrix<vcl_scalar_type>> chow_patel_ilu(
        vl_C, chow_patel_ilu_config);

    // Set up convergence tolerance to have an average value for each unknown
    double suspension_gmres_tol_per_unknown = 1e-8;
    double suspension_gmres_tol = suspension_gmres_tol_per_unknown * nLayer * ntri;
    // Set max iterations and maximum Krylov dimension before restart
    size_t suspension_gmres_max_iterations = 1000;
    size_t suspension_gmres_krylov_dimension = 30;

    // compute result and copy back to CPU device (if an accelerator was used),
    // otherwise access is slow
    viennacl::linalg::gmres_tag suspension_custom_gmres(suspension_gmres_tol, suspension_gmres_max_iterations,
                                                        suspension_gmres_krylov_dimension);
    viennacl::vector<vcl_scalar_type> vl_x = viennacl::linalg::solve(vl_C, b, suspension_custom_gmres, chow_patel_ilu);
    std::vector<vcl_scalar_type> x(vl_x.size());
    viennacl::copy(vl_x, x);

    // Log final state of the linear solve
//    LOG_DEBUG << "Suspension_GMRES # of iterations: " << suspension_custom_gmres.iters();
//    LOG_DEBUG << "Suspension_GMRES final residual : " << suspension_custom_gmres.error();

    /*
      Dump matrix to ASCII file
    */
    //    ofstream ofile;
    //    ofile.open("C.out");
    //    ofile << std::setprecision(12) << vl_C;
    //    ofile.close();
    //
    //    ofile.open("b.out");
    //    ofile << std::setprecision(12) << b;
    //    ofile.close();
    //
    //    ofile.open("x.out");
    //    ofile << std::setprecision(12) << vl_x;
    //    ofile.close();

    // Now we have the concentration, compute the suspension flux

#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {

        auto face = domain->face(i);
        auto d = face->get_module_data<data>(ID);
        double Qsusp = 0;
        double hs = d->hs;

        double Qsubl = 0;
        for (int z = 0; z < nLayer; ++z)
        {
            double c = x[ntri * z + face->cell_local_id];
            c = c < 0 || is_nan(c) ? 0 : c; // harden against some numerical issues that
                                            // occasionally come up for unknown reasons.

            double u_z = d->u_z_susp.at(z);

            Qsusp += c * u_z * v_edge_height; /// kg/m^3 ---->  kg/(m.s)

            if (debug_output)
            {
                (*face)["c" + std::to_string(z)] = c;

                // This is an approximation as it uses after transport concentrations.
                // However this will have already taken into account sublimation during the coupled transport phase
                // Eqn 20 Pomeroy 1993
                Qsubl += (*face)["csubl" + std::to_string(z)] * c *
                         v_edge_height; // kg/(m^2 *s) => per unit area of snowcover
            }
        }
        (*face)["Qsusp"_s] = Qsusp;
        if (debug_output)
        {
            (*face)["Qsubl"_s] = Qsubl;
            (*face)["Qsubl_mass"_s] = Qsubl * global_param->dt(); // kg/m^2 or mm
            d->sum_subl += (*face)["Qsubl_mass"_s];
            (*face)["sum_subl"_s] = d->sum_subl;
        }
    }

//    // Setup the matrix to be used to the solution of the gradient of the
//    // suspension flux this will give us our deposition flux
//    // get row buffer
//    unsigned int const* A_row_buffer =
//        viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(vl_A.handle1());
//    unsigned int const* A_col_buffer =
//        viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(vl_A.handle2());
//    vcl_scalar_type* A_elements =
//        viennacl::linalg::host_based::detail::extract_raw_pointer<vcl_scalar_type>(vl_A.handle());
//
//    // zero CSR vector in vl_A
//    viennacl::vector_base<vcl_scalar_type> init_temporaryA(
//        vl_A.handle(), viennacl::compressed_matrix<vcl_scalar_type>::size_type(nnz_drift + 1), 0, 1);
//    // write:
//    init_temporaryA = viennacl::zero_vector<vcl_scalar_type>(
//        viennacl::compressed_matrix<vcl_scalar_type>::size_type(nnz_drift + 1), viennacl::traits::context(vl_A));
//
//    // zero fill RHS for drift
//    bb.clear();

//#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->get_module_data<data>(ID);
        auto& m = d->m;

        double phi = (*face)["vw_dir"_s];
        Vector_2 v = -math::gis::bearing_to_cartesian(phi);

        // setup wind vector
        arma::vec uvw(3);
        uvw(0) = v.x(); // U_x
        uvw(1) = v.y(); // U_y
        uvw(2) = 0;

        double udotm[3];

        // edge lengths b/c 2d now
        double E[3] = {0, 0, 0};

        for (int j = 0; j < 3; ++j)
        {
            // just unit vectors as qsusp/qsalt flux has magnitude
            udotm[j] = arma::dot(uvw, m[j]);

            E[j] = face->edge_length(j);
        }

        double V = face->get_area();
        double qdep = 0;

        for (int j = 0; j < 3; j++)
        {
            double Qtj = 0;
            double Qsj = 0;

            if(udotm[j] > 0)
            {
                //get that upwind neighbour
                if (d->face_neigh[j])
                {
                    auto neigh = face->neighbor(j);
                    Qtj = (*neigh)["Qsusp"_s];
                    Qsj = (*neigh)["Qsalt"_s];
                }
                else
                {
                    //upwind neighbour doesn't exist, susp = 0
                    Qtj=Qsj=0;
                }
            }
            else
            {
                if (d->face_neigh[j])
                {

                    auto neigh = face->neighbor(j);
                    Qtj = (*neigh)["Qsusp"_s];
                    Qsj = (*neigh)["Qsalt"_s];

                }
                else
                {
                    //neighbour doesn't exist, susp = 0
                    Qtj=Qsj=0;
                }

            }

            qdep = qdep - (Qtj+Qsj)*udotm[j];
        }


        //take 1 fwd euler step

        double mass =  qdep/V * global_param->dt(); // kg/s*dt -> kg/m^2
        (*face)["drift_mass"_s] = mass;
        d->sum_drift += mass;

        (*face)["sum_drift"_s] = d->sum_drift;
    } // end face itr

// setup the compressed matrix on the compute device, if available
//#ifdef VIENNACL_WITH_OPENCL
//    //    viennacl::context gpu_ctx(viennacl::OPENCL_MEMORY);  <--- already
//    //    defined above
//    vl_A.switch_memory_context(gpu_ctx);
//    bb.switch_memory_context(gpu_ctx);
//#endif

//    // Solve the deposition flux --> how much drifting there is.
//
//    // configuration of preconditioner:
//    viennacl::linalg::chow_patel_tag deposition_flux_chow_patel_config;
//    deposition_flux_chow_patel_config.sweeps(3);       //  nonlinear sweeps
//    deposition_flux_chow_patel_config.jacobi_iters(2); //  Jacobi iterations per triangular 'solve' Rx=r
//    viennacl::linalg::chow_patel_icc_precond<viennacl::compressed_matrix<vcl_scalar_type>>
//        deposition_flux_chow_patel_icc(vl_A, deposition_flux_chow_patel_config);
//
//    // Set up convergence tolerance to have an average value for each unknown
//    double deposition_flux_cg_tol_per_unknown = 1e-7;
//    double deposition_flux_cg_tol = deposition_flux_cg_tol_per_unknown * ntri;
//    // Set max iterations and maximum Krylov dimension before restart
//    size_t deposition_flux_cg_max_iterations = 500;
//
//    // compute result and copy back to CPU device (if an accelerator was used),
//    // otherwise access is slow
//    viennacl::linalg::cg_tag deposition_flux_custom_cg(deposition_flux_cg_tol, deposition_flux_cg_max_iterations);
//
//    // compute result and copy back to CPU device (if an accelerator was used),
//    // otherwise access is slow
//    viennacl::vector<vcl_scalar_type> vl_dSdt =
//        viennacl::linalg::solve(vl_A, bb, deposition_flux_custom_cg, deposition_flux_chow_patel_icc);
//    // viennacl::vector<vcl_scalar_type> vl_dSdt = viennacl::linalg::solve(vl_A,
//    // bb, deposition_flux_custom_cg);
//    std::vector<vcl_scalar_type> dSdt(vl_dSdt.size());
//    viennacl::copy(vl_dSdt, dSdt);
//
//    // Log final state of the linear solve
//    LOG_DEBUG << "deposition_flux_CG # of iterations: " << deposition_flux_custom_cg.iters();
//    LOG_DEBUG << "deposition_flux_CG final residual : " << deposition_flux_custom_cg.error();
//
//    // take one FE integration step to get the total mass (SWE) that is eroded or
//    // deposited
//
//#pragma omp parallel for
//    for (size_t i = 0; i < domain->size_faces(); i++)
//    {
//
//        auto face = domain->face(i);
//        auto d = face->get_module_data<data>(ID);
//
//        double qdep = is_nan(dSdt[i]) ? 0 : dSdt[i];
//
//        double mass = 0;
//
//        mass = qdep * global_param->dt(); // kg/m^2*s *dt -> kg/m^2
//
//        (*face)["drift_mass"_s] = mass;
//        d->sum_drift += mass;
//
//        (*face)["sum_drift"_s] = d->sum_drift;
//    }
}

PBSM3D::~PBSM3D() {}
