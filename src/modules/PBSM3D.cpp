#include "PBSM3D.hpp"

//starting at row_start until row_end find offset for col
inline unsigned int offset(const unsigned int& row_start,const unsigned int& row_end,const unsigned int  * col_buffer, const unsigned int& col)
{
    for(unsigned int i=row_start; i < row_end; ++i)
    {
        if( col_buffer[i] == col)
            return i;
    }
    return -1; //wrap it and index garbage
}
PBSM3D::PBSM3D(config_file cfg)
        :module_base(parallel::domain)
{
    depends("U_2m_above_srf");
    depends("vw_dir");
    depends("swe");
    depends("t");
    depends("rh");


    depends("p_snow");
    depends("p");
//    provides("p");
//    provides("p_snow");

    optional("fetch");


    debug_output=cfg.get("debug_output",false);

//    provides("u10");
    provides("is_drifting");
//    provides("salt_limit");

    if(debug_output)
    {
        nLayer = cfg.get("nLayer", 5);
        for (int i = 0; i < nLayer; ++i)
        {
            provides("K" + std::to_string(i));
            provides("c" + std::to_string(i));
            provides("rm" + std::to_string(i));
            provides("csubl"+ std::to_string(i));
        }

        provides("Km_coeff");
        provides("Qsusp_pbsm");
        provides("inhibit_saltation");
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
        provides("LAI");

        provides("csalt");
        provides("c_salt_fetch_big");

        provides("u*_th");
        provides("u*_n");

        provides("dm/dt");
        provides("mm");
        provides("Qsubl");
    }

    provides("drift_mass"); //kg/m^2
    provides("Qsusp");
    provides("Qsalt");

    provides("sum_drift");
}

void PBSM3D::init(mesh domain)
{
    nLayer = cfg.get("nLayer",5);

    susp_depth = 5; //5m as per pomeroy
    v_edge_height = susp_depth / nLayer; //height of each vertical prism
    l__max = 40; // mixing length for diffusivity calculations

    settling_velocity = cfg.get("settling_velocity",-0.5); // m/s, Lehning, M., H. Löwe, M. Ryser, and N. Raderschall (2008), Inhomogeneous precipitation distribution and snow transport in steep terrain, Water Resour. Res., 44(7), 1–19, doi:10.1029/2007WR006545.

    if(settling_velocity > 0)
        BOOST_THROW_EXCEPTION(module_error() << errstr_info ("PBSM3D settling velocity must be negative"));


    do_sublimation = cfg.get("do_sublimation",true);
    do_lateral_diff = cfg.get("do_lateral_diff",true);
    eps = cfg.get("smooth_coeff",820);
    limit_mass= cfg.get("limit_mass",false);
    min_mass_for_trans = cfg.get("min_mass_for_trans",0);

    snow_diffusion_const = cfg.get("snow_diffusion_const",0.5); // Beta * K, this is beta and scales the eddy diffusivity
    rouault_diffusion_coeff = cfg.get("rouault_diffusion_coef",false);

    enable_veg = cfg.get("enable_veg",true);

    if(rouault_diffusion_coeff)
    {
        LOG_WARNING << "rouault_diffusion_coef overrides const snow_diffusion_const values.";
    }

    n_non_edge_tri = 0;


    //use this to build the sparsity pattern for suspension matrix
    size_t ntri = domain->number_of_faces();
    std::vector< std::map< unsigned int, vcl_scalar_type> > C(ntri * nLayer);

    //sparsity pattern for drift
    std::vector< std::map< unsigned int, vcl_scalar_type> > A(ntri);

#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<data>(ID);
        d->Tsguess = 273.0;
        d->z0Fnguess = 0.01;

        if(!face->has_vegetation() && enable_veg)
        {
            LOG_ERROR << "Vegetation is enabled, but no vegetation parameter was found.";
        }
        if(face->has_vegetation() && enable_veg)
        {
            d->CanopyHeight = face->veg_attribute("CanopyHeight");
            d->LAI = face->veg_attribute("LAI");
        } else{
            d->CanopyHeight = 0;
            d->LAI = 0;
            enable_veg = false;
        }

        auto& m = d->m;
        //edge unit normals
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

        //top
        m[3].set_size(3);
        m[3].fill(0);
        m[3](2) = 1;

        //bottom
        m[4].set_size(3);
        m[4].fill(0);
        m[4](2) = -1;

        //face areas
        for (int j = 0; j < 3; ++j)
            d->A[j] = face->edge_length(j) * v_edge_height;

        //top, bottom
        d->A[3] = d->A[4] = face->get_area();

        d->is_edge=false;
        //which faces have neighbours? Ie, are we an edge?
        for (int a = 0; a < 3; ++a)
        {
            A[i][i] = -9999;

            auto neigh = face->neighbor(a);
            if (neigh == nullptr)
            {
                d->face_neigh[a] = false;
                d->is_edge = true;
            }
            else
            {
                d->face_neigh[a] = true;

                A[i][neigh->cell_id] = -9999;
            }

        }
        if(!d->is_edge)
        {
            d->cell_id=n_non_edge_tri;
            ++n_non_edge_tri;
        }

        d->sum_drift = 0;



        // iterate over the vertical layers
        for (int z = 0; z < nLayer; ++z)
        {
            size_t idx = ntri * z + face->cell_id;
            for (int f = 0; f < 3; f++)
            {
                if (d->face_neigh[f])
                {
                    size_t nidx = ntri * z + face->neighbor(f)->cell_id;
                    C[idx][idx] = -9999;
                    C[idx][nidx] = -9999;
                } else
                {
                    C[idx][idx] = -9999;
                }
            }

            if (z == 0)
            {
                C[idx][idx] = -9999;
                C[idx][ntri * (z + 1) + face->cell_id] = -9999;
            } else if (z == nLayer - 1)
            {
                C[idx][idx] = -9999;
                C[idx][ntri * (z - 1) + face->cell_id] = -9999;

            }
            else //middle layers
            {
                C[idx][idx] = -9999;
                C[idx][ntri * (z + 1) + face->cell_id] = -9999;
                C[idx][ntri * (z - 1) + face->cell_id] = -9999;
            }
        }
    }


    
    viennacl::copy(C,vl_C); // copy C -> vl_C, sets up the sparsity pattern
    viennacl::copy(A,vl_A); // copy A -> vl_A, sets up the sparsity pattern

    b.resize(ntri * nLayer);
    bb.resize(ntri);
    nnz = vl_C.nnz();
    nnz_drift = vl_A.nnz();

}

void PBSM3D::run(mesh domain)
{
    //needed for linear system offsets
    size_t ntri = domain->number_of_faces();

    //vcl_scalar_type is defined in the main CMakeLists.txt file.
    // Some GPUs do not have double precision so the run will fail if the wrong precision is used

#ifdef VIENNACL_WITH_OPENCL
    viennacl::context host_ctx(viennacl::MAIN_MEMORY);
    vl_C.switch_memory_context(host_ctx);
    vl_A.switch_memory_context(host_ctx);

    b.switch_memory_context(host_ctx);
    bb.switch_memory_context(host_ctx);
#endif


    //zero CSR vector in vl_C
    viennacl::vector_base<unsigned int> init_temporary(vl_C.handle(), viennacl::compressed_matrix<vcl_scalar_type>::size_type(nnz+1), 0, 1);
    // write:
    init_temporary = viennacl::zero_vector<unsigned int>(viennacl::compressed_matrix<vcl_scalar_type>::size_type(nnz+1), viennacl::traits::context(vl_C));

    //zero-fill RHS
    b.clear();

    //ice density
    double rho_p = PhysConst::rho_ice;

    //get row buffer
    unsigned int const* row_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(vl_C.handle1());
    unsigned int const* col_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(vl_C.handle2());
    vcl_scalar_type*    elements   = viennacl::linalg::host_based::detail::extract_raw_pointer<vcl_scalar_type>(vl_C.handle());

#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto id = face->cell_id;
        auto d = face->get_module_data<data>(ID);
        auto& m = d->m;

        double fetch = 1000;
        if(has_optional("fetch"))
            fetch = face->face_data("fetch");

        //get wind from the face
        double phi = face->face_data("vw_dir");
        double u2 = face->face_data("U_2m_above_srf");

        double swe = face->face_data("swe"); // mm   -->    kg/m^2
        swe = is_nan(swe) ? 0 : swe; // handle the first timestep where swe won't have been updated if we override the module order

        double snow_depth = face->face_data("snowdepthavg");

        snow_depth = is_nan(snow_depth) ? 0: snow_depth;


        double u10 = Atmosphere::log_scale_wind(u2, 2, 10, 0);
        // face->set_face_data("u10",u10);

        Vector_2 vwind = -math::gis::bearing_to_cartesian(phi);

        //setup wind vector
        arma::vec uvw(3);
        uvw(0) = vwind.x(); //U_x
        uvw(1) = vwind.y(); //U_y
        uvw(2) = 0;


        // -------
        // Calculate the new value of z0 to take into account partially filled vegetation and the momentum sink

        // Li and Pomeroy 2000, eqn 5. But follow the LAI/2.0 suggestion from Raupach 1994 (DOI:10.1007/BF00709229) Section 3(a)

        double height_diff = std::max(0.0,d->CanopyHeight - snow_depth);//don't allow negative differences in height
        if(debug_output) face->set_face_data("height_diff",height_diff);
        double ustar = 1.3; //placeholder
        double lambda = 0.5 * d->LAI / d->CanopyHeight * (height_diff);
        if(debug_output) face->set_face_data("LAI",d->LAI);
        if(debug_output) face->set_face_data("lambda",lambda);

        d->saltation = true;

        // The strategy here is as follows:
        // 1) Wait for vegetation to fill up until it is within 30 cm of the top
        //     then we can apply Pomeroy and Li 2000 eqn 4 to calculate a z0. This param doesn't seem to work for exposed veg > 30 cm
        //     i.e, what you'd expect for crop stubble where it was derived.
        // 2) Within 30 cm of the top, use P&l2000 eqn 4 to calculate a z0, and effectively allow wind to blow snow out of the vegetation
        // 3) Once lambda is close to 0, just solve eqn 4 directly via lambert fn without the veg sink as this is faster than the iter sol'n



        //if lambda is -> 0, we can solve for a ustar and the z0 value directly without calculating an iterative sol'n
        // or if we have disabled veg, just use the non-veg param
        // or if we are still filling it up, calculate a 'sane' zo/ustar and just disable saltation in this triangle.
        if (height_diff < 0.05  ||  height_diff > 0.3 || !enable_veg)
        {
            ustar = -.2000000000 * u2 / gsl_sf_lambert_Wm1(-0.1107384167e-1 * u2);
            d->z0 = std::max(0.001, 0.1203 * ustar * ustar / (2.0 * 9.81));
            lambda = 0;

            if( height_diff > 0.3 && enable_veg)
            {
                //just disable saltation
                d->saltation = false;
            }
        }
        else  /*if( height_diff < 0.3)*/
        {
            auto z0Fn = [&](double z0) -> double {
                //This formulation has the following coeffs built in
                // c_2 = 1.6;
                // c_3 = 0.07519;
                // c_4 - 0.5;
                return (1.6 * 0.7519e-1) * pow(.41 * u2 / log(2.0 / z0), 2.0) / (2.0 * 9.81) +  lambda - z0;  //0.5
            };

            double min = 0.0001;
            double max = 1;
            boost::uintmax_t max_iter = 500;

            auto tol = [](double a, double b) -> bool {
                return fabs(a - b) < 1e-8;
            };

            try {
//                auto r = boost::math::tools::toms748_solve(z0Fn, min, max, tol, max_iter);
                auto r = boost::math::tools::bracket_and_solve_root(z0Fn,d->z0Fnguess,2.0,false,tol,max_iter);
                d->z0 = r.first + (r.second - r.first) / 2;
                ustar = u2 * PhysConst::kappa / log(2.0 / d->z0);
            }
            catch (...) {
                // for cases of large LAI, the above will not converge within the min/max given
                // this will trap that error, and sets that cell to have no saltation,
                // however re calculate the z0 and ustar as if there was no vegetation.
                // TODO: something with zeroplane displacement should replace this, but for now this works as it limits saltation
                // in ares of trees.
                d->saltation = false;

                ustar = -.2000000000 * u2 / gsl_sf_lambert_Wm1(-0.1107384167e-1 * u2);
                d->z0 = std::max(0.001, 0.1203 * ustar * ustar / (2.0 * 9.81));
            }
        }

        d->z0 = std::max(0.001,d->z0);
        ustar = std::max(0.01,ustar);

        if(debug_output) face->set_face_data("z0",d->z0);
        // ------------------

        //depth of saltation layer
        double hs = 0;

        //lehning
//        hs = d->z0 + 2.4025 * pow(u2, 2.) * pow(PhysConst::kappa, 2.) * pow(cos(25. * M_PI / 180.), 2.) /
//                         (pow(log(2. / d->z0), 2.) * 9.81);

        //pomeroy
        hs = 0.08436*pow(ustar,1.27);
        d->hs = hs;
        if(debug_output) face->set_face_data("hs",hs);

        //eddy diffusivity (m^2/s)
        // 0,1,2 will all be K = 0, as no horizontal diffusion process
        double K[5] = {0, 0, 0, 0, 0};

        // holds A_ * K_ / h_
        // _0 -> _2 are the horizontal sides
        // _3 -> is the top of the prism
        // _4 -> is the bottom the prism
        double alpha[5] = {0, 0, 0, 0, 0};


        double t = face->face_data("t")+273.15;

        double rho_f =  mio::Atmosphere::stdDryAirDensity(face->get_z(),t); //kg/m^3, comment in mio is wrong.1.225; // air density, fix for T dependency

        //threshold friction velocity, paragraph below eqn 3 in Saltation of Snow, Pomeroy 1990
        //double u_star_t = pow(tau_t_f/rho_f,0.5);

        //Pomeroy and Li, 2000
        // Eqn 7
        double T = face->face_data("t");
        double u_star_saltation = 0.35+(1.0/150.0)*T+(1.0/8200.0)*T*T;


        double Qsalt = 0;
        double c_salt = 0;


        double c_salt_fetch_big=0;
        if(debug_output) face->set_face_data("u*_th",u_star_saltation);

        if(debug_output) face->set_face_data("ustar",ustar);

        if(debug_output) face->set_face_data("is_drifting",0);
        if(debug_output) face->set_face_data("Qsusp_pbsm",0); //for santiy checks against pbsm

        if(!d->saltation && debug_output)
        {
             face->set_face_data("inhibit_saltation",1);
        }

        if( ustar > u_star_saltation &&
               swe > min_mass_for_trans &&
                d->saltation)
        {

            double pbsm_qsusp = pow(u10,4.13)/674100.0;
            if(debug_output) face->set_face_data("Qsusp_pbsm",pbsm_qsusp);
            if(debug_output) face->set_face_data("is_drifting",1);

            //Pomeroy and Li 2000, eqn 8
            double Beta =  202.0; //170.0;
            double ustar_n = 0;

//            if (snow_depth < d->CanopyHeight && enable_veg)
//            {
                double m = 0.16;
                //this is ustar_n / ustar
                ustar_n = sqrt(m*Beta*lambda)*pow(1.0+m* Beta*lambda,-0.5); //MacDonald 2009 eq 3;
//            }

            if(debug_output) face->set_face_data("u*_n",ustar_n);
            //Pomeroy 1992, eqn 12, see note above for ustar_n calc
            c_salt = rho_f / (3.29 * ustar) * (1.0 - ustar_n*ustar_n -  (u_star_saltation*u_star_saltation) / (ustar * ustar));

            // occasionally happens to happen at low wind speeds where the parameterization breaks.
            // hasn't happened since changed this to the threshold velocity though...
            if(c_salt < 0 || std::isnan(c_salt))
            {
                c_salt = 0;
            }

            c_salt_fetch_big = c_salt;

            //use the exp decay of Liston, eq 10
            //95% of max saltation occurs at fetch = 500m
            //Liston, G., & Sturm, M. (1998). A snow-transport model for complex terrain. Journal of Glaciology.
            if(fetch < 500)
            {
                double fetch_ref = 500;
                double mu = 3.0;
                c_salt *= 1.0-exp(-mu * fetch/fetch_ref);
            }
            //mean wind speed in the saltation layer
            double uhs = std::max(0.1,Atmosphere::log_scale_wind(u2, 2, hs, 0,d->z0)/2.);

            // kg/(m*s)
            Qsalt =  c_salt * uhs * hs; //integrate over the depth of the saltation layer, kg/(m*s)
            // face->set_face_data("saltation_mass",c_salt*hs * face->get_area());

            //calculate the surface integral of Qsalt, and ensure we aren't saltating more mass than what exists
            //in the triangle. I
//            double salt=0;
//            double udotm[3];
//            double A = face->get_area();
//            double E[3];
//            for(int j=0; j<3; ++j)
//            {
//                E[j]=face->edge_length(j);
//                udotm[j] = arma::dot(uvw, m[j]);
//                salt += -.5000000000*E[j]*Qsalt*udotm[j]/A;
//            }
//
//            //https://www.wolframalpha.com/input/?i=(m*(kg%2Fm%5E3*(m%2Fs)*m)%2Fm%5E2)*s
//            salt = std::fabs(salt) *  global_param->dt(); // ->kg/m^2
//
            // face->set_face_data("salt_limit",salt);
//
//            //Figure out the max we can transport, ie total mass in cell
//            if( limit_mass && salt > swe) // could we move more than the total mass in the cell during this timestep?
//            {
//
//                //back out what the max conc should be based on our swe
//                //units: ((kg/m^2)*m^2)/( s*m*(m/s)*m ) -> kg/m^3
//                c_salt = -2.*swe*A/(uhs*hs*(E[0]*udotm[0]+E[1]*udotm[1]+E[2]*udotm[2]));
//
//                Qsalt = c_salt * uhs * hs;
//                if(is_nan(c_salt)) //shoouldn't happen but....
//                {
//                    c_salt = 0;
//                    Qsalt = 0;
//                }
//
//
////                LOG_DEBUG << "More saltation than snow, limiting conc to " << c_salt << " triangle="<<i;
////                LOG_DEBUG << "Avail mass = " << swe << ", would have salted =  " << salt;
//            }
        }

        // can use for point scale plume testing.
//
//        if (i != 7105 )
//        {
//            Qsalt = 0;
//            c_salt = 0;
//        }


        if(debug_output) face->set_face_data("csalt", c_salt);
        if(debug_output) face->set_face_data("c_salt_fetch_big", c_salt_fetch_big);
        face->set_face_data("Qsalt", Qsalt);

        double rh = face->face_data("rh")/100.;
        double RH = rh*100.;
        double es = mio::Atmosphere::saturatedVapourPressure(t);
        double ea = rh * es / 1000.; // ea needs to be in kpa
 
        double v = 1.88*10e-5; //kinematic viscosity of air, below eqn 13 Pomeroy 1993

        //vapour pressure, Pa
        auto e = [](double T,double RH)
        {
            return RH/100.*611.0*exp(17.3*T/(237.3+T)); //Pa
        };

        // iterate over the vertical layers
        for (int z = 0; z < nLayer; ++z)
        {
            //height in the suspension layer
            double cz = z + hs+ v_edge_height/2.; //cell center height
            //compute new U_z at this height in the suspension layer
            double u_z = std::max(0.1, Atmosphere::log_scale_wind(u2, 2, cz, 0,d->z0));

            //calculate dm/dt from
            // equation 13 from Pomeroy and Li 2000
            // To do so, use equations 12 - 16 in Pomeroy et al 2000
            //Pomeroy, J. W., and L. Li (2000), Prairie and arctic areal snow cover mass balance using a blowing snow model, J. Geophys. Res., 105(D21), 26619–26634, doi:10.1029/2000JD900149. [online] Available from: http://www.agu.org/pubs/crossref/2000/2000JD900149.shtml

            //these are from
            // Pomeroy, J. W., D. M. Gray, and P. G. Landine (1993), The prairie blowing snow model: characteristics, validation, operation, J. Hydrol., 144(1–4), 165–192.
            double rm = 4.6e-5 * pow(cz,-0.258); // eqn 18, mean particle size
            if(debug_output) face->set_face_data("rm"+std::to_string(z), rm);

            double xrz = 0.005 * pow(u_z,1.36);  //eqn 16
            double omega = 1.1*10e7 * pow(rm,1.8); //eqn 15
            double Vr = omega + 3.0*xrz*cos(M_PI/4.0); //eqn 14

            double Re = 2.0*rm*Vr / v; //eqn 12

            double Nu, Sh;
            Nu = Sh = 1.79 + 0.606 * pow(Re,0.5);

            //define above, T is in C, t is in K

            // (A.6)
            double D = 2.06 * 10e-5 * pow(t/273.15,1.75); //diffusivity of water vapour in air, t in K, eqn A-7 in Liston 1998 or Harder 2013 A.6

            // (A.9)
            double lambda_t = 0.000063 * t + 0.00673; // J/(kmol K) thermal conductivity, user Harder 2013 A.9, Pomeroy's is off by an order of magnitude, this matches this https://www.engineeringtoolbox.com/air-properties-d_156.html

            // (A.10) (A.11)
            double L = 1000.0 * (2834.1 - 0.29 *T - 0.004*T*T); // T in C

            /*
             * The *1000 and /1000 are important unit conversions. Doesn't quite match the harder paper, but Phil assures me it is correct.
             */
            double mw = 0.01801528 * 1000.0; //[kg/mol]  ---> g/mol
            double R = 8.31441/1000.0; // [J mol-1 K-1]

            double rho = (mw * ea) / (R*t);

            //use Harder 2013 (A.5) Formulation, but Pa formulation for e
            auto fx = [=](double Ti)
            {
                return boost::math::make_tuple(
                        T+D*L*(rho/(1000.0)-.611*mw*exp(17.3*Ti/(237.3+Ti))/(R*(Ti+273.15)*(1000.0)))/lambda_t-Ti,
                        D*L*(-0.6110000000e-3*mw*(17.3/(237.3+Ti)-17.3*Ti/pow(237.3+Ti,2))*exp(17.3*Ti/(237.3+Ti))/(R*(Ti+273.15))+0.6110000000e-3*mw*exp(17.3*Ti/(237.3+Ti))/(R*pow(Ti+273.15,2)))/lambda_t-1);
            };


            double guess = T;
            double min = -50;
            double max = 0;
            int digits = 6;

            double Ti = boost::math::tools::newton_raphson_iterate(fx, guess, min, max, digits);
//            if(debug_output) face->set_face_data("Ti",Ti);

            double Ts = Ti+273.15; //dmdtz expects in K

            //now use equation 13 with our solved Ts to compute dm/dt(z)
            double dmdtz = 2.0 * M_PI * rm * lambda_t / L * Nu * (Ts - (t+273.15));  //eqn 13 in Pomeroy and Li 2000
            if(debug_output) face->set_face_data("dm/dt",dmdtz);

            //calculate mean mass, eqn 23, 24 in Pomeroy 1993 (PBSM)
            double mm_alpha = 4.08 + 12.6*cz; //24
            double mm = 4./3. * M_PI * rho_p * rm*rm*rm *(1.0 + 3.0/mm_alpha + 2./(mm_alpha*mm_alpha)); //mean mass, eqn 23
            if(debug_output) face->set_face_data("mm",mm);
            double csubl = dmdtz/mm; //EQN 21 POMEROY 1993 (PBSM)
            if(debug_output) face->set_face_data("csubl"+std::to_string(z),dmdtz);
            //compute alpha and K for edges
            if(do_lateral_diff)
            {
                for (int a = 0; a < 3; ++a)
                {
                    auto neigh = face->neighbor(a);
                    alpha[a] = d->A[a];

//                    //if we have a neighbour, use the distance
//                    if (neigh != nullptr)
//                    {
//                        alpha[a] /= math::gis::distance(face->center(), neigh->center());
//
//                    } else
//                    {
//                        //otherwise assume 2x the distance from the center of face to one of it's vertexes, so ghost triangle is approx the same size
//                        alpha[a] /= 5.0 * math::gis::distance(face->center(), face->vertex(0)->point());
//
//                    }

                    //do just very low horz diffusion for numerics
                     K[a] =  0.00001;    //PhysConst::kappa * cz * ustar;// std::max(ustar * l, PhysConst::kappa * 2. * ustar);
                    alpha[a] *= K[a];
                }
            }
            //Li and Pomeroy 2000
            double l = PhysConst::kappa * (cz + d->z0) / ( 1.0  + PhysConst::kappa * (cz+d->z0)/ l__max);
            if(debug_output) face->set_face_data("l",l);
            double w = 1.1*10e7*pow(rm,1.8);
            if(debug_output) face->set_face_data("w",w);

            double diffusion_coeff = snow_diffusion_const; //snow_diffusion_const is a shared param so need a seperate copy here we can overwrite
            if (rouault_diffusion_coeff)
            {
                double c2 = 1.0;
                double dc = 1.0 / (1.0 + (c2 * w * w) / (1.56 * ustar * ustar));
                diffusion_coeff = dc;     //nope, snow_diffusion_const is shared, use a new
            }
            if(debug_output) face->set_face_data("Km_coeff", diffusion_coeff);

             //snow_diffusion_const is pretty much a calibration constant. At 1 it seems to over predict transports.
             K[3] = K[4] = diffusion_coeff * ustar * l;

//            K[3] = K[4] = snow_diffusion_const * PhysConst::kappa * cz * ustar;// std::max(ustar * l, PhysConst::kappa * cz * ustar);

            if(debug_output) face->set_face_data("K"+std::to_string(z), K[3] );
            //top
            alpha[3] = d->A[3] * K[3] / v_edge_height;
            //bottom
            alpha[4] = d->A[4] * K[4] / v_edge_height;


            double length = arma::norm(uvw, 2);
            double scale = u_z / length;

            uvw *= scale;
            uvw(2) = -w; //settling_velocity;

            //holds wind dot face normal
            double udotm[5];
            for (int j = 0; j < 5; ++j)
            {
                udotm[j] = arma::dot(uvw, m[j]);
            }
            //lateral
            size_t idx = ntri*z + face->cell_id;

            //[idx][idx]
            size_t idx_idx_off = offset(row_buffer[idx],
                                row_buffer[idx+1],
                                col_buffer, idx);
            b[idx] = 0;
            double V = face->get_area()  * v_edge_height;

            //the sink term is added on for each edge check, which isn't right and ends up double counting it
            //so / by 5 for csubl and V so it's not 5x counted.
            csubl /= 5.0;
            V/=5.0;
            if(!do_sublimation)
            {
                csubl = 0.0;
            }

            for(int f = 0; f < 3; f++)
            {
                if(udotm[f] > 0)
                {

                    if (d->face_neigh[f])
                    {
                        size_t nidx = ntri * z + face->neighbor(f)->cell_id;

                        elements[ idx_idx_off ] += V*csubl-d->A[f]*udotm[f]-alpha[f];

                        size_t idx_nidx_off = offset(row_buffer[idx],
                                            row_buffer[idx+1],
                                            col_buffer, nidx);

                        elements[ idx_nidx_off ] += alpha[f];

                    }
                    else
                    {
                        //no mass in
//                       vl_C(idx,idx) += V*csubl-d->A[f]*udotm[f]-alpha[f];

                        //allow mass in
//                        vl_C(idx,idx) += V*csubl-d->A[f]*udotm[f];


                        elements[ idx_idx_off ] += -0.1e-1*alpha[f]-1.*d->A[f]*udotm[f]+csubl*V;

                    }
                } else
                {
                    if (d->face_neigh[f])
                    {
                        size_t nidx = ntri * z + face->neighbor(f)->cell_id;

                        size_t idx_nidx_off = offset(row_buffer[idx],
                                                     row_buffer[idx+1],
                                                     col_buffer, nidx);

                        elements[ idx_idx_off ] +=V*csubl-alpha[f];
                        elements[ idx_nidx_off ]  += -d->A[f]*udotm[f]+alpha[f];
                    }
                    else
                    {
                        //No mass in
//                        vl_C(idx,idx) += V*csubl-alpha[f];

                        //allow mass in
//                        vl_C(idx,idx) += -.99*d->A[f]*udotm[f]+csubl*V;


                        elements[ idx_idx_off ] += -0.1e-1*alpha[f]-.99*d->A[f]*udotm[f]+csubl*V;

                    }
                }
            }


            //vertical layers
            if (z == 0)
            {

                double alpha4 = d->A[4] * K[4] / (hs/2.0 + v_edge_height/2.0);

                //bottom face, no advection
//                vl_C(idx,idx) += V*csubl-alpha4;
//                b[idx] += -alpha4*c_salt;

                //includes advection term
                elements[ idx_idx_off ] += V*csubl-d->A[4]*udotm[4]-alpha4;

                b[idx] += -alpha4*c_salt;

                //ntri * (z + 1) + face->cell_id
                size_t idx_nidx_off = offset(row_buffer[idx],
                                             row_buffer[idx+1],
                                             col_buffer, ntri * (z + 1) + face->cell_id);

                if (udotm[3] > 0)
                {
                    elements[ idx_idx_off ] += V*csubl-d->A[3]*udotm[3]-alpha[3];
                    elements[ idx_nidx_off ] += alpha[3];
                } else
                {
                    elements[ idx_idx_off ] += V*csubl-alpha[3];
                    elements[ idx_nidx_off ] += -d->A[3]*udotm[3]+alpha[3];
                }
            } else if (z == nLayer - 1)// top z layer
            {
                //(kg/m^2/s)/(m/s)  ---->  kg/m^3
                double cprecip = 0;//face->face_data("p_snow")/global_param->dt()/w;

                // face->set_face_data("p_snow",0);
                // face->set_face_data("p",0);

                if (udotm[3] > 0)
                {
                    elements[ idx_idx_off ] += V*csubl-d->A[3]*udotm[3]-alpha[3];
                    b[idx] += -alpha[3] * cprecip;
                } else
                {
                    elements[ idx_idx_off ] += V*csubl-alpha[3];
                    b[idx] += d->A[3]*cprecip*udotm[3] - alpha[3] * cprecip;
                }

                //ntri * (z - 1) + face->cell_id
                size_t idx_nidx_off = offset(row_buffer[idx],
                                             row_buffer[idx+1],
                                             col_buffer, ntri * (z - 1) + face->cell_id);
                if (udotm[4] > 0)
                {
                    elements[ idx_idx_off ] += V*csubl-d->A[4]*udotm[4]-alpha[4];
                    elements[ idx_nidx_off ] += alpha[4];

                } else
                {
                    elements[ idx_idx_off ] += V*csubl-alpha[4];
                    elements[ idx_nidx_off ] += -d->A[4]*udotm[4]+alpha[4];
                }


            } else //middle layers
            {
                //ntri * (z + 1) + face->cell_id
                size_t idx_nidx_off = offset(row_buffer[idx],
                                             row_buffer[idx+1],
                                             col_buffer, ntri * (z + 1) + face->cell_id);

                if (udotm[3] > 0)
                {
                    elements[ idx_idx_off ] += V*csubl-d->A[3]*udotm[3]-alpha[3];
                    elements[ idx_nidx_off ] += alpha[3];
                } else
                {
                    elements[ idx_idx_off ] += V*csubl-alpha[3];
                    elements[ idx_nidx_off ] += -d->A[3]*udotm[3]+alpha[3];
                }

                idx_nidx_off = offset(row_buffer[idx],
                                      row_buffer[idx+1],
                                      col_buffer, ntri * (z - 1) + face->cell_id);
                if (udotm[4] > 0)
                {
                    elements[ idx_idx_off ] += V*csubl-d->A[4]*udotm[4]-alpha[4];
                    elements[ idx_nidx_off ]  += alpha[4];
                } else
                {
                    elements[ idx_idx_off ] += V*csubl-alpha[4];
                    elements[ idx_nidx_off ]  += -d->A[4]*udotm[4]+alpha[4];
                }
            }
        } // end z iter
    } //end face iter

    //setup the compressed matrix on the compute device, if available
#ifdef VIENNACL_WITH_OPENCL
    viennacl::context gpu_ctx(viennacl::OPENCL_MEMORY);
    vl_C.switch_memory_context(gpu_ctx);
    b.switch_memory_context(gpu_ctx);
#endif
 
    // configuration of preconditioner:
    viennacl::linalg::chow_patel_tag chow_patel_ilu_config;
    chow_patel_ilu_config.sweeps(3);       //  nonlinear sweeps
    chow_patel_ilu_config.jacobi_iters(2); //  Jacobi iterations per triangular 'solve' Rx=r
    viennacl::linalg::chow_patel_ilu_precond< viennacl::compressed_matrix<vcl_scalar_type> > chow_patel_ilu(vl_C, chow_patel_ilu_config);


    //compute result and copy back to CPU device (if an accelerator was used), otherwise access is slow
    viennacl::linalg::gmres_tag gmres_tag(1e-8, 500, 30);
    viennacl::vector<vcl_scalar_type> vl_x = viennacl::linalg::solve(vl_C, b, gmres_tag, chow_patel_ilu);
    std::vector<vcl_scalar_type> x(vl_x.size());
    viennacl::copy(vl_x,x);


#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->get_module_data<data>(ID);
        double Qsusp = 0;
        double hs = d->hs;

        double u2 = face->face_data("U_2m_above_srf");

        double total_mass = 0;

        double Qsubl=0;
        for (int z = 0; z<nLayer;++z)
        {
            double c = x[ntri * z + face->cell_id];
            c = c < 0 || is_nan(c) ? 0 : c; //harden against some numerical issues that occasionally come up for unknown reasons.

            double cz = z + hs+ v_edge_height/2.; //cell center height

            double u_z =std::max(0.1, Atmosphere::log_scale_wind(u2, 2, cz, 0,d->z0)); //compute new U_z at this height in the suspension layer
            Qsusp += c * u_z * v_edge_height; /// kg/m^3 ---->  kg/(m.s)

            total_mass += c * v_edge_height * face->get_area();

            if(debug_output) face->set_face_data("c"+std::to_string(z),c);

            if(debug_output) Qsubl+=face->face_data("csubl"+std::to_string(z))* c*v_edge_height;

        }
         face->set_face_data("Qsusp",Qsusp);
        if(debug_output) face->set_face_data("Qsubl",Qsubl);
        // face->set_face_data("suspension_mass",total_mass);


    }

    //get row buffer
    unsigned int const* A_row_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(vl_A.handle1());
    unsigned int const* A_col_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(vl_A.handle2());
    vcl_scalar_type*    A_elements   = viennacl::linalg::host_based::detail::extract_raw_pointer<vcl_scalar_type>(vl_A.handle());

    //zero CSR vector in vl_A
    viennacl::vector_base<unsigned int> init_temporaryA(vl_A.handle(), viennacl::compressed_matrix<vcl_scalar_type>::size_type(nnz_drift+1), 0, 1);
    // write:
    init_temporaryA = viennacl::zero_vector<unsigned int>(viennacl::compressed_matrix<vcl_scalar_type>::size_type(nnz_drift+1), viennacl::traits::context(vl_A));

    //zero fill RHS for drift
    bb.clear();

#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->get_module_data<data>(ID);
        auto &m = d->m;

        double phi = face->face_data("vw_dir");
        Vector_2 v = -math::gis::bearing_to_cartesian(phi);

        //setup wind vector
        arma::vec uvw(3);
        uvw(0) = v.x(); //U_x
        uvw(1) = v.y(); //U_y
        uvw(2) = 0;

        double udotm[3];

        //edge lengths b/c 2d now
        double E[3] = {0, 0, 0};

        for (int j = 0; j < 3; ++j)
        {
            //just unit vectors as qsusp/qsalt flux has magnitude
            udotm[j] = arma::dot(uvw, m[j]);

            E[j] = face->edge_length(j);
        }

        double dx[3] = {2.0, 2.0, 2.0};

        double V = face->get_area();


        //[i][i]
        size_t i_i_off = offset(A_row_buffer[i],
                                A_row_buffer[i+1],
                                A_col_buffer, i);

        for (int j = 0; j < 3; j++)
        {
            if (d->face_neigh[j])
            {
                auto neigh = face->neighbor(j);
                auto Qtj = neigh->face_data("Qsusp") + face->face_data("Qsusp");
                auto Qsj = neigh->face_data("Qsalt") + face->face_data("Qsalt");
                double Qt = Qtj / 2.0 + Qsj / 2.0;

                dx[j] = math::gis::distance(face->center(), neigh->center());

                A_elements[i_i_off] += V + eps * E[j] / dx[j];

                //[i][neigh->cell_id]
                size_t i_ni_off = offset(A_row_buffer[i],
                                        A_row_buffer[i+1],
                                        A_col_buffer, neigh->cell_id);
                A_elements[i_ni_off] += - eps * E[j] / dx[j];
                bb[i] += - E[j] * Qt * udotm[j];

            } else
            {
                auto Qtj = 2. * face->face_data("Qsusp"); // const flux across, 0 -> drifts!!!
                auto Qsj = 2. * face->face_data("Qsalt");
                double Qt = Qtj / 2.0 + Qsj / 2.0;
                A_elements[i_i_off] += V;
                bb[i] += -E[j] * Qt * udotm[j];
            }
        }



    } // end face itr

//setup the compressed matrix on the compute device, if available
#ifdef VIENNACL_WITH_OPENCL
//    viennacl::context gpu_ctx(viennacl::OPENCL_MEMORY);  <--- already defined above
    vl_A.switch_memory_context(gpu_ctx);
    bb.switch_memory_context(gpu_ctx);
#endif
    // configuration of preconditioner:
    viennacl::linalg::chow_patel_ilu_precond< viennacl::compressed_matrix<vcl_scalar_type> > chow_patel_ilu2(vl_A, chow_patel_ilu_config);

    //compute result and copy back to CPU device (if an accelerator was used), otherwise access is slow
    viennacl::vector<vcl_scalar_type> vl_dSdt = viennacl::linalg::solve(vl_A, bb, gmres_tag,chow_patel_ilu2);
    std::vector<vcl_scalar_type> dSdt(vl_dSdt.size());
    viennacl::copy(vl_dSdt,dSdt);

#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->get_module_data<data>(ID);

        double qdep = is_nan(dSdt[i]) ? 0 : dSdt[i];

        double mass = 0;

        mass = qdep * global_param->dt();// kg/m^2*s *dt -> kg/m^2

        face->set_face_data("drift_mass", mass);
        d->sum_drift += mass;

        face->set_face_data("sum_drift", d->sum_drift);


    }
}

PBSM3D::~PBSM3D()
{

}
