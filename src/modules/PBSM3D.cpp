#include "PBSM3D.hpp"


PBSM3D::PBSM3D(config_file cfg)
        :module_base(parallel::domain)
{
    depends("U_2m_above_srf");
    depends("vw_dir");
    depends("swe");
    depends("t");
    depends("rh");

    provides("u10");
    provides("is_drifting");


    provides("c0");
    provides("c1");
    provides("c2");
    provides("c3");
    provides("c4");

    provides("K0");
    provides("K1");
    provides("K2");
    provides("K3");
    provides("K4");

//    nLayer=10;
//    for(int i=0; i<nLayer;++i)
//    {
//        provides("K"+std::to_string(i));
//        provides("c"+std::to_string(i));
//    }

    provides("Qsusp_pbsm");

    provides("hs");
    provides("ustar");

    provides("csalt");


    provides("u*_th");

    provides("drift_mass"); //kg/m^2
    provides("drift_mass_no_subl");

    provides("Qsusp");
    provides("Qsubl");
    provides("Qsalt");

    provides("sum_drift");
    provides("sum_subl");
}

void PBSM3D::init(mesh domain)
{
    nLayer = 5;


    susp_depth = 5; //5m as per pomeroy
    v_edge_height = susp_depth / nLayer; //height of each vertical prism
    l__max = 40; // mixing length for diffusivity calculations

    settling_velocity = cfg.get("settling_velocity",-0.5); // m/s, Lehning, M., H. Löwe, M. Ryser, and N. Raderschall (2008), Inhomogeneous precipitation distribution and snow transport in steep terrain, Water Resour. Res., 44(7), 1–19, doi:10.1029/2007WR006545.

    if(settling_velocity > 0)
        BOOST_THROW_EXCEPTION(module_error() << errstr_info ("PBSM3D settling velocity must be negative"));

    snow_diffusion_const = cfg.get("snow_diffusion_const",0.005); // Beta * K, this is beta and scales the eddy diffusivity
    do_vertical_advection = cfg.get("vertical_advection",true);
    do_sublimation = cfg.get("do_sublimation",true);
    eps = cfg.get("smooth_coeff",1e-5);
    limit_mass= cfg.get("limit_mass",true);
    min_mass_for_trans = cfg.get("min_mass_for_trans",10);
    n_non_edge_tri = 0;
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);


        auto d = face->make_module_data<data>(ID);

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
            auto neigh = face->neighbor(a);
            if (neigh == nullptr)
            {
                d->face_neigh[a] = false;
                d->is_edge = true;
            }
            else
            {
                d->face_neigh[a] = true;
            }

        }
        if(!d->is_edge)
        {
            d->cell_id=n_non_edge_tri;
            ++n_non_edge_tri;
        }


        face->set_face_data("sum_drift",0);
        face->set_face_data("sum_subl",0);
    }
}

void PBSM3D::run(mesh domain)
{
    //needed for linear system offsets
    size_t ntri = domain->number_of_faces();

    //vcl_scalar_type is defined in the main CMakeLists.txt file.
    // Some GPUs do not have double precision so the run will fail if the wrong precision is used
    std::vector< std::map< unsigned int, vcl_scalar_type> > C(ntri * nLayer);
    std::vector<vcl_scalar_type> b(ntri * nLayer , 0.0);

    //ice density
    double rho_p = PhysConst::rho_ice;

#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->get_module_data<data>(ID);
        auto& m = d->m;

        //get wind from the face
        double phi = face->face_data("vw_dir");
        double u2 = face->face_data("U_2m_above_srf");
        double u10 = Atmosphere::log_scale_wind(u2, 2, 10, 0);
        face->set_face_data("u10",u10);

        Vector_2 v = -math::gis::bearing_to_cartesian(phi);

        //setup wind vector
        arma::vec uvw(3);
        uvw(0) = v.x(); //U_x
        uvw(1) = v.y(); //U_y
        uvw(2) = 0;

        //solve for ustar as perturbed by blowing snow
        //  not 100% sure this should be done w/o blowing snow. Might need to revisit this
        double ustar = -.2000000000*u2/gsl_sf_lambert_Wm1(-0.1107384167e-1*u2);
        ustar = std::max(0.1,ustar);
        d->z0 = std::max(0.001,0.1203 * ustar * ustar / (2.0*9.81));

        //depth of saltation layer
        double hs = 0;

        //lehning
//        hs = z0 + 2.4025 * pow(u2, 2.) * pow(PhysConst::kappa, 2.) * pow(cos(25. * M_PI / 180.), 2.) /
//                         (pow(log(2. / z0), 2.) * 9.81);

        //pomeroy
        hs = 0.08436*pow(ustar,1.27);
        d->hs = hs;
        face->set_face_data("hs",hs);

        // Assuming no horizontal diffusion of blowing snow. Thus the below section does not need to be computed
        // If we add in horizontal diffusion (not sure why), then this will have to be computed on a per-layer basis.

        //eddy diffusivity (m^2/s)
        // 0,1,2 will all be K = 0, as no horizontal diffusion process
        double K[5] = {0, 0, 0, 0, 0};

        // holds A_ * K_ / h_
        // _0 -> _2 are the horizontal sides
        // _3 -> is the top of the prism
        // _4 -> is the bottom the prism
        double alpha[5] = {0, 0, 0, 0, 0};

        //compute alpha and K for edges
//        for (int a = 0; a < 3; ++a)
//        {
//            auto neigh = face->neighbor(a);
//            alpha[a] = d->A[a];
//            double l=40;
//            //if we have a neighbour, use the distance
//            if (neigh != nullptr)
//            {
//                l = math::gis::distance(face->center(), neigh->center());
//
//                alpha[a] /= math::gis::distance(face->center(), neigh->center());
//
//            } else
//            { //otherwise assume 2x the distance from the center of face to one of it's vertexes, so ghost triangle is approx the same size
//                alpha[a] /= 2.0 * math::gis::distance(face->center(), face->vertex(0)->point());
//
//            }
//
//            //no horizontal diffusion .....
////            double l = PhysConst::kappa * (cz + z0) * l__max /
////                       (PhysConst::kappa * cz + PhysConst::kappa * z0 + l__max);
//
//            K[a] = std::max(ustar * l, PhysConst::kappa * 2. * ustar);
//            alpha[a] *= K[a];
//        }

        //saltation concentration
        double tau_t_f = 0.2; // m/s threshold shear stress for aerodynamic entrainment
        double rho_f = 1.225; // air density, fix for T dependenc
        double u_star_t = pow(tau_t_f/rho_f,0.5); //threshold friction velocity
        double Qsalt = 0;
        double c_salt = 0;

        //ensure we can have snow saltation

        double thresh_A = .18;
        double g = 9.81;
        double _d = 0.48e-3; //d in paper
        double u_star_saltation = thresh_A*sqrt((rho_p-rho_f)/(rho_f)*_d*g); // threshold for saltation to begin
        double t = face->face_data("t");

        face->set_face_data("u*_th",u_star_saltation);

        double swe = face->face_data("swe"); // mm   -->    kg/m^2
        swe = is_nan(swe) ? 0 : swe; // handle the first timestep where swe won't have been updated if we override the module order

        face->set_face_data("ustar",ustar);

        face->set_face_data("is_drifting",0);
        face->set_face_data("Qsusp_pbsm",0); //for santiy checks against pbsm

        if( ustar > u_star_saltation && swe > min_mass_for_trans)
        {

            double pbsm_qsusp = pow(u10,4.13)/674100.0;
            face->set_face_data("Qsusp_pbsm",pbsm_qsusp);
            face->set_face_data("is_drifting",1);
            //Pomeroy 1990

            //it's not really clear what the u_star_t is, but I think it's actually the threshold.
//            c_salt = rho_f / (3.29 * ustar) * (1.0 - (u_star_t*u_star_t) / (ustar * ustar));
            c_salt = rho_f / (3.29 * ustar) * (1.0 - (u_star_saltation*u_star_saltation) / (ustar * ustar));

            // occasionally happens to happen at low wind speeds where the parameterization breaks.
            // hasn't happened since changed this to the threshold velocity though...
            if(c_salt < 0 || std::isnan(c_salt))
            {
                c_salt = 0;
            }

            //mean wind speed in the saltation layer
            double uhs = std::max(0.1,Atmosphere::log_scale_wind(u2, 2, hs, 0,d->z0)/2.);

            // kg/(m*s)
            Qsalt =  c_salt * uhs * hs; //integrate over the depth of the saltation layer, kg/(m*s)

            //calculate the surface integral of Qsalt, and ensure we aren't saltating more mass than what exists
            //in the triangle. In this case, we are approximating the edge value with just Qsalt, and are not considering
            //the neighbour values.
            double salt=0;

            double udotm[3];


            for(int j=0; j<3; ++j)
            {
                udotm[j] = arma::dot(uvw, m[j]);
                salt += face->edge_length(j) * udotm[j] * Qsalt;
            }
            salt /= face->get_area(); // -> kg/m^2*s

            salt *= global_param->dt(); // ->kg/m^2



            //Figure out the max we can transport, ie total mass in cell
            if( limit_mass && salt > swe) // could we move more than the total mass in the cell during this timestep?
            {
                double el0,el1,el2;
                el0 = face->edge_length(0);
                el1 = face->edge_length(1);
                el2 = face->edge_length(2);
                double dt = global_param->dt();

                //back out what the max conc should be based on our swe
                //units: ((kg/m^2)*m^2)/( s*m*(m/s)*m ) -> kg/m^3
                c_salt = swe*face->get_area()/(dt*hs*uhs*(el0*udotm[0]+el1*udotm[1]+el2*udotm[2]));

                Qsalt = c_salt * uhs * hs;
                if(is_nan(c_salt)) //if we have no swe this happens
                {
                    c_salt = 0;
                    Qsalt = 0;
                }

//                LOG_DEBUG << "More saltation than snow, limiting conc to " << c_salt << " triangle="<<i;
//                LOG_DEBUG << "Avail mass = " << swe << ", would have salted =  " << salt;
            }
        }

        // can use for point scale plume testing.
//        c_salt = 0;
//        if (i == 12125 )
//        {
//            Qsalt = 0;
//            c_salt = 0.8;
//        }


        face->set_face_data("csalt", c_salt);
        face->set_face_data("Qsalt", Qsalt);


        // iterate over the vertical layers
        for (int z = 0; z < nLayer; ++z)
        {
            //height in the suspension layer
            double cz = z + hs+ v_edge_height/2.; //cell center height
            double l = PhysConst::kappa * (cz + d->z0) * l__max / (PhysConst::kappa * cz + PhysConst::kappa * d->z0 + l__max);

            //snow_diffusion_const is pretty much a calibration constant. At 1 it seems to over predict transports.
            K[3] = K[4] = snow_diffusion_const * std::max(ustar * l, PhysConst::kappa * cz * ustar);
            face->set_face_data("K"+std::to_string(z), K[3] );
            //top
            alpha[3] = d->A[3] * K[3] / v_edge_height;
            //bottom
            alpha[4] = d->A[4] * K[4] / v_edge_height;

            //compute new U_z at this height in the suspension layer
            double u_z = std::max(0.1, Atmosphere::log_scale_wind(u2, 2, cz, 0,d->z0));
            double length = arma::norm(uvw, 2);
            double scale = u_z / length;

            uvw *= scale;
            uvw(2) = settling_velocity; //settling velocity,

            //holds wind dot face normal
            double udotm[5];
            for (int j = 0; j < 5; ++j)
            {
                udotm[j] = arma::dot(uvw, m[j]);
            }
            //lateral
            size_t idx = ntri*z + face->cell_id;

            //I'm unclear as to what the map inits the double value to.  I assume undef'd like always, so this is called
            //to ensure that they are inited to zero so we can call the += operator w/o issue below. Only needs to be done for the diagonal elements
            //as the rest are straight assignments and not +=

            if( C[idx].find(idx) == C[idx].end() )
                C[idx][idx] = 0.0;

            for(int f = 0; f < 3; f++)
            {
                if(udotm[f] > 0)
                {
                    if (d->face_neigh[f])
                    {
                        size_t nidx = ntri * z + face->neighbor(f)->cell_id;
                        if( C[idx].find(nidx) == C[idx].end() )
                            C[idx][nidx] = 0.0;

                       C[idx][idx] += -d->A[f]*udotm[f]-alpha[f];
                       C[idx][nidx] += alpha[f];

                    }
                    else
                    {
                       C[idx][idx] += -d->A[f]*udotm[f];

                    }
                } else
                {
                    if (d->face_neigh[f])
                    {
                        size_t nidx = ntri * z + face->neighbor(f)->cell_id;
                        if( C[idx].find(nidx) == C[idx].end() )
                            C[idx][nidx] = 0.0;

                        C[idx][idx] += -alpha[f];
                        C[idx][nidx] += -d->A[f]*udotm[f]+alpha[f];
                    }
                    else
                    {
                        C[idx][idx] += -alpha[f];

                    }
                }
            }



            //init to zero the vertical component regardless of do_vertical_advection
            if( z != nLayer -1 &&
                C[idx].find(ntri * (z + 1) + face->cell_id) == C[idx].end() )
                C[idx][ntri * (z + 1) + face->cell_id] = 0.0;

            if(z!=0 &&
               C[idx].find(ntri * (z - 1) + face->cell_id) == C[idx].end() )
                C[idx][ntri * (z - 1) + face->cell_id] = 0.0;

            if(do_vertical_advection)
            {
                //vertical layers
                //this formulation includes the 3D advection term
                if (z == 0)
                {
                    //bottom face, no advection
                    C[idx][idx] += -d->A[4] * K[4];
                    b[idx] = -d->A[4] * K[4] * c_salt+0.0;

                    if (udotm[3] > 0)
                    {
                        C[idx][idx] += (-d->A[3] * udotm[3] - alpha[3]);
                        C[idx][ntri * (z + 1) + face->cell_id] += alpha[3];
                    } else
                    {
                        C[idx][idx] += -alpha[3];
                        C[idx][ntri * (z + 1) + face->cell_id] += -d->A[3] * udotm[3] + alpha[3];
                    }


                } else if (z == nLayer - 1)// top z layer
                {
                    if (udotm[3] > 0)
                    {
                        C[idx][idx] += -d->A[3] * udotm[3] - alpha[3];
                    } else
                    {
                        C[idx][idx] += -alpha[3];
                    }


                    if (udotm[4] > 0)
                    {
                        C[idx][idx] += -d->A[4] * udotm[4] - alpha[4];
                        C[idx][ntri * (z - 1) + face->cell_id] += alpha[4];
                    } else
                    {
                        C[idx][idx] += -alpha[4];
                        C[idx][ntri * (z - 1) + face->cell_id] += -d->A[4] * udotm[4] + alpha[4];
                    }
                } else //middle layers
                {
                    if (udotm[3] > 0)
                    {
                        C[idx][idx] += -d->A[3] * udotm[3] - alpha[3];
                        C[idx][ntri * (z + 1) + face->cell_id] += alpha[3];
                    } else
                    {
                        C[idx][idx] += -alpha[3];
                        C[idx][ntri * (z + 1) + face->cell_id] += -d->A[3] * udotm[3] + alpha[3];
                    }

                    if (udotm[4] > 0)
                    {
                        C[idx][idx] += -d->A[4] * udotm[4] - alpha[4];
                        C[idx][ntri * (z - 1) + face->cell_id] += alpha[4];
                    } else
                    {
                        C[idx][idx] += -alpha[4];
                        C[idx][ntri * (z - 1) + face->cell_id] += -d->A[4] * udotm[4] + alpha[4];
                    }
                }
            } else
            {
                //This section has the vertical case for diffusion only.

                if (z == 0)
                {
                    //bottom face
                   C[idx][idx] += -d->A[4]*K[4];
                   b[ntri * z + face->cell_id] = -d->A[4]*K[4]*c_salt;

                    //top face
                    C[idx][ntri * z + face->cell_id] +=  -alpha[3];
                    C[idx][ntri * (z + 1) + face->cell_id] +=  alpha[3];


                } else if (z == nLayer - 1)// top z layer
                {
                    //top
                    C[idx][ntri * z + face->cell_id] += -alpha[3]-alpha[4];
                    //bottom
                    C[idx][ntri * (z - 1) + face->cell_id] += alpha[4];


                } else // internal cell
                {

                    C[idx][idx] += -alpha[3]-alpha[4];
                    //top
                    C[idx][ntri * (z + 1) + face->cell_id] += alpha[3];
                    //bottom
                    C[idx][ntri * (z - 1) + face->cell_id] += alpha[4];

                }

            }
        } // end z iter
    } //end face iter


    //setup the compressed matrix on the compute device, if available
    viennacl::compressed_matrix<vcl_scalar_type>  vl_C(ntri * nLayer, ntri * nLayer);
    viennacl::copy(C,vl_C);
    viennacl::vector<vcl_scalar_type> rhs(ntri * nLayer);
    viennacl::copy(b,rhs);

    // configuration of preconditioner:
    viennacl::linalg::chow_patel_tag chow_patel_ilu_config;
    chow_patel_ilu_config.sweeps(3);       //  nonlinear sweeps
    chow_patel_ilu_config.jacobi_iters(2); //  Jacobi iterations per triangular 'solve' Rx=r
    viennacl::linalg::chow_patel_ilu_precond< viennacl::compressed_matrix<vcl_scalar_type> > chow_patel_ilu(vl_C, chow_patel_ilu_config);


    //compute result and copy back to CPU device (if an accelerator was used), otherwise access is slow
    viennacl::linalg::gmres_tag gmres_tag(1e-3, 500, 30);
    viennacl::vector<vcl_scalar_type> vl_x = viennacl::linalg::solve(vl_C, rhs, gmres_tag, chow_patel_ilu);
    std::vector<vcl_scalar_type> x(vl_x.size());
    viennacl::copy(vl_x,x);

//    LOG_DEBUG << "No. of iters: " << gmres_tag.iters();

#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->get_module_data<data>(ID);
        double Qsusp = 0;
        double Qsubl = 0;
        double hs = d->hs;

        double u2 = face->face_data("U_2m_above_srf");

        double rh = face->face_data("rh")/100.;
        double t = face->face_data("t")+273.15;
        double es = mio::Atmosphere::saturatedVapourPressure(t);
        double ea = rh * es / 1000.; // e in kpa
        double P = mio::Atmosphere::stdAirPressure(face->get_z())/1000.;// kpa //might be a better fun for this

        //specific humidity of the air at air temp
        double q = 0.633*ea/P;


        double v = 1.88*pow(10.,-5.); //kinematic viscosity of air

        for (int z = 0; z<nLayer;++z)
        {
            double c = x[ntri * z + face->cell_id];
            c = c < 0 || is_nan(c) ? 0 : c; //harden against some numerical issues that occasionally come up for unknown reasons.

            double cz = z + hs+ v_edge_height/2.; //cell center height

            double u_z =std::max(0.1, Atmosphere::log_scale_wind(u2, 2, cz, 0,d->z0)); //compute new U_z at this height in the suspension layer
            Qsusp += c * u_z * v_edge_height; /// kg/m^3 ---->  kg/(m.s)


            face->set_face_data("c"+std::to_string(z), c );

            //calculate dm/dt from
            // equation 13 from Pomeroy and Li 2000
            // To do so, use equations 12 - 16 in Pomeroy et al 1999

            double rm = 4.6*pow(10.0,-5.) * pow((double)cz,-0.258); // eqn 18, mean particle size, also in liston 1998, eq A-4
            double xrz = 0.005 * pow(u_z,1.36);//eqn 16
            double omega = 1.1*pow(10.,7.) * pow(rm,1.8);
            double Vr = omega + 3.0*xrz*cos(M_PI/4.0);

            double Re = 2.0*rm*Vr / v;

            double Nu, Sh;
            Nu = Sh = 1.79 + 0.606 * pow(Re,0.5);

            double D = 2.06*pow(10.,-5.)*pow(t/(273.0),1.75); //diffusivity of water vapour in air, t in K, eqn A-7 in Liston 1998

            double lambda_t = 0.000063*(t-273.15)+0.00673; // J/(kmol K) thermal conductivity looks like this is degC, not K. order of magnitude off if K, eqn 11 Pomeroy 1993
            double Ls = 2.838*pow(10.,6.); // latent heat of sublimation
            double rho_a = mio::Atmosphere::stdDryAirDensity(face->get_z(),t); //kg/m^3, comment in mio is wrong.
            auto Tsfn = [&](double Ts) -> double
            {
                double es_Ts = mio::Atmosphere::saturatedVapourPressure(Ts);
                double ea_ts = 1. * es_Ts /1000.; //kpa

                //specific humidity of the particle at the particle temp
                double qTs = 0.633*ea_ts/P;
                double result = (D*Sh*Ls*q*rho_a-D*Sh*Ls*qTs*rho_a+Nu*t*lambda_t)/(lambda_t*Nu)-Ts;
                return result;
            };

            double guess = t;
            double factor = 1;
            int digits = std::numeric_limits<double>::digits - 3;
            double min = 200;
            double max = 300;
            boost::uintmax_t max_iter=500;
            boost::math::tools::eps_tolerance<double> tol(30);
            auto r = boost::math::tools::toms748_solve(Tsfn, min, max, tol, max_iter);
            double Ts = r.first + (r.second - r.first)/2;

            //now use equation 13 with our solved Ts to compute dm/dt(z)
            double dmdtz = 2.0 * M_PI * rm * lambda_t / Ls * Nu * (Ts - t);  //eqn 13 in Pomeroy and Li 2000

            //calculate mean mass, eqn 23, 24 in Pomeroy 1993
            double alpha = 4.08 + 12.6*cz;
            double mm = 4./3. * M_PI * rho_p * rm*rm*rm *(1.0 + 3.0/alpha + 2./(alpha*alpha)); //mean mass, eqn 23

            double csubl = dmdtz/mm;

            //eqn 20 in Pomeroy 1993
            Qsubl += csubl * c * v_edge_height; //kg/(m^2*s)
        }
        face->set_face_data("Qsusp",Qsusp);
        face->set_face_data("Qsubl",Qsubl);
    }

    std::vector< std::map< unsigned int, vcl_scalar_type> > A(ntri);
    std::vector<vcl_scalar_type> bb(ntri, 0.0);

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

        if (A[i].find(i) == A[i].end())
            A[i][i] = 0.0;


        for (int j = 0; j < 3; j++)
        {
            if (d->face_neigh[j])
            {
                auto neigh = face->neighbor(j);
                auto Qtj = neigh->face_data("Qsusp") + face->face_data("Qsusp");
                auto Qsj = neigh->face_data("Qsalt") + face->face_data("Qsalt");
                double Qt = Qtj / 2.0 + Qsj / 2.0;

                if (A[i].find(neigh->cell_id) == A[i].end())
                    A[i][neigh->cell_id] = 0.0;

                dx[j] = math::gis::distance(face->center(), neigh->center());

                A[i][i] += -eps * E[j] / dx[j] - V;
                A[i][neigh->cell_id] += eps * E[j] / dx[j];
                bb[i] += E[j] * Qt * udotm[j];

            } else
            {
                auto Qtj = 2 * face->face_data("Qsusp"); // const flux across, 0 -> drifts!!!
                auto Qsj = 2 * face->face_data("Qsalt");
                double Qt = Qtj / 2.0 + Qsj / 2.0;
                A[i][i]  += -V;
                bb[i] += E[j] * Qt * udotm[j];
            }

        }



    } // end face itr


    viennacl::compressed_matrix<vcl_scalar_type>  vl_A(ntri, ntri);
    viennacl::copy(A,vl_A);
    viennacl::vector<vcl_scalar_type> rhs_bb(ntri);
    viennacl::copy(bb,rhs_bb);

    // configuration of preconditioner:
    viennacl::linalg::chow_patel_ilu_precond< viennacl::compressed_matrix<vcl_scalar_type> > chow_patel_ilu2(vl_A, chow_patel_ilu_config);

    //compute result and copy back to CPU device (if an accelerator was used), otherwise access is slow
    viennacl::vector<vcl_scalar_type> vl_dSdt = viennacl::linalg::solve(vl_A, rhs_bb, viennacl::linalg::gmres_tag(),chow_patel_ilu2);
    std::vector<vcl_scalar_type> dSdt(vl_dSdt.size());
    viennacl::copy(vl_dSdt,dSdt);


//    auto dSdt = viennacl::linalg::solve(A, bb, viennacl::linalg::bicgstab_tag());

#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->get_module_data<data>(ID);

        double subl_mass_flux = face->face_data("Qsubl");

        double qdep = is_nan(dSdt[i]) ? 0 : dSdt[i]; //d->is_edge ||

        double mass = 0;
        if(do_sublimation)
            mass = (qdep + subl_mass_flux) * global_param->dt(); // kg/m^2*s *dt -> kg/m^2
        else
            mass = qdep * global_param->dt();// kg/m^2*s *dt -> kg/m^2

        face->set_face_data("drift_mass", mass);

        double sum_drift = face->face_data("sum_drift");
        face->set_face_data("sum_drift", sum_drift + mass);

        double sum_subl = face->face_data("sum_subl");
        sum_subl += subl_mass_flux * global_param->dt();
        face->set_face_data("sum_subl", sum_subl);
    }
}

PBSM3D::~PBSM3D()
{

}
