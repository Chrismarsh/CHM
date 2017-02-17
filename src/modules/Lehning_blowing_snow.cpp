#include "Lehning_blowing_snow.hpp" 


Lehning_blowing_snow::Lehning_blowing_snow(config_file cfg)
        :module_base(parallel::domain)
{
    depends("U_2m_above_srf");
    depends("vw_dir");

    provides("c0");
    provides("c1");
    provides("c2");
    provides("c3");
    provides("c4");
    provides("c5");
    provides("c6");
    provides("c7");
    provides("c8");
    provides("c9");
    provides("c10");
    provides("c11");
    provides("c12");
    provides("c13");
    provides("c14");

    provides("hs");
    provides("ustar");


    provides("u_hs");

//    provides("ustar");
    provides("csalt");
    provides("Qsalt");

    provides("drift_depth");
    provides("drift_mass");
    provides("Qsusp");
    provides("ustar_t");
    provides("sum_Qdep");
}

void Lehning_blowing_snow::init(mesh domain)
{
    nLayer = 5;
    susp_depth = 5; //5m as per pomeroy
    v_edge_height = susp_depth / nLayer; //height of each vertical prism
    l__max = 40;

    do_vertical_advection = cfg.get("vertical_advection",true);

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

        //which faces have neighbours? Ie, are we an edge?
        for (int a = 0; a < 3; ++a)
        {
            auto neigh = face->neighbor(a);
            if (neigh == nullptr)
                d->face_neigh[a] = false;
            else
                d->face_neigh[a] = true;
        }



        face->set_face_data("sum_Qdep",0);
    }
}

void Lehning_blowing_snow::run(mesh domain)
{
    //hardcode at the moment

    //needed for linear system offsets
    size_t ntri = domain->number_of_faces();

    //we need to hold 1 matrix per thread so we get thread safety
//    std::vector<arma::sp_mat> sp_C;
//    for(int i = 0; i < omp_get_max_threads(); i++ )
//    {
//        sp_C.push_back(arma::sp_mat(ntri * nLayer, ntri * nLayer ) );
//    }

    //vcl_scalar_type is defined in the main CMakeLists.txt file.
    // Some GPUs do not have double precision so the run will fail is the wrong precision is used
    std::vector< std::map< unsigned int, vcl_scalar_type> > C(ntri * nLayer);

//    arma::sp_mat C(ntri * nLayer, ntri * nLayer ) ;
//    arma::vec b(ntri * nLayer,arma::fill::zeros );
//    arma::vec x(ntri * nLayer );

    std::vector<vcl_scalar_type> b(ntri * nLayer , 0.0);

    double z0 = 0.01; //m

    double k=0;
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
//        auto& C = sp_C.at(omp_get_thread_num());

        auto face = domain->face(i);
        auto d = face->get_module_data<data>(ID);
        auto& m = d->m;

        //get wind from the face
        double phi = face->face_data("vw_dir");
        double u2 =  face->face_data("U_2m_above_srf");

        Vector_2 v = -math::gis::bearing_to_cartesian(phi);

        //setup wind vector
        arma::vec uvw(3);
        uvw(0) = v.x(); //U_x
        uvw(1) = v.y(); //U_y
        uvw(2) = 0;

        double ustar = std::max(0.1, u2 * PhysConst::kappa / log(2. / z0));

        //depth of saltation layer
        double hs = 0;

        //lehning
//        z0 + 2.4025 * pow(u2, 2.) * pow(PhysConst::kappa, 2.) * pow(cos(25. * M_PI / 180.), 2.) /
//                         (pow(log(2. / z0), 2.) * 9.81);

        //pomeroy
        hs = 0.08436*pow(ustar,1.27);
        d->hs = hs;
        face->set_face_data("hs",hs);
        face->set_face_data("ustar",ustar);

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
        double tau_t_f = 0.2; // m/s
        double rho_f = 1.225; // air density, fix for T dependenc
        double u_star_t = pow(tau_t_f/rho_f,0.5); //threshold friction velocity
        double Qsalt = 0;
        double c_salt = 0;

        //ensure we can have snow saltation
        double rho_p = 917;
        double thresh_A = .18;
        double g = 9.81;
        double _d = 0.48e-3; //d in paper
        double u_star_saltation = thresh_A*sqrt((rho_p-rho_f)/(rho_f)*_d*g); // threshold for saltation to begin


        if( ustar > u_star_saltation)
        {
            //Pomeroy 1990
            c_salt = rho_f / (3.29 * ustar) * (1.0 - (u_star_t*u_star_t) / (ustar * ustar));

            Qsalt =  c_salt *std::max(0.1, Atmosphere::log_scale_wind(u2, 2, hs, 0))/2. * hs; //integrate over the depth of the saltation layer

            face->set_face_data("u_hs",std::max(0.1, Atmosphere::log_scale_wind(u2, 2, hs, 0)));

            if(Qsalt < 0) // this apparently happens and bad things happen bc of it
            {
                Qsalt = 0;
                c_salt = 0;
            }
        }

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
            double l = PhysConst::kappa * (cz + z0) * l__max / (PhysConst::kappa * cz + PhysConst::kappa * z0 + l__max);

            K[3] = K[4] = std::max(ustar * l, PhysConst::kappa * cz * ustar);

            //top
            alpha[3] = d->A[3] * K[3] / v_edge_height;
            //bottom
            alpha[4] = d->A[4] * K[4] / v_edge_height;

            //compute new U_z at this height in the suspension layer
            double u_z = std::max(0.1, Atmosphere::log_scale_wind(u2, 2, cz, 0));
            double length = arma::norm(uvw, 2);
            double scale = u_z / length;

            uvw *= scale;
            uvw(2) = -0.5; //settling velocity

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

            auto lol = &C[idx];

            for(int f = 0; f < 3; f++)
            {
                if(udotm[f] > 0)
                {
                    if (d->face_neigh[f])
                    {
                       C[idx][idx] += -d->A[f]*udotm[f]-alpha[f];
                       C[idx][idx] += alpha[f];
                    }
                    else
                    {

                       C[idx][idx] += -d->A[f]*udotm[f];
                    }
                } else
                {
                    if (d->face_neigh[f])
                    {
                        auto nidx = ntri * z + face->neighbor(f)->cell_id;

                        C[idx][idx] += -alpha[f];
                        C[idx][nidx] = -d->A[f]*udotm[f]+alpha[f];

                    }
                    else
                    {
                       C[idx][idx] += -d->A[f]*udotm[f]-alpha[f];
                    }
                }

            }

            //vertical layers
            //this formulation includes the 3D advection term
            if( z != nLayer -1 &&
                    C[idx].find(ntri * (z + 1) + face->cell_id) == C[idx].end() )
                C[idx][ntri * (z + 1) + face->cell_id] = 0.0;


            if(z!=0 &&
                    C[idx].find(ntri * (z - 1) + face->cell_id) == C[idx].end() )
                C[idx][ntri * (z - 1) + face->cell_id] = 0.0;

            if(do_vertical_advection)
            {
                if (z == 0)
                {

                    //bottom face, no advectiono
                    //                C[idx][idx] += -d->A[4] * K[4];
                    //                b[idx] = -d->A[4] * K[4] * c_salt+0.0;

                    if (udotm[3] > 0)
                    {
                        C[idx][idx] += (-d->A[3] * udotm[3] - alpha[3]);
                        C[idx][ntri * (z + 1) + face->cell_id] += alpha[3];
                    } else
                    {
                        C[idx][idx] += -alpha[3];
                        C[idx][ntri * (z + 1) + face->cell_id] += -d->A[3] * udotm[3] + alpha[3];
                    }


                    // just the smallest little bit of advection to keep it stable :(
                    auto uvw_0 = uvw;
                    uvw_0(2) = sng(uvw(2)) * 0.1; //match the sign of the top of this cell.

                    auto udotm_0 = arma::dot(uvw_0, m[4]);

                    if (udotm[4] > 0)
                    {
                        C[idx][idx] += -d->A[4] * K[4] - d->A[4] * udotm_0;
                        b[idx] = -d->A[4] * K[4] * c_salt;
                    } else
                    {
                        C[idx][idx] += -d->A[4] * K[4];
                        b[idx] += (d->A[4] * K[4] + d->A[4] * udotm_0) * c_salt;
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
                } else
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


//        arma::sp_mat result(ntri * nLayer, ntri * nLayer );
//    for(int i = 0; i < omp_get_max_threads(); i++ )
//    {
//        result += sp_C.at(i);
//    }

    //setup the compressed matrix on the compute device, in available
//    viennacl::compressed_matrix<vcl_scalar_type>  vl_C(ntri * nLayer, ntri * nLayer);
//    viennacl::copy(C,vl_C);
//    viennacl::vector<vcl_scalar_type> rhs(ntri * nLayer);
//    viennacl::copy(b,rhs);
////
////
////
////    // configuration of preconditioner:
//    viennacl::linalg::chow_patel_tag chow_patel_ilu_config;
//    chow_patel_ilu_config.sweeps(3);       // three nonlinear sweeps
//    chow_patel_ilu_config.jacobi_iters(2); // two Jacobi iterations per triangular 'solve' Rx=r
//    viennacl::linalg::chow_patel_ilu_precond< viennacl::compressed_matrix<vcl_scalar_type> > chow_patel_ilu(vl_C, chow_patel_ilu_config);


    LOG_DEBUG << "solving";
    //docs say 'This will use appropriate ViennaCL objects internally.'  http://viennacl.sourceforge.net/doc/iterative_8cpp-example.html
    auto x = viennacl::linalg::solve(C, b, viennacl::linalg::bicgstab_tag());
//    auto x = viennacl::linalg::solve(result, b, viennacl::linalg::bicgstab_tag());

//   auto x = viennacl::linalg::solve(vl_C, rhs, viennacl::linalg::gmres_tag(),chow_patel_ilu);


#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->get_module_data<data>(ID);
        double Qsusp = 0;
        double hs = d->hs;

        double u2 = face->face_data("U_2m_above_srf");
        for (int z = 0; z<nLayer;++z)
        {
            double c = x[ntri * z + face->cell_id];
            double cz = z + hs+ v_edge_height/2.; //cell center height

            double u_z =std::max(0.1, Atmosphere::log_scale_wind(u2, 2, cz, 0)); //compute new U_z at this height in the suspension layer
            Qsusp += c * u_z * v_edge_height; /// kg/m^3 ---->  kg/(m.s)

            if(z < 15)
                face->set_face_data("c"+std::to_string(z), c ); //<0?0:c
        }
        face->set_face_data("Qsusp",Qsusp);
    }

#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->get_module_data<data>(ID);
        auto& m = d->m;

        double u2 = face->face_data("U_2m_above_srf");
        double phi = face->face_data("vw_dir");
        Vector_2 v = -math::gis::bearing_to_cartesian(phi);

        //setup wind vector
        arma::vec uvw(3);
        uvw(0) = v.x(); //U_x
        uvw(1) = v.y(); //U_y
        uvw(2) = 0;

        double udotm[3];

        //edge lengths b/c 2d now
        double E[3]={0,0,0};

        for (int j = 0; j < 3; ++j)
        {
            auto neigh = face->neighbor(j);

            //just unit vectors as qsusp/qsalt flux has magnitude
            udotm[j] = arma::dot(uvw, m[j]);

            E[j] = face->edge_length(j);
        }

        interpolation interp(interp_alg::tpspline);

        std::vector< boost::tuple<double, double, double> > vec_qs;
        std::vector< boost::tuple<double, double, double> > vec_qt;


        //build up the interpreted values
        for(int j = 0; j < 3; j++)
        {
            if (d->face_neigh[j])
            {
                auto neigh = face->neighbor(j);
                vec_qs.push_back( boost::make_tuple(neigh->center().x(),neigh->center().y(), neigh->face_data("Qsalt")));
                vec_qt.push_back( boost::make_tuple(neigh->center().x(),neigh->center().y(), neigh->face_data("Qsusp")));
            }
        }


        double qdep = 0;
        for(int j = 0; j < 3; j++)
        {
            auto emp = face->edge_midpoint(j);
            auto query = boost::make_tuple(emp.x(), emp.y(), 0.0); //z isn't used in the interp call

            double qs  = interp(vec_qs, query);
            double qt  = interp(vec_qt, query);

            qdep += E[j]*udotm[j]*(qs+qt);
        }

        double mass = -qdep * global_param->dt();
        double depth = -qdep * global_param->dt() / 400. / face->get_area();  // 1000 kg/m^3 -> density of water
        face->set_face_data("drift_depth",depth);
        face->set_face_data("drift_mass",mass);

        face->set_face_data("sum_Qdep", face->face_data("sum_Qdep") + qdep);
    }
}

Lehning_blowing_snow::~Lehning_blowing_snow()
{

}
