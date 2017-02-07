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

    provides("hs");
    provides("ustar");


    provides("u_hs");

//    provides("ustar");
    provides("csalt");
    provides("Qsalt");
    provides("Qsalt_pbsm");
    provides("csalt_pbsm");
    provides("Qdep");
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

    std::vector< std::map< unsigned int, double> > C(ntri * nLayer);

//    arma::sp_mat C(ntri * nLayer, ntri * nLayer ) ;
//    arma::vec b(ntri * nLayer );
//    arma::vec x(ntri * nLayer );

    std::vector<double> b(ntri * nLayer ,0.0);
//    b.clear();


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
        double u2 = face->face_data("U_2m_above_srf");

        Vector_2 v = -math::gis::bearing_to_cartesian(phi);

        //setup wind vector
        arma::vec uvw(3);
        uvw(0) = v.x(); //U_x
        uvw(1) = v.y(); //U_y
        uvw(2) = 0;

        double ustar = std::max(0.1, u2 * PhysConst::kappa / log(2. / z0));

        //depth of saltation layer
        double hs = z0 + 2.4025 * pow(u2, 2.) * pow(PhysConst::kappa, 2.) * pow(cos(25. * M_PI / 180.), 2.) /
                         (pow(log(2. / z0), 2.) * 9.81);

        hs = 0.08436*pow(ustar,1.27);
        d->hs = hs;

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
//            alpha[a] = A[a];
//
//            //if we have a neighbour, use the distance
//            if (neigh != nullptr)
//            {
//                l__max = math::gis::distance(face->center(), neigh->center());
//
//                alpha[a] /= math::gis::distance(face->center(), neigh->center());
//
//            } else
//            { //otherwise assume 2x the distance from the center of face to one of it's vertexes, so ghost triangle is approx the same size
//                alpha[a] /= 2.0 * math::gis::distance(face->center(), face->vertex(0)->point());
//
//                face_neigh[a] = false;
//            }
//
//            //no horizontal diffusion .....
////            double l = PhysConst::kappa * (cz + z0) * l__max /
////                       (PhysConst::kappa * cz + PhysConst::kappa * z0 + l__max);
//
//            K[a] = 0;//std::max(ustar * l, PhysConst::kappa * 2. * ustar);
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


        face->set_face_data("ustar",ustar);
        if( ustar > u_star_saltation)
        {
            double Qsalt_pbsm = 0.68 * rho_f * u_star_t / (g *ustar) *( ustar*ustar - u_star_t*u_star_t);
            double csalt_pbsm = Qsalt_pbsm / std::max(0.1, Atmosphere::log_scale_wind(u2, 2, hs, 0))/2. / hs;
            face->set_face_data("Qsalt_pbsm", Qsalt_pbsm);
            face->set_face_data("csalt_pbsm", csalt_pbsm);


            //Pomeroy 1990
            c_salt = rho_f / (3.29 * ustar) * (1.0 - (u_star_t*u_star_t) / (ustar * ustar));
//            c_salt = csalt_pbsm;

            Qsalt =  c_salt *std::max(0.1, Atmosphere::log_scale_wind(u2, 2, hs, 0))/2. * hs; //integrate over the depth of the saltation layer


            face->set_face_data("u_hs",std::max(0.1, Atmosphere::log_scale_wind(u2, 2, hs, 0)));
            face->set_face_data("hs",hs);

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


        // iterate over the vertical layers, each 1m in height
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
            double u_z =std::max(0.1, Atmosphere::log_scale_wind(u2, 2, cz, 0));
            double length = arma::norm(uvw, 2);
            double scale = u_z / length;

            uvw *= scale;
//          uvw(2) = -.5;

            //holds wind dot face normal
            double udotm[5];
            for (int j = 0; j < 5; ++j)
            {
                udotm[j] = arma::dot(uvw, m[j]);
            }

            //I'm unclear as to what the map inits the double value to.  I assume undef'd like always, so this is called
            //to ensure that they are inited to zero so we can call the += operator w/o issue below.
            auto init_C = [&](std::vector< std::map< unsigned int, double> >& C, unsigned int i, unsigned int j){
                if( C[i].find(j) == C[i].end() )
                    C[i][j] = 0;
            };

            //lateral
            size_t idx = ntri*z + face->cell_id;
            for(int f = 0; f < 3; f++)
            {

                if(udotm[f] > 0)
                {
                    if (d->face_neigh[f])
                    {
//                        init_C(C,idx,idx);

                       C[idx][idx] += -d->A[f]*udotm[f]-alpha[f];
                       C[idx][idx] += alpha[f];
                    }
                    else
                    {
//                        init_C(C,idx,idx);

                       C[idx][idx] += -d->A[f]*udotm[f];
                    }
                } else
                {
                    if (d->face_neigh[f])
                    {
//                        init_C(C,idx,idx);
//                        init_C(C,idx,nidx);
                        auto nidx = ntri * z + face->neighbor(f)->cell_id;

                        C[idx][idx] += -alpha[f];
                        C[idx][nidx]+= -d->A[f]*udotm[f]+alpha[f];
                    }
                    else
                    {
//                        init_C(C,idx,idx);;

                       C[idx][idx] += -d->A[f]*udotm[f]-alpha[f];
                    }
                }

            }



//            for(int f = 0; f < 3; f++)
//            {
//
//
//                if (d->face_neigh[f])
//                {
//                    auto nidx = ntri * z + face->neighbor(f)->cell_id;
//                    if(udotm[f] > 0 )
//                    {
//                        if( C[idx].find(idx) == C[idx].end() )
//                           C[idx][idx] = 0;
//
//                       C[idx][idx] += -d->A[f]*udotm[f]-alpha[f];
//                       C[idx][idx] += alpha[f];
//                    } else
//                    {
//                        if( C[idx].find(nidx) == C[idx].end() )
//                           C[idx][nidx] = 0;
//
//                       C[idx][idx] += -alpha[f];
//                       C[idx][nidx]+= -d->A[f]*udotm[f]+alpha[f];
//                    }
//
//                }
//                else
//                {
//                    if( C[idx].find(idx) == C[idx].end() )
//                       C[idx][idx] = 0;
//
//                   C[idx][idx] += -d->A[f]*udotm[f]+alpha[f];
//                }
//            }


            // 1 layer special case
            if (nLayer == 1)
            {
               C[idx][idx] += -d->A[4]*K[4]+alpha[3];
                b[idx] += -d->A[4]*K[4]*c_salt;
                //top
                //0
            } else
            {
                //now do the vertical fluxes
                //we are the bottom?
                if (z == 0)
                {
                    //bottom face
                   C[idx][idx] += -d->A[4]*K[4];
                    b[ntri * z + face->cell_id] += -d->A[4]*K[4]*c_salt;

//                    init_C(C,idx,ntri * z + face->cell_id);
//                    init_C(C,idx,ntri * (z + 1) + face->cell_id);
                    //top face
                    C[idx][ntri * z + face->cell_id] +=  -alpha[3];
                    C[idx][ntri * (z + 1) + face->cell_id] =  alpha[3];


                } else if (z == nLayer - 1)// top z layer
                {
//                    init_C(C,idx,ntri * z + face->cell_id);
//                    init_C(C,idx,ntri * (z - 1) + face->cell_id);
                    //top
                    C[idx][ntri * z + face->cell_id] += -alpha[3]-alpha[4];

                    //bottom
                    C[idx][ntri * (z - 1) + face->cell_id] = alpha[4];

                } else // internal cell
                {
                   C[idx][idx] += -alpha[3]-alpha[4];

//                    init_C(C,idx,ntri * (z + 1) + face->cell_id);
//                    init_C(C,idx,ntri * (z - 1) + face->cell_id);

                    //top
                    C[idx][ntri * (z + 1) + face->cell_id]=alpha[3];

                    //bottom
                    C[idx][ntri * (z - 1) + face->cell_id] = alpha[4];
                }
            }
        }
    }

    viennacl::compressed_matrix<double>  vl_C(ntri * nLayer, ntri * nLayer);
    viennacl::copy(C,vl_C);
    viennacl::vector<double> rhs;
    viennacl::copy(b,rhs);

    LOG_DEBUG << "solving";
//    arma::sp_mat result(ntri * nLayer, ntri * nLayer );
//    for(int i = 0; i < omp_get_max_threads(); i++ )
//    {
//        result += sp_C.at(i);
//    }

    // configuration of preconditioner:
//    viennacl::linalg::chow_patel_tag chow_patel_ilu_config;
//    chow_patel_ilu_config.sweeps(3);       // three nonlinear sweeps
//    chow_patel_ilu_config.jacobi_iters(2); // two Jacobi iterations per triangular 'solve' Rx=r
//// create and compute preconditioner:
//    viennacl::linalg::chow_patel_ilu_precond< viennacl::compressed_matrix<double> > chow_patel_ilu(vl_C, chow_patel_ilu_config);

//    std::vector<double> x(ntri * nLayer ,0.0);
    viennacl::vector<double> x = viennacl::linalg::solve(vl_C, rhs, viennacl::linalg::bicgstab_tag());

    LOG_DEBUG << "done";

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
//            if (d->face_neigh[j])
//            {
//                auto neigh = face->neighbor(j);
                auto emp = face->edge_midpoint(j);
                auto query = boost::make_tuple(emp.x(), emp.y(), 0.0); //z isn't used in the interp call

                double qs  = interp(vec_qs, query);
                double qt  = interp(vec_qt, query);

                qdep += E[j]*udotm[j]*(qs+qt);
//            }
//            else
//            {
//                //use the linear average if we
//                qdep += .5000000000*E[j]*udotm[j]*face->face_data("Qsalt")+.5000000000*E[j]*udotm[j]*face->face_data("Qsusp")+
//                        .5000000000*E[j]*face->face_data("Qsalt")*udotm[j]+.5000000000*E[j]*face->face_data("Qsusp")*udotm[j];
//            }
        }



        qdep = -qdep * global_param->dt() / 400. / face->get_area();  // 1000 kg/m^3 -> density of water

        face->set_face_data("Qdep",qdep);
        face->set_face_data("sum_Qdep", face->face_data("sum_Qdep") + qdep);
    }
}

Lehning_blowing_snow::~Lehning_blowing_snow()
{

}
