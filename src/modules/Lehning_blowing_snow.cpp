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
    provides("csalt");
    provides("Qsalt");
    provides("Qdep");


    provides("K");

    provides("Qsusp");
    provides("susp_dep");
}

void Lehning_blowing_snow::init(mesh domain)
{
//domain    for (size_t i = 0; i < domain->size_faces(); i++)
//    {
//        auto face = domain->face(i);
//        face->set_face_data("sum_q_t",0);
//    }
}

void Lehning_blowing_snow::run(mesh domain)
{
    //hardcode at the moment
    double nLayer = 5;
    double susp_depth = 5; //5m as per pomeroy
    double v_edge_height = susp_depth / nLayer; //height of each vertical prism

    //needed for linear system offsets
    size_t ntri = domain->number_of_faces();

    arma::sp_mat C(ntri * nLayer, ntri * nLayer ) ;
    arma::vec b(ntri * nLayer , arma::fill::zeros);
    arma::vec x(ntri * nLayer , arma::fill::zeros);

    double z0 = 0.01; //m

    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);

        //edge unit normals
        arma::vec m[5];
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
        double A[5];
        for (int j = 0; j < 3; ++j)
            A[j] = face->edge_length(j) * v_edge_height;

        //top, bottom
        A[3] = A[4] = face->get_area();

        //get wind from the face
        double phi = face->face_data("vw_dir");
        double u2 = face->face_data("U_2m_above_srf");

        //depth of saltation layer
        double hs = z0 + 2.4025 * pow(u2, 2.) * pow(PhysConst::kappa, 2.) * pow(cos(25. * M_PI / 180.), 2.) /
                         (pow(log(2. / z0), 2.) * 9.81);
        face->set_face_data("hs", hs);

        double ustar = std::max(0.1, u2 * PhysConst::kappa / log(2. / z0));
        face->set_face_data("ustar", ustar);


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

        //which faces have neighbours? Ie, are we an edge?
        bool face_neigh[3] = {true, true, true};

        double l__max = 40; //mixing length

        for (int a = 0; a < 3; ++a)
        {
            auto neigh = face->neighbor(a);
            if (neigh == nullptr)
                face_neigh[a] = false;
        }

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
        double u_star_th = pow(tau_t_f/rho_f,0.5);
        double Qsalt = 0;
        double c_salt = 0;

        //ensure we can have snow saltation
        double rho_p = 917;
        double thresh_A = .18;
        double g = 9.81;
        double d = 0.48e-3;
        double u_star_t = thresh_A*sqrt((rho_p-rho_f)/(rho_f)*d*g); // threshold

        if( ustar > u_star_t)
        {
            Qsalt = 0.68 * rho_f * u_star_th / (9.81 *ustar) *( ustar*ustar - u_star_th*u_star_th);
            c_salt = Qsalt / std::max(0.1, Atmosphere::log_scale_wind(u2, 2, hs, 0)); // back out the concentration from flux

        }


//        if (i == 1579 )//|| i==  13065
//        {
//            Qsalt = 0.8/std::max(0.1, Atmosphere::log_scale_wind(u2, 2, hs, 0));
//            c_salt = .80;
//        } else{
//            Qsalt = 0;
//            c_salt = 0.;
//        }
//            double c_salt = nzr *exp(-1.55*(1./pow(0.05628*ustar,.544)-1./pow(cz,.544)));
        face->set_face_data("csalt", c_salt);
        face->set_face_data("Qsalt", Qsalt);

        // iterate over the vertical layers, each 1m in height
        for (int z = 0; z < nLayer; ++z)
        {
            //height in the suspension layer
            double cz = z + hs+ v_edge_height/2.; //cell center height
            double l = PhysConst::kappa * (cz + z0) * l__max / (PhysConst::kappa * cz + PhysConst::kappa * z0 + l__max);

            K[3] = K[4] = std::max(ustar * l, PhysConst::kappa * cz * ustar);

            if(z==0)
                face->set_face_data("K",K[3]);

            //top
            alpha[3] = A[3] * K[3] / v_edge_height;
            //bottom
            alpha[4] = A[4] * K[4] / v_edge_height;

//            face->set_face_data("K0", K[0]);
//            face->set_face_data("K1", K[1]);
//            face->set_face_data("K2", K[2]);
//            face->set_face_data("K3", K[3]);
//            face->set_face_data("K4", K[4]);
//
//            face->set_face_data("l_vert", l);

            Vector_2 v = -math::gis::bearing_to_cartesian(phi);

            //setup wind vector
            arma::vec uvw(3);
            uvw(0) = v.x(); //U_x
            uvw(1) = v.y(); //U_y
            uvw(2) = 0;

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

            //lateral
            size_t idx = ntri*z + face->cell_id;

            for(int f = 0; f < 3; f++)
            {
                if (face_neigh[f])
                {
                    auto nidx = ntri * z + face->neighbor(f)->cell_id;

                    if(udotm[f] > 0 )
                    {
                        C(idx,idx) += -A[f]*udotm[f]-alpha[f];
                        C(idx, nidx) += alpha[f];
                    } else
                    {
                        C(idx,idx) += -alpha[f];
                        C(idx, nidx)+= -A[f]*udotm[f]+alpha[f];
                    }

                }
                else
                {
                    C(idx, idx) += -A[f]*udotm[f]+alpha[f];
                }
            }


            // 1 layer special case
            if (nLayer == 1)
            {

                C(idx, idx) += -A[4]*K[4]+alpha[3]; // A[4]*K[4]+alpha[3];

                b(idx) += -A[4]*K[4]*c_salt; // A[4]*K[4]*c_salt;
                //top
                //0
            } else
            {
                //now do the vertical fluxes
                //we are the bottom?
                if (z == 0)
                {
                    //bottom face
                    C(idx,idx) += -A[4]*K[4];
                    b(ntri * z + face->cell_id) += -A[4]*K[4]*c_salt;

                    //top face
                    C(idx, ntri * z + face->cell_id) +=  -alpha[3];
                    C(idx, ntri * (z + 1) + face->cell_id) =  alpha[3];


                } else if (z == nLayer - 1)// top z layer
                {
                    //top
                    C(idx, ntri * z + face->cell_id) += -alpha[3]-alpha[4];

                    //bottom
                    C(idx, ntri * (z - 1) + face->cell_id) = alpha[4];

                } else // internal cell
                {
                    C(idx, idx) += -alpha[3]-alpha[4];

                    //top
                    C(idx, ntri * (z + 1) + face->cell_id) =alpha[3];

                    //bottom
                    C(idx, ntri * (z - 1) + face->cell_id) = alpha[4];
                }
            }
        }
    }

    LOG_DEBUG << "solving";

//    x = arma::solve(C,b);
//    x = arma::spsolve(C,b);

    x = viennacl::linalg::solve(C, b, viennacl::linalg::bicgstab_tag());

    LOG_DEBUG << "done";

//    arma::sp_mat Q(ntri * 2, ntri * 2 ) ;
//    arma::vec bq(ntri * 2 , arma::fill::zeros);
//    arma::vec xq(ntri * 2 , arma::fill::zeros);

    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        double Qsusp = 0;
        double hs = face->face_data("hs");
        double u2=face->face_data("U_2m_above_srf");
        for (int z = 0; z<nLayer;++z)
        {
            double c = x[ntri * z + face->cell_id];
            double cz = z + hs+ v_edge_height/2.; //cell center height

            double u_z =std::max(0.1, Atmosphere::log_scale_wind(u2, 2, cz, 0)); //compute new U_z at this height in the suspension layer
            Qsusp += c * u_z * v_edge_height;

            face->set_face_data("c"+std::to_string(z), c ); //<0?0:c
        }
        face->set_face_data("Qsusp",Qsusp);
    }

    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);

        //edge unit normals, only horz
        arma::vec m[3];
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

        double u2 = face->face_data("U_2m_above_srf");
        double phi = face->face_data("vw_dir");
        Vector_2 v = -math::gis::bearing_to_cartesian(phi);

        //setup wind vector
        arma::vec uvw(3);
        uvw(0) = v.x(); //U_x
        uvw(1) = v.y(); //U_y
        uvw(2) = 0;

        //which faces have neighbours? Ie, are we an edge?
        bool face_neigh[3] = {true, true, true};
        double udotm[3];

        //edge lengths b/c 2d now
        double E[3]={0,0,0};
        for (int j = 0; j < 3; ++j)
        {
            auto neigh = face->neighbor(j);
            if (neigh == nullptr)
                face_neigh[j] = false;

            //just unit vectors as qsusp/qsalt flux has magnitude
            udotm[j] = arma::dot(uvw, m[j]);

            E[j] = face->edge_length(j);// * v_edge_height;
        }



        double qdep = 0;
        for(int j = 0; j < 3; j++)
        {
            if (face_neigh[j])
            {
                auto neigh = face->neighbor(j);
//                if(udotm[j] > 0)
//                {
//                    qdep += face->face_data("Qsalt")*E[j]*udotm[j]+face->face_data("Qsusp")*E[j]*udotm[j];
//                }
//                else
//                {
//                    qdep += E[j]*neigh->face_data("Qsalt")*udotm[j]+E[j]*neigh->face_data("Qsusp")*udotm[j];
//                }


                qdep += .5000000000*E[j]*udotm[j]*face->face_data("Qsalt")+.5000000000*E[j]*udotm[j]*face->face_data("Qsusp")+
                        .5000000000*E[j]*neigh->face_data("Qsalt")*udotm[j]+.5000000000*E[j]*neigh->face_data("Qsusp")*udotm[j];
            }
            else
            {
                qdep += .5000000000*E[j]*udotm[j]*face->face_data("Qsalt")+.5000000000*E[j]*udotm[j]*face->face_data("Qsusp")+
                        .5000000000*E[j]*face->face_data("Qsalt")*udotm[j]+.5000000000*E[j]*face->face_data("Qsusp")*udotm[j];
            }
        }
        qdep = -qdep/ 1000. * global_param->dt();
        face->set_face_data("Qdep",qdep);

    }
}

Lehning_blowing_snow::~Lehning_blowing_snow()
{

}
