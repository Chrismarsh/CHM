#include "Lehning_blowing_snow.hpp" 


Lehning_blowing_snow::Lehning_blowing_snow(config_file cfg)
        :module_base(parallel::domain)
{
    depends("U_2m_above_srf");
    depends("vw_dir");
    depends("t");

    provides("c0");
    provides("q_t");
    provides("sum_q_t");
    provides("u_z");

}

void Lehning_blowing_snow::init(mesh domain)
{
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        face->set_face_data("sum_q_t",0);
    }
}

void Lehning_blowing_snow::run(mesh domain)
{
    //hardcode at the moment
    int nLayer = 2;
    //needed for linear system offsets
    size_t ntri = domain->number_of_faces();

    arma::sp_mat C(ntri * nLayer , ntri * nLayer ) ;
//    arma::mat C(ntri * nLayer , ntri * nLayer, arma::fill::zeros) ;
    arma::vec b(ntri * nLayer , arma::fill::zeros);

    //concentration of snow at reference height in saltation layer
    //Pomeroy et al 1993, p.169
    double c_salt = 0.8; // kg/m^3

    double z0 = 0.01; //m

    double susp_depth = 4; //5m as per pomeroy

    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);

        double T = face->face_data("t");
        double phi = face->face_data("vw_dir");
        double u2 = face->face_data("U_2m_above_srf");

        arma::vec u(2);
        Vector_2 v = math::gis::bearing_to_cartesian(phi);
        u(0) = v.x(); //U_x
        u(1) = v.y(); //U_y

        //edge unit normal
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

        //face areas, todo: needs to be area's b/c 3D
        double A[5];
        for(int j=0; j < 3; ++j)
            A[j] = face->edge_length(j); //assuming 1m height

        //top, bottom
        A[3] = A[4] = face->get_area();


        //depth of saltation layer
        double hs = z0 + 2.4025 * pow(u2,2.)* pow(PhysConst::kappa,2.) * pow(cos(25. * M_PI/180.),2.) / ( pow(log(2./z0),2.) *9.81);

        double ustar = std::max(0.1,u2*PhysConst::kappa/log(2./z0));

        double l_max = 1;



        // iterate over the 5 vertical layers, each 1m in height, up to 5m
        // this is the height Pomeroy found for suspension layer to matter
        for(int z = 0; z < nLayer;++z)
        {

            //height of the suspension layer
            double cz = z + hs ; //cell height

            //compute new U_z
            double u_z = std::max(0.1,Atmosphere::log_scale_wind(u2, 2, cz, 0));

            double K[5];

            //top and bottom = 5m
            K[3] = K[4] = 40.;


            // holds A_ * K_ / h_
            // _0 -> _2 are the horizontal sides
            // _3 -> is the top of the prism
            // _4 -> is the bottom the prism
            double alpha[5];

            //do calculation of alpha for the sides
            for(int a = 0; a < 3; ++a)
            {
                auto neigh = face->neighbor(a);
                alpha[a] = A[a];

                double l__max = 1; // 1m if we are edge cell

                //if we have a neighbour, use the distance
                if (neigh != nullptr)
                {
                    l__max =  math::gis::distance(face->center(), neigh->center());

                    alpha[a] /= math::gis::distance(face->center(), neigh->center());
                }
                else
                { //otherwise assume 2x the distance from the center of face to one of it's vertexes, so ghost triangle is approx the same size
                    alpha[a] /= 2.0*math::gis::distance(face->center(), face->vertex(0)->point());
                }


                double l = PhysConst::kappa*(z+z0)*l__max/(PhysConst::kappa*z+PhysConst::kappa*z0+l__max);
                K[a] = std::max(ustar * l, PhysConst::kappa * 2. * ustar);

                alpha[a] *= K[a];
            }

            //top
            alpha[3] = A[3] * K[3] / (susp_depth / nLayer);
            //bottom
            alpha[4] = A[4] * K[4] / (susp_depth / nLayer);

            arma::vec uvw(3);
            uvw(0) = u(0);
            uvw(1) = u(1);
            uvw(2) = 0;

            double length = arma::norm(uvw,2);
            double scale = u_z/length;

            //calculate u_z at this height in the suspension layer
            uvw *= scale;

//            uvw(2) = -.5;


            //holds wind dot face normal
            double udotm[5];
            for(int j=0; j < 5; ++j)
            {
                udotm[j] = arma::dot(uvw, m[j]);
            }

            double temp = 0;
            //we are the bottom?
            if(z == 0)
            {
                face->set_face_data("u_z", u_z);
                temp = -1.*alpha[0]-1.*alpha[1]-1.*alpha[2]-1.*alpha[3]+A[4]*K[4]/(.50000*hs+.50000)-.50000*A[0]*udotm[0]-.50000*A[1]*udotm[1]-.50000*A[2]*udotm[2]-.50000*A[3]*udotm[3];
                C(i, ntri * z + face->cell_id) = temp;

                b(ntri * z + face->cell_id) = A[4]*K[4]*c_salt/(.50000*hs+.50000);


                //do edges
                for(int j = 0; j<3; ++j)
                {
                    auto neigh = face->neighbor(0);
                    if (neigh != nullptr)
                    {
                        temp = alpha[j]-.50000*A[0]*udotm[0];
                        C(i,  ntri * z + neigh->cell_id) = temp;
                    }
                    //0 flux BC if not
                }

                //top
                temp = alpha[3]-.50000*A[3]*udotm[3]; //top
                C(i,  ntri * (z+1) + face->cell_id) = temp;

            }
            else if(z == nLayer -1)// top z layer
            {
                temp = -1.*A[3]*K[3]-4.-.50000*A[0]*udotm[0]-.50000*A[1]*udotm[1]-.50000*A[2]*udotm[2]-.50000*A[3]*udotm[3]-.50000*A[4]*udotm[4];
                C(i, ntri * z + face->cell_id) = temp;
                //do edges
                    for(int j = 0; j<3; ++j)
                    {
                        auto neigh = face->neighbor(0);
                        if (neigh != nullptr)
                        {
                           temp =  alpha[j]-.50000*A[0]*udotm[0];
                            C(i,  ntri * z + neigh->cell_id) = temp;
                        }
                        //0 flux BC if not
                    }

                //bottom
                temp = alpha[4]-.50000*A[4]*udotm[4];
                C(i,  ntri * (z-1) + face->cell_id) = temp;

                //0 flux top
                // nothing here
            }
            else // internal cell
            {

                //value in our cell
                temp = -1.*alpha[0]-1.*alpha[1]-1.*alpha[2]-1.*alpha[3]-1.*alpha[4]-.50000*A[0]*udotm[0]-.50000*A[1]*udotm[1]-.50000*A[2]*udotm[2]-.50000*A[3]*udotm[3]-.50000*A[4]*udotm[4];

                C(i, ntri * z + face->cell_id) = temp;
                //do edges first
                for(int j = 0; j<3; ++j)
                {
                    auto neigh = face->neighbor(0);

                    if (neigh != nullptr)
                    {
                        temp = alpha[j] - 0.5 * A[j]*udotm[j];
                        C(i,  ntri * z + neigh->cell_id) = temp;
                    }
                    //else{} 0 flux BC if not
                }

                //top
                temp = alpha[3] - 0.5 * A[3]*udotm[3];
                C(i,  ntri * (z+1) + face->cell_id) = temp;

                //bottom
                temp = alpha[4] - 0.5 * A[4]*udotm[4];
                C(i,  ntri * (z-1) + face->cell_id) = temp;
            }

        }
    }

//    LOG_DEBUG << "Cond: " << arma::cond(C);
//    arma::mat L, U, P;
//    arma::lu(L, U, P, C);
//
//    LOG_DEBUG << "U diag min/max ratio: " << std::max(U.diag() / U.diag().min() );

//    LOG_DEBUG << "saving";
////    C.save("A.dat", arma::raw_ascii);
//    C.save("A.dat");
//    LOG_DEBUG << "done";

//    BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("lol"));
//    L.reset();
//    U.reset();
//    P.reset();

    arma::vec x = arma::spsolve(C,b);
//    arma::vec x = arma::solve(C,b);
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        face->set_face_data("c0", x[i]);
//        face->set_face_data("q_t", x(i)*global_param->dt());
//        face->set_face_data("sum_q_t",  face->face_data("sum_q_t") + x(i)*global_param->dt());
    }
}

Lehning_blowing_snow::~Lehning_blowing_snow()
{

}
