#include "Lehning_blowing_snow.hpp" 


Lehning_blowing_snow::Lehning_blowing_snow(config_file cfg)
        :module_base(parallel::domain)
{
    depends("U_2m_above_srf");
    depends("vw_dir");
    depends("t");

    provides("Lehning_blowing_snow_pq");
    provides("q_t");
    provides("sum_q_t");

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
    int nLayer = 5;
    //needed for linear system offsets
    size_t ntri = domain->number_of_faces();

//    arma::sp_mat C(ntri * nLayer + ntri, ntri * nLayer + ntri) ;
    arma::mat C(ntri * nLayer + ntri, ntri * nLayer + ntri, arma::fill::zeros) ;
    arma::vec b(ntri * nLayer + ntri, arma::fill::zeros);

    //concentration of snow at reference height in saltation layer
    //Pomeroy et al 1993, p.169
    double c_salt = 0.8; // kg/m^3

    double z0 = 0.01; //m

    double susp_depth = 5; //5m as per pomeroy

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

//        auto fx = [=](double z)
//        {
//            //calculate wind at the height z
//            //TODO: should this take into account snowdepth?
//            double u_z = u2;// Atmosphere::log_scale_wind(u2, 2, z, 0);
//            double g= 9.81;
//            double alpha = 25 * M_PI/180.;
//            double k = PhysConst::kappa * PhysConst::kappa;
//            double new_z= z0 + 2.4025 * pow(u_z,2.)* k * pow(cos(alpha),2.) / ( pow(log(z/z0),2.) *g);
//
//            return new_z;
//        };
//
//        double guess = 1;
//        double min = 0.1;
//        double max = 5;
//        size_t max_iter = 50;
//        auto tol = [=](double a, double b)
//        {
//            return fabs(a-b) < 0.01;
//        };
//        //depth of saltation layer
//        auto z =  boost::math::tools::toms748_solve(fx,min,max,tol,max_iter);
//        double hs = (z.first + z.second)/2.;
//
        //depth of saltation layer
        double hs = z0 + 2.4025 * pow(u2,2.)* pow(PhysConst::kappa,2.) * pow(cos(25 * M_PI/180.),2.) / ( pow(log(2/z0),2.) *9.81);



        // iterate over the 5 vertical layers, each 1m in height, up to 5m
        // this is the height Pomeroy found for suspension layer to matter
        for(int z = 0; z < nLayer;++z)
        {

            //height of the suspension layer
            double cz = z + hs ; //cell height

            //compute new U_z
            double u_z = Atmosphere::log_scale_wind(u2, 2, cz, 0);
            double ustar = u_z*PhysConst::kappa/log(cz/z0);
            double l_max = 1;
            double l = (PhysConst::kappa * (cz+z0) * l_max)/(PhysConst::kappa * cz+ PhysConst::kappa * z0+l_max);

            double K = ustar * l;

            // holds A_ * K_ / h_
            // _0 -> _2 are the horizontal sides
            // _3 -> is the top of the prism
            // _4 -> is the bottom the prism
            double alpha[5];

            //do calculation of alpha for the sides
            for(int a = 0; a < 3; ++a)
            {
                alpha[a] = A[a] * K;
                auto neigh = face->neighbor(a);

                //if we have a neighbour, use the distance, otherwise assumed h = 1 (this reasonable?)
                if (neigh != nullptr)
                {
                    alpha[a] /= math::gis::distance(face->center(), neigh->center());
                }
            }

            //top
            alpha[3] = A[3] * K / (susp_depth / nLayer);
            //bottom
            alpha[4] = A[4] * K / (susp_depth / nLayer);



            arma::vec uvw(3);
            uvw(0) = u(0);
            uvw(1) = u(1);
            uvw(2) = 0;
            //calculate u_z at this height in the suspension layer
            uvw *= Atmosphere::log_scale_wind(u2, 2, cz, 0);

            //we are the bottom?
            if(z == 0)
            {
                double temp =
                        -4-2*A[4]*K/hs-(1/2)*A[0]*arma::dot(uvw, m[0])-(1/2)*A[1]*arma::dot(uvw, m[1])-(1/2)*A[2]*arma::dot(uvw, m[2])-(1/2)*A[3]*arma::dot(uvw, m[3])-(1/2)*A[4]*arma::dot(uvw, m[4]);
                C(i, ntri * z + face->cell_id) = temp;
                //value for csalt, robin BC, bottom boundary
                //ntri * layer gets us in the offsets for the robin BC
                C(i,  ntri * nLayer + face->cell_id) = 2*A[4]*K/hs;

                //do edges
                for(int j = 0; j<3; ++j)
                {
                    auto neigh = face->neighbor(0);
                    if (neigh != nullptr)
                    {
                        C(i,  ntri * z + neigh->cell_id) = alpha[j] - 0.5 * A[j]*arma::dot(uvw, m[j]);
                    }
                    //0 flux BC if not
                }

                //top
                C(i,  ntri * (z+1) + face->cell_id) = alpha[3] - 0.5 * A[3]*arma::dot(uvw, m[3]); //top

            }
            else if(z == nLayer -1)// are we the top?
            {
                C(i, ntri * z + face->cell_id) =
                        -A[3]*K-2-(1/2)*A[0]*arma::dot(uvw, m[0])-(1/2)*A[1]*arma::dot(uvw, m[1])-(1/2)*A[2]*arma::dot(uvw, m[2])-(1/2)*A[3]*arma::dot(uvw, m[3])-(1/2)*A[4]*arma::dot(uvw, m[4]);

                //do edges
                    for(int j = 0; j<3; ++j)
                    {
                        auto neigh = face->neighbor(0);
                        if (neigh != nullptr)
                        {
                            C(i,  ntri * z + neigh->cell_id) = alpha[j] - 0.5 * A[j]*arma::dot(uvw, m[j]);
                        }
                        //0 flux BC if not
                    }

                //bottom
                C(i,  ntri * (z-1) + face->cell_id) = alpha[4] - 0.5 * A[4]*arma::dot(uvw, m[4]);

                //0 flux top
                // nothing here //
            }
            else // we're somehwere in the middle
            {
                C(i, ntri * z + face->cell_id) =
                        -5-(1/2)*A[0]*arma::dot(uvw, m[0])-(1/2)*A[1]*arma::dot(uvw, m[1])-(1/2)*A[2]*arma::dot(uvw, m[2])-(1/2)*A[3]*arma::dot(uvw, m[3])-(1/2)*A[4]*arma::dot(uvw, m[4]);
                //do edges first
                for(int j = 0; j<3; ++j)
                {
                    auto neigh = face->neighbor(0);
                    if (neigh != nullptr)
                    {
                        C(i,  ntri * z + neigh->cell_id) = alpha[j] - 0.5 * A[j]*arma::dot(uvw, m[j]);
                    }
                    //0 flux BC if not
                }

                //top
                C(i,  ntri * (z+1) + face->cell_id) = alpha[3] - 0.5 * A[3]*arma::dot(uvw, m[3]);

                //bottom
                C(i,  ntri * (z-1) + face->cell_id) = alpha[4] - 0.5 * A[4]*arma::dot(uvw, m[4]);
            }

        }
    }

//    arma::vec x = arma::spsolve(C,b);
    arma::vec x = arma::solve(C,b);
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
