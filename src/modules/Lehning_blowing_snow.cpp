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
    arma::sp_mat A(domain->size_faces(), domain->size_faces());
    arma::vec b(domain->size_faces(), arma::fill::zeros);

    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);

        double T = face->face_data("t");
        double phi = face->face_data("vw_dir");
        double u2 = face->face_data("U_2m_above_srf");

        double alpha_i = F/(6*face->get_area());

        arma::vec u(2);
        Vector_2 v = math::gis::bearing_to_cartesian(phi);
        u(0) = v.x(); //U_x
        u(1) = v.y(); //U_y

        arma::vec m1(2);
        m1(0) = face->edge_unit_normal(0).x();
        m1(1) = face->edge_unit_normal(0).y();


        arma::vec m2(2);
        m2(0) = face->edge_unit_normal(1).x();
        m2(1) = face->edge_unit_normal(1).y();


        arma::vec m3(2);
        m3(0) = face->edge_unit_normal(2).x();
        m3(1) = face->edge_unit_normal(2).y();


        //edge lengths
        double E1 = face->edge_length(0);
        double E2 = face->edge_length(1);
        double E3 = face->edge_length(2);

        A(i,face->cell_id) =  alpha_i*(E1*arma::dot(u,m1) + E2*arma::dot(u,m2)+E3*arma::dot(u,m3) +1/alpha_i );

        //0 flux BC
        auto neigh = face->neighbor(0);
        if(neigh != nullptr)
        {
            A(i, neigh->cell_id) = alpha_i * E1 * arma::dot(u, m1);
        }

        neigh = face->neighbor(1);
        if(neigh != nullptr)
        {
            A(i, neigh->cell_id) = alpha_i * E2 * arma::dot(u, m2);
        }

        neigh = face->neighbor(2);
        if(neigh != nullptr)
        {
            A(i, neigh->cell_id) = alpha_i * E3 * arma::dot(u, m3);
        }

        double pq = p10_dry(T, u2) * Qt(T, u2);
        b[i] = pq < 1e-6? 0 : -pq;

        face->set_face_data("Lehning_blowing_snow_pq", b[i] );
    }


    arma::vec x = arma::spsolve(A,b);

    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);

        face->set_face_data("q_t", x(i)*global_param->dt());
        face->set_face_data("sum_q_t",  face->face_data("sum_q_t") + x(i)*global_param->dt());
    }
}

Lehning_blowing_snow::~Lehning_blowing_snow()
{

}
