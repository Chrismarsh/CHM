#include "deform_mesh.hpp"

deform_mesh::deform_mesh(config_file cfg)
        :module_base(parallel::domain)
{

}

deform_mesh::~deform_mesh()
{

}

void deform_mesh::run(mesh domain)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_vertex(); i++)
    {
        Point_3 p;

        auto vert = domain->vertex(i);
        double z = vert->point().z();

        if(z > domain->min_z())
        {
           z -= (z-domain->min_z()) * 0.25;
        }


        p = Point_3(vert->point().x(), vert->point().y(), z);
        vert->set_point(p);
    }

    domain->_terrain_deformed = true;
}
