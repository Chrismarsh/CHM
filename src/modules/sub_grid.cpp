#include "sub_grid.hpp"

sub_grid::sub_grid(config_file cfg)
        :module_base(parallel::data)

{
    depends("snowdepthavg");

    provides("snowcoverfraction");

    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

sub_grid::~sub_grid()
{


}

void sub_grid::run(mesh_elem &face) {

    double snowdepthavg = face->face_data("snowdepthavg");
    double snowcoverfraction;

    if(snowdepthavg>0) {
        snowcoverfraction = 1;
    } else {
        snowcoverfraction = 0;
    }

    face->set_face_data("snowcoverfraction",snowcoverfraction);

}