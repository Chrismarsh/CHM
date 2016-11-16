#include "crop_rotation.hpp"

crop_rotation::crop_rotation(config_file cfg)
        :module_base(parallel::data)
{

}

crop_rotation::~crop_rotation()
{

}

void crop_rotation::run(mesh_elem& face)
{

    int year = global_param->year();

    if(year == 2010)
        face->set_parameter("crop", face->get_parameter("annual_crop_inventory_2010"));
    else if( year == 2011)
        face->set_parameter("crop", face->get_parameter("annual_crop_inventory_2011"));

}
