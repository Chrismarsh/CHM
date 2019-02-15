//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//

#include "Winstral_index.hpp"
REGISTER_MODULE_CPP(Winstral_index);

Winstral_index::Winstral_index(config_file cfg)
        : module_base("Winstral_index", parallel::domain, cfg)

{

    depends("vw_dir");
    depends("snowdepthavg");

    provides("Sx");

    //max distance to search
    dmax = cfg.get("dmax",300.0);

    //size of the step to take
    size_of_step = cfg.get("size_of_step",30);

    //number of steps along the search vector to check for a higher point
    //steps = cfg.get("steps",10);
    steps = dmax / size_of_step;

    //height parameter to accound for instrument height or the impact of small terrain perturbation on Sx 
    // see Winstral et al. (2013) for me details
    height_param = cfg.get("height_param",0.0);

    //separation distance to compute drift separators
    //sepdist = cfg.get("dmax",60.0);

    angular_window = cfg.get("angular_window",30.);
    delta_angle = cfg.get("delta_angle",5.);
    nangle = angular_window / delta_angle+1;

    // Logical to include vegetation in the computation of Sx
    incl_veg = cfg.get("incl_veg",false);

    // Logical to include snow depth in the computation of Sx
    incl_veg = cfg.get("incl_snw",false);

    // Option to compute the elevation of the point considered to compute Sx
    use_subgridz = cfg.get("use_subgridz",true);


    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

void Winstral_index::run(mesh domain)
{


  #pragma omp parallel for
  for (size_t ii = 0; ii < domain->size_faces(); ii++)
  {
    auto face = domain->face(ii);

    // Reference point: center of the triangle
    Point_3 me = face->center();

    // Reference height: elevation of the center of the triangle
    double Z_loc = face->center().z()+height_param;

    if (incl_veg && face->has_vegetation())
    {
         Z_loc = Z_loc + face->veg_attribute("CanopyHeight");
    }
    if (incl_snw)
    {
         Z_loc = Z_loc + face->face_data("snowdepthavg");
    }

    double sx_mean  = 0.;

    // Extract Wind direction
    double wind_dir = face->face_data("vw_dir") ;

    for (int i = 1; i <= nangle; ++i)
    {
     
        double sx = -9999.0;
        double max_tan_sx = 0.;

        //direction it is from,i need upwind fetch
        double wdir = wind_dir- angular_window/2+ (i-1)*delta_angle;

       // search along wdir azimuth in j step increments
        for (int j = 1; j <= steps; ++j)
        {
           double distance = j * size_of_step;

           // Select point along the line
           Point_2 pref =  math::gis::point_from_bearing(me,wind_dir,distance);
           // Find corresponding triangle
           auto f = domain->find_closest_face (pref );

           double Z_dist = 0.;
           if(use_subgridz)
           {
              Z_dist = f->get_subgrid_z(pref);
           }
           else
           {
              Z_dist = f->center().z();
           } 
            
           if (incl_veg && f->has_vegetation())
           {
               Z_dist = Z_dist + f->veg_attribute("CanopyHeight");
            }

           if (incl_snw)
           {
               Z_dist = Z_dist+ f->face_data("snowdepthavg");
           }

           double tan_sx = (Z_dist-Z_loc) / distance;
           if(std::abs(tan_sx) > std::abs(max_tan_sx))
           {
              max_tan_sx = tan_sx;
           }
        }

        sx = atan(max_tan_sx);
        sx_mean = sx_mean+sx; 

     }

    // Derive Sx averaged over the angular windows 
    sx_mean = sx_mean/nangle;
   
    face->set_face_data("Sx", sx_mean);    
   }

}
Winstral_index::~Winstral_index()
{

}




