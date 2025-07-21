#include "Soil.h"
#include <iostream>

namespace Soil
{
    double _soils_base::lookup(const mymap& map, const std::string key) const
    {
        auto it = map.find(key);
        if (it != map.end())
        {
           return it->second;
        }
        else
        {

           CHM_THROW_EXCEPTION(module_error, "Soil Type does not exist in map"); 
        } 
    }

    soils_na::soils_na()
    {
        _make_hash();
    }

    soils_na::~soils_na()
    {

    }

    double soils_na::porosity(std::string soil_type) const
    {
        return lookup(_porosity,soil_type);
    }

    double soils_na::pore_size_dist(std::string soil_type) const
    {
        return lookup(_pore_size_dist,soil_type);//_pore_size_dist[soil_type];
    }

    double soils_na::wilt_point(std::string soil_type) const
    {
        return lookup(_wilt_point,soil_type); //_wilt_point[soil_type];
    }

    double soils_na::air_entry_tension(std::string soil_type) const
    {
        return lookup(_air_entry_tension,soil_type); //_air_entry_tension[soil_type];
    }

    double soils_na::capillary_suction(std::string soil_type) const
    {
        return lookup(_capillary_suction,soil_type);//_capillary_suction[soil_type];
    }

    double soils_na::saturated_conductivity(std::string soil_type) const
    {
        return lookup(_saturated_conductivity,soil_type);
    }
    
    double soils_na::ayers_texture(std::string texture, std::string ground_cover) const
    {
        auto it = _ayers_texture.find(texture);
        if (it != _ayers_texture.end())
        {
            return lookup(it->second, ground_cover);
        }
        else
            CHM_THROW_EXCEPTION(module_error, "Soil Type does not exist in map"); 

    }

    void soils_na::check_map()
    {
        print_ayers_texture(_ayers_texture);
    }

    void soils_na::print_ayers_texture(const auto& texture_map) {
        for (const auto& [texture, inner_map] : texture_map) {
            for (const auto& [ground_cover, value] : inner_map) {
                std::cout << texture << " -> " 
                          << ground_cover << " = " 
                          << value << "\n";
            }
        }
    }

    void soils_na::_make_hash()
    {
        // illegal to not include this step for the dense_hash_map
        _porosity.set_empty_key("");
        _pore_size_dist.set_empty_key("");
        _wilt_point.set_empty_key("");
        _air_entry_tension.set_empty_key("");
        _capillary_suction.set_empty_key("");
        _saturated_conductivity.set_empty_key("");
        _ayers_texture.set_empty_key("");

        std::string t1 = "sand";
        std::string t2 = "loam";
        std::string t3 = "clay";
        std::string t4 = "loamsand";
        std::string t5 = "sandloam";
        std::string t6 = "siltloam";
        std::string t7 = "saclloam";
        std::string t8 = "clayloam";
        std::string t9 = "siclloam";
        std::string t10 = "sandclay";
        std::string t11 = "siltclay"; 

        _porosity[t1] = 0.395;
        _porosity[t2] = 0.451;
        _porosity[t3] = 0.482;
        _porosity[t4] = 0.410;
        _porosity[t5] = 0.435;
        _porosity[t6] = 0.485;
        _porosity[t7] = 0.420;
        _porosity[t8] = 0.476;
        _porosity[t9] = 0.477;
        _porosity[t10] = 0.426;
        _porosity[t11] = 0.492;

        _pore_size_dist[t1] = 4.05;
        _pore_size_dist[t2] = 5.39;
        _pore_size_dist[t3] = 11.4;
        _pore_size_dist[t4] = 4.38;
        _pore_size_dist[t5] = 4.9;
        _pore_size_dist[t6] = 5.3;
        _pore_size_dist[t7] = 7.12;
        _pore_size_dist[t8] = 8.52;
        _pore_size_dist[t9] = 7.75;
        _pore_size_dist[t10] = 10.4;
        _pore_size_dist[t11] = 11.4;

        _wilt_point[t1] = 0.040;
        _wilt_point[t2] = 0.110;
        _wilt_point[t3] = 0.215;
        _wilt_point[t4] = 0.060;
        _wilt_point[t5] = 0.100;
        _wilt_point[t6] = 0.130;
        _wilt_point[t7] = 0.140;
        _wilt_point[t8] = 0.150;
        _wilt_point[t9] = 0.190;
        _wilt_point[t10] = 0.200;
        _wilt_point[t11] = 0.210;

        _air_entry_tension[t1] = 0.121;
        _air_entry_tension[t2] = 0.478;
        _air_entry_tension[t3] = 0.405;
        _air_entry_tension[t4] = 0.090;
        _air_entry_tension[t5] = 0.218;
        _air_entry_tension[t6] = 0.786;
        _air_entry_tension[t7] = 0.299;
        _air_entry_tension[t8] = 0.630;
        _air_entry_tension[t9] = 0.356;
        _air_entry_tension[t10] = 0.153;
        _air_entry_tension[t11] = 0.490;

        _capillary_suction[t1] = 0.4950;
        _capillary_suction[t2] = 0.8890;
        _capillary_suction[t3] = 0.3163;
        _capillary_suction[t4] = 0.6130;
        _capillary_suction[t5] = 0.1101;
        _capillary_suction[t6] = 0.1668;
        _capillary_suction[t7] = 0.2185;
        _capillary_suction[t8] = 0.2088;
        _capillary_suction[t9] = 0.2733;
        _capillary_suction[t10] = 0.2390;
        _capillary_suction[t11] = 0.2922;

        _saturated_conductivity[t1] = 117.8;
        _saturated_conductivity[t2] = 3.4;
        _saturated_conductivity[t3] = 0.3;
        _saturated_conductivity[t4] = 29.9;
        _saturated_conductivity[t5] = 10.9;
        _saturated_conductivity[t6] = 6.5;
        _saturated_conductivity[t7] = 1.5;
        _saturated_conductivity[t8] = 1.0;
        _saturated_conductivity[t9] = 1.0;
        _saturated_conductivity[t10] = 0.6;
        _saturated_conductivity[t11] = 0.5;


        std::string c1 = "bare_soil";
        std::string c2 = "row_crop";
        std::string c3 = "poor_pasture";
        std::string c4 = "small_grains";
        std::string c5 = "good_pasture";
        std::string c6 = "forested";
        google::dense_hash_map<std::string, double> inner_map;
        inner_map.set_empty_key(""); 


        // coarse_over_coarse
        inner_map[c1] = 7.6;
        inner_map[c2] = 12.7;
        inner_map[c3] = 15.2;
        inner_map[c4] = 17.8;
        inner_map[c5] = 25.4;
        inner_map[c6] = 76.2;

        _ayers_texture["coarse_over_coarse"] = inner_map;

        // medium over medium
        inner_map[c1] = 2.5;
        inner_map[c2] = 5.1;
        inner_map[c3] = 7.6;
        inner_map[c4] = 10.2;
        inner_map[c5] = 12.7;
        inner_map[c6] = 15.2;

        _ayers_texture["medium_over_medium"] = inner_map;

        // medium/fine over fine
        inner_map[c1] = 1.3;
        inner_map[c2] = 1.8;
        inner_map[c3] = 2.5;
        inner_map[c4] = 3.8;
        inner_map[c5] = 5.1;
        inner_map[c6] = 6.4;

        _ayers_texture["medium_over_fine"] = inner_map;
        _ayers_texture["fine_over_fine"] = inner_map;


        // soil over shallow bedrock
        inner_map[c1] = 0.5;
        inner_map[c2] = 0.5;
        inner_map[c3] = 0.5;
        inner_map[c4] = 0.5;
        inner_map[c5] = 0.5;
        inner_map[c6] = 0.5;

        _ayers_texture["soil_over_bedrock"] = inner_map;
    }

}
