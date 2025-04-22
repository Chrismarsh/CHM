#pragma once

template<class soil_data, auto max_infil>
Ayers<soil_data,max_infil>::Ayers(const double& _rainfall,
        const double& _snowmelt, const std::string& _texture,
        const std::string& _ground_cover, const soil_data& _soil_data) : 
        rainfall(_rainfall), 
        snowmelt(_snowmelt), 
        texture(_texture), 
        ground_cover(_ground_cover), 
        soils(_soil_data) {};

template<class soil_data, auto max_infil>
void Ayers<soil_data,max_infil>::run()
{   
    if (rainfall == 0.0 && snowmelt == 0.0) 
        return;
    
    if (rainfall > 0.0)
    {
        double maxinfil = (soils.*max_infil)(texture,ground_cover);
        inf = std::min(maxinfil,rainfall);
        runoff = rainfall - inf;
        if (runoff < 1e-12)
            runoff = 0.0;

    }

    if (snowmelt > 0.0)
    {
        inf += snowmelt;
        snow_inf = snowmelt; 
    } 
    
}
