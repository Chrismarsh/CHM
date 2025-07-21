#include "submodule_base.hpp"

template<class soil_data, auto max_infil>
class Ayers : public submodule_base
{
public:
    Ayers(const double& _rainfall,const double& _snowmelt,
            const std::string& _texture,const std::string& _ground_cover, 
            const soil_data& _soil_data); 
 
    ~Ayers() {};

    virtual void run() override;
   
    const double& rainfall;
    const double& snowmelt;
    const std::string& texture;
    const std::string& ground_cover;
    const soil_data& soils; 

    double runoff = 0.0;
    double get_runoff() { return runoff; };
    double inf = 0.0;
    double get_inf() { return inf; };
    double snow_inf = 0.0;
    double get_snow_inf() { return snow_inf; };

};

#include "Ayers.ipp"
