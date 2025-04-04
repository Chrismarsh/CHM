
#include "submodule_base.hpp"

class Crack : submodule_base
{
public:
    struct info;
    Crack(double& _major,double& _min_swe_to_freeze,unsigned int& _infDays,bool& _AllowPriorInf,double& _lenstemp,double& _steps_per_day,info& _d);
    ~Crack() {};

    virtual void run() override;
    void init_inputs(double _snowmelt,double _rainfall, double _swe, double _soil_storage_at_freeze, double _airtemp, bool _newday);
   
    double runoff = 0.0;
    double get_runoff() { return runoff; };
    double melt_runoff = 0.0;
    double get_melt_runoff() { return melt_runoff; };
    double inf = 0.0;
    double get_inf() { return inf; };
    double snow_inf = 0.0;
    double get_snow_inf() { return snow_inf; };
    double rain_on_snow = 0.0;
    double get_rain_on_snow() { return rain_on_snow; };

    double snowmelt;
    double rainfall;
    double swe;
    double soil_storage_at_freeze;
    double airtemp;
    bool is_newday;
        
    const double& major;
    const double& min_swe_to_freeze;
    const unsigned int& infDays;
    const bool& AllowPriorInf;
    const double& lenstemp;
    const double& steps_per_day;
    
    struct info
    {
        bool frozen;
        unsigned int major_melt_count;
        double index;
        double max_major_per_melt;
        double init_SWE;
        double daily_melt_total;
        double daily_rain_total;
        bool current_day_is_major;
        double tmax;
        double current_inf;
        double current_snow_inf;
        double current_runoff;
        double current_melt_runoff;
        double yesterday_melt;
        
        void init()
        {
            frozen = false;
            major_melt_count = 0;
            index = 0.0;
            max_major_per_melt = 0.0;
            init_SWE = 0.0;
            daily_melt_total = 0.0;
            current_day_is_major = false;
            tmax = 0.0;
            current_inf = 0.0;
            current_snow_inf = 0.0;
            current_runoff = 0.0;
            current_melt_runoff = 0.0;
            yesterday_melt = 0.0;
        };

        void begin_freeze()
        {
            frozen = true;
            index = 0.0;
            max_major_per_melt = 0.0;
            init_SWE = 0.0;
            current_inf = 0.0;
            current_runoff = 0.0;
            yesterday_melt = 0.0;
        };

        void end_freeze()
        {
            frozen = false;
            major_melt_count = 0;
        };

        void increment_daily_melt(const double& snowmelt);
        void increment_daily_rain(const double& rain);

    };
    info& d;

    void melt_to_infil();

    // Crack Functions
    void Calc_Index();
    void Calc_Actual_Inf();
    void Check_for_ice_lens(); 
    void increment_major_count();
    bool is_first_major();
    bool is_limited_phase();
    bool is_prior_first_major();
    bool is_major_melt();
    void daily_melt_increment();
    void update_t_max();


};


