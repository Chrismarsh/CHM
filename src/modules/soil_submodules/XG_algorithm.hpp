// This code was written as an exact copy of the CRHM module ClassXG
// Therefore, any bugs in ClassXG must also be present here
// Tests were written by examining the original code to 'guess' at behaviour
#include "I_freeze_thaw_depths.hpp"
#include <cmath>
#include <vector>

class XG_algorithm : public I_freeze_thaw_depths
{
public:
    class state;
    class params;
    XG_algorithm(const double& t, const double& _moist, 
            const double& _rechr, state& _S, const params& _P);
    
    ~XG_algorithm() {};

    virtual void run() override;
	void init_freezethaw_degreedays(const double& Zdf_init, const double& Zdt_init, const double& Zpf_init);
    
    static constexpr double ko = 0.21;  // W/(m K) organic material
    static constexpr double km = 2.50;  // W/(m K) mineral
    static constexpr double ka = 0.025; // W/(m K) air
    static constexpr double ki = 2.24;  // W/(m K) ice  2.24
    static constexpr double kw = 0.57;  // W/(m K) water   0.57

private:
    const double& surface_temp;
    const double& soil_moist;
    const double& soil_rechr;
    state& S;
    const params& P;
	std::vector<double> depths;    
    
    void freeze(void);
    void thaw(void);
    double stefan_equation(double& surface_index, 
            double& thermal_conductivity, size_t& layer);
    double Interpolated_ttc(double Za, size_t layer);
    double Interpolated_ftc(double Za, size_t layer);
    void push_front(double D);
    void find_thaw_D(double dt);
    void find_freeze_D(double df);

public:
    class state
    {
    public:
        // Depth-related
        double Zdf;                // (m) depth of freezing front  
        double Zdt;                // (m) depth of thawing front  
        std::vector<double> Zd_front; // (m) depths of all freezing/thawing fronts (thaw: +, freeze: -)  
        size_t Th_low;            // lowest thawed layer  
        size_t Fz_low;            // lowest frozen layer  
        size_t nfront;                // number of freezing/thawing fronts  

        // Degree-day counters  
        double Bfr;                // (ºC*d) freeze degree days  
        double Bth;                // (ºC*d) thaw degree days  

        // Layer ratios  
        std::vector<double> pf;    // soil layers freezing ratios  
        std::vector<double> pt;    // soil layers thawing ratios  

        // Thermal properties  
        std::vector<double> ttc;   // (W/(m*K)) thawing thermal conductivity  
        std::vector<double> ftc;   // (W/(m*K)) freezing thermal conductivity  
        std::vector<double> tc_composite;  // (W/(m*K)) composite ftc/ttc  
        std::vector<double> tc_composite2; // (W/(m*K)) composite2 ftc/ttc  

        // Moisture/water content  
        std::vector<double> theta;           // (m³/m³) layer theta  
        std::vector<double> layer_h2o;       // (kg/m³) layer water content  
        double XG_moist_d;      // (m) layer depth of soil moisture in XG  
        double XG_rechr_d;      // (m) layer depth of soil recharge in XG  
        std::vector<double> XG_max;          // (mm) layer max soil moisture  
        std::vector<double> XG_moist;        // (mm) layer moisture content  
        double check_XG_moist;     // (mm) sum of XG soil_moist  

        // Local/temporary  
        double B;                  // (ºC*d) interval degree-day sum  
        double TrigAcc;            // (ºC*d) freeze/thaw cycle detection  
        int TrigState;             // 1/0/-1 → thaw/idle/freeze  
        std::vector<size_t> ttc_contents;          // 0/1 → thaw/freeze  
        std::vector<size_t> ftc_contents;          // 0/1 → freeze/thaw  
        double t_trend;            // (°C) temperature long-term trend  

        bool is_newday = false;
        bool last_step_new_day = false;
        
		// Fractions  
        std::vector<double> rechr_fract;   // fraction of layer (soil_rechr_max)  
        std::vector<double> moist_fract;   // fraction of layer (soil_moist_max)  
        std::vector<double> default_fract; // fraction of layer (theta_default)  
        // Constructor that initializes vector sizes
        //
       
        // TODO: refactor P as a unique pointer which is created by the module, then ownership is moved to B 
        const params& P;

        explicit state(size_t N,const params& _P) :
            Zd_front(N,0.0),
            pf(N,0.0),
            pt(N,0.0),
            ttc(N,0.0),
            ftc(N,0.0),
            tc_composite(N,0.0),
            tc_composite2(N,0.0),
            theta(N,0.0),
            layer_h2o(N,0.0),
            XG_max(N,0.0),
            XG_moist(N,0.0),
            ttc_contents(N,0.0),
            ftc_contents(N,0.0),
            rechr_fract(N,0.0),
            moist_fract(N,0.0),
            default_fract(N,0.0),
            P(_P)
        {
            assert(N > 0 && "State requires a positive number of layers");
            // Initialize other members
            Zdf = 0.0;
            Zdt = 0.0;
            Th_low = 1;
            Fz_low = 1;
            nfront = 0;
            Bfr = 0.0;
            Bth = 0.0;
            XG_moist_d = 0.0;
            XG_rechr_d = 0.0;
            check_XG_moist = 0.0;
            B = 0.0;
            TrigAcc = 0.0;
            TrigState = 0;
            t_trend = 0.0;
        }
        double get_ftc(size_t layer, const params& P);
        double get_ttc(size_t layer, const params& P);

        //void set_XG_max(params& P);
        //void set_theta(params& P);
        state& set_layer_moisture_maximums(const params& P);
        state& set_thermal_conductivities(const params& P, 
                double const& soil_moist, double const& soil_rechr);
        state& set_freezethaw_ratios(const params& P);
		 
        bool freezing();
        bool thawing();
        bool net_negative_degree_days();
        void reset_degree_day_counter();
        void determine_freeze_thaw_idle();
        bool thaw_front_overtaken();
        bool freeze_front_overtaken();
        void merge_freeze_fronts(XG_algorithm* XG);
        void merge_thaw_fronts(XG_algorithm* XG);
        bool exist_excess_fronts();
        void push_front(double D);
        void accumulate_degree_days(const double& surface_temp);
        double last_front(void);
        double pop_front(void);
        
    };

    class params
    {
	private:
        std::vector<double> depths;  // (m) soil layer thicknesses  
    public:
    // Core parameters  
        double Trigthrhld;           // (ºC*d) Trigger reference level 
		const std::vector<double>& getdepths() const 
		{ return depths; };
        std::vector<double> por;     // soil porosity  
        size_t N_Soil_layers;           // number of soil layers (≤ nlay)  
        std::vector<double> theta_default;        // (m³/m³) default theta  
        double theta_min;            // (m³/m³) minimum theta  
        std::vector<double> soil_solid_km;        // (W/(m*K)) dry soil conductivity  
        std::vector<double> soil_solid_km_ki;     // (W/(m*K)) saturated frozen conductivity  
        std::vector<double> soil_solid_km_kw;     // (W/(m*K)) saturated unfrozen conductivity  
        double SWE_k;                // (W/(m*K)) snow thermal conductivity UNUSED IN CRHM 
        //Next 3 are initial conditions
        //const double Zdf_init;             // (m) initial freezing front depth  
        //const double Zdt_init;             // (m) initial thawing front depth  
        const double Zpf_init;             // (m) initial permafrost depth  
        bool freeze_kw_ki_update;     // update kw→ki behind freeze front  
        bool thaw_ki_kw_update;       // update ki→kw behind thaw front  
        size_t k_update;                // 0=never, 1=post-layer, 2=continuous  
        double soil_rechr_max;       // (mm) max recharge zone capacity  
        double soil_moist_max;       // (mm) max rooting zone capacity  
        double time_step_per_day;
        bool calc_conductivity;
        bool is_crhm_test = false;
        ~params() {}; 
        params(
            const std::vector<double> d, const double& t, const std::vector<double> p,
            const size_t& n, const std::vector<double> td, const double& tm, const std::vector<double> skm,
            const std::vector<double> ski, const std::vector<double> skw, const double swk, const double zpf, 
            const size_t& fku, const size_t& tku, const size_t& ku, const double& srm, const double& smm, 
            const double& tspd, const bool& cc
        ) : 
            depths(d), 
            Trigthrhld(t), 
            por(p), 
            N_Soil_layers(n), 
            theta_default(td),
            theta_min(tm), 
            soil_solid_km(skm), 
            soil_solid_km_ki(ski),
            soil_solid_km_kw(skw), 
            SWE_k(swk), 
            Zpf_init(zpf), 
            freeze_kw_ki_update(fku), 
            thaw_ki_kw_update(tku), 
            k_update(ku), 
            soil_rechr_max(srm), 
            soil_moist_max(smm), 
            time_step_per_day(tspd), 
            calc_conductivity(cc)
        {
            //CHM
        };
    };
    // Variation #1 parameters  
    //const double& n_factor_a;           // surface-to-air temp ratio  
    //const double& n_factor_b;           // surface-to-air temp ratio  
    //const double& n_factor_c;           // surface-to-air temp ratio  
};

//class StateBuilder {
//public:
//    explicit StateBuilder(int num_layers) {
//        state_ = std::make_unique<state>();
//        initialize_vectors(num_layers);
//    }
//
//    StateBuilder& set_initial_freezethaw_depths(const double Zdf, const double Zdt)
//    {
//        state_->Zdf = Zdf;
//        state_->Zdt = Zdt;
//
//        return *this;
//    };

//XG_algorithm::state& XG_algorithm::state::check_initial_Zdf_depths(const std::vector<double>& depths)
//{
//    double sum = 0.0;
//    for (double val : depths)
//    {
//        sum += val;
//    }
//
//    if (sum < this->Zdf || sum < this->Zdt)
//    {
//        // TODO add exception
//    }
//
//    return *this;
//};

//XG_algorithm::state& XG_algorithm::state::set_XG_max(std::vector<double> por,std::vector<double> depths) 
//{
//
//    for (auto&& [m,p,d] : std::views::zip(state_->XG_max,por,depths))
//        m = p * d * 1000.0;
//
//    return *this;
//
//};
//
//XG_algorithm::state& XG_algorithm::state::set_theta(std::vector<double> theta_default) 
//{
//    state_->theta = theta_default;
//
//    return *this;
//};




