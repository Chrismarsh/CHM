#include <gtest/gtest.h>
#include "data_base.hpp"
#include "triangulation.hpp"
#include <limits>
// Mock Cache for Testing
struct MockCache : public cache_base {
    double value = default_value<double>();
};

// Mock data class
class data : public data_base<MockCache>
{
public:
    data(mesh_elem face, boost::shared_ptr<global> param, pt::ptree& cfg)
        : data_base<MockCache>(face,param,cfg,true) {};
    ~data() {};

    // update_value is private to data and therefore must be access in a function
    // rather than use directly in the test.
    template<typename V,typename L>
    void update(V&& v, L&& l)
    {
        update_value(v, l);
    };
    
    // Same as update_value
    template<typename V,typename T>
    void output(V&& output, const T& t)
    {
        set_output(output,t);
    };

    // Again, to access private member
    double& get_value()
    {
        return cache_->value;
    };
};

// Test Fixture
class DataBaseTest : public ::testing::Test {
protected:
    void SetUp() override {
        mock_global = boost::make_shared<global>();
    }


    // mock components for the data constructor
    mesh_elem mock_face;
    boost::shared_ptr<global> mock_global;
    pt::ptree mock_cfg;

    // Lambdas for testing, standing in for face access
    static inline auto A = []() -> auto& { static double a = 42.0; return a;};
    static inline auto B = []() -> auto& { static double b = 123.0; return b;};
    

};


TEST_F(DataBaseTest, CacheInitializesOnFirstUpdate) {
    data db(mock_face, mock_global, mock_cfg);

    db.update([&]() -> auto& { return db.get_value();},A);

    EXPECT_FALSE(std::isnan(db.get_value()));
    EXPECT_EQ(db.get_value(), A());
}

TEST_F(DataBaseTest, CacheRespectsTimestepChanges) {
    data db(mock_face, mock_global, mock_cfg);
    
    // First call at timestep 0
    mock_global->timestep_counter = 0;
    db.update([&]() -> auto& {return db.get_value();} ,B);

    // Second call at same timestep - should use cached value
    db.update([&]() -> auto& {return db.get_value();} ,A);
    EXPECT_EQ(db.get_value(), B());  // value remains B. Should not update from B to A 
                                     // at the same time step

    // Force cache reset by changing timestep
    mock_global->timestep_counter = 1;
    db.update([&]() -> auto& {return db.get_value();} ,A);
    EXPECT_EQ(db.get_value(), A()); // Updates
}

TEST_F(DataBaseTest, ManualCacheResetWorks) {
    data db(mock_face, mock_global, mock_cfg);

    db.update([&]() -> auto& { return db.get_value(); },A);
    db.reset_cache();

    db.update([&]() -> auto& {return db.get_value(); },B);
    EXPECT_EQ(db.get_value(), B());
}

TEST_F(DataBaseTest, SetOutputUpdatesValue) {
    data db(mock_face, mock_global, mock_cfg);
    double output = 0.0;
    static constexpr double pi = 3.14; 
    EXPECT_FALSE(db.get_cache());
    db.output([&]() -> auto& {return output;}, pi);
    EXPECT_EQ(output, pi);
    EXPECT_TRUE(db.get_cache());
}

TEST_F(DataBaseTest, OnlyUpdatesNaNValues) {
    data db(mock_face, mock_global, mock_cfg);
    double test_value = 10.0;  // Not NaN
    auto L = [&]() -> auto& {return test_value;};
    db.update(L, A);
    EXPECT_EQ(L(), 10.0);  // Should remain unchanged
}

TEST_F(DataBaseTest, ChecksDefaultValue)
{
    double double_var = MockCache::default_value<double>();
    float float_var = MockCache::default_value<float>();
    int int_var = MockCache::default_value<int>();
    size_t size_t_var = MockCache::default_value<size_t>();
    
    EXPECT_TRUE(std::isnan(double_var));
    EXPECT_TRUE(std::isnan(float_var));
    EXPECT_TRUE(int_var == std::numeric_limits<int>::min());
    EXPECT_TRUE(size_t_var == std::numeric_limits<size_t>::min());
};
