#include <concepts>
#include <gtest/gtest.h>
#include "data_base.hpp"
#include "base_step.hpp"
#include "triangulation.hpp"

template<class T>

concept SubmoduleData =  requires(T& t)
{
    {t.input1()} -> std::floating_point;
    {t.input2()} -> std::floating_point;

    {t.output(std::declval<const double>())} -> std::same_as<void>;
};

template<SubmoduleData data>
class submodule : public base_step<submodule<data>,data>
{
public:
    submodule() {};
    ~submodule() {};
    
    void execute_impl(data& d) const
    {
        auto input1 = d.input1();
        auto input2 = d.input2();
        
        auto output = input1 * input2;
        d.output(output);
    }
};

struct Cache : public cache_base
{
    double input1 = default_value<double>();
    double input2 = default_value<double>();

    double output = 0.0;
};

class data2 : public data_base<Cache>
{
public:
    data2(mesh_elem face_in,boost::shared_ptr<global> param,pt::ptree& cfg)
        : data_base(face_in,param,cfg,true) {};

    double input1()
    {
        update_value(
                [this]() -> auto& { return cache_->input1;},
                []() { return 4.0; }
                );

        return cache_->input1;
    };
    
    double input2()
    {
        update_value(
                [this]() -> auto& { return cache_->input2;},
                []() { return 2.5; }
                );

        return cache_->input1;
    };

    void output(const double t)
    {
        set_output([this]() -> auto& { return cache_->output; },t);
    };
};

class TestModule : public ::testing::Test
{
protected:
    TestModule() : d(face,param,cfg) {};
    mesh_elem face{nullptr};
    boost::shared_ptr<global> param;
    pt::ptree cfg;

    data2 d;
    submodule<data2> submodule_instance;    
    void run() 
    { 
        submodule_instance.execute(d);
    };
};

TEST_F(TestModule,SimpleSubmoduleTest)
{
    const auto& cache = d.get_cache();
    EXPECT_FALSE(cache);
    run();
    EXPECT_TRUE(cache);
    EXPECT_EQ(cache->input1,4.0);
    EXPECT_EQ(cache->input2,2.5);
    EXPECT_DOUBLE_EQ(cache->output,cache->input1 * cache->input2);
};
