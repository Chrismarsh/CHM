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



#include "triangulation.hpp"
#include "gtest/gtest.h"
#include "readjson.hpp"
#include <boost/property_tree/ptree.hpp>

struct test_module_data : face_info
{
    double x;
    double y;
    std::string ID;
};

class TriangulationTest : public testing::Test
{
  protected:

    virtual void SetUp()
    {
        logging::core::get()->set_logging_enabled(false);
        mesh_json = read_json("meshes/granger1m.mesh");
        param_json = read_json("meshes/granger1m.param");
        ic_json = read_json("meshes/granger1m.ic");

        for(auto& ktr : param_json)
        {
            //use put to ensure there are no duplciate parameters...
            std::string key = ktr.first.data();
            mesh_json.put_child( "parameters." + key ,ktr.second);
        }

        for(auto& ktr : ic_json)
        {
            //use put to ensure there are no duplciate parameters...
            std::string key = ktr.first.data();
            mesh_json.put_child( "initial_conditions." + key ,ktr.second);
        }
        mesh.from_json(mesh_json);
        mesh.init_timeseries(variables);

    }

    pt::ptree mesh_json;
    pt::ptree param_json;
    pt::ptree ic_json;

    std::set<std::string> variables = {"t","rh","u"};
    std::set<std::string> modules = {"module_1","module_42"};

    triangulation mesh;
};


//tests for proper default init

TEST_F(TriangulationTest, DefaultInit)
{
    triangulation mesh;
    ASSERT_NO_THROW(mesh.from_json(mesh_json));


}
TEST_F(TriangulationTest, DefaultInitWithVars)
{
    triangulation mesh;
    ASSERT_NO_THROW(mesh.from_json(mesh_json));
    ASSERT_NO_THROW(mesh.init_timeseries(variables));

}


TEST_F(TriangulationTest, VarReadValue)
{
    auto f = mesh.face(0);
    ASSERT_EQ((*f)["t"],-9999.0);
}

TEST_F(TriangulationTest, VarReadWriteValue)
{
    triangulation mesh;
    ASSERT_NO_THROW(mesh.from_json(mesh_json));
    ASSERT_NO_THROW(mesh.init_timeseries(variables));

    auto f = mesh.face(0);
    ASSERT_EQ((*f)["t"],-9999.0);

    (*f)["t"] = -10;
    (*f)["rh"] = 80;
    (*f)["u"] = 15.0;

    ASSERT_EQ((*f)["t"],-10.0);
    ASSERT_EQ((*f)["rh"],80.0);
    ASSERT_EQ((*f)["u"],15.0);
}

TEST_F(TriangulationTest, VarReadWriteValue_s)
{
    triangulation mesh;
    ASSERT_NO_THROW(mesh.from_json(mesh_json));
    ASSERT_NO_THROW(mesh.init_timeseries(variables));

    auto f = mesh.face(0);
    ASSERT_EQ((*f)["t"_s],-9999.0);

    (*f)["t"_s] = -10;
    (*f)["rh"_s] = 80;
    (*f)["u"_s] = 15.0;

    ASSERT_EQ((*f)["t"_s],-10.0);
    ASSERT_EQ((*f)["rh"_s],80.0);
    ASSERT_EQ((*f)["u"_s],15.0);
}

TEST_F(TriangulationTest, ParamReadValue)
{
    auto f = mesh.face(0);
    ASSERT_DOUBLE_EQ(f->parameter("MS0"),0.972731475402661);

    auto f2 = mesh.face(1);
    ASSERT_DOUBLE_EQ(f2->parameter("MS0"),0.984907757406954);
}

TEST_F(TriangulationTest, ParamReadWriteValue)
{
    triangulation mesh;
    ASSERT_NO_THROW(mesh.from_json(mesh_json));
    ASSERT_NO_THROW(mesh.init_timeseries(variables));

    auto f = mesh.face(0);
    ASSERT_DOUBLE_EQ(f->parameter("MS0"),0.972731475402661);


    auto f2 = mesh.face(1);
    ASSERT_DOUBLE_EQ(f2->parameter("MS0"),0.984907757406954);
}

TEST_F(TriangulationTest, ModuleData)
{
    triangulation mesh;
    ASSERT_NO_THROW(mesh.from_json(mesh_json));
    ASSERT_NO_THROW(mesh.init_timeseries(variables));
    ASSERT_NO_THROW(mesh.init_module_data(modules));

    auto f = mesh.face(0);

    auto& d = f->make_module_data<test_module_data>("module_42");
//    ASSERT_NE(d,nullptr);

    d.x = -100;
    ASSERT_DOUBLE_EQ(d.x,-100.0);

    // should be able to call make data again and just get back what we have
    auto& d2 = f->make_module_data<test_module_data>("module_42");
    ASSERT_DOUBLE_EQ(d2.x,-100.0);

    //retrieve the data
    auto& d3 = f->get_module_data<test_module_data>("module_42");
    ASSERT_DOUBLE_EQ(d3.x,-100.0);



}