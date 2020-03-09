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
    auto f = mesh.face(0);
    ASSERT_EQ((*f)["t"],-9999.0);

    (*f)["t"] = -10;
    (*f)["rh"] = 80;
    (*f)["u"] = 15.0;

    ASSERT_EQ((*f)["t"],-10.0);
    ASSERT_EQ((*f)["rh"],80.0);
    ASSERT_EQ((*f)["u"],15.0);
}