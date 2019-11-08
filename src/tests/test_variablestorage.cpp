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

#include "variablestorage.hpp"
#include "gtest/gtest.h"

class VariableStorageTest : public testing::Test
{
  protected:

    virtual void SetUp()
    {
        logging::core::get()->set_logging_enabled(false);

        // some test variables
        variables.insert("t");
        variables.insert("rh");
        variables.insert("vw");
        variables.insert("p");
    }

    std::set< std::string> variables;
};

//basic default init sanity checks
TEST_F(VariableStorageTest, DefaultInit)
{
    variablestorage v;
    ASSERT_EQ(v.size() , 0);
    ASSERT_NO_THROW(v.variables());
    ASSERT_EQ(v.variables().size() , 0);
    ASSERT_FALSE(v.has("t"));
}

// check if the ctor correctly inits the storage pool
TEST_F(VariableStorageTest, ctorInit)
{
    variablestorage v(variables);

    // check no odd throw on var access
    ASSERT_NO_THROW(v.variables());

    // check size is correct
    ASSERT_EQ(v.size() , 4);
    ASSERT_EQ(v.variables().size() , 4);

    // check they got defaulted correctly
    ASSERT_EQ(v["t"] , -9999);
    ASSERT_EQ(v["rh"] , -9999);
    ASSERT_EQ(v["vw"] , -9999);
    ASSERT_EQ(v["p"] , -9999);

    // check they got defaulted correctly
    ASSERT_EQ(v["t"_s] ,-9999);
    ASSERT_EQ(v["rh"_s] , -9999);
    ASSERT_EQ(v["vw"_s] , -9999);
    ASSERT_EQ(v["p"_s] , -9999);
}

TEST_F(VariableStorageTest, valueAccess)
{
    variablestorage v (variables);

    v["t"] = 1;
    v["rh"] = 2;
    v["vw"] = 3;
    v["p"] = 4;

    ASSERT_EQ(v["t"] , 1);
    ASSERT_EQ(v["rh"] , 2);
    ASSERT_EQ(v["vw"] , 3);
    ASSERT_EQ(v["p"] , 4);

}

TEST_F(VariableStorageTest, valueAccess_s)
{
    variablestorage v (variables);

    v["t"_s] = 1;
    v["rh"_s] = 2;
    v["vw"_s] = 3;
    v["p"_s] = 4;

    ASSERT_EQ(v["t"_s] , 1);
    ASSERT_EQ(v["rh"_s] , 2);
    ASSERT_EQ(v["vw"_s] , 3);
    ASSERT_EQ(v["p"_s] , 4);
}

TEST_F(VariableStorageTest, has)
{
    variablestorage v (variables);

    v["t"] = 1;
    v["rh"] = 2;
    v["vw"] = 3;
    v["p"] = 4;

    ASSERT_TRUE(v.has("t"));
    ASSERT_TRUE(v.has("rh"));
    ASSERT_TRUE(v.has("vw"));
    ASSERT_TRUE(v.has("p"));

    ASSERT_FALSE(v.has("tttt"));
}

TEST_F(VariableStorageTest, size)
{
    variablestorage v (variables);

    v["t"] = 1;
    v["rh"] = 2;
    v["vw"] = 3;
    v["p"] = 4;

    ASSERT_EQ(v.size(),4);

}

TEST_F(VariableStorageTest, uninit)
{
    variablestorage v;
    ASSERT_ANY_THROW(v["t"] = 1);
}