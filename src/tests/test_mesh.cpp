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

class MeshTest : public testing::Test
{
protected:

    virtual void SetUp()
    {
        logging::core::get()->set_logging_enabled(false);
        EXPECT_NO_THROW(t.from_file("tin_30mdem_30mtol_nodes.csv"));
      
    }

    triangulation t;
};

TEST_F(MeshTest,ThrowsOnInvalidFile)
{
    triangulation t0;
    EXPECT_ANY_THROW(t0.from_file("invalid_file.txt"));
}

TEST_F(MeshTest,ToFileVTU)
{
    EXPECT_NO_THROW(t.mesh_to_vtu("test.vtu"));
}