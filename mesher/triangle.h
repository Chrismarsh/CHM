#pragma once

#include <utility>
#include <boost/shared_ptr.hpp>

#include <ogr_spatialref.h>
#include <ogr_geometry.h>
#include <gdal_alg.h>
#include <ogrsf_frmts.h>

#include <cmath>

#include "raster.h"

typedef double vertex[3];

class triangle
{
public:
    triangle();
    void make_rasterized(vertex v0_in, vertex v1_in, vertex v2_in, const raster& r);


    boost::shared_ptr<raster> rasterized_triangle;
    virtual ~triangle();

    bool is_nan;

    //vertexes transformed to rasterized pixel offsets
    //x,y,z
    double v0[3];
    double v1[3];
    double v2[3];
};

