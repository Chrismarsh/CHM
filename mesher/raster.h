#pragma once
#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()

#include <iostream>
#include <string>
#include <utility>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

class raster
{
public:
    raster();

    virtual ~raster();

    void open(std::string path);
    void create_memory_raster(double xsize, double ysize, GDALDataType type);

    double getXY(double X, double Y);
    double getpXpY(int px, int py);
    GDALRasterBand *getBand() const;

    GDALDataset *getDs() const;
    void setDs(GDALDataset* ds);
    std::pair<int,int> xy_to_pxpy(double pX, double pY);
    void setMask(float *mask);
    void setBand(float *data, int xsize, int ysize);
    double *getGt() const;

private:
    double* gt;
    GDALDataset* ds;
    GDALRasterBand* band;
    float* mask;
    float* data;



};
