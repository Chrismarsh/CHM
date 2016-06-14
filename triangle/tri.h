#pragma  once
#include "triangle_defines.h"
#include "masked_array.h"

#include <gdal.h>
#include "ogr_api.h"
#include <gdal_alg.h>
#include "ogr_srs_api.h"
#include "cpl_conv.h" /* for CPLMalloc() */
#include <cpl_string.h>
#include <cpl_error.h>
struct tri
{
    //holds the rasterized triangle along with the masked array.
    struct masked_array* rasterize_triangle;


    //vertexes transformed to rasterized pixel offsets
    //x,y,z
    double v0[3];
    double v1[3];
    double v2[3];

    //if 1, then just ignore this triangle due to edge effects
    int is_nan;

};

void xyToPixel(double x, double y, int* px, int* py, const double *gt, const void *raster);
double getRasterCell(const double *gt, const void *raster, double x, double y);
struct tri* createTriangle(vertex triorg, vertex tridest, vertex triapex, GDALDatasetH raster);
void destory_triangle(struct tri* t);