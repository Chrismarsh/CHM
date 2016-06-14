#pragma  once
#include "triangle_defines.h"
#include "masked_array.h"

#include <gdal.h>
#include "ogr_api.h"
#include <gdal_alg.h>
#include "ogr_srs_api.h"
#include "cpl_conv.h" /* for CPLMalloc() */
#include <cpl_string.h>
struct tri
{
    //holds the rasterized triangle along with the masked array.
    struct masked_array* rasterize_triangle;


    //vertexes transformed to rasterized pixel offsets
    //x,y,z
    double v0[3];
    double v1[3];
    double v2[3];

};

void rasterizeTriangle(vertex triorg, vertex tridest, vertex triapex, OGRSpatialReferenceH srs, const double *gt,
                       GDALDatasetH* rvds, int* x1, int* y1, int* xsize, int* ysize, double** new_gt);
void xyToPixel(double x, double y, int* px, int* py, const double *gt, const void *raster);
double getRasterCell(const double *gt, const void *raster, double x, double y);
double getRasterizedPixel(double x, double y, int* px, int* py, int x1, int y1, const double* gt, const void *raster);
struct tri* createTriangle(vertex triorg, vertex tridest, vertex triapex, GDALDatasetH raster);