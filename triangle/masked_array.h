#pragma once

#include <gdal.h>
#include <math.h>
struct masked_array
{
    float* data;
    float* mask;
    int xsize;
    int ysize;
    int x1;
    int y1;
    double* gt;
} ;


double get_MA(const struct masked_array* array, int row, int col);
struct masked_array* malloc_MA(int xsize, int ysize);