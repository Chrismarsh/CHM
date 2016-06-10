#pragma once

#include <gdal.h>
#include <math.h>
struct masked_array
{
    float* data;
    float* mask;
    int xsize;
    int ysize;
} ;

double get_MA(const struct masked_array* array, int i, int j);
void malloc_MA(struct masked_array* array, int xsize, int ysize);