//
// Created by Chris Marsh on 2016-06-09.
//

#include "masked_array.h"

double get_MA(const struct masked_array* array, int i, int j)
{

    if (array->mask[i*array->xsize+j] == 0)
    {
        return nan("");
    }
    else
    {
        return array->data[i*array->xsize+j];
    }

}

void malloc_MA(struct masked_array* array, int xsize, int ysize)
{
//    array->mask = malloc(sizeof(double*) * ysize);
//    array->data = malloc(sizeof(double*) * ysize);
//    for(int i =0; i < ysize; ++i )
//    {
//        array->mask[i] = (double*)malloc(sizeof(double) * xsize);
//        array->data[i] = (double*)malloc(sizeof(double) * xsize);
//    }

    array->mask = malloc(sizeof(float) * ysize*xsize);
    array->data = malloc(sizeof(float) * ysize*xsize);

    array->xsize = xsize;
    array->ysize = ysize;

}