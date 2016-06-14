//
// Created by Chris Marsh on 2016-06-09.
//

#include "masked_array.h"

void destroy_MA(struct masked_array** array)
{
    free((*array)->data);
    free((*array)->mask);

    (*array)->data = NULL;
    (*array)->mask = NULL;

    free( (*array)->gt);
    (*array)->gt = NULL;
}
double get_MA(const struct masked_array* array, int row, int col)
{

    //col + NCOLS * row
    int idx = col+array->xsize*row;
    double value;
    if (array->mask[idx] == 0)
    {
        value = nan("");
    }
    else
    {
        value = array->data[idx];
    }

    value = value == -9999.?nan(""):value;
    return value;

}

double get_MA_data(const struct masked_array* array, int row, int col)
{
    int idx = col+array->xsize*row;
    return array->data[idx];
}
struct masked_array* malloc_MA(int xsize, int ysize)
{

    struct masked_array*  array = malloc(sizeof(struct masked_array));
    array->mask = malloc(sizeof(float) * ysize*xsize);
    array->data = malloc(sizeof(float) * ysize*xsize);

    array->xsize = xsize;
    array->ysize = ysize;

    return array;

}