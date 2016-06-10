#include "triangle_defines.h"
#include "masked_array.h"

double getRasterCell(const double *gt, const void *raster, double x, double y);

//returns 1->if invalid
int is_invalid_mean_elevation_diff(vertex triorg, vertex tridest, vertex triapex, double maxtolerance, const double *gt, const struct masked_array* masked_dem, const GDALRasterBandH* dem)
{
    double triangle_z_mean = 0;
    double tri_count = 0;

    double mx = triorg[0];
    double my = triorg[1];
    double value = getRasterCell(gt, dem, mx, my);
    if (value != -9999.)
    {
        triangle_z_mean += value;
        ++tri_count;
    }


    mx = triapex[0];
    my = triapex[1];
    value = getRasterCell(gt, dem, mx, my);
    if (value != -9999.)
    {
        triangle_z_mean += value;
        ++tri_count;
    }


    mx = tridest[0];
    my = tridest[1];
    value = getRasterCell(gt, dem, mx, my);
    if (value != -9999.)
    {
        triangle_z_mean += value;
        ++tri_count;
    }

    //planar mean elevation
    triangle_z_mean /= tri_count;



    if (tri_count == 0)
        return 0;

    double sum = 0;
    double count = 0;

    for (int i = 0; i < masked_dem->ysize; i++)
    {
        for (int j = 0; j < masked_dem->xsize; j++)
        {
            double value = get_MA(masked_dem,i,j);

            if (!isnan(value) && value != -9999.)
            {
                sum += value;
                ++count;
            }
        }
    }


    if (count == 0)
        return 0;

    double mean = sum / count;
//    printf("DEM mean = %f, tri mean = %f\n",mean,triangle_z_mean);

    int is_invalid = 0;
    double diff = fabs(triangle_z_mean - mean);
    if (diff > maxtolerance)
    {
//        printf("Tolerance exceeded!\n");
        is_invalid = 1;
    }
    return is_invalid;
}
