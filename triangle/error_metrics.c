#include "error_metrics.h"
int is_invalid_mean_elevation_diff(struct tri* t, double maxtolerance)
{
    double triangle_z_mean = 0;
    double tri_count = 0;

    if ( !isnan(t->v0[2]))
    {
        triangle_z_mean += t->v0[2];
        ++tri_count;
    }

    if ( !isnan(t->v1[2]))
    {
        triangle_z_mean += t->v1[2];
        ++tri_count;
    }

    if ( !isnan(t->v2[2]))
    {
        triangle_z_mean += t->v2[2];
        ++tri_count;
    }


    //planar mean elevation
    triangle_z_mean /= tri_count;


    if (tri_count == 0)
        return 0;

    double sum = 0;
    double count = 0;

    for (int i = 0; i < t->rasterize_triangle->ysize; i++)
    {
        for (int j = 0; j < t->rasterize_triangle->xsize; j++)
        {
            double value = get_MA(t->rasterize_triangle,i,j);

            if (!isnan(value))
            {
                sum += value;
                ++count;
            }
        }
    }


    if (count == 0)
        return 0;

    double mean = sum / count;


    int is_invalid = 0;
    double diff = fabs(triangle_z_mean - mean);
    if (diff > maxtolerance)
    {
        printf("DEM mean = %f, tri mean = %f\n",mean,triangle_z_mean);
        printf("Tolerance exceeded!\n");
        is_invalid = 1;
//        exit(1);

    }
    return is_invalid;

}

int is_invalid_tolerance(struct tri* t, double maxtolerance)
{


    return 0;

}