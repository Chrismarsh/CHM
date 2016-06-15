#include "error_metrics.h"
int is_invalid_mean_elevation_diff(struct tri* t, double maxtolerance)
{
    // Initialize triangle mean elevation (m)
    double triangle_z_mean = 0;
    // Initialize count of trianlge vertices found not-nan
    double tri_count = 0;

    // Check if elevations of vertex is not-nan
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

    // If no non-nan vertices were found, return zero
    if (tri_count == 0)
        return 0;

    // Initialize sum and count of grid cell elevations
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

    // If none were found return zero
    if (count == 0)
        return 0;

    // Take mean of all grid cell elevations found
    double mean = sum / count;

    // Take difference of means
    int is_invalid = 0;
    double diff = fabs(triangle_z_mean - mean);
    if (diff > maxtolerance)
    {
//        printf("DEM mean = %f, tri mean = %f\n",mean,triangle_z_mean);
//        printf("Tolerance exceeded!\n");
        is_invalid = 1;
//        exit(1);

    }
    return is_invalid;

}

int is_invalid_tolerance(struct tri* t, double maxtolerance)
{


    return 0;

}
