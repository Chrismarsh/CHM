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

    //create the vectors veco0 and
    double u1,u2,u3;
    double v1,v2,v3;
    double o1,o2,o3; //origin of tri

    o1 = t->v0[0];
    o2 = t->v0[1];
    o3 = t->v0[2];

    //if we have only nan z values, bail.
    if(isnan(t->v0[2]) || isnan(t->v1[2]) || isnan(t->v2[2]))
    {
        t->is_nan;
        return 0;
    }

    //following http://www.had2know.com/academics/equation-plane-through-3-points.html
    //create the two vectors
    u1 = t->v1[0] - t->v0[0];
    u2 = t->v1[1] - t->v0[1];
    u3 = t->v1[2] - t->v0[2];

    v1 = t->v2[0] - t->v0[0];
    v2 = t->v2[1] - t->v0[1];
    v3 = t->v2[2] - t->v0[2];

    //calculate the normal vector via cross product
    double a, b, c;
    a = u2 * v3 - v2 * u3;
    b = v1 * u3 - u1 * v3;
    c = u1 * v2 - v1 * u2;


    //solve for d
    double d =a*o1 + b*o2 + c*o3;

    double rmse = 0;
    double n = 0;
    for (int y = 0; y < t->rasterize_triangle->ysize; y++)
    {
        for (int x = 0; x < t->rasterize_triangle->xsize; x++)
        {
            double value = get_MA(t->rasterize_triangle,y,x);

            if (!isnan(value))
            {
                //value predicted by the triangle
                double z = -(a*x+b*y-d)/c; //plane eqn solved for z. allows us to predict z values via x,y coords
                double diff = z - value;

                rmse += diff * diff;
                n++;

            }
        }
    }

    //bail
    if (n == 0.)
        return 0; // don't bother checking again

    rmse /= n;

    rmse = sqrt(rmse);

    if(rmse > maxtolerance)
        return 1;

    return 0;

}
