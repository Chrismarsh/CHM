#include "triangle_defines.h"
#include "error_metrics.h"
#include "tri.h"

//#include <gdal.h>
//#include "ogr_api.h"
//#include <gdal_alg.h>
//#include "ogr_srs_api.h"
//#include "cpl_conv.h" /* for CPLMalloc() */
//#include <cpl_string.h>



int triunsuitable(vertex triorg, vertex tridest, vertex triapex, REAL area, struct behavior *b)
{
//there is no need to test the area here. This is because -aX comes first in the if chain, prior to -u
//so we only need to check if we are violating other constraints as we don't even get this far if area > b->max area

    struct tri* t = createTriangle(triorg, tridest, triapex, b->hDataset);

    int errormetric = 1;
    int is_invalid = 0;

    if(t)
    {
        switch(errormetric)
        {
            case 1:
                is_invalid = is_invalid_mean_elevation_diff(t, b->maxtolerance);
                break;
            case 2:
                is_invalid = is_invalid_tolerance(t, b->maxtolerance);
                break;
            default:
                printf("Invalid error metric chosen");
                exit(1);
        }

        destory_triangle(t);
        free(t);
        t = NULL;

    }
    else
    {
        is_invalid = 0;
    }


    return is_invalid;

}

