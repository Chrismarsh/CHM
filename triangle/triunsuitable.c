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

    struct tri* t=NULL;

//    const char *wkt = GDALGetProjectionRef(b->hDataset);
//    OGRSpatialReferenceH srs = OSRNewSpatialReference(wkt);
//    int isGeographic = OSRIsGeographic(srs);
////    printf("DEM is geographic = %d\n",b->isGeographic);
//
//    //if we are lat/long geographic, we have to project the triangle, which requires a call to create triangle
//    if(isGeographic)
//    {
//        t = createTriangle(triorg, tridest, triapex, b->hDataset);
//        area = t->area; //override with our projection calculated area
//    }

    //if we are already projected, we can try to skip the createTriangle call.
    if(area >  b->maxarea)
        return 1;

    if(area <= b->minarea)
        return 0;

//    if(!isGeographic)
        t = createTriangle(triorg, tridest, triapex, b->hDataset);


    int errormetric = b->tol_method; // 1 is the default
    int is_invalid = 0;

    if(!t->is_nan)
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
    }
    else
    {
        is_invalid = 0;
    }


    destory_triangle(t);
    free(t);
    t = NULL;


    return is_invalid;

}

