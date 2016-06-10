#include "triangle_defines.h"
#include "masked_array.h"
#include "error_metrics.h"

#include <gdal.h>
#include "ogr_api.h"
#include <gdal_alg.h>
#include "ogr_srs_api.h"
#include "cpl_conv.h" /* for CPLMalloc() */
#include <cpl_string.h>

void rasterizeTriangle(vertex triorg, vertex tridest, vertex triapex, OGRSpatialReferenceH srs, const double *gt,
                       GDALDatasetH* rvds, int* x1, int* y1, int* xsize, int* ysize);

double getRasterCell(const double *gt, const void *raster, double x, double y);


int triunsuitable(vertex triorg, vertex tridest, vertex triapex, REAL area, struct behavior *b)
{

    //get spatial reference data from our raster
    double *gt = malloc(6 * sizeof(double));
    GDALGetGeoTransform(b->hDataset, gt);

    const char *wkt = GDALGetProjectionRef(b->hDataset);
    OGRSpatialReferenceH srs = OSRNewSpatialReference(wkt);

    //extents of the rasterized triangle
    int x1,y1,xsize,ysize;
    GDALDatasetH rvds=NULL; //output rasterized triangle
    rasterizeTriangle(triorg, tridest, triapex, srs, gt,
                                      &rvds, &x1, &y1, &xsize, &ysize);

    OSRDestroySpatialReference(srs);

    GDALRasterBandH dem = GDALGetRasterBand(b->hDataset, 1);
    GDALRasterBandH masked = GDALGetRasterBand(rvds, 1);

    //check for trivial out of bounds
    if (y1 + ysize >= GDALGetRasterBandYSize(dem))
    {
        ysize = GDALGetRasterBandYSize(dem) - y1;

    }
    if (x1 + xsize >= GDALGetRasterBandXSize(dem))
    {
        xsize = GDALGetRasterBandXSize(dem) - x1;
    }

    struct masked_array array;
    malloc_MA(&array,xsize,ysize);

    GDALRasterIO(dem, GF_Read, x1, y1, xsize, ysize,
                 array.data, xsize, ysize,GDT_Float32,
                 0, 0);

    GDALRasterIO(masked, GF_Read, 0, 0, xsize, ysize,
                 array.mask , xsize, ysize, GDT_Float32,
                 0, 0);


    int errormetric = 1;
    int is_invalid = 0;

    switch(errormetric)
    {
        case 1:
            is_invalid = is_invalid_mean_elevation_diff(triorg, tridest, triapex,  b->maxtolerance,  gt,  &array,  dem);
            break;

        default:
            printf("Invalid error metric chosen");
            exit(1);
    }

//
    free(gt);
    free(array.data);
    free(array.mask);

//    free(pafScanline);
//    free(masked_pafScanline);

//
//    GDALDestroyDriver(driver);
//    GDALDestroyDriver(mem_drv);

//    OGR_DS_Destroy(DS);

//    OGR_DS_Destroy(layer);
//    OGR_G_DestroyGeometry(ring);
//    OGR_G_DestroyGeometry(poly);
//
//    OGR_DS_Destroy(rvds);



    return is_invalid;

}


//there is no need to test the area here. This is because -aX comes first in the if chain, prior to -u
//so we only need to check if we are violating other constraints as we don't even get this far if area > b->max area

/**
 * Rasterizes a triangle as given by 3 verticies.
 * @param triorg
 * @param tridest
 * @param triapex
 * @param srs Spatial reference to use
 * @param gt Geotransform of the raster within which we are operating
 * @return
 */
void rasterizeTriangle(vertex triorg, vertex tridest, vertex triapex, OGRSpatialReferenceH srs, const double *gt,
        GDALDatasetH* rvds, int* x1, int* y1, int* xsize, int* ysize)
{
    GDALDriverH driver = GDALGetDriverByName("Memory"); //MEM ESRI Shapefile
    GDALDatasetH DS = GDALCreate(driver, "", 0, 0, 0, GDT_Unknown, NULL);

    GDALDriverH driver = GDALGetDriverByName("Memory"); //MEM ESRI Shapefile
    GDALDatasetH DS   = GDALCreate	(driver,"",0,0,0,GDT_Unknown,NULL);



    OGRLayerH layer = GDALDatasetCreateLayer(DS, "poly", srs, wkbPolygon, NULL);

    if (layer == NULL)
    {
        printf("Failed to create layer");
        exit(1);
    }
    OGRGeometryH ring = OGR_G_CreateGeometry(wkbLinearRing);
    if (ring == NULL)
    {
        printf("Failed to create ring");
        exit(1);
    }
    OGR_G_AddPoint_2D(ring, triorg[0], triorg[1]);
    OGR_G_AddPoint_2D(ring, tridest[0], tridest[1]);
    OGR_G_AddPoint_2D(ring, triapex[0], triapex[1]);
    OGR_G_AddPoint_2D(ring, triorg[0], triorg[1]);

    OGRGeometryH poly = OGR_G_CreateGeometry(wkbPolygon);
    if (poly == NULL)
    {
        printf("Failed to create poly");
        exit(1);
    }
    OGR_G_AddGeometry(poly, ring);

    OGRFeatureH feature = OGR_F_Create(OGR_L_GetLayerDefn(layer));
    OGR_F_SetGeometry(feature, poly);

    OGR_L_CreateFeature(layer, feature);

    //we need to rasterize our triangle into a small raster

    double originX = gt[0];
    double originY = gt[3];

    double pixel_width = gt[1];
    double pixel_height = gt[5];

    OGREnvelope bbox;
    OGR_G_GetEnvelope(poly, &bbox);


    *x1 = (int) ((bbox.MinX - originX) / pixel_width);
    int x2 = (int) ((bbox.MaxX - originX) / pixel_width) + 1;
    *y1 = (int) ((bbox.MaxY - originY) / pixel_height);
    int y2 = (int) ((bbox.MinY - originY) / pixel_height) + 1;

    *xsize = x2 - *x1;
    *ysize = y2 - *y1;

    double *new_gt = malloc(sizeof(double) * 6);
    new_gt[0] = gt[0] + (*x1 * gt[1]);
    new_gt[1] = gt[1];
    new_gt[2] = 0.0;
    new_gt[3] = gt[3] + (*y1 * gt[5]);
    new_gt[4] = 0.0;
    new_gt[5] = gt[5];

    GDALDriverH mem_drv = GDALGetDriverByName("MEM");
    *rvds = GDALCreate(mem_drv, "", *xsize, *ysize, 1, GDT_Byte, NULL);
    GDALSetGeoTransform(*rvds, new_gt);

    int nBandCount = 1;
    int *panBandList = (int *) malloc(sizeof(int) * nBandCount);
    panBandList[0] = 1; //raster band to use
    int nLayerCount = 1;
    double burnValue = 1;

    OGRLayerH *pahLayers = (OGRLayerH *) malloc(sizeof(OGRLayerH) * nLayerCount);
    pahLayers[0] = layer;

    char **options = NULL;

    options = CSLSetNameValue(options, "ALL_TOUCHED", "TRUE");
    GDALRasterizeLayers(*rvds, nBandCount, panBandList, nLayerCount, pahLayers, NULL, NULL, &burnValue, options, NULL,
                        NULL);

//    free(pahLayers);

}

/**
 * Returns the value given by projected coordinates from a raster
 * @param gt
 * @param raster
 * @param x
 * @param y
 * @return
 */
double getRasterCell(const double *gt, const void *raster, double x, double y)
{

    float pafScanline;
//    pafScanline = (float *) CPLMalloc(sizeof(float)*1);
    int px = (int) ((x - gt[0]) / gt[1]);  // x pixel
    int py = (int) ((y - gt[3]) / gt[5]);  // y pixel

    //trivial out of bound issue
    if (px == GDALGetRasterBandXSize(raster))
        px = px - 1;
    if (py == GDALGetRasterBandYSize(raster))
        py = py - 1;


    GDALRasterIO(raster, GF_Read, px, py, 1, 1,
                 &pafScanline, 1, 1, GDT_Float32,
                 0, 0);
    return pafScanline;
}