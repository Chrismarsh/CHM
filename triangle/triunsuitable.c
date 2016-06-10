#include "triangle_defines.h"

#include <gdal.h>
#include "ogr_api.h"
#include <gdal_alg.h>
#include "ogr_srs_api.h"
#include "cpl_conv.h" /* for CPLMalloc() */
#include <cpl_string.h>


double getRasterCell(const double *gt, const void *raster, double x, double y)
{//extract the elevation from the 3 vertexes
    float pafScanline;
//    pafScanline = (float *) CPLMalloc(sizeof(float)*1);
    int px = (int)((x - gt[0]) / gt[1]);  // x pixel
    int py = (int)((y - gt[3]) / gt[5]);  // y pixel
    if (px == GDALGetRasterBandXSize(raster))
        px = px -1;
    if (py == GDALGetRasterBandYSize(raster))
        py = py - 1;


    GDALRasterIO( raster, GF_Read, px, py, 1, 1,
                  &pafScanline, 1, 1, GDT_Float32,
                  0, 0 );
    return pafScanline;
}

int triunsuitable(vertex triorg, vertex tridest, vertex triapex, REAL area, struct behavior *b)
{

    //there is no need to test the area here. This is because -aX comes first in the if chain, prior to -u
    //so we only need to check if we are violating other constraints as we don't even get this far if area > b->max area

    GDALDriverH driver = GDALGetDriverByName("ESRI Shapefile"); //MEM ESRI Shapefile
    GDALDatasetH DS   = GDALCreate	(driver,"lol.shp",0,0,0,GDT_Unknown,NULL);

//    GDALDriverH driver = GDALGetDriverByName("MEM"); //MEM ESRI Shapefile
//    GDALDatasetH DS   = GDALCreate	(driver,"out",0,0,0,GDT_Unknown,NULL);

    //get spatial reference from our raster
    const char* wkt = GDALGetProjectionRef(b->hDataset);

    OGRSpatialReferenceH srs = OSRNewSpatialReference(wkt);

    OGRLayerH layer = GDALDatasetCreateLayer(DS,"poly",srs,wkbPolygon,NULL);

    if(layer == NULL)
    {
        printf("Failed to create layer");
        exit(1);
    }
    OGRGeometryH ring = OGR_G_CreateGeometry(wkbLinearRing);
    if(ring == NULL)
    {
        printf("Failed to create ring");
        exit(1);
    }
    OGR_G_AddPoint_2D(ring,triorg[0],triorg[1]);
    OGR_G_AddPoint_2D(ring,tridest[0],tridest[1]);
    OGR_G_AddPoint_2D(ring,triapex[0],triapex[1]);
    OGR_G_AddPoint_2D(ring,triorg[0],triorg[1]);

    OGRGeometryH poly = OGR_G_CreateGeometry(wkbPolygon);
    if(poly == NULL)
    {
        printf("Failed to create poly");
        exit(1);
    }
    OGR_G_AddGeometry(poly, ring);

    OGRFeatureH feature = OGR_F_Create( OGR_L_GetLayerDefn(layer));
    OGR_F_SetGeometry(feature,poly);

    OGR_L_CreateFeature(layer,feature);

    //we need to rasterize our triangle into a small raster, so we can find out which cells are under it

    double* gt = malloc( 6 * sizeof(double));
    GDALGetGeoTransform( b->hDataset,gt);

    double originX = gt[0];
    double originY = gt[3];

    double pixel_width = gt[1];
    double pixel_height = gt[5];

    OGREnvelope bbox;
    OGR_G_GetEnvelope(poly,&bbox);


    int x1 = (int)((bbox.MinX - originX) / pixel_width);
    int x2 = (int)((bbox.MaxX- originX) / pixel_width) + 1;
    int y1 = (int)((bbox.MaxY - originY) / pixel_height);
    int y2 = (int)((bbox.MinY - originY) / pixel_height) + 1;

    int xsize = x2 - x1;
    int ysize = y2 - y1;

    double* new_gt = malloc( sizeof(double) * 6);
    new_gt[0] = gt[0] + (x1 * gt[1]);
    new_gt[1] = gt[1];
    new_gt[2] = 0.0;
    new_gt[3] = gt[3] + (y1 * gt[5]);
    new_gt[4] = 0.0;
    new_gt[5] = gt[5];

    GDALDriverH mem_drv = GDALGetDriverByName("MEM");
    GDALDatasetH rvds = GDALCreate(mem_drv,"",xsize,ysize,1,GDT_Byte,NULL);
    GDALSetGeoTransform(rvds,new_gt);

    int nBandCount = 1;
    int* panBandList = (int *) malloc(sizeof(int)*nBandCount);
    panBandList[0] = 1; //raster band to use
    int nLayerCount = 1;
    double burnValue = 1;

    OGRLayerH* pahLayers = (OGRLayerH *) malloc(sizeof(OGRLayerH)*nLayerCount);
    pahLayers[0] = layer;

    char** options = NULL;

    options = CSLSetNameValue(options, "ALL_TOUCHED", "TRUE");
    GDALRasterizeLayers(rvds,nBandCount,panBandList,nLayerCount,pahLayers,NULL,NULL,&burnValue,options,NULL,NULL);

    //now read the raster band that is underneath our triangle



    GDALRasterBandH dem_band = GDALGetRasterBand( b->hDataset, 1 );


    double triangle_z_mean = 0;

    double mx = triorg[0];
    double my = triorg[1];
    double tri_count = 0;
    int px;
    int py;

    double value = getRasterCell(gt, dem_band, mx, my);

    if(value != -9999.)
    {
        triangle_z_mean +=value;
        ++tri_count;
    }



    mx = triapex[0];
    my = triapex[1];
    value = getRasterCell(gt, dem_band, mx, my);
    if(value != -9999.)
    {
        triangle_z_mean +=value;
        ++tri_count;
    }


    mx = tridest[0];
    my = tridest[1];
    value = getRasterCell(gt, dem_band, mx, my);
    if(value != -9999.)
    {
        triangle_z_mean += value;
        ++tri_count;
    }

    //planar mean elevation
    triangle_z_mean /= tri_count;

    if (tri_count ==0)
        return 0;

    GDALRasterBandH masked = GDALGetRasterBand( rvds, 1 );

//    if (GDALGetRasterBandXSize( dem_band ) != GDALGetRasterBandXSize( masked ) ||
//            GDALGetRasterBandYSize(dem_band)!= GDALGetRasterBandYSize(masked)  )
//    {
//        printf("Error: Masked size != array size.");
//        exit(1);
//    }


    int   nXSize, nYSize,dem_x,dem_y;
    nXSize = GDALGetRasterBandXSize( masked );
    nYSize = GDALGetRasterBandYSize(masked);

    dem_x = GDALGetRasterBandXSize( dem_band );
    dem_y = GDALGetRasterBandYSize(dem_band);

    int masked_ysize =0;
    int masked_xsize =0;

    double sum = 0;
    double count = 0;

    if (y1+ysize >=  GDALGetRasterBandYSize(dem_band))
        ysize = GDALGetRasterBandYSize(dem_band) - y1;
    if (x1+xsize >=  GDALGetRasterBandXSize(dem_band))
        xsize = GDALGetRasterBandXSize(dem_band) - x1;


    float* pafScanline = (float *) malloc(sizeof(float)*xsize);
    float* masked_pafScanline = (float *) malloc(sizeof(float)*nXSize);
    for (int i=0;i<ysize;i++)
    {
        //read the DEM via offsets


        GDALRasterIO( dem_band, GF_Read, x1, y1+i, xsize, 1,
                      pafScanline, xsize, 1, GDT_Float32,
                      0, 0 );

        //read the masked

        GDALRasterIO( masked, GF_Read, 0, i, nXSize, 1,
                      masked_pafScanline, nXSize, 1, GDT_Float32,
                      0, 0 );


        for (int j=0;j<xsize;j++)
        {

            double m = masked_pafScanline[j];
            double value = pafScanline[j];
//            printf("\t(%f,%f)",m,value);
            if (m == 1 && value != -9999)
            {

                sum += value;
                ++count;
            }
        }
//        printf("\n");
    }


    if(count == 0)
        return 0;

    double mean = sum/ count;
//    printf("DEM mean = %f, tri mean = %f\n",mean,triangle_z_mean);

    int is_invalid=0;
    double diff = fabs(triangle_z_mean - mean);
    if(diff  > b->maxtolerance )
    {
//        printf("Tolerance exceeded!\n");
        is_invalid=1;
    }



    free(pahLayers);
    free(gt);
    free(new_gt);
    free(pafScanline);
    free(masked_pafScanline);

//
//    GDALDestroyDriver(driver);
//    GDALDestroyDriver(mem_drv);

    OGR_DS_Destroy(DS);

//    OGR_DS_Destroy(layer);
    OGR_G_DestroyGeometry(ring);
    OGR_G_DestroyGeometry(poly);

    OGR_DS_Destroy(rvds);

    OSRDestroySpatialReference(srs);

    return is_invalid;

}