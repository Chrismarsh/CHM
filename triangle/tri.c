#include <ogr_core.h>
#include "tri.h"


struct tri* createTriangle(vertex triorg, vertex tridest, vertex triapex, GDALDatasetH raster)
{
    struct tri* t = malloc(sizeof(struct tri));
    t->is_nan = 0;

    //get spatial reference data from our raster
    double *gt = malloc(6 * sizeof(double));
    GDALGetGeoTransform(raster, gt);

    const char *wkt = GDALGetProjectionRef(raster);
    OGRSpatialReferenceH srs = OSRNewSpatialReference(wkt);

    //determine if we need to enable the projection specific codepath for dealing with lat-long
    int is_geographic = OSRIsGeographic(srs);
//    printf("Is geographic = %d",is_geographic);

    GDALDriverH driver = GDALGetDriverByName("Memory"); //MEM ESRI Shapefile
    GDALDatasetH DS = NULL;
    DS = GDALCreate(driver, "", 0, 0, 0, GDT_Unknown, NULL);

//    GDALDriverH driver = GDALGetDriverByName("ESRI Shapefile"); //MEM ESRI Shapefile
//    GDALDatasetH DS = GDALCreate(driver, "tri.shp", 0, 0, 0, GDT_Unknown, NULL);

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
//    OGR_G_AssignSpatialReference(ring,srs);
    OGR_G_AddPoint_2D(ring, triorg[0], triorg[1]);
    OGR_G_AddPoint_2D(ring, tridest[0], tridest[1]);
    OGR_G_AddPoint_2D(ring, triapex[0], triapex[1]);
    OGR_G_AddPoint_2D(ring, triorg[0], triorg[1]);

//    printf("%f %f\n", triorg[0], triorg[1]);
//    printf("%f %f\n",  tridest[0], tridest[1]);
//    printf("%f %f\n",  triapex[0], triapex[1]);


    OGRGeometryH poly = OGR_G_CreateGeometry(wkbPolygon);
    if (poly == NULL)
    {
        printf("Failed to create poly");
        exit(1);
    }
//    OGR_G_AssignSpatialReference(poly,srs);
    OGR_G_AddGeometry(poly, ring);

    //triangle area in m^2
    double area=0;
/**
 * Uses a Equal area projection to determine the area in m^2 of a triangle
 */
    if(is_geographic)
    {
        const char *wkt_out = "PROJCS[\"North_America_Albers_Equal_Area_Conic\",     GEOGCS[\"GCS_North_American_1983\",         DATUM[\"North_American_Datum_1983\",             SPHEROID[\"GRS_1980\",6378137,298.257222101]],         PRIMEM[\"Greenwich\",0],         UNIT[\"Degree\",0.017453292519943295]],     PROJECTION[\"Albers_Conic_Equal_Area\"],     PARAMETER[\"False_Easting\",0],     PARAMETER[\"False_Northing\",0],     PARAMETER[\"longitude_of_center\",-96],     PARAMETER[\"Standard_Parallel_1\",20],     PARAMETER[\"Standard_Parallel_2\",60],     PARAMETER[\"latitude_of_center\",40],     UNIT[\"Meter\",1],     AUTHORITY[\"EPSG\",\"102008\"]]";

        OGRSpatialReferenceH srs_out = OSRNewSpatialReference(wkt_out);
        if (!srs_out)
        {
            printf("error");
            exit(1);
        }
        OGRSpatialReferenceH trans = OCTNewCoordinateTransformation(srs, srs_out);

        OGRGeometryH poly_prj = OGR_G_Clone(poly);
        OGR_G_Transform(poly_prj, trans);

        area = OGR_G_Area(poly_prj);   //    OGR_G_GetArea() is deprecated

        OGR_G_DestroyGeometry(poly_prj);
        OCTDestroyCoordinateTransformation(trans);
    }
    else
    {
        area = OGR_G_Area(poly); //use a non projected triangle
    }
/**
 * Done transformation
 */

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

    //extents of the rasterized triangle
    double* new_gt = malloc(sizeof(double) * 6);

    int x1 = (int) ((bbox.MinX - originX) / pixel_width);
    int x2 = (int) ((bbox.MaxX - originX) / pixel_width) + 1;
    int y1 = (int) ((bbox.MaxY - originY) / pixel_height);
    int y2 = (int) ((bbox.MinY - originY) / pixel_height) + 1;

    int xsize = x2 - x1;
    int ysize = y2 - y1;

    //make a window into the main extent
    new_gt[0] = gt[0] + (x1 * gt[1]);
    new_gt[1] = gt[1];
    new_gt[2] = 0.0;
    new_gt[3] = gt[3] + (y1 * gt[5]);
    new_gt[4] = 0.0;
    new_gt[5] = gt[5];


    GDALDatasetH rvds=NULL; //output rasterized triangle

    GDALDriverH mem_drv = GDALGetDriverByName("MEM");
//    GDALDriverH mem_drv = GDALGetDriverByName("GTiff");
    rvds = GDALCreate(mem_drv, "", xsize, ysize, 1, GDT_Byte , NULL); //GDT_Byte
//    rvds = GDALCreate(mem_drv, "triraster.tiff", xsize, ysize, 1, GDT_Byte , NULL); //GDT_Byte
    GDALSetGeoTransform(rvds, new_gt);

    char *pszSRS_WKT = NULL;
    OSRExportToWkt( srs, &pszSRS_WKT );
    GDALSetProjection(rvds, pszSRS_WKT);

    int nBandCount = 1;
    int *panBandList = (int *) malloc(sizeof(int) * nBandCount);
    panBandList[0] = 1; //raster band to use
    int nLayerCount = 1;
    double* burnValue = malloc(sizeof(double));
    *burnValue = 1;

    OGRLayerH *pahLayers = (OGRLayerH *) malloc(sizeof(OGRLayerH) * nLayerCount);
    pahLayers[0] = layer;

    char **options = NULL;

    options = CSLSetNameValue(options, "ALL_TOUCHED", "TRUE");
    CPLErr err = GDALRasterizeLayers(rvds, nBandCount, panBandList, nLayerCount, pahLayers, NULL, NULL, burnValue, options, NULL,
                                     NULL);

    CPLFree(*options);
    CPLFree(options);
    if(err != 0)
    {
        printf("Return err was %d\n",err);
        exit(1);
    }


    OSRDestroySpatialReference(srs);

    GDALRasterBandH dem = GDALGetRasterBand(raster, 1);
    GDALRasterBandH masked = GDALGetRasterBand(rvds, 1);


    //we can end up 1 pixel out of bounds of the main DEM.
    //if this happens, the calls to fill the masked array will fail. So if we end up off by one
    //shrink the domain to make sure it's inline with the main DEM bounds

    int y_raster_size = GDALGetRasterBandYSize(dem);
    int x_raster_size = GDALGetRasterBandXSize(dem);
    if (y1 + ysize >= y_raster_size)
    {
        ysize = y_raster_size - y1;
    }
    if (x1 + xsize >=x_raster_size)
    {
        xsize = x_raster_size - x1;
    }

    struct masked_array* array=  malloc_MA(xsize,ysize);
    array->x1 = x1;
    array->y1 = y1;
    array->gt = new_gt; //store our new reference

    GDALRasterIO(dem, GF_Read, x1, y1, xsize, ysize,
                 array->data, xsize, ysize, GDT_Float32,
                 0, 0);

    GDALRasterIO(masked, GF_Read, 0, 0, xsize, ysize,
                 array->mask , xsize, ysize, GDT_Float32,
                 0, 0);

    int px;
    int py;
    double pz;

    // because of how the rasterization works we MIGHT have a cell that lies out side of our actual dem domain.
    //for these edge effects we might have to average out the triangle vertexes
    int is_nan[3]={0,0,0};
    /*
     * VERTEX 0
     */
    xyToPixel(triorg[0],triorg[1],&px,&py,new_gt,xsize,ysize);


    pz = get_MA_data(array,py,px);

    t->v0[0] = px;
    t->v0[1] = py;
    t->v0[2] = pz;
    if(isnan(pz))
    {
        is_nan[0]=1;
    }
    /*
     * VERTEX 1
     */
    xyToPixel(tridest[0],tridest[1],&px,&py,new_gt,xsize,ysize);
    pz = get_MA_data(array,py,px);

    t->v1[0] = px;
    t->v1[1] = py;
    t->v1[2] = pz;
    if(isnan(pz))
    {
        is_nan[1]=1;
    }

    /*
     * VERTEX 2
     */
    xyToPixel(triapex[0],triapex[1],&px,&py,new_gt,xsize,ysize);
    pz = get_MA_data(array,py,px);

    t->v2[0] = px;
    t->v2[1] = py;
    t->v2[2] = pz;
    if(isnan(pz))
    {
        is_nan[2]=1;
    }

    if(is_nan[0] && is_nan[1] && is_nan[2])
    {
//        printf("Triangle vertexes are all nan\n");
        t->is_nan = 1;

    }


    t->rasterize_triangle = array;
    t->area = area;

//    GDALDestroyDriver(driver);
//    GDALDestroyDriver(mem_drv);

//   CPLFree(wkt);

    GDALClose(DS);
    GDALClose(rvds);



//    OGR_F_Destroy(layer);
    OGR_G_DestroyGeometry(poly);
    OGR_G_DestroyGeometry(ring);
    OGR_F_Destroy(feature);


    free(gt);
    free(pszSRS_WKT);
    free(panBandList);
    free(burnValue);
    free(pahLayers);

    return t;
}


/**
 * Converts a raster projected x,y coordinate into a pixel coordinate
 * @param x
 * @param y
 * @param px
 * @param py
 * @param gt
 * @param raster
 */
void xyToPixel(double x, double y, int* px, int* py, const double *gt, int max_x, int max_y)
{
//    adfGeoTransform[0] /* top left x */
//    adfGeoTransform[1] /* w-e pixel resolution */
//    adfGeoTransform[2] /* 0 */
//    adfGeoTransform[3] /* top left y */
//    adfGeoTransform[4] /* 0 */
//    adfGeoTransform[5] /* n-s pixel resolution (negative value) */
//
    *px = (int) ((x - gt[0]) / gt[1]);  // x pixel
    *py = (int) ((y - gt[3]) / gt[5]);  // y pixel


    //out of bound issue when we are off by one
    if (*px == max_x)
        *px = max_x - 1;

    if (*py == max_y)
        *py = max_y - 1;

    if(*py == -1)
        *py = 0;

    if(*px == -1)
        *px = 0;

}

/**
 * Returns the value given by projected coordinates from a raster
 * @param gt
 * @param raster
 * @param x
 * @param y
 * @return
 */
//double getRasterCell(const double *gt, const void *raster, double x, double y)
//{
//
//    int px,py;
//    xyToPixel(x,y,&px,&py,gt,max_x,max_y);
//    float pafScanline;
//
//    GDALRasterIO(raster, GF_Read, px, py, 1, 1,
//                 &pafScanline, 1, 1, GDT_Float32,
//                 0, 0);
//
//    return pafScanline;
//}

void destory_triangle(struct tri* t)
{
    free(t->rasterize_triangle->data);
    free(t->rasterize_triangle->mask);
    free(t->rasterize_triangle->gt);

    t->rasterize_triangle->data = NULL;
    t->rasterize_triangle->mask = NULL;
    t->rasterize_triangle->gt = NULL;

    free(t->rasterize_triangle);
    t->rasterize_triangle = NULL;
//    destroy_MA( &((*t)->rasterize_triangle));
//    free( (*t)->rasterize_triangle);
//    (*t) = NULL;
//

}