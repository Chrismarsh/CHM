#include "tri.h"


struct tri* createTriangle(vertex triorg, vertex tridest, vertex triapex, GDALDatasetH raster)
{
    struct tri* t = malloc(sizeof(struct tri));

    //get spatial reference data from our raster
    double *gt = malloc(6 * sizeof(double));
    GDALGetGeoTransform(raster, gt);

    const char *wkt = GDALGetProjectionRef(raster);
    OGRSpatialReferenceH srs = OSRNewSpatialReference(wkt);


//    rasterizeTriangle(triorg, tridest, triapex, srs, gt,
//                      &rvds, &x1, &y1, &xsize, &ysize, &new_gt);
        GDALDriverH driver = GDALGetDriverByName("Memory"); //MEM ESRI Shapefile
    GDALDatasetH DS = GDALCreate(driver, "", 0, 0, 0, GDT_Unknown, NULL);
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

    //extents of the rasterized triangle



    double* new_gt = malloc(sizeof(double) * 6);

    int x1 = (int) ((bbox.MinX - originX) / pixel_width);
    int x2 = (int) ((bbox.MaxX - originX) / pixel_width) + 1;
    int y1 = (int) ((bbox.MaxY - originY) / pixel_height);
    int y2 = (int) ((bbox.MinY - originY) / pixel_height) + 1;

    int xsize = x2 - x1;
    int ysize = y2 - y1;

    new_gt[0] = gt[0] + (x1 * gt[1]);
    new_gt[1] = gt[1];
    new_gt[2] = 0.0;
    new_gt[3] = gt[3] + (y1 * gt[5]);
    new_gt[4] = 0.0;
    new_gt[5] = gt[5];

//    GDALDriverH mem_drv = GDALGetDriverByName("MEM");
//    *rvds = GDALCreate(mem_drv, "", *xsize, *ysize, 1, GDT_Byte, NULL);

    GDALDriverH mem_drv = GDALGetDriverByName("GTiff");
    GDALDatasetH rvds=NULL; //output rasterized triangle
    rvds = GDALCreate(mem_drv, "rasterized_tri.tif", xsize, ysize, 1, GDT_Byte, NULL);
    GDALSetGeoTransform(rvds, new_gt);

    char *pszSRS_WKT = NULL;
    OSRExportToWkt( srs, &pszSRS_WKT );
    GDALSetProjection(rvds, pszSRS_WKT);

    int nBandCount = 1;
    int *panBandList = (int *) malloc(sizeof(int) * nBandCount);
    panBandList[0] = 1; //raster band to use
    int nLayerCount = 1;
    double burnValue = 1;

    OGRLayerH *pahLayers = (OGRLayerH *) malloc(sizeof(OGRLayerH) * nLayerCount);
    pahLayers[0] = layer;

    char **options = NULL;

    options = CSLSetNameValue(options, "ALL_TOUCHED", "TRUE");
    CPLErr err = GDALRasterizeLayers(rvds, nBandCount, panBandList, nLayerCount, pahLayers, NULL, NULL, &burnValue, options, NULL,
                                     NULL);
    printf("Return err was %d\n",err);

    OSRDestroySpatialReference(srs);

    GDALRasterBandH dem = GDALGetRasterBand(raster, 1);
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

    struct masked_array* array=  malloc_MA(xsize,ysize);
    array->x1 = x1;
    array->y1 = y1;
    array->gt = new_gt; //store our new reference

    GDALRasterIO(dem, GF_Read, x1, y1, xsize, ysize,
                 array->data, xsize, ysize,GDT_Float32,
                 0, 0);

    GDALRasterIO(masked, GF_Read, 0, 0, xsize, ysize,
                 array->mask , xsize, ysize, GDT_Float32,
                 0, 0);

    int px;
    int py;
    double pz;

    /*
     * VERTEX 0
     */
    xyToPixel(triorg[0],triorg[1],&px,&py,new_gt,rvds);
    pz = get_MA(array,py,px);

    t->v0[0] = px;
    t->v0[1] = py;
    t->v0[2] = pz;

    /*
     * VERTEX 1
     */
    xyToPixel(tridest[0],tridest[1],&px,&py,new_gt,raster);
    pz = get_MA(array,py,px);

    t->v0[0] = px;
    t->v0[1] = py;
    t->v0[2] = pz;

    /*
     * VERTEX 2
     */
    xyToPixel(triapex[0],triapex[1],&px,&py,new_gt,raster);
    pz = get_MA(array,py,px);

    t->v0[0] = px;
    t->v0[1] = py;
    t->v0[2] = pz;

    t->rasterize_triangle = array;

}
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
                       GDALDatasetH* rvds, int* x1, int* y1, int* xsize, int* ysize, double** new_gt)
{
//    GDALDriverH driver = GDALGetDriverByName("Memory"); //MEM ESRI Shapefile
//    GDALDatasetH DS = GDALCreate(driver, "", 0, 0, 0, GDT_Unknown, NULL);
    GDALDriverH driver = GDALGetDriverByName("ESRI Shapefile"); //MEM ESRI Shapefile
    GDALDatasetH DS = GDALCreate(driver, "tri.shp", 0, 0, 0, GDT_Unknown, NULL);

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



//    *new_gt = malloc(sizeof(double) * 6);

    (*new_gt)[0] = gt[0] + (*x1 * gt[1]);
    (*new_gt)[1] = gt[1];
    (*new_gt)[2] = 0.0;
    (*new_gt)[3] = gt[3] + (*y1 * gt[5]);
    (*new_gt)[4] = 0.0;
    (*new_gt)[5] = gt[5];

//    GDALDriverH mem_drv = GDALGetDriverByName("MEM");
//    *rvds = GDALCreate(mem_drv, "", *xsize, *ysize, 1, GDT_Byte, NULL);

    GDALDriverH mem_drv = GDALGetDriverByName("GTiff");
    *rvds = GDALCreate(mem_drv, "rasterized_tri.tif", *xsize, *ysize, 1, GDT_Byte, NULL);
    GDALSetGeoTransform(*rvds, *new_gt);

    char *pszSRS_WKT = NULL;
    OSRExportToWkt( srs, &pszSRS_WKT );

    GDALSetProjection(*rvds, pszSRS_WKT);

    int nBandCount = 1;
    int *panBandList = (int *) malloc(sizeof(int) * nBandCount);
    panBandList[0] = 1; //raster band to use
    int nLayerCount = 1;
    double burnValue = 1;

    OGRLayerH *pahLayers = (OGRLayerH *) malloc(sizeof(OGRLayerH) * nLayerCount);
    pahLayers[0] = layer;

    char **options = NULL;

    options = CSLSetNameValue(options, "ALL_TOUCHED", "TRUE");
    CPLErr err = GDALRasterizeLayers(*rvds, nBandCount, panBandList, nLayerCount, pahLayers, NULL, NULL, &burnValue, options, NULL,
                        NULL);
    printf("Return was %d\n",err);
//    free(pahLayers);

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
void xyToPixel(double x, double y, int* px, int* py, const double *gt, const void *raster)
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

    //trivial out of bound issue
    if (*px == GDALGetRasterBandXSize(raster))
        *px = *px - 1;
    if (*py == GDALGetRasterBandYSize(raster))
        *py = *py - 1;
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

    int px,py;
    xyToPixel(x,y,&px,&py,gt,raster);
    float pafScanline;

    GDALRasterIO(raster, GF_Read, px, py, 1, 1,
                 &pafScanline, 1, 1, GDT_Float32,
                 0, 0);

    return pafScanline;
}
/**
 * Converts the global x,y coordinate into the pixel x,y of a rasterized triangle
 * @param x
 * @param y
 * @param px
 * @param py
 * @param x1
 * @param y1
 * @param gt
 * @param raster
 * @return
 */
double getRasterizedPixel(double x, double y, int* px, int* py, int x1, int y1, const double* gt, const void *raster)
{
    int _px, _py;
    //get px and py have x,y in terms of the full raster
    xyToPixel(x,y,&_px,&_py,gt,raster);

    //uper x and y coords in pixel
    double Ux = gt[0] + (x1 * gt[1]); //upper x
    double Uy =  gt[3] + (y1 * gt[5]); //upper y

    int pUx, pUy;
    xyToPixel(Ux,Uy,&pUx,&pUy,gt,raster);

    //subtract the (0,0) pixel coordinates
    *px = _px - pUx;
    *py = _py - pUy;

}
//    OGR_DS_Destroy(DS);

