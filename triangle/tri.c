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

//    printf("%f %f\n", triorg[0], triorg[1]);
//    printf("%f %f\n",  tridest[0], tridest[1]);
//    printf("%f %f\n",  triapex[0], triapex[1]);


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

//    bbox.MinX -= 10;
//    bbox.MaxX += 10;
//    bbox.MinY -= 10;
//    bbox.MaxY += 10;

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

    GDALDriverH mem_drv = GDALGetDriverByName("MEM"); //MEM
    GDALDatasetH rvds=NULL; //output rasterized triangle
    rvds = GDALCreate(mem_drv, "", xsize, ysize, 1, GDT_Byte , NULL); //GDT_Byte
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

    if(err != 0)
    {
        printf("Return err was %d\n",err);
        exit(1);
    }
//

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

    // because of how the rasterization works we MIGHT have a cell that lies out side of our actual dem domain.
    //for these edge effects we might have to average out the triangle vertexes
    int is_nan[3]={0,0,0};
    /*
     * VERTEX 0
     */
    xyToPixel(triorg[0],triorg[1],&px,&py,new_gt,rvds);


    pz = get_MA(array,py,px);

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
    xyToPixel(tridest[0],tridest[1],&px,&py,new_gt,raster);
    pz = get_MA(array,py,px);

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
    xyToPixel(tridest[0],triapex[1],&px,&py,new_gt,raster);
    pz = get_MA(array,py,px);

    t->v2[0] = px;
    t->v2[1] = py;
    t->v2[2] = pz;
    if(isnan(pz))
    {
        is_nan[2]=1;
    }

    if(is_nan[0] && is_nan[1] && is_nan[2])
    {
        printf("Triangle vertexes are all nan\n");
        return NULL;
//        exit(-1);
    }

    t->rasterize_triangle = array;


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

void destory_triangle(struct tri* t)
{
    free( t->rasterize_triangle->data);
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