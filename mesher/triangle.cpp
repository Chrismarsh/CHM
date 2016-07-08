//
// Created by chris on 04/07/16.
//

#include "triangle.h"

triangle::triangle()
{
    GDALAllRegister();
    is_nan = false;
}

triangle::~triangle()
{

}

void triangle::make_rasterized(vertex v0, vertex v1, vertex v2, const raster& r)
{
    auto gt = r.getGt();
    auto wkt = r.getDs()->GetProjectionRef();
    auto srs = OGRSpatialReference(wkt);

    //memory shape file
    auto driver = GetGDALDriverManager()->GetDriverByName("Memory");
    auto shp = driver->Create("",0,0,0,GDT_Unknown,NULL);

    auto layer = shp->CreateLayer("poly",&srs,wkbPolygon,NULL);

//    OGRSpatialReferenceRLinearRing* ring = (OGRLinearRing*) OGRGeometryFactory().createGeometry(wkbLinearRing);
    OGRLinearRing ring;
    ring.addPoint(v0[0],v0[1]);
    ring.addPoint(v1[0],v1[1]);
    ring.addPoint(v2[0],v2[1]);
    ring.addPoint(v0[0],v0[1]); //close it

    OGRPolygon poly;
    poly.addRing(&ring);

    auto feature = OGRFeature::CreateFeature( layer->GetLayerDefn());
    feature->SetGeometry(&poly);
    layer->CreateFeature(feature);

    //we need to rasterize our triangle into a small raster

    double originX = gt[0];
    double originY = gt[3];

    double pixel_width = gt[1];
    double pixel_height = gt[5];

    OGREnvelope bbox;
    poly.getEnvelope(&bbox);

    //extents of the rasterized triangle
    double new_gt[6];

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

    int y_raster_size = r.getDs()->GetRasterYSize();
    int x_raster_size = r.getDs()->GetRasterXSize();
    if (y1 + ysize >= y_raster_size)
    {
        ysize = y_raster_size - y1;
    }
    if (x1 + xsize >=x_raster_size)
    {
        xsize = x_raster_size - x1;
    }

//    auto mem_drv = GetGDALDriverManager()->GetDriverByName("GTiff");
//    auto rvds = mem_drv->Create("rtri.tiff",xsize,ysize,1,GDT_Float32,NULL);

    auto mem_drv = GetGDALDriverManager()->GetDriverByName("MEM");
    auto rvds = mem_drv->Create("",xsize,ysize,1,GDT_Float32,NULL);
    rvds->SetGeoTransform(new_gt);

    char *pszSRS_WKT = NULL;
    srs.exportToWkt(&pszSRS_WKT);
    rvds->SetProjection(pszSRS_WKT);
    CPLFree(pszSRS_WKT);

    int nBandCount = 1;
    int *panBandList = (int *) malloc(sizeof(int) * nBandCount);
    panBandList[0] = 1; //raster band to use
    int nLayerCount = 1;
    double* burnValue =  new double;
    *burnValue = 1;

    OGRLayerH *pahLayers = (OGRLayerH *) malloc(sizeof(OGRLayerH) * nLayerCount);
    pahLayers[0] = layer;

    char **options = NULL;
    options = CSLSetNameValue(options, "ALL_TOUCHED", "TRUE");
    CPLErr err = GDALRasterizeLayers(rvds, nBandCount, panBandList, nLayerCount, pahLayers, NULL, NULL, burnValue, options, NULL,
                                     NULL);

//    GDALClose(rvds);
//    rvds = (GDALDataset*) GDALOpen("rtri.tiff",GA_Update);
//    exit(1);
    delete burnValue;
    CPLFree(*options);
    CPLFree(options);
//    delete *options;
//    delete options;
    free(panBandList);
    free(pahLayers);



    rasterized_triangle = boost::make_shared<raster>();
    rasterized_triangle->setDs(rvds);

    float* mask = new float[xsize*ysize];
    //copy the mask into the float array for the mask
    rasterized_triangle->getBand()->RasterIO(GF_Read, 0, 0, xsize, ysize,
                 mask, xsize, ysize, GDT_Float32,
                 0, 0);


    rasterized_triangle->setMask(mask);

    //overwrite the band data (that has the mask) to be the dem data
    float* dem = new float[xsize*ysize];

    r.getBand()->RasterIO(GF_Read, x1, y1, xsize, ysize,
                                             dem, xsize, ysize, GDT_Float32,
                                             0, 0);
    rasterized_triangle->setBand(dem,xsize,ysize);


    delete[] dem;
    OGRFeature::DestroyFeature(feature);
    GDALClose(shp);

    int px;
    int py;
    double pz;

    // because of how the rasterization works we MIGHT have a cell that lies out side of our actual dem domain.
    //for these edge effects we might have to average out the triangle vertexes
    int is_nan_v[3]={0,0,0};


    pz = rasterized_triangle->getXY(v0[0],v0[1]);
    this->v0[0] = v0[0];
    this->v0[1] = v0[1];
    this->v0[2] = pz;

    if(isnan(pz))
    {
        is_nan_v[0]=1;
    }


    pz = rasterized_triangle->getXY(v1[0],v1[1]);
    this->v1[0] = v1[0];
    this->v1[1] = v1[1];
    this->v1[2] = pz;
    if(isnan(pz))
    {
        is_nan_v[1]=1;
    }


    pz = rasterized_triangle->getXY(v2[0],v2[1]);
    this->v2[0] = v2[0];
    this->v2[1] = v2[1];
    this->v2[2] = pz;
    if(isnan(pz))
    {
        is_nan_v[2]=1;
    }
    if(is_nan_v[0] && is_nan_v[1] && is_nan_v[2])
    {
       is_nan = true;
    }

}
