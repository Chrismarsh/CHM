//
// Created by chris on 04/07/16.
//

#include <csetjmp>
#include <csignal>
jmp_buf env;
void on_sigabrt (int signum)
{
    longjmp (env, 1);
}

#include "triangle.h"

triangle::triangle()
{

    is_nan = false;
}

triangle::~triangle()
{

}

void triangle::make_rasterized(vertex v0_in, vertex v1_in, vertex v2_in, const raster& r)
{
    auto gt = r.getGt();
    auto wkt = r.getDs()->GetProjectionRef();
    auto srs = OGRSpatialReference(wkt);

    OGRLinearRing ring;
    ring.addPoint(v0_in[0],v0_in[1]);
    ring.addPoint(v1_in[0],v1_in[1]);
    ring.addPoint(v2_in[0],v2_in[1]);
    ring.addPoint(v0_in[0],v0_in[1]); //close it

    OGRPolygon poly;
    poly.addRing(&ring);

    //memory shape file
//    auto driver = GetGDALDriverManager()->GetDriverByName("Memory");
//    auto shp = driver->Create("",0,0,0,GDT_Unknown,NULL);
//    auto layer = shp->CreateLayer("poly",&srs,wkbPolygon,NULL);

//    auto feature = OGRFeature::CreateFeature( layer->GetLayerDefn());
//    feature->SetGeometry(&poly);
//    layer->CreateFeature(feature);

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

    if(!rvds)
    {
        std::cout << "Yikes!" << std::endl;
    }
    rvds->SetGeoTransform(new_gt);


    rvds->SetProjection(wkt);

    std::vector<int> bandlist;
    bandlist.push_back(1);

    double burnValue = 1;

    std::vector<OGRPolygon*> geoms;
    geoms.push_back(&poly);

    char **options = nullptr;
    options = CSLSetNameValue(options, "ALL_TOUCHED", "TRUE");

    //in som very rare, and still not totally sure why cases, this will fail with SIGABRT. So catch it, mark the triangle as isnan, and try to continue
    if (setjmp (env) == 0) {
        signal(SIGABRT, &on_sigabrt);
        CPLErr err = GDALRasterizeGeometries(rvds, 1, &bandlist[0], 1, (OGRGeometryH*)&geoms[0], NULL, NULL, &burnValue, options, NULL, NULL);
    }
    else {
        //all kinds of problems at this point, bail out
        this->is_nan = true;
        //        std::cout << "disaster!" << std::endl;
        exit(1);
        return;



//        //testing code to write out the triangle that is causing issues.
//        auto driver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
//        auto shp = driver->Create("broken_tri.shp",0,0,0,GDT_Unknown,NULL);
//        auto layer = shp->CreateLayer("poly",&srs,wkbPolygon,NULL);
//
//        auto feature = OGRFeature::CreateFeature( layer->GetLayerDefn());
//        feature->SetGeometry(&poly);
//        layer->CreateFeature(feature);
//        GDALClose(shp);
//
//
//        auto drv = GetGDALDriverManager()->GetDriverByName("GTiff");
//        auto rvds = drv->Create("rtri.tiff",xsize,ysize,1,GDT_Float32,NULL);
//
//        char *pszSRS_WKT = NULL;
//        srs.exportToWkt(&pszSRS_WKT);
//        rvds->SetProjection(pszSRS_WKT);
//        CPLFree(pszSRS_WKT);
//
//        GDALClose(rvds);

//        exit(1);
//
//        CPLErr err = GDALRasterizeGeometries(rvds, 1, &bandlist[0], 1, (OGRGeometryH*)&geoms[0], NULL, NULL, &burnValue, options, NULL, NULL);


    }

    CSLDestroy(options);

    rasterized_triangle = boost::make_shared<raster>();
    rasterized_triangle->setDs(rvds); //rvds owned now

    float* mask = new float[xsize*ysize]; //mask will be owned
    //copy the mask into the float array for the mask
    rasterized_triangle->getBand()->RasterIO(GF_Read, 0, 0, xsize, ysize,
                 mask, xsize, ysize, GDT_Float32,
                 0, 0);

    rasterized_triangle->setMask(mask);

    //overwrite the band data (that has the mask) to be the dem data
    float* dem = new float[xsize*ysize]; // triangle owns this now

    r.getBand()->RasterIO(GF_Read, x1, y1, xsize, ysize,
                                             dem, xsize, ysize, GDT_Float32,
                                             0, 0);
    rasterized_triangle->setBand(dem,xsize,ysize);


    double pz;

    // because of how the rasterization works we MIGHT have a cell that lies out side of our actual dem domain.
    //for these edge effects we might have to average out the triangle vertexes
    int is_nan_v[3]={0,0,0};


    pz = rasterized_triangle->getXY(v0_in[0],v0_in[1]);
    this->v0[0] = v0_in[0];
    this->v0[1] = v0_in[1];
    this->v0[2] = pz;

    if(std::isnan(pz))
    {
        is_nan_v[0]=1;
    }


    pz = rasterized_triangle->getXY(v1_in[0],v1_in[1]);
    this->v1[0] = v1_in[0];
    this->v1[1] = v1_in[1];
    this->v1[2] = pz;
    if(std::isnan(pz))
    {
        is_nan_v[1]=1;
    }


    pz = rasterized_triangle->getXY(v2_in[0],v2_in[1]);
    this->v2[0] = v2_in[0];
    this->v2[1] = v2_in[1];
    this->v2[2] = pz;
    if(std::isnan(pz))
    {
        is_nan_v[2]=1;
    }
    if(is_nan_v[0] && is_nan_v[1] && is_nan_v[2])
    {
       is_nan = true;
    }

}
