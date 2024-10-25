#include "utility/gis.hpp"

namespace gis
{

    void xy2shp(std::vector<std::tuple<float, float>> xy, std::string shp_path, std::string proj4)
    {
        GDALAllRegister();
        const char *pszDriverName = "ESRI Shapefile";
        GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName);
        if (poDriver == nullptr)
        {
            CHM_THROW_EXCEPTION(chm_error, "Shapefile driver not available");
        }

        GDALDataset *poDS = poDriver->Create(shp_path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
        if (poDS == nullptr)
        {
            CHM_THROW_EXCEPTION(chm_error, "Creation of output file failed");
        }

        OGRLayer *poLayer = poDS->CreateLayer("layer", nullptr, wkbPoint, nullptr);
        if (poLayer == nullptr)
        {
            CHM_THROW_EXCEPTION(chm_error, "Layer creation failed");
        }

//        OGRSpatialReference oSRS;

//        if (oSRS.importFromProj4(proj4.c_str()) != OGRERR_NONE)
//        {
//            CHM_THROW_EXCEPTION(chm_error, "Failed to import PROJ.4 string");
//        }
//        if (poDS->SetProjection(proj4.c_str()) != CE_None)
//        {
//            CHM_THROW_EXCEPTION(chm_error, "Failed to set spatial reference");
//        }

        OGRFieldDefn oFieldX("X", OFTReal);
        oFieldX.SetWidth(32);
        poLayer->CreateField(&oFieldX);


        OGRFieldDefn oFieldY("Y", OFTReal);
        oFieldY.SetWidth(32);
        poLayer->CreateField(&oFieldY);

//        double xCoordinates[] = {10.0, 20.0, 30.0}; // Example x coordinates
//        double yCoordinates[] = {40.0, 50.0, 60.0}; // Example y coordinates
        int numPoints = xy.size();

        for (auto itr : xy)
        {
            float xCoordinates = std::get<0>(itr);
            float yCoordinates = std::get<1>(itr);

            OGRFeature *poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());

            poFeature->SetField("X", xCoordinates);
            poFeature->SetField("Y", yCoordinates);

            OGRPoint pt;
            pt.setX(xCoordinates);
            pt.setY(yCoordinates);
            poFeature->SetGeometry(&pt);

            if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
            {
                CHM_THROW_EXCEPTION(chm_error, "Failed to create feature in shapefile");
            }

            OGRFeature::DestroyFeature(poFeature);
        }
        GDALClose(poDS);
    }

}