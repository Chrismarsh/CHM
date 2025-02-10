#include "utility/gis.hpp"

namespace gis
{

    void bbox2geojson(double xmin, double ymin, double xmax, double ymax, const std::string& filepath)
    {
        GDALAllRegister();
        const char *pszDriverName = "GeoJSON";
        GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName);
        if (poDriver == nullptr)
        {
            CHM_THROW_EXCEPTION(chm_error, "GeoJSON driver not available");
        }

        GDALDataset *poDS = poDriver->Create(filepath.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
        if (poDS == nullptr)
        {
            CHM_THROW_EXCEPTION(chm_error, "Creation of output file failed");
        }

        OGRLayer *poLayer = poDS->CreateLayer("bbox", nullptr, wkbPolygon, nullptr);
        if (poLayer == nullptr)
        {
            CHM_THROW_EXCEPTION(chm_error, "Layer creation failed");
        }

        // Create a polygon from the bounding box
        OGRPolygon polygon;
        OGRLinearRing ring;

        // Add points in clockwise order
        ring.addPoint(xmin, ymin); // bottom left
        ring.addPoint(xmax, ymin); // bottom right
        ring.addPoint(xmax, ymax); // top right
        ring.addPoint(xmin, ymax); // top left
        ring.addPoint(xmin, ymin); // close the ring

        polygon.addRing(&ring);

        // Create and set the feature
        OGRFeature *poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
        poFeature->SetGeometry(&polygon);

        if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
        {
            OGRFeature::DestroyFeature(poFeature);
            GDALClose(poDS);
            CHM_THROW_EXCEPTION(chm_error, "Failed to create feature in GeoJSON");
        }

        OGRFeature::DestroyFeature(poFeature);
        GDALClose(poDS);
    }

    void xy2geojson(std::vector<std::tuple<float, float>> xy, std::string filepath, std::string proj4)
    {
        GDALAllRegister();
        const char *pszDriverName = "GeoJSON";
        GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName);
        if (poDriver == nullptr)
        {
            CHM_THROW_EXCEPTION(chm_error, "GeoJSON driver not available");
        }

        GDALDataset *poDS = poDriver->Create(filepath.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
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
                CHM_THROW_EXCEPTION(chm_error, "Failed to create feature in GeoJSON");
            }

            OGRFeature::DestroyFeature(poFeature);
        }
        GDALClose(poDS);
    }

}