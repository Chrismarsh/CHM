#include <vtkXMLUnstructuredGridReader.h>
#include <vtkSmartPointer.h>
#include <vtkDataSetMapper.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <iostream>

#include <ogrsf_frmts.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
namespace pt = boost::property_tree;
int main ( int argc, char *argv[] )
{

  std::string filename = "granger0.vtu";//argv[1];
 
  //read all the data from the file
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  vtkSmartPointer<vtkUnstructuredGrid> mesh = vtkSmartPointer<vtkUnstructuredGrid>(reader->GetOutput());

  // std::cout << mesh->GetNumberOfCells() << std::endl;

 vtkCellData *cd = mesh->GetCellData();
      if (cd)
        {
        std::cout << " contains cell data with "
             << cd->GetNumberOfArrays()
             << " arrays." << std::endl;
        for (int i = 0; i < cd->GetNumberOfArrays(); i++)
          {
          std::cout << "\tArray " << i
               << " is named "
               << (cd->GetArrayName(i) ? cd->GetArrayName(i) : "NULL")
               << std::endl;
          }
            auto* z = cd->GetArray("Elevation");
          auto*  xy = mesh->GetPoint(2);
            std::cout << xy[0] <<","<<xy[1]<<","<<*z->GetTuple(0) << std::endl;


        }


    const char *pszDriverName = "ESRI Shapefile";
    OGRSFDriver *poDriver;

    OGRRegisterAll();

    poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(
            pszDriverName );


    OGRDataSource *poDS;

    poDS = poDriver->CreateDataSource( "granger.shp", NULL );
    if( poDS == NULL )
    {
        printf( "Creation of output file failed.\n" );
        exit( 1 );
    }

    OGRLayer *poLayer;

    poLayer = poDS->CreateLayer( "grangerpoly", NULL, wkbPolygon, NULL );
    if( poLayer == NULL )
    {
        printf( "Layer creation failed.\n" );
        exit( 1 );
    }

    OGRFieldDefn oField( "elevation", OFTString );
    oField.SetWidth(32);

    if( poLayer->CreateField( &oField ) != OGRERR_NONE )
    {
        printf( "Creating Name field failed.\n" );
        exit( 1 );
    }



    for(int i=0; i < mesh->GetNumberOfCells(); i++)
    {


//        auto* xy = mesh->GetPoint(i);


//        if (xy[0] != 0 && xy[1]!=0)
//        {
//            auto *z = cd->GetArray("Elevation")->GetTuple(i);


            OGRLinearRing *ring =  (OGRLinearRing*) OGRGeometryFactory::createGeometry(wkbLinearRing); ;

            auto v0 = mesh->GetCell(i)->GetPointId(0);
            auto v1 = mesh->GetCell(i)->GetPointId(1);
            auto v2 = mesh->GetCell(i)->GetPointId(2);

            ring->addPoint( mesh->GetPoint(v0)[0],mesh->GetPoint(v0)[1]);
            ring->addPoint( mesh->GetPoint(v1)[0],mesh->GetPoint(v1)[1]);
            ring->addPoint( mesh->GetPoint(v2)[0],mesh->GetPoint(v2)[1]);
            ring->addPoint( mesh->GetPoint(v0)[0],mesh->GetPoint(v0)[1]);
            ring->closeRings();

        std::cout << i << std::endl;

            OGRPolygon *poly = (OGRPolygon*) OGRGeometryFactory::createGeometry(wkbPolygon);
        poly->addRing(ring);
//            poly->addRing(ring);


            OGRFeature *poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
//            poFeature->SetField("elevation", *z);
            poFeature->SetGeometry(poly);

            poLayer->CreateFeature(poFeature);

//            if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
//            {
//                printf("Failed to create feature in shapefile.\n");
//                exit(1);
//            }

            OGRFeature::DestroyFeature(poFeature);
    }
//    }

    OGRDataSource::DestroyDataSource( poDS );

  return EXIT_SUCCESS;
}