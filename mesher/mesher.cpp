#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <fstream>

#include "mesh.h"
#include "version.h"

#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "ogrsf_frmts.h"

#include "raster.h"

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
    std::string version = "mesher 0.1 " GIT_BRANCH "/" GIT_COMMIT_HASH ;
    std::string poly_file;
    double max_area = 0;
    double min_area = 1;

    std::string error_metric = "rmse"; //default of RMSE
    //holds all rasters we perform tolerance checking on
    std::vector< std::pair< boost::shared_ptr<raster>,double> > rasters;

    //these are category based rasters, e.g., landcover, so use a fractional % to figure out if we should split on the tri.
    std::vector< std::pair< boost::shared_ptr<raster> ,double> > category_rasters;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "This message")
            ("version,v", "Version number")
            ("poly-file,p", po::value<std::string>(&poly_file),
             "PLGS file to use to bound triangulation. Same format as Triangle.")
            ("raster,r", po::value<std::vector<std::string>>(), "If tolerance checking is used,"
                    "this provides a list of the rasters to provide the tolerance checking against. "
                    "Order needs to match the order the tolerances are given in.")
            ("tolerance,t", po::value<std::vector<double>>(), "Tolerances, same units as the method"
                    " Must be give in the same order as the rasters.")

            ("category-raster,R", po::value<std::vector<std::string>>(), "Optional landcover raster to conform mesh to.")
            ("category-frac,T",  po::value<std::vector<double>>(), "Franctional percent of continous landcover required to not-split a triangle.")

            ("area,a", po::value<double>(&max_area), "Maximum area a triangle can be. Square unit.")
            ("min-area,m", po::value<double>(&min_area), "Minimum area a triangle can be. Square unit.")
            ("error-metric,M", po::value<std::string>(&error_metric), "Error metric. One of: rmse, mean_tol, max_tol."
                                                         "mean_tol compares the mean triangle vertex value to the mean raster value. "
                                                        "max_tol mimics the ArcGIS TIN tolerance, and is the maximum difference between the triangle and any single raster cell.");

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        exit(1);
    }
    if (vm.count("version"))
    {
        std::cout << version << std::endl;
        exit(1);
    }


    if ( (vm.count("raster") && ! vm.count("tolerance")) ||
            (!vm.count("raster") &&  vm.count("tolerance")) )
    {
        std::cout << "Both raster and tolerance must be specified!" << std::endl;
        exit(1);
    }

    if ( (vm.count("category-raster") && ! vm.count("category-frac")) ||
         (!vm.count("category-raster") &&  vm.count("category-frac")) )
    {
        std::cout << "Both category rasters and fractions must be specified!" << std::endl;
        exit(1);
    }

    if(max_area <= 0)
    {
        std::cout << "Area must be positive" << std::endl;
        exit(1);
    }

    if(min_area <= 0 )
    {
        std::cout << "Min area must be greater than zero" << std::endl;
        exit(1);
    }


    if( vm.count("raster") && vm.count("tolerance"))
    {
        auto files = vm["raster"].as<std::vector<std::string>>();
        auto tols = vm["tolerance"].as<std::vector<double>>();

        if(files.size() != tols.size())
        {
            std::cout << "Mismatched lengths in rasters and tolerances. Must be equale."<<std::endl;
            exit(1);
        }
        for(int i = 0 ; i< tols.size();++i)
        {
            auto r = boost::make_shared<raster>();
            r->open(files.at(i));
            rasters.push_back(std::make_pair( r,tols.at(i) ));
        }
    }

    if( vm.count("category-raster") && vm.count("category-frac"))
    {
        auto files = vm["category-raster"].as<std::vector<std::string>>();
        auto tols = vm["category-frac"].as<std::vector<double>>();

        if(files.size() != tols.size())
        {
            std::cout << "Mismatched lengths in category-raster and category-frac. Must be equale."<<std::endl;
            exit(1);
        }
        for(int i = 0 ; i< tols.size();++i)
        {
            auto r = boost::make_shared<raster>();
            r->open(files.at(i));

            if(tols.at(i) > 1)
            {
                std::cout << "Fractional percentage required for category-frac!"<<std::endl;
                exit(1);
            }
            category_rasters.push_back(std::make_pair( r,tols.at(i) ));
        }
    }

    CDT cdt;


    std::ifstream infile(poly_file);
    if(!infile)
    {
        std::cout << "Failed to open poly file" << std::endl;
        exit(1);
    }

    boost::filesystem::path p(poly_file);
    auto path = p.parent_path();

    std::string line;
    //header
    std::getline(infile, line);
    std::istringstream iss(line);
    int row,col;
    iss >> row >> col;
    std::cout << std::setprecision(20);
    std::vector<Vertex_handle> vertex;
    if (col != 2)
    {
        std::cout << "Wrong number of columns" <<std::endl;
        exit(1);
    }
    for(int i = 0 ; i < row -1 ; i++)
    {
        std::getline(infile, line);
        std::istringstream iss(line);
        double id,x,y;
        if (!(iss >> id >> x >> y))
        {
            std::cout << "Error on line " << i << std::endl;
            exit(1);
        }

        Vertex_handle vh = cdt.insert(Point(x,y));
        vh->info()=i;
        vertex.push_back(vh);
    }

    //skip last line of vertexes. this is a duplicate, so don't need it (Triangle does though)
    std::getline(infile, line);

    //empty line
    std::getline(infile, line);

    //header
    std::getline(infile, line);
    std::istringstream header2(line);
    double special;
    header2 >> row >> special;

    //ignoring the last 2 items as they are for Triangle, but we aren't using it here
    for(int i=0;i<row-2; i++)
    {
        std::getline(infile, line);

        int rowid,i0,i1;
        std::istringstream iss(line);
        if(! (iss >> rowid >> i0 >> i1))
        {
            std::cout << "Error in section 2 on line " << rowid << std::endl;
            exit(1);
        }
        //1 indexed,fix
        --i0;
        --i1;
        Vertex_handle v0,v1;
        v0 = vertex.at(i0);
        v1 = vertex.at(i1);

        cdt.insert_constraint(v0,v1);
    }

    //and lastly, close the loop, ignore
    Vertex_handle v0,v1;
    v0 = vertex.at(row-2);
    v1 = vertex.at(0);
    cdt.insert_constraint(v0,v1);

    std::cout << "Number of input PLGS vertices: " << cdt.number_of_vertices() << std::endl;
    std::cout << "Meshing the triangulation..." << std::endl;
    CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125,max_area,min_area,rasters,category_rasters,error_metric));


    auto nodefilepath = path; //eg PLGSwolf_lidar1.1.node
    nodefilepath /= p.filename().replace_extension(".1.node");
    std::ofstream nodefile( nodefilepath.string());

    //header
    nodefile << cdt.number_of_vertices() << " 2 0 1" << std::endl;
    nodefile << std::setprecision(15);

    size_t mesh_vertex_i=1;

    for(auto itr = cdt.finite_vertices_begin(); itr != cdt.finite_vertices_end(); ++itr)
    {

        itr->info() = mesh_vertex_i;

        nodefile << mesh_vertex_i << "   " << itr->point().x() << "   " << itr->point().y() << "   0" << std::endl;
        ++mesh_vertex_i;
    }


    nodefile.close();


    size_t elem_i=1; //1 indexing

    for(auto itr = cdt.finite_faces_begin(); itr != cdt.finite_faces_end(); itr++ )
    {
        if(itr->is_in_domain())
        {
            itr->id = elem_i;
            ++elem_i;
        }

    }

    //i has total number of faces that are in domain
    std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
    std::cout << "Number of triangles: " << elem_i << std::endl;

    //niehgbour file
    auto neighfilepath = path; //eg PLGSwolf_lidar1.1.neigh
    neighfilepath /= p.filename().replace_extension(".1.neigh");
    std::ofstream neighfile(neighfilepath.string());
    neighfile << elem_i-1 << " 3" << std::endl;

       //element file, defines the triangles
    auto elefilepath = path; //eg "PLGSwolf_lidar1.1.ele"
    elefilepath /= p.filename().replace_extension(".1.ele");
    std::ofstream elemfile(elefilepath.string());
    elemfile << elem_i-1 << " 3 0" << std::endl;


    int i=1;

    for(auto itr = cdt.finite_faces_begin(); itr != cdt.finite_faces_end(); itr++ )
    {
        if(itr->is_in_domain())
        {
            size_t v0 = itr->vertex(0)->info();
            size_t v1 = itr->vertex(1)->info();
            size_t v2 = itr->vertex(2)->info();

            elemfile << i << "    " << v0 << "    " << v1 << "    "<< v2 << std::endl;

            auto n0 = itr->neighbor(0); cdt.is_infinite(n0) ? NULL : n0; //only want finite nieghbours
            auto n1 = itr->neighbor(1); cdt.is_infinite(n1) ? NULL : n1;
            auto n2 = itr->neighbor(2); cdt.is_infinite(n2) ? NULL : n2;

            neighfile << i <<  "  " << (n0 != NULL ? n0->id : -1) << "  " << (n1 != NULL ? n1->id : -1) <<"  "<< (n1 != NULL ? n1->id : -1) << std::endl;
            ++i;
        }

    }

    elemfile.close();
    neighfile.close();
    return 0;
}

