#include "mesh.hpp"
#include "exception.hpp"

mesh::mesh(boost::shared_ptr<maw::matlab_engine> engine)
{
    _engine = engine;    
    _gfx = boost::make_shared<maw::graphics>(_engine.get());
}


mesh::~mesh()
{
    
    
}
void mesh::plot(std::string ID)
{
    mesh_hash::const_accessor a;
    bool r = _meshes.find(a,ID);
    if(!r)
    {
        BOOST_THROW_EXCEPTION( mesh_error() << errstr_info("Specified mesh not found: " + ID));
    }
    
    LOG_DEBUG << "Sending triangulation to matlab...";

    _engine->put_double_matrix("tri",a->second->matlab_tri_matrix());
    
    _engine->evaluate("ff=figure; set(gcf,'units','normalized','outerposition',[0 0 1 1]);");
    _engine->evaluate("set(ff,'Renderer','OpenGL')");
    double handle = _gfx->plot_patch("[mxDomain(:,1) mxDomain(:,2) mxDomain(:,3)]","tri","mxDomain(:,3)");
    
    _gfx->spin_until_close(handle);
    
    
}

void mesh::add_mesh(std::string file, std::string ID)
{
    _engine->evaluate(std::string("load ") + file);
    maw::d_mat xyz = _engine->get_double_matrix(file.substr(0,file.length()-4));
    
    //_engine->evaluate( std::string("clear ") + file.substr(0,file.length()-4) );
    
    LOG_DEBUG << "Creating triangulation for " + ID;
    
    boost::shared_ptr<triangulation> tri = boost::make_shared<triangulation>(_engine.get());
    
    auto x = xyz->unsafe_col(0);
    auto y = xyz->unsafe_col(1);
    auto z = xyz->unsafe_col(2);
    tri->create_delaunay(&x,&y,&z);
    
    LOG_DEBUG << "Sending triangulation to matlab...";

    _engine->put_double_matrix("tri",tri->matlab_tri_matrix());
                
    LOG_DEBUG << "Sending domain data to matlab...";
    _engine->put_double_matrix("mxDomain",xyz);

//    LOG_DEBUG << "Creating 3D bounding box...";
//    _engine->evaluate("[~,cornerpoints,~,~,~] = minboundbox(mxDomain(:,1),mxDomain(:,2),mxDomain(:,3))");
//    maw::d_mat cornerpoints = _engine->get_double_matrix("cornerpoints");
//    _engine->evaluate("clear cornerpoints");

//    _engine->evaluate("mxDomain=mxDomain(:,1:3)"); //because we are reading in the skyview data along with the dem, future calls assume a nx3 matrix.
//    std::cout << "Creating face normals...";
//    tri->compute_face_normals();
   
    mesh_hash::accessor a;
    bool r = _meshes.insert(a,ID);
    if(!r)
    {
        BOOST_THROW_EXCEPTION( mesh_error()
                            << errstr_info("Failed to insert mesh " + ID) 
                            );
    }
    a->second = tri;
    a.release();
    
//    _meshes.push_back(tri);
    

    
    

}

size_t mesh::hash_compare::hash( const std::string& x )
{
        boost::crc_32_type crc32;
        std::string xlower = boost::algorithm::to_lower_copy<std::string>(x);
        crc32.process_bytes(xlower.c_str(),xlower.length());

        return crc32.checksum();
}

bool  mesh::hash_compare::equal( const std::string& s1, const std::string& s2 )
{
        std::string ss1 = boost::algorithm::to_lower_copy<std::string>(s1);
        std::string ss2 = boost::algorithm::to_lower_copy<std::string>(s2);

        return ss1 == ss2;
}
