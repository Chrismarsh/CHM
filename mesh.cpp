#include "mesh.hpp"

mesh::mesh()
{
    _engine = boost::make_shared<maw::matlab_engine>();
    _gfx = boost::make_shared<maw::graphics>();
    
    _engine->start();
    _engine->set_working_dir();
   
    LOG_DEBUG << "Matlab engine started";
    
}

mesh::~mesh()
{
    
    
}

void mesh::read_mesh(std::string file)
{
    std::ifstream is( file );

    //holds the top level object { ... }
    json_spirit::Value value;

    json_spirit::read( is, value );
 
    //get the top-level 
    const json_spirit::Object& top_level = value.get_obj();

    
    //loop over the top-level mesh list
    for(auto& itr : top_level )
    {
        const json_spirit::Pair& pair = itr;
        const std::string& name  = pair.name_;
        const json_spirit::Value&  value = pair.value_;
	
	std::cout << name <<  std::endl;
	
	//deal with the mesh enumeration
	if (name == "meshes")
	{
	  for(auto& jtr : value.get_obj())
	  {
	    _mesh_config m;
	    const json_spirit::Pair& pair = jtr;
   
	    if( pair.name_ == "")
	    {
	     std::cout << "empty mesh name!" << std::endl;
	     return;
	    }
	    std::cout <<  "Found mesh " << pair.name_ << std::endl;
	    m.name =  pair.name_;
	    
	    //will only be key-value pairs now
	    for(auto& ktr : pair.value_.get_obj())
	    {
		const json_spirit::Pair& pair = ktr;
		if(pair.name_ == "file")
		{
		  m.file = pair.value_.get_str();
		}
		std::cout << "\t" << pair.name_ << ":" <<pair.value_.get_str()  << std::endl;
	    }
	    _meshes.push_back(m);
	  }
	  
	}
    }
    

}
