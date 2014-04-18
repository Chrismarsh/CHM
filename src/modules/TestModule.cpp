
#include "TestModule.hpp"

TestModule::TestModule( std::string ID)
{
    LOG_DEBUG << "Successfully instatiated module " << ID;
}
void TestModule::run(mesh_elem& elem, boost::shared_ptr<global> global_param)
{
    
    elem.add_face_data("rand",std::rand());
}

TestModule::~TestModule()
{
    
    
}