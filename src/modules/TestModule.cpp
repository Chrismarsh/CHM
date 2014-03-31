
#include "TestModule.hpp"

TestModule::TestModule( std::string ID)
{
    sev_logger& lg = logger::get();
    BOOST_LOG_SEV(lg,log_level::debug) << "Successfully instatiated module " << ID;
}
void TestModule::run()
{
    sev_logger& lg = logger::get();
    BOOST_LOG_SEV(lg,log_level::debug) << "Hello world from Test Module";
}

TestModule::~TestModule()
{
    
    
}