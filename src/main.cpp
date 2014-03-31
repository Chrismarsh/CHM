#include <iostream>
#include "mesh.h"
#include "core.h"



int main (int argc, char *argv[])
{

    try
    {
        //setup logging first
        
        core kernel;

        kernel.read_config_file("CHM.config") ;
    }
    catch( boost::exception& e)
    {
        std::cout << boost::diagnostic_information(e);
    }


	
// 	mesh terrain;
// 	terrain.read_mesh("mesh.config");


	return 0;
}
