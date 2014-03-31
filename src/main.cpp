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
       LOG_ERROR << boost::diagnostic_information(e);
    }



	return 0;
}
