#include <iostream>
#include <string>
#include "core.h"




int main (int argc, char *argv[])
{

    try
    {

        core kernel;

        kernel.init(argc, argv) ;
        
        kernel.run();
    }
    catch( boost::exception& e)
    {
       LOG_ERROR << boost::diagnostic_information(e);
    }



	return 0;
}
