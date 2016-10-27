#include <iostream>
#include <string>

#define BOOST_SPIRIT_THREADSAFE

#include "core.h"




int main (int argc, char *argv[])
{
    core kernel;

    try
    {
        kernel.init(argc, argv) ;
        
        kernel.run();

        kernel.end();
    }
    catch( boost::exception& e)
    {
        kernel.end(); //make sure endwin() is called so ncurses cleans up

        LOG_ERROR << boost::diagnostic_information(e);

        return -1;
    }

	return 0;
}
