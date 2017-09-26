//
// Created by chris on 18/11/15.
//

#include "macdonald_undercatch.h"


macdonald_undercatch::macdonald_undercatch()
{

}
macdonald_undercatch::~macdonald_undercatch()
{

}
void macdonald_undercatch::process(boost::shared_ptr<station> station)
{

    //look at the config data to determine what we are modifying
    std::string var = cfg.get<std::string>("variable");

    do{
        double data = station->now().get(var);
        //trap missing data, just ignore it.
        if( !is_nan(data))
        {
            double u = station->now().get("u");
            data = data * (1.010 * exp(-0.09*u));
        }

        station->now().set(var,data);
    }while(station->next());
}