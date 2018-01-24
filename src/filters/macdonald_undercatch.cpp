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
void macdonald_undercatch::init(boost::shared_ptr<station> station)
{
    //look at the config data to determine what we are modifying
    var = cfg.get<std::string>("variable");
}
void macdonald_undercatch::process(boost::shared_ptr<station> station)
{

    double data = station->now().get(var);
    double u = station->now().get("u");
    //trap missing data, just ignore it.
    if( !is_nan(data) && !is_nan(u))
    {
        data = data * (1.010 * exp(-0.09*u));
    } else
    {
        data = -9999;
    }

    station->now().set(var,data);

}