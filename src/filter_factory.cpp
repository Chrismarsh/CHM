//
// Created by chris on 18/11/15.
//

#include "filter_factory.h"


filter_base* filter_factory::get(std::string ID, pt::ptree config)
{
    LOG_VERBOSE << "Filter ID=" << ID;

    filter_base* filter = nullptr;

    if (ID == "macdonald_undercatch")
        filter = new macdonald_undercatch();
    else if (ID == "scale_wind_speed")
        filter = new scale_wind_speed();

    if(filter == nullptr)
    {
        BOOST_THROW_EXCEPTION(module_not_found()
                              << errstr_info( std::string("Filter not found ") + ID)
        );
    }

    filter->ID = ID;
    filter->cfg = config;

    return filter;
}