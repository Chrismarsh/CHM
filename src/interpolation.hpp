#pragma once


#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>

#include <math.h>
#include <iostream>



class interpolation_base
{
public:
       interpolation_base();
       virtual ~interpolation_base();
       void operator()=0;
};

class interpolation
{
public:
    interpolation();
    ~interpolation();
    boost::shared_ptr<interpolation_base> operator()(std::string method);
};

class IDW : interpolation_base
{
public:
    IDW();
    ~IDW();
    void operator()();
           
};

class spline : public interpolation_base
{
public:
    void operator()();
       
};



