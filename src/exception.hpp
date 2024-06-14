//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//

#pragma once

#include "logger.hpp"
#include <boost/exception/all.hpp>
#include <boost/throw_exception.hpp>
#include <mutex>
#include <string>

struct exception_base: virtual std::exception, virtual boost::exception { };

// IO errors
struct io_error: virtual exception_base { };
struct file_read_error: virtual io_error { };
struct file_write_error: virtual io_error { };

// Generic CHM errors
struct chm_error: virtual exception_base { };
struct chm_done: virtual exception_base { };
struct model_init_error: virtual chm_error { };
struct not_impl: virtual chm_error { };

// Module errors
struct module_error : virtual exception_base { };
struct module_not_found : virtual module_error { };
struct no_modules_defined : virtual module_error {};
struct module_data_error : virtual module_error {};

//regex tokenizer errors
struct tokenizer_error : virtual exception_base{};
struct invalid_rexp : virtual tokenizer_error {} ;
struct bad_lexical_cast : virtual tokenizer_error {};

//forcing data errors
struct forcing_error : virtual exception_base{};
struct forcing_insertion_error : virtual forcing_error{};
struct forcing_lookup_error : virtual forcing_error{};
struct forcing_badcast : virtual forcing_error{};
struct forcing_no_regexmatch : virtual forcing_error{};
struct forcing_timestep_mismatch : virtual forcing_error{};
struct forcing_timestep_notfound : virtual forcing_error{};
struct forcing_no_stations : virtual forcing_error{};

//interpolation errors
struct interp_error : virtual exception_base{};
struct interp_unknown_type : virtual interp_error {};
struct interpolation_error : virtual interp_error{};

//configuration errors
struct config_error : virtual exception_base{};
struct missing_value_error : virtual exception_base{};

typedef boost::error_info<struct errstr_info_,std::string> errstr_info;

//mesh errors
struct mesh_error : virtual exception_base{};
struct mesh_insertion_error : virtual mesh_error{};
struct mesh_lookup_error : virtual mesh_error{};



//http://stackoverflow.com/questions/11828539/elegant-exceptionhandling-in-openmp
/**
    ompException e;

    #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        e.Run([&]{
                  // code that might throw
                  // ...
              });
    }
    e.Rethrow()
 */
class ompException
{
    std::exception_ptr Ptr;
    std::mutex Lock;
public:
    ompException() : Ptr(nullptr)
    {}

    ~ompException()
    {

    }

    void Rethrow()
    {
        if (this->Ptr)
            std::rethrow_exception(this->Ptr);
    }

    template <typename Function, typename... Parameters>
    void Run(Function f, Parameters... params)
    {
        try
        {
            f(params...);
        }
        catch (...)
        {
            std::unique_lock<std::mutex> guard(this->Lock);
            this->Ptr = std::current_exception();
        }
    }
};

// Convenience macro for exception throwing
#define CHM_THROW_EXCEPTION(exception_type,message) \
    BOOST_THROW_EXCEPTION( exception_type() << errstr_info( message ) )

