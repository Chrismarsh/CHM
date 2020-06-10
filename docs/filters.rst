Filters
========
Filters are a mechanism whereby the input forcing data can be modified
in some way prior to the model run. For example, this could be use to
apply a gauge undercatch to precip. Filters modify the data of a virtual station
*in situ*.

.. warning::

   Filters run in the order defined in the configuration file.

Implementation
---------------

Filters inherent from the base ``filter_base`` class.

It must impliment the ``process`` function which recives as input a
station.

.. code:: cpp

   void process(boost::shared_ptr<station> station);

init()
-------

``init`` can be used to determine what variable names should be used and accessed. For example,

.. code:: cpp

    var = cfg.get<std::string>("variable");
    fac = cfg.get<double>("factor");    


process
--------

The filter is given the current timestep to modify

.. code:: cpp

   double data = (*station)[var];
    if(!is_nan(data))
    {
         data = data + fac;
    }
    
    (*station)[var]=data;