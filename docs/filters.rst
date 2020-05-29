Filters
========
Filters are a mechanism whereby the input forcing data can be modified
in some way prior to the model run. For example, this could be use to
apply a gauge undercatch to precip. Filters modify the data of a station
*in situ*.

Note! Filters run in the order defined in the configuration file.

Implementation
==============

Filters inherent from the base ``filter_base`` class.

It must impliment the ``process`` function which recives as input a
station.

.. code:: cpp

   void process(boost::shared_ptr<station> station);

process
=======

The timeseries data have a ``next()`` iterator so can be traversed in a
do…while loop.

An example for wind undercatch is as follows:

.. code:: cpp

       std::string var = cfg.get<std::string>("variable");
       do{
           double data = station->now().get(var);
           double u = station->now().get("u");
           data = data * (1.010 * exp(-0.09*u));
           station->now().set(var,data);
       }while(station->next());
