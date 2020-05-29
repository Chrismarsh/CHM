Modules
========

Modules are the short-hand for a process representation. A principal
design goal of a module is that it may depend either upon some set of
variables produced by either other modules or on input forcing data.
Modules define a set of variables which it provides globally to other
modules. A module may not ever write any other variable global variable
which it does declare. It should also not overwrite the variables of
another module.

There are two types of modules: + Forcing data interpolant + Standard
module

Forcing data interpolants (``modules/interp_met``) depend upon
point-scale input forcing data variables and interpolate these data onto
every domain element. Standard modules (``modules/*``) depend only upon
the output of the ``interp_met`` modules as well as other module
outputs.

All modules may either be data parallel or domain parallel. Data
parallel modules are point-scale models that are applied to every
triangle. Domain parallel modules are modules that require knowledge of
surrounding mesh points.

Implementation Details
======================

Modules inherent from ``module_base.hpp`` and implement a standard
interface. In the most simple case, a module must have a constructor
which defines all variable dependencies and a run function.

Data parallel
-------------

Data parallel modules implement a ``run`` function that takes as input a
single mesh element to operate upon. These modules do not need to
implement any special-sauce to be parallel. The main model loop
automatically schedules modules to execute in parallel. Domain parallel
modules may access the ``elem`` variable directly to get access to the
triangle element.

The constructor is used to set a module to be the correct parallel type.

.. code:: cpp

   class data_parallel_example : public module_base
   {
   public:
       data_parallel_example(config_file cfg);
       ~data_parallel_example();
       void run(mesh_elem& face, boost::shared_ptr<global> global_param);
   }; 

   data_parallel_example::data_parallel_example(config_file cfg)
           :module_base(parallel::data)
   {
   ...
   }

Domain parallel
---------------

Domain parallel modules implement a run function that takes the entire
mesh domain. The module must iterate over the faces of the domain to
gain access to each element. This may be done in parallel but must be
explicitly done by the module.

.. code:: cpp

   class domain_parallel_example : public module_base
   {
   public:
       domain_parallel_example(config_file cfg);
       ~domain_parallel_example();
       void run(mesh domain, boost::shared_ptr<global> global_param);
   }; 

   domain_parallel_example::domain_parallel_example(config_file cfg)
           :module_base(parallel::domain)
   {
   ...
   }

   void run(mesh domain, boost::shared_ptr<global> global_param)
   {
    #pragma omp parallel for
       for (size_t i = 0; i < domain->size_faces(); i++)
       {
           auto face = domain->face(i);
          /** do stuff with face **/
       }
   }

init()
------

In all cases a module may implement the ``init`` function.

.. code:: cpp

   void example_module::init(mesh domain, boost::shared_ptr<global> global_param)

Regardless of if the module is data or domain parallel, this function
receives the entire mesh. ``init`` is called exactly once, after all
other model setup has occurred, but prior to the main model execution
loop. It is responsible for any initialization required by the model. An
example of a complicated ``init`` is found in
`Liston_wind <https://github.com/Chrismarsh/CHM/blob/master/src/modules/interp_met/Liston_wind.cpp>`__
where the ``init`` function is used to pre-compute the wind curvature
function.

In some cases, a module may be able to work in either a domain parallel
or a data parallel mode with little modification. To avoid duplicating
code, a module may provide two ``run`` methods, one for each. Then, in
the ``init`` function, it can change the type of parallelism that is
declared. This is the only place where this change can be safely done.
To do so, both run interfaces are exposed:

::

       virtual void run(mesh domain);
       virtual void run(mesh_elem &face);

and then in ``init``, the module can query ``global`` as if CHM is in
point-mode. If not, it can safely switch to domain parallel. E.g.:

::

       if(!global_param->is_point_mode())
           _parallel_type =  parallel::domain;

``scale_wind_vert.cpp`` is an example of this.

Variables
---------

Dependencies
~~~~~~~~~~~~

In the constructor, a module declares itself to ``provides`` a set of
variables and optionally ``depends`` upon other variables. Lastly, it
may ``optionally`` depend upon a variable. If the the variable is not
present, module dependency checks will still succeed, but the module
*must* check prior to access to avoid a segfault. If a

.. code:: cpp

   # from another modules
   depends("ilwr");

   #optionally depend on another modules output
   optional("snow_albedo");

   #provide for another module.
   provides("dQ");

Variable access
~~~~~~~~~~~~~~~

Modules read from a variable stored on the mesh element via

.. code:: cpp

   auto albedo = elem->face_data("snow_albedo");

Modules may *only* write to variables they provide via

.. code:: cpp

   elem->set_face_data("dQ", 100.0);

If ``optional`` has been used, a module can test for existance via

.. code:: cppp

    if(has_optional("snow_albedo"))
       {
          #do stuff
       }
       else
       {
          #default behaviour?
       }

Data storage
------------

Frequently, the module must maintain a set of data that is separate from
the variables that are exposed to other modules with the
``set_face_data`` function. These data can be stored in two ways: a) as
a member variable in the module class; b) in a per-triangle data store.
If the data is stored as a member variable, this is global to every call
of the module and shared across the entire mesh. Remember, there is only
1 instance of a module class. To achieve per-triangle data storage, a
module should create a sub-class that inherants from ``face_info``

.. code:: cpp

   class test : public module_base
   {
    struct data : public face_info
       {
          double my_data;
       }
   };

This sub-class then should be initialized on each element using
``make_module_data``. As the class’ member variable ``ID`` is passed to
the call to create and access the data, other modules’ data is
technically available for access. *Don’t do this*.

.. code:: cpp

   auto d = face->make_module_data<test::data>(ID);  #returns the instance just created
   d->my_data = 5;

   #access later
   auto d = face->get_module_data<test::data>(ID); 

interp_met modules
------------------

Meterological interpolation functions are slightly different than the
above. They should all declare an interpolant in their per-face data
store. This must be on a per-element basis to ensure parallelism is
possible. If this is not done, large wait-locks must be used to prevent
the internal consistency of the linear systems. The other benefit of
this design is the interpolant is on a per-module basis, allowing each
module to use a different interpolant.

.. code:: cpp

       struct data : public face_info
       {
           interpolation interp;
       };

The ``interpolation`` object is a combo functor and factory. Its goal is
to abstract the creation of different types of spatial interpolators.
Currently Inverse-Distance-Weighting (IDW) and Thin Plate Spline with
Tension (TPSwT) are implemented (although only TPSwT is usable). The
interpolation method is choosable via the ``interp_alg`` enum. This is
passed to the constructor

.. code:: cpp

   interpolation::init(interp_alg ia, size_t size)

For performance reasons, it is best to initialize the interpolator prior
to the module running. The size parameter should be used to denote the
number of locations to be used in the interpolation.

The interpolation is performed by calling operator () on the
interpolation instance

.. code:: cpp

   operator()(std::vector< boost::tuple<double,double,double> >& sample_points, boost::tuple<double,double,double>& query_point)

where ``sample_points`` is a vector of (x,y,value) location tuples of
each input data. ``query_point`` is then the (x,y,z) location we wish to
interpolate. Frequently values cannot be interpolate directly and
requires lowering to a common reference level. An example of what this
looks like for constant temperature lapse rate is shown.

.. code:: cpp

       double lapse_rate = 0.0065;
       //lower all the station values to sea level prior to the interpolation
       std::vector< boost::tuple<double, double, double> > lowered_values;
       for (auto& s : global_param->stations)
       {
           if( is_nan(s->get("t")))
               continue;
           double v = s->get("t") - lapse_rate * (0.0 - s->z());
           lowered_values.push_back( boost::make_tuple(s->x(), s->y(), v ) );
       }

       auto query = boost::make_tuple(elem->get_x(), elem->get_y(), elem->get_z());
       double value = elem->get_module_data<data>(ID)->interp(lowered_values, query);

       //raise value back up to the face's elevation from sea level
       value =  value + lapse_rate * (0.0 - elem->get_z());

If the interpolant requires knowledge of the number of stations (e.g.,
TPSwT), and less stations are input (e.g., a NaN value is present – see
`timeseries <Timeseries>`__, the the interpolant will on-the-fly
reinitialize itself with the new size.

#Inter-module dependencies Intermodule dependencies, and thus the order
to run modules, is resolved during to run time. The order of module
execution is not dependent upon the order listed in the configuration
file (see `Configuration <Configuration>`__). The interpolation modules
always come prior to the process modules.

Inter-module variable dependencies is determined via the ``provides``
and ``depends`` declarations in the constructor. A module’s dependencies
are *every* other module that provides that output. This connectivity is
represented internally with a graph. Thus, the linear sequential
execution of the modules is determined via a topological sort.

If `Graphviz <http://www.graphviz.org/>`__ is installed, ``modules.pdf``
is generated which contains the graph of the inter-module dependencies.
|image0|

Once the linear order is determined, the modules are chunked into
execution groups depending on their ``parallel::`` flag. For example,
consider the following set of modules, sorted via the topological sort:

::

   mod_A (parallel::data)
   mod_B (parallel::data)
   mod_C (parallel::data)
   mod_D (parallel::domain)
   mod_E (parallel::data)

These are then chunked into 3 sub groups:

::

   mod_A (parallel::data)
   mod_B (parallel::data)
   mod_C (parallel::data)

::

   mod_D (parallel::domain)

::

   mod_E (parallel::data)

In the first data parallel subgroup, ``mod_A``, ``mod_B``, ``mod_C`` are
executed sequentially on each triangle, but each triangle is done in
parallel. Then subgroup 2 is run over the entire domain. Then subgroup 3
runs in parallel.

This purpose of this chunking is to attempt to schedule as many modules
as possible, to avoid the increase in overhead of running M modules over
N mesh points.

Registration with module factory
================================

Once the module has been written, it needs to be registered with the
module factory.

1. In the ``hpp`` file, within the class definition add
   ``REGISTER_MODULE_HPP(module_name);`` where ``module_name`` exactly
   matches the class name
2. In the ``cpp`` file, outside of all the other definitions add
   ``REGISTER_MODULE_CPP(module_name);``
3. All configuration options, use of the module, etc will be refered to
   as ``module_name`` in the config file.

.. |image0| image:: https://github.com/Chrismarsh/CHM/blob/master/modules_readme.png
