
Modules
---------

Modules contain the physical process representation. A principal
design goal of a module is that it may depend either upon some set of
variables produced by either other modules or on input forcing data.
Modules define a set of variables which it provides globally to other
modules.


1. All modules have pre-/post-conditions;

   Pre condition

   -  input forcing data or post-conditions from other modules;

   Post condition

   -  see pre condition;

   Variables

   -  provide global variables to other modules, but these
      variables do not change in other modules.

2. There are two types of modules:

   Forcing data interpolant

   -  depends upon point-scale input forcing data variables and
      interpolate these data onto every domain element;

   Standard process module

   -  depends only upon the output of interpolation modules or
      other modules’ output.

3. Parallelizations are offered in two ways, each module belongs to one
   of them:

   Data parallel

   -  point-scale models that are applied to every triangle;

   Domain parallel

   -  requires knowledge of surrounding mesh points.

   Parallelization process group modules with same parallel ty


A module may not ever write any other variable global variable
which it does declare. It should also not overwrite the variables of
another module.

There are two types of modules: 

- Forcing data interpolant 
- Standard process module

Forcing data interpolants (``src/modules/interp_met``) depend upon
point-scale input forcing data variables and interpolate these data onto
every domain element. Standard modules (``src/modules/*``) depend only upon
the output of the ``interp_met`` modules as well as other module
outputs.

All modules may either be data parallel or domain parallel. Data
parallel modules are point-scale models that are applied to every
triangle, e.g., snowmodel. Domain parallel modules are modules that require knowledge of
surrounding mesh points, e.g., blowing snow transport model.


Modules inherent from ``module_base.hpp`` and implement a standard
interface. In the most simple case, a module must have a constructor
which defines all variable dependencies and a run function.

Data parallel
~~~~~~~~~~~~~~

Data parallel modules implement a ``run`` method that takes as input a
single mesh element to operate upon. These modules are automatically made parallel by CHM. The main model loop
automatically schedules modules to execute in parallel. Domain parallel
modules may access the ``elem`` variable directly to get access to the
triangle element.

The constructor is used to set a module to be the correct parallel type: ``parallel::data``.

.. code:: cpp

   class data_parallel_example : public module_base
   {
   public:
       data_parallel_example(config_file cfg);
       ~data_parallel_example();
       void run(mesh_elem& face);
   }; 

   data_parallel_example::data_parallel_example(config_file cfg)
           :module_base("data_parallel_example", parallel::data, cfg)
   {
   ...
   }

Domain parallel
~~~~~~~~~~~~~~~~

Domain parallel modules implement a ``run`` method that takes the entire
mesh domain. The module must iterate over the faces of the domain to
gain access to each element. This may be done in parallel but must be
explicitly done by the module. The constructor specifies the paralle type: ``parallel::domain``.

.. code:: cpp

   class domain_parallel_example : public module_base
   {
   public:
       domain_parallel_example(config_file cfg);
       ~domain_parallel_example();
       void run(mesh& domain);
   }; 

   domain_parallel_example::domain_parallel_example(config_file cfg)
           :module_base("domain_parallel_example", parallel::domain, cfg)
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

Iteration over mesh
~~~~~~~~~~~~~~~~~~~~

Because the triangle iterators provided by CGAL have a non-deterministic
order, as well as being incompatible with OpenMP, the way to access the
i-th triangle is via

.. code:: cpp

   #pragma omp parallel for
       for (size_t i = 0; i < domain->size_faces(); i++)
       {
           auto elem = domain->face(i);
        ...
       }


init()
~~~~~~~~

In all cases a module may implement the ``init`` method.

.. code:: cpp

   void example_module::init(mesh& domain)

Regardless of if the module is data or domain parallel, this function
receives the entire mesh. ``init`` is called exactly once, after all
other model setup has occurred, but prior to the main model execution
loop. It is responsible for any initialization required by the model. 

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



Dependencies
~~~~~~~~~~~~

In the constructor, a module declares itself to ``provides`` a set of
variables and optionally ``depends`` upon other variables. Lastly, it
may ``optionally`` depend upon a variable. If the the variable is not
present, module dependency checks will still succeed, but the module
*must* check prior to access to avoid a segfault. 

.. code:: cpp

   # from another modules
   depends("ilwr");

   #optionally depend on another modules output
   optional("snow_albedo");

   #provide for another module.
   provides("dQ");

Conflicts
~~~~~~~~~~

Sometimes two modules absolutely should not be used together. The ``conflicts`` allows for specifying the name of a module to conflict against. 
When a conflict is detected, the setup stops. This should be used sparingly.

.. code:: cpp

   conflicts("snow_slide"); 



Variable access
~~~~~~~~~~~~~~~

Modules read from a variable stored on the mesh element via

.. code:: cpp

   auto albedo = (*elem)["snow_albedo"];

Modules may *only* write to variables they provide via

.. code:: cpp

   (*elem)["dQ"] = 100.0;

If ``optional`` has been used, a module can test for existance via

.. code:: cpp

    if(has_optional("snow_albedo"))
    {
       #do stuff
    }

Variable names
~~~~~~~~~~~~~~~

Variable access via the above variable access incurs some computational cost to convert the string to a hash for lookup in the underlying data-structure. If possible, suffix a variable name string with ``_s``. For example ``(*elem)["snow_albedo"_s]``. This will replace the string with a compile-time hash value, making the runtime lookup significantly faster. This can be done as long as the variable is known at compile time. For example if diagnostic output is done for *n* layers at run time

.. code:: cpp

   for(int i = 0; i < n; ++i)
   {
      provides("my_debug_layer_"+i);
   }

then these are ineligible for the ``_s`` suffix and speedup.

Registration with module factory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once the module has been written, it needs to be registered with the
module factory.

1. In the ``hpp`` file, within the class definition add
   ``REGISTER_MODULE_HPP(module_name);`` where ``module_name`` exactly
   matches the class name
2. In the ``cpp`` file, outside of all the other definitions add
   ``REGISTER_MODULE_CPP(module_name);``
3. All configuration options, use of the module, etc will be refered to
   as ``module_name`` in the config file.

Data storage
~~~~~~~~~~~~~

Frequently, the module must maintain a set of data that is separate from
the variables that are exposed to other modules (i.e., via ``provides``). These data can be stored in two ways: a) as
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

The ``make_module_data`` should be called in the ``init`` setup method.


interp_met modules
------------------

Meteorological interpolation functions are slightly different than the
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

The ``interpolation`` object abstracts the creation of different types of spatial interpolators.
Currently Inverse-Distance-Weighting (IDW) and Thin Plate Spline with
Tension (TPSwT) are implemented. The
interpolation method is chosen via the ``interp_alg`` enum. This is
passed to the constructor

.. code:: cpp

   interpolation::init(interp_alg ia, size_t size)

For performance reasons, it is best to initialize the interpolator in the ``init`` method. The size parameter should be used to denote the
number of locations to be used in the interpolation.

.. code:: cpp

   void test_module::init(mesh& domain)
   {

       #pragma omp parallel for
       for (size_t i = 0; i < domain->size_faces(); i++)
       {
           auto face = domain->face(i);
           auto d = face->make_module_data<const_llra_ta::data>(ID);
           d->interp.init(global_param->interp_algorithm,face->stations().size() );
       }
       LOG_DEBUG << "Successfully init module " << this->ID;

   }


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
   for (auto& s : face->stations())
   {
     if( is_nan((*s)["t"_s]))
         continue;
     double v = (*s)["t"_s] - lapse_rate * (0.0 - s->z());
     lowered_values.push_back( boost::make_tuple(s->x(), s->y(), v ) );
   }


   auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
   double value = face->get_module_data<data>(ID)->interp(lowered_values, query);

   //raise value back up to the face's elevation from sea level
   value =  value + lapse_rate * (0.0 - face->get_z());

   (*face)["t"_s]=value;


If the interpolant requires knowledge of the number of stations (e.g.,
TPSwT), and less stations are input (e.g., a NaN value is present), the the interpolant will on-the-fly
reinitialize itself with the new size.

Execution order
--------------------------

Inter-module dependencies, and thus the order
to run modules, is resolved during to run time. The order of module
execution is not dependent upon the order listed in the configuration
file. The interpolation modules always come prior to the process modules.

Inter-module variable dependencies is determined via the ``provides``
and ``depends`` declarations in the constructor. A module’s dependencies
are *every* other module that provides that output. This connectivity is
represented internally with a directed acyclic graph. Thus, the linear sequential
execution of the modules is determined via a topological sort.

If `Graphviz <http://www.graphviz.org/>`__ is installed, ``modules.pdf``
is generated which contains the graph of the inter-module dependencies.
|image0|

Once the executed order is determined, the modules are chunked into
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



.. |image0| image:: https://github.com/Chrismarsh/CHM/blob/master/modules_readme.png
