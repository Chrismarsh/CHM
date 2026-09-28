How To: Submodule Framework
============================

The CHM modules, for good reason, allow the user to design their module however they'd like as long as they couple the calculations to the CHM coupler that coordinates the modules and the data passing between them. Modules can be subdivided into two categories:

1. External models coupled to CHM.
2. Models coded in CHM directly.

Being able to accommodate both is one of CHMs biggest strengths (the same can be said for most modular models). Modules of the first type are only possible when the models are reasonably portable and adaptable to any framework. It can be readily seen by reading these modules that a lot of work is required in order to make it mesh (pun intended) with CHM. And ultimately that is not surprising. It was not designed with CHM in mind. However, given that we have full control over CHM modules, it makes sense to design a consistent framework for writing them. And that framework should go beyond implementing the virtual functions ``void module_base::run(mesh_elem&)`` and ``void module_base::init(mesh&)``, and using a class derived from ``face_info`` to store per-triangle data.

If no framework - or *encouraged* - design exists, then modules will either become so completely different that learning each is equivalent to learning a new model each time, or they will be coded in such a manner that testability, extensibility, flexibility, and readability fall to the back burner.

In this document, I'll present basic instructions and guidelines for designing submodules and then how to write (or modify) a module that uses them.

Designing Submodules
--------------------

Avoiding the dirty details of the why, here is the how.

Step 1
^^^^^^

Open a new file in ``src/modules/submodules`` with a descriptive name. For this we will use ``submodule_example.hpp``. In this file, immediately create a namespace with the same name as the file and nest inside of it a class called ``Model`` derived from ``base_step<Model<Data>,Data>`` where ``Data`` is a template parameter. More on this later...

The ``Model`` class should have one public function: ``execute_impl(Data& d)``. It will look something like this:

.. code-block:: cpp

   #include "base_step.hpp"

   namespace submodule_example
   {
       template<typename Data>
       class Model : public base_step<Model<Data>,Data>
       {
       public:
           void execute_impl(Data& d);
       };
   };

Step 2
^^^^^^

Take stock of your calculation and write out a list of inputs it needs and outputs it provides (similar to CHM). Then write a concept to enforce that ``Data`` has getters and setters for both of these.

.. code-block:: cpp

   #include "base_step.hpp"
   #include <concepts>

   namespace submodule_example
   {
       template<typename T>
       concept submodule_data = requires(T& t)
       {
           { t.input1(); } -> std::floating_point;
           { t.input2(); } -> std::floating_point;

           { t.output1(std::declval<double>()); } -> std::same_as<void>;
           { t.output2(std::declval<double>()); } -> std::same_as<void>;
       };
       
       template<submodule_data Data>
       class Model : public base_step<Model<Data>,Data>
       {
       public:
           void execute_impl(Data& d);
       };
   };

I've left out implementing ``execute_impl(Data&)`` on purpose for now.

Step 3 (optional)
^^^^^^^^^^^^^^^^^

Immediately write tests. Navigate to ``src/tests/submoduletests/`` and create a file with a good name like ``test_submodule_example.cpp``. The body will look something like this:

.. code-block:: cpp

   #include <gtest/gtest.h>
   #include "submodule_example.hpp"

   class data_for_test
   {
       // Implement example data class to satisfy the concept above  
   };

   class SubmoduleExampleTest : public ::testing::Test
   {
   protected: // Must be protected!!!
       submodule_example<data> submodule;
   };

   TEST_F(SubmoduleExampleTest,TestZerosReturnsZeros)
   {
       // TEST_F is a macro to create a class derived from SubmoduleExampleTest. So the `submodule` object already exists.
       // Operate on that are use macros like EXPECT_EQ(expected,produced); to check if the submodule is working as expected.
       // Test BEHAVIOURS.
   }

If you are having trouble designing a test as the result of designing increasingly contrived situations to test obscure edge cases that will actually compile, or you feel yourself wanting to test private methods, this is a `code smell <https://en.wikipedia.org/wiki/Code_smell>`_. I would recommend then further encapsulating these calculations into another class, defined within the ``submodule_example`` namespace, with a public, testable interface. Then test that. If you had an ``std::vector`` private member, you wouldn't test the ``std::vector`` through the public interface, you'd assume (rightly) it was tested elsewhere.

Why write tests first
^^^^^^^^^^^^^^^^^^^^^

It's called Test Driven Development. It exists to speed up development by spending time early on tests and later benefiting from the pre-planning. A common expression I've read is "If you have a bug, that means there is a missing test".

Step 4
^^^^^^

Return to ``submodule_example.hpp`` and iterate until it passes your tests. Tada! Now you know its working. Don't forget to add ``test_submodule_example.cpp`` to your ``CMakeLists.txt`` and enable the testing flag. I prefer to compile once, then open ``CMakeCache.txt`` in my build directory and change ``BUILD_TESTS`` to ``ON``, then recompile. It saves having to deal with changing the ``CMakeLists.txt`` every time while trying to commit changes to the actual source later.

**You're done!** You've successfully written a submodule. The next step is to include it in a CHM module.

Designing Modules with Submodules
---------------------------------

For simplicity, we will assume that we are starting from scratch to build this module, but these instructions can easily be applied to an existing module, either to change a calculation to a submodule or add a submodule as a separate calculation. Write the header as normal with a nice an descriptive name, here we will use ``module_example.hpp``. The way we do this is make use of the class ``data_base`` for automatic caching and some other tools.

Step 1: Write the Cache
^^^^^^^^^^^^^^^^^^^^^^^

First, you need a struct (classes are OK) to act as the cache. The short of it is that you require a member for outputs and a member for inputs. Therefore our Cache is

.. code-block:: cpp

   struct Cache : public cache_base
   {
       double input1 = default_value<double>();
       double input2 = default_value<double>();

       double output1 = 0.0;
       double output2 = 0.0;
   };

``default_value<T>()`` comes from the ``cache_base`` parent class. The details don't matter too much here, but it initializes the input as a sentinel value so that each member can separately be lazy-initialized.

Step 2: Write the data class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``data`` class is the class that holds the per-triangle data for a module between time steps. This is the class that will be passed as a template parameter to the submodule. Therefore, it must satisfy the concept we defined in ``submodule_example.hpp``. In order to use the cache, it also must be derived from ``data_base<Cache>``. The Cache as a template parameter sets the type of the cache used behind the scenes.

Lazy initialization of the input variables on the cache is only possible if ``data`` has a reference to the ``mesh_elem`` object it corresponds to, the ``config_file`` object, and the ``global`` shared pointer. So ``data_base<T>`` has a constructor to pass these references. therefore, the data class is as follows:

.. code-block:: cpp

   class data : data_base<Cache>
   {
   public:
       double input1() const;
       double input2() const;

       double output1(const double out_val) const;
       double output2(const double out_val) const;

       using data_base<Cache>::data_base;
       // Alternatively, you can write out the constructor as follows:
       // data(mesh_elem& face_in, std::shared_ptr<global> param, config_file& cfg) 
       // :  data_base<Cache>(face_in,param,cfg);
   };

Step 3: Write the Module header
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Put it all together. I recommend placing the ``data`` and ``Cache`` classes inside the body of the module. Like so:

.. code-block:: cpp

   // Essential CHM includes here
   #include "submodule_example.hpp"

   class module_example : public module_base
   {
   public:
       void run(mesh_elem&) override;
       void init(mesh&) override;

       struct Cache : public cache_base
       {
           double input1 = default_value<double>();
           double input2 = default_value<double>();

           double output1 = 0.0;
           double output2 = 0.0;
       };
       
       class data : data_base<Cache>
       {
       public:
           double input1() const;
           double input2() const;

           double output1(const double out_val) const;
           double output2(const double out_val) const;

           using data_base<Cache>::data_base;
       };
   private:
       submodule_example<data> submodule_obj;
   }

Step 4: Implement the ``data`` member functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Switching to ``module_example.cpp``, we implement the auto-caching using the protected member function from the ``data_base`` class for the inputs and another protected member for outputs. These functions automatically lazy initialize the cache. I will admit here that this is possibly the ugliest part of the implementation and could use improvements. Nonetheless, here is how you do it. You call the function ``update_value`` which expects two arguments, each argument must be a callable. I intend for these to be lambda functions but it may work with function pointers, but I have not tested that. The first lambda returns a reference to the cache member to be accessed, and the second lambda returns a copy of the dereferenced ``face`` object. ``update_value`` calls the second lambda if the cache has not be initialized yet or if the input variable is equal to ``default_value<double>()`` that we initialized the inputs to in the ``module_example::Cache``.

.. code-block:: cpp

   double module_example::data::input1()
   {
       update_value( [this]() -> auto& { return cache_->input1; },
                     [this]() { return (*face)["input1"_s];}
                   );

       return cache_->input1;
   }

   double module_example::data::input2()
   {
       update_value( [this]() -> auto& { return cache_->input2; },
                     [this]() { return (*face)["input2"_s];}
                   );

       return cache_->input2;
   }

Note that in all four lambdas, ``this`` is captured. ``This`` allows the use of ``face`` and ``cache_`` objects. I means that both lambdas have a copy of the ``this`` pointer (not the whole class) so the overhead is minimal and only stored during the call to ``update_value`` and freed immeidately after the function call.

``update_value`` is a templated function and therefore the compiler generates a unique function for each unique pair of lambdas. In addition, the arguments are r-value references so these temporary lambdas are passes correctly. The first lambda **must** return a reference. Luckily, concepts hidden in ``data_base`` ensure this is the case, so it won't compile if you forget. This allows the first lambda to be assigned.

Now for the output functions, its very similar and uses a function called ``set_output`` that likewise accepts two arguments but in this case, only the first argument is a lambda reference to the cache member, and the second is simply forwarding the function argument.

.. code-block:: cpp

   double module_example::data::output1(const double out_val)
   {
       set_output( [this]() -> auto& { return cache_->output1; }, out_val);  
   };

   double module_example::data::output2(const double out_val)
   {
       set_output( [this]() -> auto& { return cache_->output2; }, out_val);  
   };

Other kinds of inputs
^^^^^^^^^^^^^^^^^^^^^

Not all inputs are created equal. Some are variables, accessed by a dereference from ``face``, but others might be per-triangle or domain-wide parameters. You can approach this however you'd like. But remember that ``cfg`` and ``global`` and ``face`` are all available. So... make use of them.

For domain-wide parameters, I recommend adding a static variable to the data class, and assigning it in the ``module_example::init`` function, like so:

.. code-block:: cpp

   // .hpp

   data : public data_base<Cache>
   {
   public:
       // other stuff
       static double param1;
   };

   // .cpp

   void module_example::init(mesh& domain)
   {
       data::param1 = cfg.get<double>("param1");
   };

   static double module_example::data::param1() // Doesn't modify the state and only uses static members, so this can be static.
   {
       return param1;
   };

   // Another option
   double module_example::data::param2()
   {
       //thread safe
       //initialized once and then never again
       //techincally has a slight overhead to check if its been set yet.
       // Doesn't pollute the public interface of the data class.
       static double param2 = cfg.get<double>("param2");
       return param2;
   }

For triangle specific parameters, I would proceed exactly the same way as ``param1`` without the static keyword.

.. code-block:: cpp

   // .hpp

   data : public data_base<Cache>
   {
   public:
       // other stuff
       double param1;
   };

   // .cpp

   void module_example::init(mesh& domain)
   {
       // loop over every face
       // Get the local instance of data via make_module_data, skipped here because I've modified it.
       param1 = face->veg_parameter("param1"_s);
   };

   double module_example::data::param1() 
   {
       return param1;
   };

   // Another option
   // Static variagble in a functin cannot work because 
   // it is per CLASS not per instance and you need a unique value per triangle
   // However, you could do the following:

   module_example::data::data(mesh_elem& face_in,std::shared_ptr<global> param,config_file& cfg)
       : data_base<Cache>(face_in,param,cfg)
   {
       param1 = face->veg_parameter("param1"_s);
       param2 = cfg.get("param2"_s); // here param2 is a static member of data.
   };

In fact, that last example of setting them in the constructor is a bit genius and I actually haven't done it yet myself... welp guess I'll start doing that.

The final type is something accessed via the ``global`` instance ``param``. These will often have to be called dynamically because things change. If it is a time step, set it in the constructor for ``data``, but if it is the time of day then one must instead set it at each call. If one wants to avoid calling ``global`` more than once per time step, feel free to add it as a cached member of ``Cache``.

.. code-block:: cpp

   double module_example::data::hour()
   {
       // or whatever other function
       return param->hour();
   }

Step 5: Implement the module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Most proceeds as normal. Constructor sets the depends/provides, destructor is empty unless you were evil and used ``new``, be sure to ``delete`` those.

If we use the method to set parameters in the constructor, the ``.cpp`` file will look like this:

.. code-block:: cpp

   #include "module_example.hpp"
   REGISTER_MODULE_CPP(module_example);

   module_example::module_example(config_file cfg) : module_base("module_example", parallel::data, cfg)
   {
       depends("input1");
       depends("input2");

       provides("output1");
       provides("output2");
   }

   //empty deconstructor

   void module_example::init(mesh& domain)
   {
       // Recommend setting these outside the loop and not in the constructor.
       // You create it so many times, it is unwise to set it each time, but 
       // overhead might be minimal
       data::static_parameter = cfg.get<double>("static_parameter");
       for (size_t i = 0; i < domain->size_local_faces(); ++i)
       {
           auto face = domain->face(i);
           face->make_module_data<module_example::data>(ID,face,global_param,cfg);
           // Everything is set in the constructor... so no need do actually store a reference
           // there may be some situations where a parameter determines how something is set
           // In that case, store a reference to data and make whatever assignments you'd like.
       }
   }

   void module_example::run(mesh_elem& face)
   {
       auto& d = face->get_module_data<module_example::data>(ID);

       submodule_example(d);

       d.set_module_outputs();
   }

   void module_example::data::set_module_outputs()
   {
       auto& c = get_cache();
       if (c)
       {
           (*face)["output1"_s] = (c->output1 == default_value<double>()) ? c.output1 : 0.0;
           (*face)["output2"_s] = (c->output2 == default_value<double>()) ? c.output2 : 0.0;
       }
       else
       {
           (*face)["output1"_s] = 0.0;
           (*face)["output2"_s] = 0.0;
       }
   }

And now your module using submodules is complete! Notice that a function called ``make_module_data`` is called in the module virtual function ``init`` and its constructor now accepts more than just the module ``ID`` as an argument. I've rewritten just the ``make_module_data`` function to have a template of variadic parameters to pass to the constructor of the ``data`` object. 

Final Notes
-----------

There are a few things that could be improved:

1. The ``update_value`` and ``set_output`` protected ``data_base`` member functions. They are a bit ugly.

   One idea is to set these lambdas in the constructor and store them as std::functions. Not ideal really but doable.

2. The style of setting input cache members is not ideal. ``double input1 = default_value<double>();`` is redundant and ugly.

   One idea here is to create an ``input`` class that takes a parameter to a primitive type and auto sets it to ``default_value<double>()``. Like so:

   .. code-block:: cpp

      class Cache : public cache_base
      {
          input<double> input1;
          input<double> input2;

          output<double> output1; // initialized to 0.0
          output<double> output2;
      }

   Perhaps ``input`` and ``output`` classes are where we could store lambdas for note 1 above?

3. When you have many state variables, I recommend defining a struct like:

   .. code-block:: cpp

      namespace submodule_example
      {
          struct State
          {
              //state members here
          };

          // modify concept to expect a get_state function as a reference.
          // Store a State object privately on data.
      }

Then, if during the testing you decide to break-up your submodule into several subclasses you can design them to operate with the ``State`` object, rather than on the ``data`` object. This nifty trick means you can test these classes (defined in the ``submodule_example`` namespace) without having to worry about the concept to enforce on the template parameter of ``submodule_example``. This also means that these classes can be compiled in a ``submodule_example.cpp`` file. 

4. Consider a situation with N modules. Now you'll have a very busy ``data`` class. Depending on what functions are enforced by their respective concepts, you could have multiple functions for air temperature: ``air_temp()``, ``air_temperature()``, ``temp()``, and so on. For these, implement one and have the others simply pass through to the implemented version.

   Likewise, only store a single cache member for all the air temperatures! Unit conversions can be done in the respective calls, if one module expects Celsius but the other Kelvin!
