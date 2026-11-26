Modules
=========

Overview
------------
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

3. A module may not ever write any other variable global variable
   which it does declare. It should also not overwrite the variables of
   another module.

Modules are specified to run in the :ref:`configuration:modules` section of the configuration file and configured as detailed in
the :ref:`configuration:config` section. For more details on how modules work, please see :ref:`development:modules`.

.. note::
   Unlike filters, modules are executed in an order to resolve the inter-module variable dependencies.

Available modules
---------------------
.. groups::
   :group: modules


By provides variable
------------------------

.. note::
   These provides are automatically produced from the module code. They do not include provides that are dynamically generated
    i.e.,  :code:`provides("t_"+some_var);`

.. provides::

Crocus snowpack
---------------

Select ``Crocus`` as the snow module in the CHM module list. Its lower boundary is a prescribed temperature in kelvin, constant across faces and time. The interface heat flux still responds to the temperature and properties of the basal snow layer. For example:

.. code-block:: json

   {
     "config": {
       "Crocus": {
         "ground_temperature_k": 263.0,
         "ground_conductivity": 1.5,
         "ground_layer_thickness_m": 0.1,
         "ground_albedo": 0.2,
         "max_layers": 20
       }
     }
   }

These are the defaults. Conductivity is in W m⁻¹ K⁻¹. Crocus allocates 28 bytes per configured layer and owned face for layer state. With 20 layers this is 560 bytes per face plus diagnostics. Forcing uses 2 m wind, air temperature and humidity, incoming shortwave and longwave, and separate rain and snow increments. Where a canopy module provides a subcanopy value it is used for that forcing independently.

Crocus supplies ``swe`` (mm), ``snowdepthavg`` (m perpendicular to the slope), ``snowdepthavg_vert`` (m vertical), ``Nsnow``, ``snow_albedo``, and the usual FSM energy and mass flux fields. ``H``, ``E``, and ``LWout`` / ``ilwr_out`` are W m⁻² positive upward. ``crocus_ground_heat_flux`` is positive toward the ground. ``melt_rate``, ``roff``, and ``subl`` are kg m⁻² s⁻¹, as in the FSM implementation. ``sum_snowpack_subl`` is kg m⁻². Outflow includes rain draining through the pack or onto bare ground. Trace snowfall counts toward SWE immediately and forms a layer once it exceeds the native threshold. The first three ``Tsnow`` and ``Sliq`` fields describe the surface layers. All four compatibility ``Tsoil`` fields report the prescribed temperature.

Start the first Crocus run with a bare column. Subsequent runs can resume a Crocus checkpoint with the same layer capacity and prescribed temperature. Checkpoints retain every layer's water equivalent, density, enthalpy, optical diameter, sphericity, grain history and age, plus unlayered snowfall, albedo, fluxes and cumulative sublimation.
