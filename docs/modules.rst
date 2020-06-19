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