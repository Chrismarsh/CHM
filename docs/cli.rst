Command Line
=============

Command line options come in two flavours: short hand (e.g. -f) and long
(e.g., –config-file).

Options: 

   - ``--help``
   - ``--version``, ``-v``
   - ``--config-file``, ``-f``
   - ``--config``, ``-c``
   - ``--remove``, ``-r``
   - ``--remove-module``, ``-d``
   - ``--add-module``, ``-m``


In addition to specifying the configuration file to run with, it can be used to specific configuration options. Any configuration that is configurable via configuration files can be specified on the command line. This is done so that configuration files do not need to be written
during, for example, uncertainty analysis.

.. note::

    Values specified on the commandline over-ride existing values in the configuration file 

help
*****

Prints the available options.

version
********
Prints the version of CHM along with the git commit.

config-file
**************

Specifies a configuration file to use. This is required.

``./CHM --config-file CHM_test.json``

``./CHM -f CHM_test.json``


config
*******

This enables changing a configuration value from the command line. The value is specified with a fully qualified configuation path. 
This can change any aspect of the configuration file, but an example with a module's configuration is given below.

For example, to change the module ``moduleA`` value of ``b`` the configuration file would be

.. code:: json

   {

      "config":
      {
         "moduleA":
         {
            "b": 3.14
         }
      }

   }


Therefore, the commandline variant would be:

::

   -c config.moduleA.b:3.14

Multiple configuration values can be specified via multiple ``-c`` arguments:

::

   -c config.moduleA.b:3.14 -c option.debug_level:"error"

.. note::
   This cannot be used to modify the list values (``[ ]``) for the ``modules`` section. Please use the add/remove module feature to modify those.


remove
*******

Removes a configuration parameter. Removals are processed after
any ``-–config`` parameters are parsed, so ``-r`` will override ``-c``. A fully qualified path to the configuration is required. For example this configuration

.. code:: json

   {

      "option":
      {
         "option_a": true,
         "option_b": 1234,
      }

   }

could have ``option_a`` removed 

::

   -r option.option_b

Because the removal happens after any modification with ``-c``, the following to add ``option_c`` has no effect:

::

     -c option.option_c:2 -r option.option_c


.. note::
   This cannot be used to modify the list values (``[ ]``) for the ``modules`` section. Please use the add/remove module feature to modify those.


remove-module
***************
Removes a module

::

   -d Marsh_shading_iswr

add-module
**********

Adds a module to the list. Adding configuration options with ``-c`` can be done before or after this call.

::

   -m snobal -m Marsh_shading_iswr



