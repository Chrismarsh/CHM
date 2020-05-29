CHM can be controlled on the command line, and anything that is
configurable via configuration files can be specified on the command
line. This is done so that configuration files do not need to be written
during, for example, uncertainty analysis.

Command line options come in two flavours: short hand (e.g. -f) and long
(e.g., –config-file).

Options: - help - version, v - config-file, f - config, c - remove, r -
remove-module, d - add-module, m

config-file
===========

``./CHM --config-file CHM_test.json``

``./CHM -f CHM_test.json``

#config Specifies a configuration parameter. This can over-ride existing
values. The value is specified with a fully qualified config path. Does
not support list values []. For example:

::

   -c config.Harder_precip_phase.const.b:1.5 -c config.debug.debug_level:"error"

#remove Removes a configuration parameter. Removals are processed after
any –config parameters are parsed, so -r will override -c. Does not
support remove of individual list [] items. For example:

::

     -c nproc:2 -r nproc

will result in nproc being removed despite being added via ``-c``.

#remove-module Removes a module. Removals are processed after any
–config parameters are parsed, so –remove will override -c. 

::

   -d Marsh_shading_iswr

#add-module Adds a module to the list.

::

   -m snobal -m Marsh_shading_iswr
