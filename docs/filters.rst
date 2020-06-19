Filters
========
Overview
-----------
Filters are a mechanism whereby the input forcing data can be modified
in some way prior to the model run. For example, this could be use to
apply a gauge undercatch to precip. Filters modify the data of a virtual station
*in situ*.

.. warning::

   Filters run in the order defined in the configuration file.

Please see :ref:`config:forcing` for how to attach a filter to an input forcing file.

Available filters
------------------
.. groups::
   :group: filters


