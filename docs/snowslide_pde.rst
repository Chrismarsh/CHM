SnowSlide on the triangular mesh
================================

The ``snowslide_pde`` module releases snow above a slope-dependent holding capacity and
routes the mobile mass over CHM's triangular faces. A face can receive snow, retain
part of it, and pass the rest farther downslope during the same CHM timestep.

Configuration
-------------

Replace ``snow_slide`` with ``snowslide_pde`` in the modules list. Use the same
coupling order as the previous run when comparing modules. For example, these
entries in a Snobal configuration schedule redistribution before Snobal. The event
reads the previous snow-column state, and Snobal consumes its increment once.

.. code-block:: json

   {
       "remove_depency": {"snowslide_pde": "snobal"},
       "config": {
           "snowslide_pde": {
               "release_profile": "bernhardt_schulz_2010",
               "mobility_length_m": 150
           }
       }
   }

Use the registered name ``FSM`` in place of ``snobal`` for FSM. On restart the
module restores both the pending increments and the cumulative volume changes.

With Crocus and PBSM3D, put ``"snowslide_pde": "Crocus"`` in the existing
``remove_depency`` section and set ``config.snowslide_pde.reserve_drift_erosion``
to ``true``. PBSM3D runs first. The PDE deducts wind erosion from each donor's
previous SWE and reduces its normal and vertical depths in proportion before
evaluating the release threshold. Crocus then applies the wind and avalanche
increments once. Wind-deposited snow becomes eligible for release on the next
CHM timestep. The reservation defaults to ``false`` and requires a provider
of ``drift_mass`` when enabled.

Crocus includes unlayered snowfall in ``swe`` and publishes that amount as
``crocus_pending_snow`` (mm SWE). A face with zero snow depth and only pending
snow can receive avalanche deposits. Its trace snow stays in the column.

Holding depth
-------------

* ``bernhardt_schulz_2010`` (default) releases only above a terrain slope of
  25 degrees. The topographic holding depth is linearly interpolated between
  measurements of the plotted curve in Bernhardt and Schulz (2010), Fig. 2b.
  The curve is interpreted as vertical depth. At 25, 30, 35, 40, 45, 50, 55,
  60, 65, 70 and 75 degrees the digitized depths are 5.0, 3.4, 2.4, 1.8, 1.4,
  1.1, 0.9, 0.8, 0.6, 0.6 and 0.5 m. One pixel on the depth axis spans about
  0.17 m. At and above 75 degrees the topographic capacity is 50 mm SWE.
  A minimum of 50 mm SWE applies throughout. Use the larger of the vegetation
  and topographic capacities. The host pack's SWE and vertical depth determine
  the release conversion.
* The paper profile takes vegetation holding depth in metres vertically. Supply
  a face parameter ``snow_holding_depth`` or a ``snow_holding_depth`` value in
  each CHM land-cover class (including open terrain) to specify the site's
  calibration. The fallback is CHM canopy height where supplied, otherwise zero.
* ``chm_vionnet_2021`` uses the same slope-normal holding depth calculation
  as CHM's ``snow_slide``, including the 10-degree regression floor, 20 m cap,
  canopy protection and legacy options ``avalache_mult`` (3178.4) and
  ``avalache_pow`` (-1.998). Use this setting to compare redistribution schemes
  while preserving the historical CHM release threshold.

Transport and deposition
------------------------

Let :math:`m` be mobile SWE in kg per square metre of plan area and :math:`s`
be dimensionless event time. The mobile layer satisfies

.. math::

   \partial_s m + \nabla\cdot(-\ell m\nabla\eta)=-d,\qquad
   \eta = z + H_{\rm stationary,vertical} + m/\rho_{\rm mobile}.

CHM supplies the original stationary depths. In each event substep a backward
Euler finite-volume solve transports the mobile inventory with an upwind flux.
The unknown at each Newton iteration is the end-of-substep mobile mass on
every face. Each trial updates the snow surface and the limited, nonorthogonal
edge slopes. The stationary surface stays fixed during this solve. The
approximate Jacobian uses PBSM's one-layer nearest-neighbour graph, Belos GMRES
and Ifpack2 ILUT. A line search requires a decreasing full residual and
nonnegative trial masses.

After transport a fraction :math:`1-\exp(-\Delta s)` of the mobile layer can
settle. A release-eligible face with a downslope route retains up to its remaining
holding capacity. Low-angle faces under the paper profile and faces with no
outward route can retain the whole settling fraction. Settled snow updates the
stationary surface for the next substep. A closed depression can fill and spill.

The new routing parameter ``mobility_length_m`` has an experimental default
of 150 m. On a uniform grade :math:`g` with negligible snow-thickness gradients
and unrestricted settling, :math:`\ell g` is the e-folding travel distance.
The SnowSlide papers constrain release. Field runout or deposit footprints are
needed to estimate this routing parameter. ``mobile_density`` (300 kg/m³)
converts SWE to vertical mobile thickness. ``deposit_density`` (500 kg/m³)
converts gross deposited SWE to slope-normal thickness, matching FSM's default
avalanche deposition density.

Event timing and boundaries
---------------------------

One CHM call runs a released event to completion. Event time does not advance
the forcing clock. The initial and maximum trial interval is ``event_step``
(0.25). With ``adaptive`` enabled (default), compare one full step against two
half steps from the same initial state. For mobile SWE, the step's deposition
and the step's export per unit area, take the maximum across faces and fields of

.. math::

   E = \frac{|u_{\rm fine}-u_{\rm coarse}|}
            {a+r\max(|u_{\rm fine}|,|u_{\rm coarse}|)}.

Defaults are ``event_abs_tol`` = 0.1 mm SWE and ``event_rel_tol`` = 0.01.
Accept the two-half-step result when :math:`E\leq 1`. Halve and retry on error
or failed nonlinear convergence. Failure at ``min_event_step`` (default
:math:`10^{-5}`) aborts the event before committing any coupling increments.
Setting ``adaptive`` to false disables the step comparison. Solver failure still
triggers subdivision. When global mobile mass falls below a fraction
:math:`10^{-8}` of release, settle the residual in place and check the mass budget.

Domain edges and interfaces with water faces are closed. A positive face
parameter ``snowslide_outlet`` opens its exterior edges to outward transport.
The boundary flux uses the reconstructed interior slope. Geographic meshes
need plan area in square metres in the ``area`` parameter.

Outputs
-------

* ``avalanche_release_swe`` and ``avalanche_deposit_swe`` are gross event
  release and deposition in mm SWE per face.
* ``avalanche_export_volume`` is water-equivalent volume exported at exterior
  edges, in m³ per face per CHM step.
* ``delta_avalanche_mass`` is net deposited minus released water-equivalent
  volume in m³ per face. ``delta_avalanche_snowdepth`` is the corresponding
  change in slope-normal snow depth multiplied by plan area, in m³ per face.
  Their ``_sum`` fields accumulate over CHM steps.
* ``maxDepth`` is the local holding threshold expressed as slope-normal depth
  in metres. Release eligibility is evaluated separately.

Face variables prefixed with ``snowslide_pde_`` are working storage for MPI
halo exchange during the solve.

References
----------

Bernhardt and Schulz (2010), *SnowSlide: A simple routine for calculating
gravitational snow transport*, https://doi.org/10.1029/2010GL043086.

Vionnet et al. (2021), *Multi-scale snowdrift-permitting modelling of mountain
snowpack*, https://doi.org/10.5194/tc-15-743-2021.
