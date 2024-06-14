Checkpointing
==================

CHM can save the current state of all modules to disk and then later resume from this point. The provides a form of
checkpointing. How to configure this for a simulation is detailed in the :ref:'configuration section <target to checkpoint>'.

The states are stored as json metadata and netcdf files located in ``output_dir/checkpoint/``. Within this folder are json meta file named in the form
``checkpoint_YYYYMMddTHHMMSS.npRANKS.json``. This file provides the list of the netcdf files that need to be loaded for each MPI rank,
along with sanity checks for how many MPI ranks were used in the checkpoint run and what time to restart from. If a checkpoint file was saved
with ``n`` ranks, then it must be loaded with ``n`` ranks. The netcdf files are saved in a sub-directory ``YYYYMMddTHHMMSS``.

.. code:: json

    {
        "ranks": "8",
        "restart_time_sec": "1577462400",
        "startdate": "20191227T160000",
        "files": [
            "20191227T160000\/chkp20191227T160000_0.nc",
            "20191227T160000\/chkp20191227T160000_1.nc",
            "20191227T160000\/chkp20191227T160000_2.nc",
            "20191227T160000\/chkp20191227T160000_3.nc",
            "20191227T160000\/chkp20191227T160000_4.nc",
            "20191227T160000\/chkp20191227T160000_5.nc",
            "20191227T160000\/chkp20191227T160000_6.nc",
            "20191227T160000\/chkp20191227T160000_7.nc"
        ]
    }



There is one netcdf file per rank, and each row in the netcdf file corresponds to a single triangle ``global_id``, indexed with ``global_id``. A simple
pointmode checkpoint file example is given below:

 .. code:: json

     netcdf chkp20200125T140000_0 {
     dimensions:
        tri_id = 1 ;
     variables:
        double Harder_precip_phase\:hours_since_snowfall(tri_id) ;
        double Harder_precip_phase\:acc_rain(tri_id) ;
        double Harder_precip_phase\:acc_snow(tri_id) ;
        double fsm\:snw(tri_id) ;
        double fsm\:snd(tri_id) ;
        double fsm\:sum_snowpack_subl(tri_id) ;
        double fsm\:albs(tri_id) ;
        double fsm\:Tsrf(tri_id) ;
        double fsm\:Dsnw\[0\](tri_id) ;
        double fsm\:Dsnw\[1\](tri_id) ;
        double fsm\:Dsnw\[2\](tri_id) ;
        double fsm\:Nsnow(tri_id) ;
        double fsm\:Qcan\[0\](tri_id) ;
        double fsm\:Qcan\[1\](tri_id) ;
        double fsm\:Sice\[0\](tri_id) ;
        double fsm\:Sice\[1\](tri_id) ;
        double fsm\:Sice\[2\](tri_id) ;
        double fsm\:Sliq\[0\](tri_id) ;
        double fsm\:Sliq\[1\](tri_id) ;
        double fsm\:Sliq\[2\](tri_id) ;
        double fsm\:Sveg\[0\](tri_id) ;
        double fsm\:Sveg\[1\](tri_id) ;
        double fsm\:Tcan\[0\](tri_id) ;
        double fsm\:Tcan\[1\](tri_id) ;
        double fsm\:Tsnow\[0\](tri_id) ;
        double fsm\:Tsnow\[1\](tri_id) ;
        double fsm\:Tsnow\[2\](tri_id) ;
        double fsm\:Tsoil\[0\](tri_id) ;
        double fsm\:Tsoil\[1\](tri_id) ;
        double fsm\:Tsoil\[2\](tri_id) ;
        double fsm\:Tsoil\[3\](tri_id) ;
        double fsm\:Tveg\[0\](tri_id) ;
        double fsm\:Tveg\[1\](tri_id) ;
        double fsm\:Vsmc\[0\](tri_id) ;
        double fsm\:Vsmc\[1\](tri_id) ;
        double global_id(tri_id) ;

     // global attributes:
            :restart_time = "2020-Jan-25 14:00:00" ;
            :restart_time_sec = 1579960800ULL ;
     }

The user does not need to manually set the restart time, CHM will override any config or user specified start times with
that specified in the checkpoint json file.

The savestate occurs at the end of the timestep, so the resume time will be one timestep into the future, i.e., the
next timestep.

The make use of checkpointing, a module must implement the ``checkpoint(mesh& domain,  netcdf& chkpt)`` and
``load_checkpoint(mesh& domain, netcdf& chkpt)`` methods.