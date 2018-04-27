![](https://github.com/Chrismarsh/CHM/blob/master/wiki/mesh.png)

# The Canadian Hydrological Model 
The Canadian Hydrological Model (CHM) is a novel modular unstructured mesh based approach for hydrological modelling. It can move between spatial scale, temporal scale, and spatial extents. It is designed for developing and testing process representations for hydrological models. 

<!-- MarkdownTOC autolink="true" -->

- [Usage](#usage)
- [Motivation](#motivation)
- [Design goals](#design-goals)
- [Publications](#publications)
- [Features](#features)
    - [Spatial Scales](#spatial-scales)
    - [Visualization](#visualization)
    - [netCDF support](#netcdf-support)
    - [Process representations](#process-representations)
    - [Unstructured mesh](#unstructured-mesh)
    - [Parallel computing](#parallel-computing)
    - [Uncertainty analysis](#uncertainty-analysis)
- [Demonstration](#demonstration)
    - [SnowCast](#snowcast)
    - [Large extent](#large-extent)
    - [Point scale](#point-scale)
    - [Blowing snow](#blowing-snow)

<!-- /MarkdownTOC -->

# Usage
Details on how to use CHM, as well as more implimentation details, can be found in the [wiki](https://github.com/Chrismarsh/CHM/wiki)


# Motivation
Modelling of hydrological processes at any scale is hampered by large uncertainties in parameters and forcing data, incomplete process representations (the scientific conceptualization of a phenomena codified numerically), and arbitrary process representation selections and linkages (collectively ‘model structure’). There is also consistent difficulty or an inability to easily test and estimate the uncertainty due to variations in model structure, parameter values, number of parameters, forcing data requirements, and spatial discretization requirements (collectively ‘model complexity’). 

In this work, a new distributed model framework is presented that can examine a variety of process representations, process linkages and levels of model complexity. Algorithms can be easily interchanged, removed, and decoupled while preserving the underlying model framework. Thus, uncertainty propagation and subsequent feedbacks within the model structure can be quantified. Unstructured meshes represent the spatial heterogeneity of surface and sub-surface features in a computationally efficient manner and also decreases number of parameters and initial conditions. The parallel architecture allows for efficient uncertainty testing of parameter ranges. By utilizing unstructured meshes, fewer than 5% of the computational elements of high-resolution structured (raster) grids are usually necessary.  This preserves surface and sub-surface heterogeneity but results in fewer parameters and initial conditions.

# Design goals
* Multi-scale, multi-physics, variable complexity and domain model
* Assessment of model structural, parameter, and data uncertainty
* Easily test multiple hypotheses, avoid rigid model structures
* Incorporate existing code
* Contribute to decision support systems

# Publications
The manuscript describing CHM is still in preperation, and is anticipated to be submited summer 2018.

# Features
## Spatial Scales
CHM is applicable to multiple scales from the basin scale, to the provincial/state scale and beyond. It may also be applied at a single point-scale.
![](https://github.com/Chrismarsh/CHM/blob/master/wiki/scale.png)

## Visualization
Output is in the vtu file format, allowing for visualization, analysis, and timeseries animation in ![ParaView](https://www.paraview.org/). Date-time support has been added to ParaView via an filter ![vtk-paraview-datetimefilter](https://github.com/Chrismarsh/vtk-paraview-datetimefilter).

![](https://github.com/Chrismarsh/CHM/blob/master/wiki/paraview.png)

## netCDF support
Input meterology may be either in a standard ASCII file, or as a netCDF file allowing for ease of use when using climate model outputs. 

The below figure shows virtual stations that correspond to the center of the 2.5 km GEM numerical weather prediction output in netCDF format.

![](https://github.com/Chrismarsh/CHM/blob/master/wiki/netcdf.png)
## Process representations 

 Process represetenation will be extented to include the entirety of the hydrological cycle. However, current representation includes mostly surface and cold regions processes

| Process | Module |
|---------|--------|
|Canopy  |Open/forest (exp/log) (Pomeroy et al., 1998; Ellis et al., 2010)|
|Snowpack  |   2-layer Snobal (Marks et al, 1999); Multi-layer Snowpack (Lehning et al., 1999); Various albedo e.g., CLASS (Verseghy 1991) |
|Soil  |   Frozen soil infiltration (Gray et al., 2001) |
|Mass redistribution | PBSM3D (Marsh et al, 2018 in review); Snowslide (Bernhardt 2010) |

Input meterology is spatially interpolated and down-scaled from the input station or virtual-station (e.g., from numerical weather prediction) to produce a spatially distributed driving dataset. There are a number of ways to downscale these meterology.

|Variable |     Type|
| ------- | ------  |
|Air temperature | Linear lapse rates (measured, seasonal, constant, neutral stability) (Kunkel, 1989, Dodson et al., 1997)|
|Relative humidity |   Linear lapse rates (measured, seasonal, constant) (Kunkel, 1989)|
| Horizontal wind | Topographic curvature (Liston, et al., 2006); Mason-Sykes (Mason and Sykes, 1979); uniform wind |
|Precipitation |   Elevation based lapse (Thornton, 1997) |
| Precipitation Phase | Linear; Psychometric (Harder and Pomeroy, 2013); Threshold |
| Solar radiation | Terrain shadows (Marsh et al., 2011, Dozier and Frew, 1990); Clear sky transmittance (Burridge, 1975); Transmittance from observations; Cloud fraction estimates (Walcek, 1994); Direct/diffuse splitting (Iqbal, 19xx) |
| Longwave |    T, RH based (Sicart et al., 2006); Constant (Marty et al., 2002) |

## Unstructured mesh
CHM uses an unstructured triangular mesh to representent the terrain. This mesh is generated by ![Mesher](https://github.com/Chrismarsh/mesher), a novel multi-objective unstructured mesh generation software that allows mesh generation to be generated from an arbitrary number of hydrologically important features while maintaining a variable spatial resolution. Triangle quality is guaranteed as well as a smooth graduation from small to large triangles. Including these additional features resulted in a better representation of spatial heterogeneity versus classic topography-only mesh generation while significantly reducing the total number of computational elements.

![](https://raw.githubusercontent.com/Chrismarsh/mesher/master/images/mesh.png)


## Parallel computing

In CHM, parallelism is currently implemented via the shared memory API OpenMP. As described above, modules may either be point-scale models that are applied to each triangle independently or require knowledge of the surrounding triangles. Mixing these two types of parallelism complicates the implementation of parallel code. To provide as much seamless parallelism as possible to the modules, each module declares the type of algorithm it is: data parallel or domain parallel. Data parallel modules are point-scale models that are applied to every triangle. Domain parallel modules are modules that require knowledge of surrounding mesh points. Thus, after the topological sort is performed to determine module execution order, the modules are scheduled together into groups that share a parallelism type

## Uncertainty analysis
A key feature of CHM is the ability to, on the command line, change any value specified by a configuration parameter. CHM provides a seamless mechanism to easily allow modules to obtain parameter data from configuration files.

```python
import subprocess
import shutil


prj_path = "CHM.config"

cf1 = "-c output.VistaView.file:vv_dodson.txt"
cf2 = "-c output.UpperClearing.file:uc_dodson.txt"
cf3 = "-c output.FiserraRidge.file:fr_dodson.txt"
cf4 = "--add-module Dodson_NSA_ta"
subprocess.check_call(['./CHM %s %s %s %s %s' % (prj_path, cf1, cf2, cf3,cf4)], shell=True)
```


# Demonstration
## SnowCast
![SnowCast](https://www.snowcast.ca) is an experimental, daily data product that uses the Global Environmental Multiscale (GEM) model forecasts from Environment and Climate Change Canada (ECCC) to drive the Canadian Hydrological Model (CHM). Estimates of snowpack are provided over the a Bow River Basin, centered over Banff, Canada.

SnowCast is developed as part of ![Global Water Futures](https://gwf.usask.ca/) and the ![Centre for Hydrology](https://www.usask.ca/hydrology/), University of Saskatchewan. 

## Large extent
Hourly solar radiation modelling for the territory of Yukon, Canada.
![](https://github.com/Chrismarsh/CHM/blob/master/wiki/yk_solar.gif)

## Point scale
Comparison of CHM driving Snobal and Snowpack at the Upper Clearing site at Marmot Creek Research Basin in Alberta, Canada
![](https://github.com/Chrismarsh/CHM/blob/master/wiki/CHM_crhm_v_chm_v_obs_swe.png)

## Blowing snow 
Comparison of blowing snow (left) versus no blowing snow (right) for a small sub-basin of Wolf Creek Reserach Basin, located in the Yukon, Canada.
![](https://github.com/Chrismarsh/CHM/blob/master/wiki/output_small.gif)







