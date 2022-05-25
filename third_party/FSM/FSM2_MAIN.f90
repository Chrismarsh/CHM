!----------------------------------------------------------------------!
! Flexible Snow Model (FSM version 2.1.0)                              !
! Example driver for CHM implementation                                !
!                                                                      !
! Richard Essery                                                       !
! School of GeoSciences                                                !
! University of Edinburgh                                              !
!----------------------------------------------------------------------!
program FSM2

use CONSTANTS, only: &
  eps,               &! Ratio of molecular weights of water and dry air
  e0                  ! Saturation vapour pressure at Tm (Pa)

use LAYERS, only: &
  Dzsnow,            &! Minimum snow layer thicknesses (m)
  Dzsoil,            &! Soil layer thicknesses (m)
  fvg1,              &! Fraction of vegetation in upper canopy layer
  Ncnpy,             &! Number of canopy layers
  Nsmax,             &! Maximum number of snow layers
  Nsoil,             &! Number of soil layers
  zsub                ! Subcanopy wind speed diagnostic height (m)

use PARAMETERS, only: &
  rgr0                ! Fresh snow grain radius (m)

use SOILPROPS, only: &
  b,                 &! Clapp-Hornberger exponent
  hcap_soil,         &! Volumetric heat capacity of dry soil (J/K/m^3)
  hcon_soil,         &! Thermal conductivity of dry soil (W/m/K)
  sathh,             &! Saturated soil water pressure (m)
  Vcrit,             &! Volumetric soil moisture at critical point
  Vsat                ! Volumetric soil moisture at saturation

implicit none

! Meteorological driving variables
integer :: &
  year,              &! Year
  month,             &! Month of year
  day                 ! Day of month
real :: &
  dt,                &! Timestep (s)
  elev,              &! Solar elevation (radians)
  hour,              &! Hour of day
  LW,                &! Incoming longwave radiation (W/m^2)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  Qs,                &! Specific humidity at saturation (kg/kg)
  Rf,                &! Rainfall rate (kg/m^2/s)
  RH,                &! Relative humidity (%)
  Sdif,              &! Diffuse shortwave radiation (W/m^2)
  Sdir,              &! Direct-beam shortwave radiation (W/m^2)
  SW,                &! Global shortwave radiation (W/m^2)
  Sf,                &! Snowfall rate (kg/m^2/s)
  Ta,                &! Air temperature (K)
  Tc,                &! Air temperature (C)
  trans,             &! Wind-blown snow transport rate (kg/m^2/s)
  Ua,                &! Wind speed (m/s)
  zT,                &! Temperature and humidity measurement height (m)
  zU                  ! Wind speed measurement height (m)

! Model state variables  
integer :: &
  Nsnow               ! Number of snow layers
real :: &
  albs,              &! Snow albedo
  Tsrf,              &! Snow/ground surface temperature (K)
  Dsnw(3),           &! Snow layer thicknesses (m)
  Qcan(2),           &! Canopy air space humidities
  Rgrn(3),           &! Snow layer grain radii (m)
  Sice(3),           &! Ice content of snow layers (kg/m^2)
  Sliq(3),           &! Liquid content of snow layers (kg/m^2)
  Sveg(2),           &! Snow mass on vegetation layers (kg/m^2)
  Tcan(2),           &! Canopy air space temperatures (K)
  Tsnow(3),          &! Snow layer temperatures (K)
  Tsoil(4),          &! Soil layer temperatures (K)
  Tveg(2),           &! Vegetation layer temperatures (K)
  Vsmc(4)             ! Volumetric moisture content of soil layers

! Diagnostics
real :: &
  H,                 &! Sensible heat flux to the atmosphere (W/m^2)
  LE,                &! Latent heat flux to the atmosphere (W/m^2)
  LWout,             &! Outgoing LW radiation (W/m^2)
  LWsub,             &! Subcanopy downward LW radiation (W/m^2)
  Melt,              &! Surface melt rate (kg/m^2/s)
  Roff,              &! Runoff from snow (kg/m^2/s)
  snd,               &! Snow depth (m)
  snw,               &! Total snow mass on ground (kg/m^2) 
  subl,              &! Sublimation rate (kg/m^2/s)
  svg,               &! Total snow mass on vegetation (kg/m^2)
  SWout,             &! Outgoing SW radiation (W/m^2)
  SWsub,             &! Subcanopy downward SW radiation (W/m^2)
  Usub,              &! Subcanopy wind speed (m/s)
  Wflx(3)             ! Water flux into snow layer (kg/m^2/s)

! Vegetation characteristics
real :: &
  alb0,              &! Snow-free ground albedo
  vegh,              &! Vegetation height (m)
  VAI                 ! Vegetation area index
integer :: &
  Ntyp                ! Vegetation type

! Grid dimensions
Ncnpy = 2
Nsmax = 3
Nsoil = 4

! Canopy, snow and soil layers
fvg1 = 0.5
zsub = 1.5
allocate(Dzsnow(Nsmax))
allocate(Dzsoil(Nsoil))
Dzsnow = (/0.1, 0.2, 0.4/)
Dzsoil = (/0.1, 0.2, 0.4, 0.8/)

! Meteorological driving data
open(8,file='met_CdP_0506.txt')
dt = 3600
zT = 2
zU = 10

! Vegetation characteristics
alb0 = 0.2
vegh = 0
VAI  = 0
Ntyp = 1

! Soil properties
b = 7.63
hcap_soil = 2.3e6
hcon_soil = 0.11
sathh = 0.41
Vcrit = 0.26
Vsat = 0.27

! Initialize state variables with no snow
albs = 0.8
Dsnw = 0
Nsnow = 0
Qcan = 0
Rgrn = rgr0
Sice = 0
Sliq = 0
Sveg = 0
Tcan = 285
Tsnow = 273
Tsoil = 285
Tsrf = 285
Tveg = 285
Vsmc = 0.5*Vsat

open(11,file='out.txt')

do

! Meteorological driving data
  read(8,*,end=1) year,month,day,hour,SW,LW,Sf,Rf,Ta,RH,Ua,Ps
! Convert relative to specific humidity
  Tc = Ta - 273.15
  Qs = eps*(e0/Ps)*exp(17.5043*Tc/(241.3 + Tc))
  Qa = (RH/100)*Qs
! Lower limit required on windspeed
  Ua = max(Ua, 0.1)
! All SW radiation assumed to be diffuse
  elev = 0
  Sdif = SW
  Sdir = 0

! Replace with snow transport rate
  trans = 0

  call FSM2_TIMESTEP(                                                  &
                     ! Driving variables                               &
                     dt,elev,zT,zU,                                    &
                     LW,Ps,Qa,Rf,Sdif,Sdir,Sf,Ta,trans,Ua,             &
                     ! Vegetation characteristics                      &
                     Ntyp,alb0,vegh,VAI,                               &
                     ! State variables                                 &
                     albs,Tsrf,Dsnw,Nsnow,Qcan,Rgrn,Sice,              &
                     Sliq,Sveg,Tcan,Tsnow,Tsoil,Tveg,Vsmc,             &
                     ! Diagnostics                                     &
                     H,LE,LWout,LWsub,Melt,Roff,snd,snw,subl,svg,      &
                     SWout,SWsub,Usub,Wflx                             )

  write(11,100) year,month,day,hour,snd,snw

end do
1 continue

100 format(3(i4),3(f12.3))

end program FSM2
