
!-----------------------------------------------------------------------
! Physical constants
!-----------------------------------------------------------------------
module CONSTANTS
real, protected :: &
  cp = 1005,         &! Specific heat capacity of air (J/K/kg)
  eps = 0.622,       &! Ratio of molecular weights of water and dry air
  e0 = 611.213,      &! Saturation vapour pressure at Tm (Pa)
  grav = 9.81,       &! Acceleration due to gravity (m/s^2)
  hcap_ice = 2100,   &! Specific heat capacity of ice (J/K/kg)
  hcap_wat = 4180,   &! Specific heat capacity of water (J/K/kg)
  hcon_air = 0.025,  &! Thermal conductivity of air (W/m/K)
  hcon_ice = 2.24,   &! Thermal conducivity of ice (W/m/K)
  hcon_wat = 0.56,   &! Thermal conductivity of water (W/m/K)
  I0 = 1367,         &! Solar constant (W/m^2)
  Lf = 0.334e6,      &! Latent heat of fusion (J/kg)
  Lv = 2.501e6,      &! Latent heat of vapourisation (J/kg)
  Ls = 2835000 ,      &! Lf + Lv Latent heat of sublimation (J/kg)
  pi = 3.14159,      &! pi
  Rair = 287,        &! Gas constant for air (J/K/kg)
  Rwat = 462,        &! Gas constant for water vapour (J/K/kg)
  rho_ice = 917,     &! Density of ice (kg/m^3)
  rho_wat = 1000,    &! Density of water (kg/m^3)
  sb = 5.67e-8,      &! Stefan-Boltzmann constant (W/m^2/K^4)
  Tm = 273.15,       &! Melting point (K)
  vkman = 0.4         ! von Karman constant
end module CONSTANTS

!-----------------------------------------------------------------------
! Canopy, snow and soil layers
!-----------------------------------------------------------------------
module LAYERS
integer :: &
  Ncnpy,             &! Number of canopy layers
  Nsmax,             &! Maximum number of snow layers
  Nsoil               ! Number of soil layers
real, allocatable :: &
  Dzsnow(:),         &! Minimum snow layer thicknesses (m)
  Dzsoil(:)           ! Soil layer thicknesses (m)
end module LAYERS

!-----------------------------------------------------------------------
! Parameters
!-----------------------------------------------------------------------
module PARAMETERS

integer :: &
  Nitr = 4            ! Iterations in energy balance calculation

! Vegetation parameters
real, protected :: &
  avg0 = 0.1,        &! Snow-free vegetation albedo
  avgs = 0.4,        &! Snow-covered vegetation albedo
  cvai = 2500,       &! Vegetation heat capacity per unit VAI (J/K/m^2)
  fvg1 = 0.5,        &! Fraction of vegetation in upper layer
  gsnf = 0.01,       &! Snow-free vegetation moisture conductance (m/s)
  kext = 0.5,        &! Vegetation light extinction coefficient
  rveg = 20,         &! Leaf boundary resistance (s/m)^(1/2)
  svai = 4.4,        &! Intercepted snow capacity per unit VAI (kg/m^2)
  tunc = 240*3600,   &! Unloading time scale for cold snow (s)
  tunm = 2.4*3600,   &! Unloading time scale for melting snow (s)
  vint = 0.7,        &! Vegetation snow interception efficiency
  wcan = 2.5,        &! Canopy wind decay coefficient
  zsub = 1.5          ! Subcanopy wind speed diagnostic height (m)

! Snow parameters
real, protected :: &
  asmn = 0.5,        &! Minimum albedo for melting snow
  asmx = 0.8,        &! Maximum albedo for fresh snow
  hfsn = 0.1,        &! Snowcover fraction depth scale (m)
  kfix = 0.24,       &! Fixed thermal conductivity of snow (W/m/K)
  rcld = 300,        &! Maximum density for cold snow (kg/m^3)
  rfix = 300,        &! Fixed snow density (kg/m^3)
  rgr0 = 5e-5,       &! Fresh snow grain radius (m)
  rhof = 100,        &! Fresh snow density (kg/m^3)
  rmlt = 500,        &! Maximum density for melting snow (kg/m^3)
  Salb = 10,         &! Snowfall to refresh albedo (kg/m^2)
  Talb = -2,         &! Snow albedo decay temperature threshold (C)
  tcld = 3.6e6,      &! Cold snow albedo decay time scale (s)
  tmlt = 3.6e5,      &! Melting snow albedo decay time scale (s)
  trho = 200*3600,   &! Snow compaction timescale (s)
  Wirr = 0.03,       &! Irreducible liquid water content of snow
  z0sn = 0.01         ! Snow roughness length (m)

! Ground surface and soil parameters
real, protected :: &
  b = 7.63,          &! Clapp-Hornberger exponent
  gsat = 0.01,       &! Surface conductance for saturated soil (m/s)
  hcap_soil = 2.3e6, &! Volumetric heat capacity of dry soil (J/K/m^3)
  hcon_soil = 0.11,  &! Thermal conductivity of dry soil (W/m/K)
  sathh = 0.41,      &! Saturated soil water pressure (m)
  Vcrit = 0.26,      &! Volumetric soil moisture at critical point
  Vsat = 0.27,       &! Volumetric soil moisture at saturation
  z0sf = 0.1          ! Snow-free surface roughness length (m)

end module PARAMETERS

subroutine allocate()
  use LAYERS, only: &
          Nsmax,             &! Maximum number of snow layers
          Nsoil,             &  ! Number of soil layers
          Dzsnow,&
          Dzsoil

  ! Snow and soil layersdo
  allocate(Dzsnow(Nsmax))
  allocate(Dzsoil(Nsoil))
  Dzsnow = (/0.1, 0.2, 0.4/)
  Dzsoil = (/0.1, 0.2, 0.4, 0.8/)

end subroutine allocate

!-----------------------------------------------------------------------
! Call FSM2 physics subroutines for one timestep at one point
!-----------------------------------------------------------------------
subroutine FSM2_TIMESTEP(dt,elev,zT,zU,LW,Ps,Qa,Rf,Sdif,Sdir,Sf,Ta,Ua, &
                         alb0,hveg,VAI,                                &
                         albs,Tsrf,Dsnw,Nsnow,Qcan,Rgrn,Sice,Sliq,     &
                         Sveg,Tcan,Tsnow,Tsoil,Tveg,Vsmc,              &
                         H,LE,LWout,LWsub,Melt,Roff,snd,Svg,SWE,       &
                         SWout,SWsub,Usub)
use iso_c_binding

use LAYERS, only: &
  Ncnpy,             &! Number of canopy layers
  Nsmax,             &! Maximum number of snow layers
  Nsoil,             &  ! Number of soil layers
  Dzsnow,&
  Dzsoil


use PARAMETERS, only: &
  Nitr                ! Iterations in energy balance calculation

implicit none

! Site characteristics
real, intent(in) :: &!
  dt,                &! Timestep (s)
  zT,                &! Temperature and humidity measurement height (m)
  zU                  ! Wind speed measurement height (m)

! Meteorological variables
real, intent(in) :: &
  elev,              &! Solar elevation (radians)
  LW,                &! Incoming longwave radiation (W/m2)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  Rf,                &! Rainfall rate (kg/m2/s)
  Sdif,              &! Diffuse shortwave radiation (W/m^2)
  Sdir,              &! Direct-beam shortwave radiation (W/m^2)
  Ta,                &! Air temperature (K)
  Ua                  ! Wind speed (m/s)
real, intent(inout) :: &
  Sf                  ! Snowfall rate (kg/m2/s)

! Vegetation characteristics
real, intent(in) :: &
  alb0,              &! Snow-free ground albedo
  hveg,              &! Vegetation height (m)
  VAI                 ! Vegetation area index

! State variables
integer, intent(inout) :: &
  Nsnow               ! Number of snow layers
real, intent(inout) :: &
  albs,              &! Snow albedo
  Tsrf,              &! Snow/ground surface temperature (K)
  Dsnw(Nsmax),       &! Snow layer thicknesses (m)
  Qcan(Ncnpy),       &! Canopy air space humidities
  Rgrn(Nsmax),       &! Snow layer grain radii (m)
  Sice(Nsmax),       &! Ice content of snow layers (kg/m^2)
  Sliq(Nsmax),       &! Liquid content of snow layers (kg/m^2)
  Sveg(Ncnpy),       &! Snow mass on vegetation layers (kg/m^2)
  Tcan(Ncnpy),       &! Canopy air space temperatures (K)
  Tsnow(Nsmax),      &! Snow layer temperatures (K)
  Tsoil(Nsoil),      &! Soil layer temperatures (K)
  Tveg(Ncnpy),       &! Vegetation layer temperatures (K)
  Vsmc(Nsoil)         ! Volumetric moisture content of soil layers

! Diagnostics
real, intent(out) :: &
  H,                 &! Sensible heat flux to the atmosphere (W/m^2)
  LE,                &! Latent heat flux to the atmosphere (W/m^2)
  LWout,             &! Outgoing LW radiation (W/m^2)
  LWsub,             &! Subcanopy downward LW radiation (W/m^2)
  Melt,              &! Surface melt (kg/m^2)
  Roff,              &! Runoff from snow (kg/m^2)
  snd,               &! Snow depth (m)
  Svg,               &! Total snow mass on vegetation (kg/m^2)
  SWE,               &! Total snow mass on ground (kg/m^2) 
  SWout,             &! Outgoing SW radiation (W/m^2)
  SWsub,             &! Subcanopy downward SW radiation (W/m^2)
  Usub                ! Subcanopy wind speed (m/s)

! Vegetation properties
real :: &
  aveg(Ncnpy),       &! Vegetation area indices in layers
  cveg(Ncnpy),       &! Vegetation layer heat capacities (J/K/m^2)
  fvegs(Ncnpy),      &! Vegetation layer snowcover fractions
  Scap(Ncnpy),       &! Vegetation layer snow capacities (kg/m^2)
  tdif(Ncnpy),       &! Veg layer transmissivities for diffuse radiation
  Tveg0(Ncnpy)        ! Veg layer temperatures at start of timestep (K)

! Snow properties
real :: &
  fsnow,             &! Ground snowcover fraction
  ksnow(Nsmax)        ! Thermal conductivities of snow layers (W/m/K)

! Surface properties
real :: &
  Ds1,               &! Surface layer thickness (m)
  gs1,               &! Surface moisture conductance (m/s)
  ks1,               &! Surface layer thermal conductivity (W/m/K)
  Ts1                 ! Surface layer temperature (K)

! Soil properties
real :: &
  csoil(Nsoil),      &! Areal heat capacity of soil layers (J/K/m^2)
  ksoil(Nsoil)        ! Thermal conductivity of soil layers (W/m/K)

! Aerodynamic conductances
real :: &
  ga,                &! Aerodynamic conductance to the atmosphere (m/s)
  gc,                &! Conductance within canopy air space (m/s)
  gs,                &! Surface to canopy air space conductance (m/s)
  gv(Ncnpy)           ! Vegetation to canopy air space conductance (m/s)

! Fluxes
real :: &
  Esrf,              &! Moisture flux from the surface (kg/m^2/s)
  Gsrf,              &! Heat flux into snow/ground surface (W/m^2)
  Gsoil,             &! Heat flux into soil (W/m^2)
  SWsrf,             &! SW absorbed by snow/ground surface (W/m^2)
  unload,            &! Snow mass unloaded from vegetation (kg/m^2)
  Eveg(Ncnpy),       &! Moisture flux from vegetation layers (kg/m^2/s)
  SWveg(Ncnpy)        ! SW absorbed by vegetation layers (W/m^2)

integer :: n          ! Iteration counter



call CANOPY(Sveg,Tveg,VAI,aveg,cveg,fvegs,Scap,Tveg0)

call SWRAD(alb0,aveg,Dsnw,dt,elev,fvegs,Sdif,Sdir,Sf,Tsrf,albs,fsnow,  &
           SWout,SWsrf,SWsub,SWveg,tdif)

call THERMAL(Dsnw,Nsnow,Sice,Sliq,Tsnow,Tsoil,Vsmc,csoil,Ds1,          &
             gs1,ksnow,ksoil,ks1,Ts1)

do n = 1, Nitr

  call SRFEXCH(aveg,fsnow,hveg,Ua,VAI,zT,zU,ga,gc,gs,gv,Usub)
 
  call SRFEBAL(cveg,Ds1,dt,fvegs,ga,gc,gs,gs1,gv,ks1,LW,Ps,Qa,SWsrf,   &
               SWveg,Ta,tdif,Ts1,Tveg0,VAI,Tsrf,Qcan,Sice,Tcan,Tveg,   &
               Esrf,Eveg,Gsrf,H,LE,LWout,LWsub,Melt)
end do

call INTERCEPT(dt,Eveg,Scap,Tveg,Sf,Sveg,Svg,unload)

call SNOW(dt,Esrf,Gsrf,ksnow,ksoil,Melt,Rf,Sf,Ta,Tsrf,unload,Nsnow,    &
          Dsnw,Rgrn,Sice,Sliq,Tsnow,Tsoil,Gsoil,Roff,snd,SWE)

call SOIL(csoil,dt,Gsoil,ksoil,Tsoil)

end subroutine FSM2_TIMESTEP

!-----------------------------------------------------------------------
! Properties of vegetation canopy layers
!-----------------------------------------------------------------------
subroutine CANOPY(Sveg,Tveg,VAI,aveg,cveg,fvegs,Scap,Tveg0)

use LAYERS, only: &
  Ncnpy               ! Number of canopy layers

use PARAMETERS, only: &
  cvai,              &! Vegetation heat capacity per unit VAI (J/K/m^2)
  fvg1,              &! Fraction of vegetation in upper layer
  svai                ! Intercepted snow capacity per unit VAI (kg/m^2)

implicit none

real, intent(in) :: &
  Sveg(Ncnpy),       &! Snow mass on vegetation layers (kg/m^2)
  Tveg(Ncnpy),       &! Vegetation layer temperatures (K)
  VAI                 ! Vegetation area index

real, intent(out) :: &
  aveg(Ncnpy),       &! Vegetation area indices in layers
  cveg(Ncnpy),       &! Vegetation layer heat capacities (J/K/m^2)
  fvegs(Ncnpy),      &! Vegetation layer snowcover fractions
  Scap(Ncnpy),       &! Vegetation layer snow capacities (kg/m^2)
  Tveg0(Ncnpy)        ! Vegetation temperatures at start of timestep (K)

aveg(1) = fvg1*VAI
aveg(2) = (1 - fvg1)*VAI

cveg = cvai*aveg
Scap = svai*aveg
fvegs = 0
if (VAI > epsilon(VAI)) fvegs = Sveg / Scap
Tveg0 = Tveg

end subroutine CANOPY

!-----------------------------------------------------------------------
! Mass balance of snow intercepted by vegetation
!-----------------------------------------------------------------------
subroutine INTERCEPT(dt,Eveg,Scap,Tveg,Sf,Sveg,Svg,unload)

use CONSTANTS, only: &
  Tm                  ! Melting point (K)

use LAYERS, only: &
  Ncnpy               ! Number of canopy layers

use PARAMETERS, only: &
  tunc,              &! Unloading time scale for cold snow (s)
  tunm,              &! Unloading time scale for melting snow (s)
  vint                ! Vegetation snow interception efficiency

implicit none

real, intent(in) :: &
  dt,                &! Timestep (s)
  Eveg(Ncnpy),       &! Moisture flux from vegetation layers (kg/m^2/s)
  Scap(Ncnpy),       &! Vegetation layer snow capacities (kg/m^2)
  Tveg(Ncnpy)         ! Vegetation layer temperatures (K)

real, intent(inout) :: &
  Sf,                &! Snowfall rate (kg/m2/s)
  Sveg(Ncnpy)         ! Snow mass on vegetation layers (kg/m^2)

real, intent(out) :: &
  Svg,               &! Total snow mass on vegetation (kg/m^2)
  unload              ! Snow mass unloaded from vegetation (kg/m^2)

integer :: k          ! Vegetation layer counter

real :: &
  intcpt              ! Vegetation layer snow interception (kg/m^2)

Svg = 0
unload = 0
if (Scap(1) > 0) then
  do k = 1, Ncnpy

    ! Interception of falling snow
    intcpt = vint*(Scap(k) - Sveg(k))*(1 - exp(-Sf*dt/Scap(k)))
    Sveg(k) = Sveg(k) + intcpt
    Sf = Sf - intcpt/dt

    ! Sublimation of canopy snow
    if (Sveg(k) > 0 .or. Tveg(k) < Tm) Sveg(k) = Sveg(k) - Eveg(k)*dt
    Sveg(k) = max(Sveg(k), 0.)

    ! Unloading of cold or melting canopy snow
    if (Tveg(k) < Tm) then
      unload = unload + Sveg(k)*dt/tunc
      Sveg(k) = (1 - dt/tunc)*Sveg(k)
    else
      unload = unload + Sveg(k)*dt/tunm
      Sveg(k) = (1 - dt/tunm)*Sveg(k)
    end if

  end do
end if
Svg = sum(Sveg)

end subroutine INTERCEPT

!-----------------------------------------------------------------------
! Solve matrix equation Ax = b for x by LU decomposition
!-----------------------------------------------------------------------
subroutine LUDCMP(N,A,b,x)

implicit none

integer, intent(in) :: &
  N                   ! Number of equations to solve

real, intent(in) :: &
 A(N,N),            & ! Matrix
 b(N)                 ! RHS of matrix equation

real, intent(out) :: &
 x(N)                 ! Solution of matrix equation

integer :: i,ii,imax,j,k,ll,indx(N)

real :: Acp(N,N),aamax,dum,sum,vv(N)

Acp(:,:) = A(:,:)
x(:) = b(:)

do i = 1, N
  aamax = 0
  do j = 1, N
    if (abs(Acp(i,j)) > aamax) aamax = abs(Acp(i,j))
  end do
  vv(i) = 1/aamax
end do

do j = 1, N
  do i = 1, j - 1
    sum = Acp(i,j)
    if (i > 1) then
      do k = 1, i - 1
        sum = sum - Acp(i,k)*Acp(k,j)
      end do
      Acp(i,j) = sum
    end if
  end do
  aamax = 0
  do i = j, N
    sum = Acp(i,j)
    do k = 1, j - 1
      sum = sum - Acp(i,k)*Acp(k,j)
    end do
    Acp(i,j) = sum
    dum = vv(i)*abs(sum)
    if (dum >= aamax) then
      imax = i
      aamax = dum
    end if
  end do
  if (j /= imax)then
    do k = 1, N
      dum = Acp(imax,k)
      Acp(imax,k) = Acp(j,k)
      Acp(j,k) = dum
    end do
    vv(imax) = vv(j)
  end if
  indx(j) = imax
  if (Acp(j,j) == 0) Acp(j,j) = 1e-20
  if (j /= N) then
    dum = 1/Acp(j,j)
    do i = j + 1, N
      Acp(i,j) = Acp(i,j)*dum
    end do
  end if
end do

ii = 0
do i = 1, N
  ll = indx(i)
  sum = x(ll)
  x(ll) = x(i)
  if (ii /= 0)then
    do j = ii, i - 1
      sum = sum - Acp(i,j)*x(j)
    end do
  else if (sum /= 0) then
    ii=i
  end if
  x(i) = sum
end do

do i = N, 1, -1
  sum = x(i)
  do j = i + 1, N
    sum = sum - Acp(i,j)*x(j)
  end do
  x(i) = sum/Acp(i,i)
end do
    
end subroutine LUDCMP

!-----------------------------------------------------------------------
! Saturation specific humidity
!-----------------------------------------------------------------------
subroutine QSAT(P,T,Qs)

use CONSTANTS, only: &
  eps,               &! Ratio of molecular weights of water and dry air
  e0,                &! Saturation vapour pressure at Tm (Pa)
  Tm                  ! Melting point (K)

implicit none

real, intent(in) :: &
  P,                 &! Air pressure (Pa)
  T                   ! Temperature (K)

real, intent(out) :: &
  Qs                  ! Saturation specific humidity

real :: &
  Tc,                &! Temperature (C)
  es                  ! Saturation vapour pressure (Pa)

Tc = T - Tm
if (Tc > 0) then
  es = e0*exp(17.5043*Tc / (241.3 + Tc))
else
  es = e0*exp(22.4422*Tc / (272.186 + Tc))
end if
Qs = eps*es / P

end subroutine QSAT

!-----------------------------------------------------------------------
! Snow thermodynamics and hydrology
!-----------------------------------------------------------------------
subroutine SNOW(dt,Esrf,Gsrf,ksnow,ksoil,Melt,Rf,Sf,Ta,Tsrf,unload,    &
                Nsnow,Dsnw,Rgrn,Sice,Sliq,Tsnow,Tsoil,Gsoil,Roff,snd,SWE)

use CONSTANTS, only: &
  hcap_ice,          &! Specific heat capacity of ice (J/K/kg)
  hcap_wat,          &! Specific heat capacity of water (J/K/kg)
  Lf,                &! Latent heat of fusion (J/kg)
  Ls,                &! Latent heat of sublimation (J/kg)
  rho_ice,           &! Density of ice (kg/m^3)
  rho_wat,           &! Density of water (kg/m^3)
  Tm                  ! Melting point (K)

use LAYERS, only: &
  Dzsnow,            &! Minimum snow layer thicknesses (m)
  Dzsoil,            &! Soil layer thicknesses (m)
  Nsmax,             &! Maximum number of snow layers
  Nsoil               ! Number of soil layers

use PARAMETERS, only: &
  rcld,              &! Maximum density for cold snow (kg/m^3)
  rfix,              &! Fixed snow density (kg/m^3)
  rgr0,              &! Fresh snow grain radius (m)
  rhof,              &! Fresh snow density (kg/m^3)
  rmlt,              &! Maximum density for melting snow (kg/m^3)
  trho,              &! Snow compaction timescale (s)
  Wirr                ! Irreducible liquid water content of snow

implicit none

real, intent(in) :: &
  dt,                &! Timestep (s)
  Esrf,              &! Moisture flux from the surface (kg/m^2/s)
  Gsrf,              &! Heat flux into snow/ground surface (W/m^2)
  Melt,              &! Surface melt (kg/m^2)
  Rf,                &! Rainfall rate (kg/m^2/s)
  Sf,                &! Snowfall rate (kg/m^2/s)
  Ta,                &! Air temperature (K)
  Tsrf,              &! Snow/ground surface temperature (K)
  unload,            &! Snow mass unloaded from vegetation (kg/m^2)
  ksnow(Nsmax),      &! Thermal conductivity of snow layers (W/m/K)
  ksoil(Nsoil)        ! Thermal conductivity of soil layers (W/m/K)

integer, intent(inout) :: &
  Nsnow               ! Number of snow layers

real, intent(inout) :: &
  Dsnw(Nsmax),       &! Snow layer thicknesses (m)
  Rgrn(Nsmax),       &! Snow layer grain radius (m)
  Sice(Nsmax),       &! Ice content of snow layers (kg/m^2)
  Sliq(Nsmax),       &! Liquid content of snow layers (kg/m^2)
  Tsnow(Nsmax),      &! Snow layer temperatures (K)
  Tsoil(Nsoil)        ! Soil layer temperatures (K)

real, intent(out) :: &
  Gsoil,             &! Heat flux into soil (W/m^2)
  Roff,              &! Runoff from snow (kg/m^2)
  snd,               &! Snow depth (m)
  SWE                 ! Total snow mass on ground (kg/m^2) 

integer :: &
  k,                 &! Snow layer counter
  knew,              &! New snow layer pointer
  kold,              &! Old snow layer pointer
  Nold                ! Previous number of snow layers

real :: &
  coldcont,          &! Layer cold content (J/m^2)
  dnew,              &! New snow layer thickness (m)
  dSice,             &! Change in layer ice content (kg/m^2)
  Esnow,             &! Snow sublimation rate (kg/m^2/s)
  ggr,               &! Grain area growth rate (m^2/s)
  mass,              &! Mass of overlying snow (kg/m^2)
  phi,               &! Porosity
  rhos,              &! Density of snow layer (kg/m^3)
  SliqMax,           &! Maximum liquid content for layer (kg/m^2)
  wt                  ! Layer weighting

real :: &
  a(Nsmax),          &! Below-diagonal matrix elements
  b(Nsmax),          &! Diagonal matrix elements
  c(Nsmax),          &! Above-diagonal matrix elements
  csnow(Nsmax),      &! Areal heat capacity of snow (J/K/m^2)
  dTs(Nsmax),        &! Temperature increments (k)
  D(Nsmax),          &! Layer thickness before adjustment (m)
  E(Nsmax),          &! Energy contents before adjustment (J/m^2)
  Gs(Nsmax),         &! Thermal conductivity between layers (W/m^2/k)
  rhs(Nsmax),        &! Matrix equation rhs
  R(Nsmax),          &! Snow grain radii before adjustment (kg/m^2)
  S(Nsmax),          &! Ice contents before adjustment (kg/m^2)
  U(Nsmax),          &! Layer internal energy contents (J/m^2)
  W(Nsmax)            ! Liquid contents before adjustment (kg/m^2)

Gsoil = Gsrf
Roff = Rf*dt

! Existing snowpack
if (Nsnow > 0) then  

  ! Heat conduction
  do k = 1, Nsnow
    csnow(k) = Sice(k)*hcap_ice + Sliq(k)*hcap_wat
  end do
  if (Nsnow == 1) then
    Gs(1) = 2 / (Dsnw(1)/ksnow(1) + Dzsoil(1)/ksoil(1))
    dTs(1) = (Gsrf + Gs(1)*(Tsoil(1) - Tsnow(1)))*dt /  &
             (csnow(1) + Gs(1)*dt)
  else
    do k = 1, Nsnow - 1
      Gs(k) = 2 / (Dsnw(k)/ksnow(k) + Dsnw(k+1)/ksnow(k+1))
    end do
    a(1) = 0
    b(1) = csnow(1) + Gs(1)*dt
    c(1) = - Gs(1)*dt
    rhs(1) = (Gsrf - Gs(1)*(Tsnow(1) - Tsnow(2)))*dt
    do k = 2, Nsnow - 1
      a(k) = c(k-1)
      b(k) = csnow(k) + (Gs(k-1) + Gs(k))*dt
      c(k) = - Gs(k)*dt
      rhs(k) = Gs(k-1)*(Tsnow(k-1) - Tsnow(k))*dt  &
               + Gs(k)*(Tsnow(k+1) - Tsnow(k))*dt 
    end do
    k = Nsnow
    Gs(k) = 2 / (Dsnw(k)/ksnow(k) + Dzsoil(1)/ksoil(1))
    a(k) = c(k-1)
    b(k) = csnow(k) + (Gs(k-1) + Gs(k))*dt
    c(k) = 0
    rhs(k) = Gs(k-1)*(Tsnow(k-1) - Tsnow(k))*dt  &
             + Gs(k)*(Tsoil(1) - Tsnow(k))*dt
    call TRIDIAG(Nsnow,Nsmax,a,b,c,rhs,dTs)
  end if 
  do k = 1, Nsnow
    Tsnow(k) = Tsnow(k) + dTs(k)
  end do
  k = Nsnow
  Gsoil = Gs(k)*(Tsnow(k) - Tsoil(1))

  ! Convert melting ice to liquid water
  dSice = Melt
  do k = 1, Nsnow
    coldcont = csnow(k)*(Tm - Tsnow(k))
    if (coldcont < 0) then
      dSice = dSice - coldcont/Lf
      Tsnow(k) = Tm
    end if
    if (dSice > 0) then
      if (dSice > Sice(k)) then  ! Layer melts completely
        dSice = dSice - Sice(k)
        Dsnw(k) = 0
        Sliq(k) = Sliq(k) + Sice(k)
        Sice(k) = 0
      else                       ! Layer melts partially
        Dsnw(k) = (1 - dSice/Sice(k))*Dsnw(k)
        Sice(k) = Sice(k) - dSice
        Sliq(k) = Sliq(k) + dSice
        dSice = 0                ! Melt exhausted
      end if
    end if
  end do

  ! Remove snow by sublimation 
  dSice = max(Esrf, 0.)*dt
  if (dSice > 0) then
    do k = 1, Nsnow
      if (dSice > Sice(k)) then  ! Layer sublimates completely
        dSice = dSice - Sice(k)
        Dsnw(k) = 0
        Sice(k) = 0
      else                       ! Layer sublimates partially
        Dsnw(k) = (1 - dSice/Sice(k))*Dsnw(k)
        Sice(k) = Sice(k) - dSice
        dSice = 0                ! Sublimation exhausted
      end if
    end do
  end if

  ! Bucket storage 
  do k = 1, Nsnow
    phi = 0
    if (Dsnw(k) > epsilon(Dsnw)) phi = 1 - Sice(k)/(rho_ice*Dsnw(k))
    SliqMax = rho_wat*Dsnw(k)*phi*Wirr
    Sliq(k) = Sliq(k) + Roff
    Roff = 0
    if (Sliq(k) > SliqMax) then       ! Liquid capacity exceeded
      Roff = Sliq(k) - SliqMax        ! so drainage to next layer
      Sliq(k) = SliqMax
    end if
    coldcont = csnow(k)*(Tm - Tsnow(k))
    if (coldcont > 0) then            ! Liquid can freeze
      dSice = min(Sliq(k), coldcont/Lf)
      Sliq(k) = Sliq(k) - dSice
      Sice(k) = Sice(k) + dSice
      Tsnow(k) = Tsnow(k) + Lf*dSice/csnow(k)
    end if
  end do

  ! Snow compaction with age
  do k = 1, Nsnow
    if (Dsnw(k) > epsilon(Dsnw)) then
      rhos = (Sice(k) + Sliq(k)) / Dsnw(k)
      if (Tsnow(k) >= Tm) then
          if (rhos < rmlt) rhos = rmlt + (rhos - rmlt)*exp(-dt/trho)
      else
          if (rhos < rcld) rhos = rcld + (rhos - rcld)*exp(-dt/trho)
      end if
      Dsnw(k) = (Sice(k) + Sliq(k)) / rhos
    end if
  end do

  ! Snow grain growth
  do k = 1, Nsnow
    ggr = 2e-13
    if (Tsnow(k) < Tm) then
      if (Rgrn(k) < 1.50e-4) then
        ggr = 2e-14
      else
        ggr = 7.3e-8*exp(-4600/Tsnow(k))
      end if
    end if
    Rgrn(k) = Rgrn(k) + dt*ggr/Rgrn(k)
  end do

end if  ! Existing snowpack

! Add snowfall and frost to layer 1 with fresh snow density and grain size
Esnow = 0
if (Esrf < 0 .and. Tsrf < Tm) Esnow = Esrf
dSice = (Sf - Esnow)*dt
Dsnw(1) = Dsnw(1) + dSice / rhof
if (Sice(1) + dSice > epsilon(Sice)) then
  Rgrn(1) = (Sice(1)*Rgrn(1) + dSice*rgr0) / (Sice(1) + dSice)
end if
Sice(1) = Sice(1) + dSice

! Add canopy unloading to layer 1 with bulk snow density and grain size
rhos = rhof
SWE = sum(Sice(:)) + sum(Sliq(:))
snd = sum(Dsnw(:))
if (snd > epsilon(snd)) rhos = SWE / snd
Dsnw(1) = Dsnw(1) + unload / rhos
if (Sice(1) + unload > epsilon(Sice)) then
  Rgrn(1) = (Sice(1)*Rgrn(1) + unload*rgr0) / (Sice(1) + unload)
end if
Sice(1) = Sice(1) + unload

! New snowpack
if (Nsnow == 0 .and. Sice(1) > 0) then
  Nsnow = 1
  Rgrn(1) = rgr0
  Tsnow(1) = min(Ta, Tm)
end if

! Store state of old layers
D(:) = Dsnw(:)
R(:) = Rgrn(:)
S(:) = Sice(:)
W(:) = Sliq(:)
do k = 1, Nsnow
  csnow(k) = Sice(k)*hcap_ice + Sliq(k)*hcap_wat
  E(k) = csnow(k)*(Tsnow(k) - Tm)
end do
Nold = Nsnow
snd = sum(Dsnw(:))

! Initialise new layers
Dsnw(:) = 0
Rgrn(:) = 0
Sice(:) = 0
Sliq(:) = 0
Tsnow(:) = Tm
U(:) = 0
Nsnow = 0

if (snd > 0) then  ! Existing or new snowpack

  ! Re-assign and count snow layers
  dnew = snd
  Dsnw(1) = dnew
  k = 1
  if (Dsnw(1) > Dzsnow(1)) then 
    do k = 1, Nsmax
      Dsnw(k) = Dzsnow(k)
      dnew = dnew - Dzsnow(k)
      if (dnew <= Dzsnow(k) .or. k == Nsmax) then
        Dsnw(k) = Dsnw(k) + dnew
        exit
      end if
    end do
  end if
  Nsnow = k

  ! Fill new layers from the top downwards
  knew = 1
  dnew = Dsnw(1)
  do kold = 1, Nold
    do
      if (D(kold) < dnew) then
        ! All snow from old layer partially fills new layer
        Rgrn(knew) = Rgrn(knew) + S(kold)*R(kold)
        Sice(knew) = Sice(knew) + S(kold)
        Sliq(knew) = Sliq(knew) + W(kold)
        U(knew) = U(knew) + E(kold)
        dnew = dnew - D(kold)
        exit
      else
        ! Some snow from old layer fills new layer
        wt = dnew / D(kold)
        Rgrn(knew) = Rgrn(knew) + wt*S(kold)*R(kold)
        Sice(knew) = Sice(knew) + wt*S(kold) 
        Sliq(knew) = Sliq(knew) + wt*W(kold)
        U(knew) = U(knew) + wt*E(kold)
        D(kold) = (1 - wt)*D(kold)
        E(kold) = (1 - wt)*E(kold)
        S(kold) = (1 - wt)*S(kold)
        W(kold) = (1 - wt)*W(kold)
        knew = knew + 1
        if (knew > Nsnow) exit
        dnew = Dsnw(knew)
      end if
    end do
  end do

  ! Diagnose snow layer temperatures
  do k = 1, Nsnow
   csnow(k) = Sice(k)*hcap_ice + Sliq(k)*hcap_wat
   Tsnow(k) = Tm + U(k) / csnow(k)
   Rgrn(k) = Rgrn(k) / Sice(k)
  end do

end if ! Existing or new snowpack

end subroutine SNOW

!-----------------------------------------------------------------------
! Update soil temperatures
!-----------------------------------------------------------------------
subroutine SOIL(csoil,dt,Gsoil,ksoil,Tsoil)

use LAYERS, only: &
  Dzsoil,            &! Soil layer thicknesses (m)
  Nsoil               ! Number of soil layers

implicit none

real, intent(in) :: &
  dt,                &! Timestep (s)
  Gsoil,             &! Heat flux into soil (W/m^2)
  csoil(Nsoil),      &! Areal heat capacity of soil layers (J/K/m^2)
  ksoil(Nsoil)        ! Thermal conductivity of soil layers (W/m/K)

real, intent(inout) :: &
  Tsoil(Nsoil)        ! Soil layer temperatures (K)

integer :: k          ! Soil layer counter

real :: &
  a(Nsoil),          &! Below-diagonal matrix elements
  b(Nsoil),          &! Diagonal matrix elements
  c(Nsoil),          &! Above-diagonal matrix elements
  dTs(Nsoil),        &! Temperature increments (k)
  gs(Nsoil),         &! Thermal conductivity between layers (W/m^2/k)
  rhs(Nsoil)          ! Matrix equation rhs

do k = 1, Nsoil - 1
  gs(k) = 2 / (Dzsoil(k)/ksoil(k) + Dzsoil(k+1)/ksoil(k+1))
end do
a(1) = 0
b(1) = csoil(1) + gs(1)*dt
c(1) = - gs(1)*dt
rhs(1) = (Gsoil - gs(1)*(Tsoil(1) - Tsoil(2)))*dt
do k = 2, Nsoil - 1
  a(k) = c(k-1)
  b(k) = csoil(k) + (gs(k-1) + gs(k))*dt
  c(k) = - gs(k)*dt
  rhs(k) = gs(k-1)*(Tsoil(k-1) - Tsoil(k))*dt  &
           + gs(k)*(Tsoil(k+1) - Tsoil(k))*dt 
end do
k = Nsoil
gs(k) = ksoil(k)/Dzsoil(k)
a(k) = c(k-1)
b(k) = csoil(k) + (gs(k-1) + gs(k))*dt
c(k) = 0
rhs(k) = gs(k-1)*(Tsoil(k-1) - Tsoil(k))*dt
call TRIDIAG(Nsoil,Nsoil,a,b,c,rhs,dTs)
do k = 1, Nsoil
  Tsoil(k) = Tsoil(k) + dTs(k)
end do

end subroutine SOIL

!-----------------------------------------------------------------------
! Surface energy balance
!-----------------------------------------------------------------------
subroutine SRFEBAL(cveg,Ds1,dt,fvegs,ga,gc,gs,gs1,gv,ks1,LW,Ps,Qa,     &
                   SWsrf,SWveg,Ta,tdif,Ts1,Tveg0,VAI,                  &
                   Tsrf,Qcan,Sice,Tcan,Tveg,                           &
                   Esrf,Eveg,Gsrf,H,LE,LWout,LWsub,Melt)

use CONSTANTS, only : &
  cp,                &! Specific heat capacity of air (J/K/kg)
  Lf,                &! Latent heat of fusion (J/kg)
  Ls,                &! Latent heat of sublimation (J/kg)
  Lv,                &! Latent heat of vapourisation (J/kg)
  Rair,              &! Gas constant for air (J/K/kg)
  Rwat,              &! Gas constant for water vapour (J/K/kg)
  sb,                &! Stefan-Boltzmann constant (W/m^2/K^4)
  Tm                  ! Melting point (K)

use LAYERS, only: &
  Ncnpy,             &! Number of canopy layers
  Nsmax               ! Maximum number of snow layers

use PARAMETERS, only: &
  gsnf                ! Snow-free vegetation moisture conductance (m/s)

implicit none

real, intent(in) :: &
  Ds1,               &! Surface layer thickness (m)
  dt,                &! Timestep (s)
  ga,                &! Aerodynamic conductance to the atmosphere (m/s)
  gc,                &! Conductance within canopy air space (m/s)
  gs,                &! Surface to canopy air space conductance (m/s)
  gs1,               &! Surface moisture conductance (m/s)
  ks1,               &! Surface layer thermal conductivity (W/m/K)
  LW,                &! Incoming longwave radiation (W/m2)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  SWsrf,             &! SW absorbed by snow/ground surface (W/m^2)
  Ta,                &! Air temperature (K)
  Ts1,               &! Surface layer temperature (K)
  VAI,               &! Vegetation area index
  cveg(Ncnpy),       &! Vegetation layer heat capacities (J/K/m^2)
  fvegs(Ncnpy),      &! Vegetation layer snowcover fractions
  gv(Ncnpy),         &! Vegetation to canopy air space conductance (m/s)
  SWveg(Ncnpy),      &! SW absorbed by vegetation layers (W/m^2)
  tdif(Ncnpy),       &! Veg layer transmissivities for diffuse radiation
  Tveg0(Ncnpy)        ! Veg layer temperatures at start of timestep (K)

real, intent(inout) :: &
  Tsrf,              &! Snow/ground surface temperature (K)
  Qcan(Ncnpy),       &! Canopy air space humidities
  Tcan(Ncnpy),       &! Canopy air space temperatures (K)
  Tveg(Ncnpy),       &! Vegetation layer temperatures (K)
  Sice(Nsmax)         ! Ice content of snow layers (kg/m^2)

real, intent(out) :: &
  Esrf,              &! Moisture flux from the surface (kg/m^2/s)
  Gsrf,              &! Heat flux into snow/ground surface (W/m^2)
  H,                 &! Sensible heat flux to the atmosphere (W/m^2)
  LE,                &! Latent heat flux to the atmosphere (W/m^2)
  LWout,             &! Outgoing LW radiation (W/m^2)
  LWsub,             &! Subcanopy downward LW radiation (W/m^2)
  Melt,              &! Surface melt (kg/m^2)
  Eveg(Ncnpy)         ! Moisture flux from vegetation layers (kg/m^2/s)

real :: &
  A(3*Ncnpy+1,3*Ncnpy+1),&! Jacobian of energy and mass balance equations
  b(3*Ncnpy+1),      &! Residuals of energy and mass balance equations
  x(3*Ncnpy+1)        ! Temperature and humidity increments

real :: &
  Dsrf,              &! dQsat/dT at ground surface temperature (1/K)
  dEs,               &! Change in surface moisture flux (kg/m^2/s)
  dGs,               &! Change in surface heat flux (kg/m^2/s)
  dHs,               &! Change in surface sensible heat flux (kg/m^2/s)
  dTs,               &! Change in surface temperature (K)
  E,                 &! Moisture flux to the atmosphere (kg/m^2/s)
  Ecan,              &! Within-canopy moisture flux (kg/m^2/s)
  Hcan,              &! Within-canopy sensible heat flux (W/m^2)
  Hsrf,              &! Sensible heat flux from the surface (W/m^2)
  Lsrf,              &! Latent heat for phase change on ground (J/kg)
  Qsrf,              &! Saturation humidity at surface temperature
  rho,               &! Air density (kg/m^3)
  Rsrf,              &! Net radiation absorbed by the surface (W/m^2)
  Ssub,              &! Mass of snow available for sublimation (kg/m^2)
  wsrf                ! Surface water availability factor

real :: &
  dEv(Ncnpy),        &! Change in vegetation moisture flux (kg/m^2/s)
  dHv(Ncnpy),        &! Change in veg sensible heat flux (kg/m^2/s)
  dQc(Ncnpy),        &! Change in canopy air humidity (kg/kg)
  dTv(Ncnpy),        &! Change in vegetation temperature (K)
  dTc(Ncnpy),        &! Change in canopy air temperature (K)
  Dveg(Ncnpy),       &! dQsat/dT at vegetation layer temperature (1/K)
  Hveg(Ncnpy),       &! Sensible heat flux from vegetation (W/m^2)
  Lveg(Ncnpy),       &! Latent heat for phase change on veg (J/kg)
  Qveg(Ncnpy),       &! Saturation humidity at vegetation temperature
  Rveg(Ncnpy),       &! Net radiation absorbed by vegetation (W/m^2)
  wveg(Ncnpy)         ! Vegetation water availability factor

integer :: k          ! Vegetation layer counter

! Saturation humidity and air density
call QSAT(Ps,Tsrf,Qsrf)
Lsrf = Ls
if (Tsrf > Tm) Lsrf = Lv
Dsrf = Lsrf*Qsrf/(Rwat*Tsrf**2)
rho = Ps/(Rair*Ta)

Gsrf = 2*ks1*(Tsrf - Ts1)/Ds1

if (VAI == 0) then  ! open

  ! Surface water availability
  if (Sice(1) > 0 .or. Qa > Qsrf) then
    wsrf = 1
  else
    wsrf = gs1/(gs1 + ga)
  end if

  ! Explicit fluxes
  Esrf = rho*wsrf*ga*(Qsrf - Qa)
  Eveg = 0
  Hsrf = cp*rho*ga*(Tsrf - Ta)
  Hveg = 0
  Melt = 0
  Rsrf = SWsrf + LW - sb*Tsrf**4

  ! Surface energy balance increments without melt
  dTs = (Rsrf - Gsrf - Hsrf - Lsrf*Esrf) /  &
        (4*sb*Tsrf**3 + 2*ks1/Ds1 + rho*(cp + Lsrf*Dsrf*wsrf)*ga)
  dEs = rho*wsrf*ga*Dsrf*dTs
  dGs = 2*ks1*dTs/Ds1 
  dHs = cp*rho*ga*dTs

  ! Surface melting
  if (Tsrf + dTs > Tm .and. Sice(1) > 0) then
    Melt = sum(Sice)
    dTs = (Rsrf - Gsrf - Hsrf - Lsrf*Esrf - Lf*Melt/dt) /  &
          (4*sb*Tsrf**3 + 2*ks1/Ds1 + rho*(cp + Ls*Dsrf*wsrf)*ga)
    dEs = rho*wsrf*ga*Dsrf*dTs
    dGs = 2*ks1*dTs/Ds1
    dHs = cp*rho*ga*dTs
    if (Tsrf + dTs < Tm) then
      call QSAT(Ps,Tm,Qsrf)
      Esrf = rho*wsrf*ga*(Qsrf - Qa)  
      Gsrf = 2*ks1*(Tm - Ts1)/Ds1
      Hsrf = cp*rho*ga*(Tm - Ta)
      Rsrf = SWsrf + LW - sb*Tm**4 
      Melt = (Rsrf - Gsrf - Hsrf - Lsrf*Esrf)*dt/Lf
      Melt = max(Melt, 0.)
      dEs = 0
      dGs = 0
      dHs = 0
      dTs = Tm - Tsrf
    end if
  end if

  LWout = sb*Tsrf**4
  LWsub = LW

else ! forest

  ! Saturation humidity
  do k = 1, Ncnpy
    call QSAT(Ps,Tveg(k),Qveg(k))
    Lveg(k) = Ls
    if (Tveg(k) > Tm) Lveg(k) = Lv
  end do
  Dveg(:) = Lveg(:)*Qveg(:)/(Rwat*Tveg(:)**2)

  ! Water availability
  if (Sice(1) > 0 .or. Qcan(Ncnpy) > Qsrf) then
    wsrf = 1
  else
    wsrf = gs1/(gs1 + gs)
  end if
  do k = 1, Ncnpy
    if (Qcan(k) > Qveg(k)) then
      wveg(k) = 1
    else
      wveg(k) = fvegs(k) + (1 - fvegs(k))*gsnf/(gsnf + gv(k))
    end if
  end do

! 2-layer canopy model

  ! Explicit fluxes
  E = rho*ga*(Qcan(1) - Qa)
  Ecan = rho*gc*(Qcan(2) - Qcan(1))
  Esrf = rho*wsrf*gs*(Qsrf - Qcan(2))
  Eveg(1) = rho*wveg(1)*gv(1)*(Qveg(1) - Qcan(1))
  Eveg(2) = rho*wveg(2)*gv(2)*(Qveg(2) - Qcan(2))
  H = rho*cp*ga*(Tcan(1) - Ta)
  Hcan = rho*cp*gs*(Tcan(2) - Tcan(1))
  Hsrf = rho*cp*gs*(Tsrf - Tcan(2))
  Hveg(1) = rho*cp*gv(1)*(Tveg(1) - Tcan(1))
  Hveg(2) = rho*cp*gv(2)*(Tveg(2) - Tcan(2))
  Melt = 0
  Rsrf = SWsrf + tdif(1)*tdif(2)*LW +              &
         (1 - tdif(1))*tdif(2)*sb*Tveg(1)**4 +     &
         (1 - tdif(2))*sb*Tveg(1)**4 - sb*Tsrf**4 
  Rveg(1) = SWveg(1) + (1 - tdif(1))*(LW - 2*sb*Tveg(1)**4    & 
            + (1 - tdif(2))*sb*Tveg(2)**4 + tdif(2)*sb*Tsrf**4) 
  Rveg(2) = SWveg(2) +                                               & 
            (1 - tdif(2))*(tdif(1)*LW + (1 - tdif(1))*sb*Tveg(1)**4  & 
               - 2*sb*Tveg(2)**4 + sb*Tsrf**4) 

! Surface energy balance increments without melt
  A(1,1) = - rho*gs*(cp + Lsrf*Dsrf*wsrf) - 4*sb*Tsrf**3 - 2*ks1/Ds1
  A(1,2) = 0
  A(1,3) = 0
  A(1,4) = 4*(1 - tdif(1))*tdif(2)*sb*Tveg(1)**3
  A(1,5) = Lsrf*rho*wsrf*gs
  A(1,6) = rho*cp*gs
  A(1,7) = 4*(1 - tdif(2))*sb*Tveg(2)**3
  A(2,1) = 4*(1 - tdif(1))*tdif(2)*sb*Tveg(2)**3 
  A(2,2) = Lveg(1)*rho*wveg(1)*gv(1)
  A(2,3) = rho*cp*gv(1)
  A(2,4) = - rho*gv(1)*(cp + Lveg(1)*Dveg(1)*wveg(1))    & 
           - 8*(1 - tdif(1))*sb*Tveg(1)**3 - cveg(1)/dt
  A(2,5) = 0
  A(2,6) = 0
  A(2,7) = 4*(1 - tdif(1))*(1 - tdif(2))*sb*Tveg(2)**3
  A(3,1) = 4*(1 - tdif(2))*sb*Tsrf**3 
  A(3,2) = 0
  A(3,3) = 0
  A(3,4) = 4*(1 - tdif(1))*(1 - tdif(2))*sb*Tveg(1)**3
  A(3,5) = Lveg(2)*rho*wveg(2)*gv(2)
  A(3,6) = rho*cp*gv(2)
  A(3,7) = - rho*gv(2)*(cp + Lveg(2)*Dveg(2)*wveg(2))    &
           - 8*(1 - tdif(2))*sb*Tveg(2)**3 - cveg(2)/dt
  A(4,1) = 0
  A(4,2) = 0
  A(4,3) = ga + gc + gv(1)
  A(4,4) = -gv(1)
  A(4,5) = 0
  A(4,6) = -gc
  A(4,7) = 0
  A(5,1) = -gs
  A(5,2) = 0
  A(5,3) = -gc
  A(5,4) = 0
  A(5,5) = 0
  A(5,6) = gc + gs + gv(2)
  A(5,7) = -gv(2)
  A(6,1) = 0
  A(6,2) = ga + gc + wveg(1)*gv(1)
  A(6,3) = 0
  A(6,4) = -Dveg(1)*wveg(1)*gv(1)
  A(6,5) = -gc
  A(6,6) = 0
  A(6,7) = 0
  A(7,1) = -Dsrf*wsrf*gs
  A(7,2) = -gc
  A(7,3) = 0
  A(7,4) = 0
  A(7,5) = gc + wsrf*gs + wveg(2)*gv(2)
  A(7,6) = 0
  A(7,7) = -Dveg(2)*wveg(2)*gv(2)
  b(1)   = -(Rsrf - Gsrf - Hsrf - Lsrf*Esrf)
  b(2)   = -(Rveg(1) - Hveg(1) - Lveg(1)*Eveg(1) -  & 
             cveg(1)*(Tveg(1) - Tveg0(1))/dt)
  b(3)   = -(Rveg(2) - Hveg(2) - Lveg(2)*Eveg(2)    &
             - cveg(2)*(Tveg(2) - Tveg0(2))/dt)
  b(4)   = -(H - Hcan - Hveg(1))/(rho*cp)
  b(5)   = -(Hcan - Hsrf - Hveg(2))/(rho*cp)
  b(6)   = -(E - Ecan - Eveg(1))/rho
  b(7)   = -(Ecan - Esrf - Eveg(2))/rho
  call LUDCMP(7,A,b,x)
  dTs    = x(1)
  dQc(1) = x(2)
  dTc(1) = x(3)
  dTv(1) = x(4)
  dQc(2) = x(5)
  dTc(2) = x(6)
  dTv(2) = x(7)
  dEs = rho*wsrf*gs*(Dsrf*dTs - dQc(2))
  dEv(:) = rho*wveg(:)*gv(:)*(Dveg(:)*dTv(:) - dQc(:))
  dGs = 2*ks1*dTs/Ds1
  dHs = rho*cp*gs*(dTs - dTc(2))
  dHv(:) = rho*cp*gv(:)*(dTv(:) - dTc(:))

  ! Surface melting
  if (Tsrf + dTs > Tm .and. Sice(1) > 0) then
    Melt = sum(Sice)
    b(1) = b(1) + Lf*Melt/dt
    call LUDCMP(7,A,b,x)
    dTs    = x(1)
    dQc(1) = x(2)
    dTc(1) = x(3)
    dTv(1) = x(4)
    dQc(2) = x(5)
    dTc(2) = x(6)
    dTv(2) = x(7)
    dEs = rho*wsrf*gs*(Dsrf*dTs - dQc(2))
    dEv(:) = rho*wveg(:)*gv(:)*(Dveg(:)*dTv(:) - dQc(:))
    dGs = 2*ks1*dTs/Ds1
    dHs = rho*cp*gs*(dTs - dTc(2))
    dHv(:) = rho*cp*gv(:)*(dTv(:) - dTc(:))
    if (Tsrf + dTs < Tm) then
      call QSAT(Ps,Tm,Qsrf)
      Esrf = rho*wsrf*gs*(Qsrf - Qcan(2))
      Hsrf = rho*cp*gs*(Tm - Tcan(2))
      Rsrf = Rsrf + sb*Tsrf**4 - sb*Tm**4
      Rveg(1) = Rveg(1) + (1 - tdif(1))*tdif(2)*sb*(Tm**4 - Tsrf**4) 
      Rveg(2) = Rveg(2) + (1 - tdif(2))*sb*(Tm**4 - Tsrf**4)
      A(1,1) = -1
      A(2,1) = 0
      A(3,1) = 0
      A(4,1) = 0
      A(5,1) = 0
      A(6,1) = 0
      A(7,1) = 0
      b(1) = -(Rsrf - Gsrf - Hsrf - Lsrf*Esrf)
      b(2) = -(Rveg(1) - Hveg(1) - Lveg(1)*Eveg(1)    & 
               - cveg(1)*(Tveg(1) - Tveg0(1))/dt)
      b(3) = -(Rveg(2) - Hveg(2) - Lveg(2)*Eveg(2)    &
               - cveg(2)*(Tveg(2) - Tveg0(2))/dt)
      b(4) = -(H - Hcan - Hveg(1))/(rho*cp)
      b(5) = -(Hcan - Hsrf - Hveg(2))/(rho*cp)
      b(6) = -(E - Ecan - Eveg(1))/rho
      b(7) = -(Ecan - Esrf - Eveg(2))/rho
      call LUDCMP(7,A,b,x)
      Melt = x(1)*dt/Lf
      dQc(1) = x(2)
      dTc(1) = x(3)
      dTv(1) = x(4)
      dQc(2) = x(5)
      dTc(2) = x(6)
      dTv(2) = x(7)
      dTs = Tm - Tsrf
      dEs = 0
      dEv(:) = rho*wveg(:)*gv(:)*(Dveg(:)*dTv(:) - dQc(:))
      dGs = 0
      dHs = 0
      dHv(:) = rho*cp*gv(:)*(dTv(:) - dTc(:))
    end if
  end if
  LWout = (1 - tdif(1))*sb*Tveg(1)**4 +          &
          (1 - tdif(2))*tdif(1)*sb*Tveg(1)**4 +  &
          tdif(1)*tdif(2)*sb*Tsrf**4
  LWsub = tdif(1)*tdif(2)*LW +                   &
          (1 - tdif(1))*tdif(2)*sb*Tveg(1)**4 +  &
          (1 - tdif(2))*sb*Tveg(2)**4

  ! Update vegetation temperatures and fluxes
  Eveg(:) = Eveg(:) + dEv(:)
  Hveg(:) = Hveg(:) + dHv(:)
  Qcan(:) = Qcan(:) + dQc(:)
  Tcan(:) = Tcan(:) + dTc(:)
  Tveg(:) = Tveg(:) + dTv(:)

end if  ! forest

! Update surface temperature and fluxes
Esrf = Esrf + dEs
Gsrf = Gsrf + dGs
Hsrf = Hsrf + dHs
Tsrf = Tsrf + dTs

! Sublimation limited by amount of snow after melt
Ssub = sum(Sice(:)) - Melt
if (Ssub > 0 .and. Esrf*dt > Ssub) then
  Esrf = Ssub / dt
  Hsrf = Rsrf - Gsrf - Ls*Esrf - Lf*Melt/dt
end if

! Fluxes to the atmosphere
E = Esrf + sum(Eveg(:))
H = Hsrf + sum(Hveg(:))
LE = Lsrf*Esrf + sum(Lveg(:)*Eveg(:))

end subroutine SRFEBAL


!-----------------------------------------------------------------------
! Surface exchange coefficients
!-----------------------------------------------------------------------
subroutine SRFEXCH(aveg,fsnow,hveg,Ua,VAI,zT,zU,ga,gc,gs,gv,Usub)

use CONSTANTS, only: &
  vkman               ! Von Karman constant

use LAYERS, only: &
  Ncnpy               ! Number of canopy layers

use PARAMETERS, only: &
  rveg,              &! Leaf boundary resistance (s/m)^(1/2)
  fvg1,              &! Fraction of vegetation in upper layer
  wcan,              &! Canopy wind decay coefficient
  zsub,              &! Subcanopy wind speed diagnostic height (m)
  z0sf,              &! Snow-free surface roughness length (m)
  z0sn                ! Snow roughness length (m)

implicit none

real, intent(in) :: &
  aveg(Ncnpy),       &! Vegetation area indices in layers
  fsnow,             &! Ground snowcover fraction
  hveg,              &! Vegetation height (m)
  Ua,                &! Wind speed (m/s)
  VAI,               &! Vegetation area index
  zT,                &! Temperature and humidity measurement height (m)
  zU                  ! Wind speed measurement height (m)

real, intent(out) :: &
  ga,                &! Aerodynamic conductance to the atmosphere (m/s)
  gc,                &! Conductance within canopy air space (m/s)
  gs,                &! Surface to canopy air space conductance (m/s)
  gv(Ncnpy),         &! Vegetation to canopy air space conductance (m/s)
  Usub                ! Subcanopy wind speed (m/s)

integer :: k          ! Vegetation layer counter

real :: &
  d,                 &! Displacement height (m)
  fveg,              &! Vegetation weighting
  Kh,                &! Eddy diffusivity at canopy top (m2/s)
  rd,                &! Dense vegetation aerodynamic resistance (s/m)
  ro,                &! Open aerodynamic resistance (s/m)
  Uc,                &! Within-canopy wind speed (m/s)
  Uh,                &! Wind speed at canopy top (m/s)
  ustar,             &! Friction velocity (m/s)
  zT1,               &! Temperature measurement height with offset (m)
  zU1,               &! Wind measurement height with offset (m)
  z0,                &! Roughness length for momentum (m)
  z0g,               &! Snow/ground surface roughness length (m)
  z0h,               &! Roughness length for heat (m)
  zh(Ncnpy)           ! Dimensionless vegetation layer heights

! Heights specified above ground
zU1 = zU
zT1 = zT

z0g = (z0sn**fsnow) * (z0sf**(1 - fsnow))
Usub = Ua*log(zsub/z0g)/log(zU1/z0g)
if (VAI == 0) then
  z0h = 0.1*z0g
  ustar = vkman*Ua/log(zU1/z0g)
  ga = vkman*ustar/log(zT1/z0h)
else
  d = 0.67*hveg
  z0 = 0.05*hveg
  fveg = 1 - exp(-0.5*VAI)
  ustar = vkman*Ua/log((zU1 - d)/z0)
  Uh = (ustar/vkman)*log((hveg - d)/z0)
  Kh = vkman*ustar*(hveg - d)
  zh(1) = 1 - fvg1/2
  zh(2) = (1 - fvg1)/2
  rd = log((zT1 - d)/(hveg - d))/(vkman*ustar) +  & 
       hveg*(exp(wcan*(1 - zh(1))) - 1)/(wcan*Kh)
  ro = log(zT1/zh(1))/(vkman*ustar)
  ga = fveg/rd + (1 - fveg)/ro
  do k = 1, Ncnpy
    Uc = exp(wcan*(zh(k) - 1))*Uh
    gv(k) = sqrt(Uc)*aveg(k)/rveg
  end do
  rd = hveg*exp(wcan)*(exp(-wcan*zh(2)) - exp(-wcan*zh(1)))/(wcan*Kh)
  ro = log(zh(1)/zh(2))/(vkman*ustar)
  gc = fveg/rd + (1 - fveg)/ro
  rd = hveg*exp(wcan)*(exp(-wcan*z0g/hveg) - exp(-wcan*zh(2)))/(wcan*Kh)
  ro = log(zh(2)*hveg/z0g)/(vkman*ustar)
  gs = fveg/rd + (1 - fveg)/ro
  Usub = fveg*exp(wcan*(zsub/hveg - 1))*Uh + (1 - fveg)*Usub
end if

end subroutine SRFEXCH

!-----------------------------------------------------------------------
! Surface and vegetation net shortwave radiation
!-----------------------------------------------------------------------
subroutine SWRAD(alb0,aveg,Dsnw,dt,elev,fvegs,Sdif,Sdir,Sf,Tsrf,       &
                 albs,fsnow,SWout,SWsrf,SWsub,SWveg,tdif)

use CONSTANTS, only: &
  Tm                  ! Melting point (K)

use LAYERS, only: &
  Ncnpy,             &! Number of canopy layers
  Nsmax               ! Maximum number of snow layers

use PARAMETERS, only: &
  asmx,              &! Maximum albedo for fresh snow
  asmn,              &! Minimum albedo for melting snow
  avg0,              &! Snow-free vegetation albedo
  avgs,              &! Snow-covered vegetation albedo
  hfsn,              &! Snowcover fraction depth scale (m)
  kext,              &! Vegetation light extinction coefficient
  Salb,              &! Snowfall to refresh albedo (kg/m^2)
  Talb,              &! Snow albedo decay temperature threshold (C)
  tcld,              &! Cold snow albedo decay time scale (s)
  tmlt                ! Melting snow albedo decay time scale (s)
 
implicit none

real, intent(in) :: &
  alb0,              &! Snow-free ground albedo
  dt,                &! Timestep (s)
  elev,              &! Solar elevation (radians)
  Sdif,              &! Diffuse shortwave radiation (W/m^2)
  Sdir,              &! Direct-beam shortwave radiation (W/m^2)
  Sf,                &! Snowfall rate (kg/m2/s)
  Tsrf,              &! Snow/ground surface temperature (K)
  aveg(Ncnpy),       &! Vegetation area indices
  Dsnw(Nsmax),       &! Snow layer thicknesses (m)
  fvegs(Ncnpy)        ! Vegetation layer snowcover fractions

real, intent(inout) :: &
  albs                ! Snow albedo

real, intent(out) :: &
  fsnow,             &! Ground snowcover fraction
  SWout,             &! Outgoing SW radiation (W/m^2)
  SWsrf,             &! SW absorbed by snow/ground surface (W/m^2)
  SWsub,             &! Subcanopy downward SW radiation (W/m^2)
  SWveg(Ncnpy),      &! SW absorbed by vegetation layers (W/m^2)
  tdif(Ncnpy)         ! Veg layer transmissivities for diffuse radiation

integer :: k          ! Vegetation layer counter

real :: &
  alim,              &! Limiting snow albedo
  asrf,              &! Snow/ground surface albedo
  snwd,              &! Snow depth (m)
  tdec                ! Snow albedo decay time scale (s)

real :: &
  A(2*Ncnpy+1,2*Ncnpy+1),   &! Vegetation radiative transfer matrix
  b(2*Ncnpy+1),      &! Vegetation layer boundary SW fluxes (W/m^2)
  x(2*Ncnpy+1),      &! Vegetation SW sources (W/m^2)
  acan(Ncnpy),       &! Dense canopy albedo
  tdir(Ncnpy)         ! Veg layer transmissivities for direct radiation

! Prognostic snow albedo
tdec = tcld
if (Tsrf >= Tm) tdec = tmlt
alim = (asmn/tdec + asmx*Sf/Salb)/(1/tdec + Sf/Salb)
albs = alim + (albs - alim)*exp(-(1/tdec + Sf/Salb)*dt)

albs = max(min(albs,asmx),asmn)

! Partial snowcover on ground
snwd = sum(Dsnw(:))

fsnow = tanh(snwd/hfsn)

! Surface and vegetation net shortwave radiation
asrf = (1 - fsnow)*alb0 + fsnow*albs
SWsrf = (1 - asrf)*(Sdif + Sdir)
SWveg(:) = 0
SWout = asrf*(Sdif + Sdir)
SWsub = Sdif + Sdir
tdif(:) = 0
tdir(:) = 0
if (aveg(1) > 0) then
  acan(:) = (1 - fvegs(:))*avg0 + fvegs(:)*avgs
  tdif(:) = exp(-1.6*kext*aveg(:))
  if (elev > 0) tdir(:) = exp(-kext*aveg(:)/sin(elev))
  A(:,:) = 0
  do k = 1, 2*Ncnpy + 1
    A(k,k) = 1
  end do
  A(1,4) = -(1 - tdif(1))*acan(1)
  A(2,1) = -tdif(2)
  A(2,3) = -(1 - tdif(2))*acan(2)
  A(3,2) = -asrf
  A(4,1) = -(1 - tdif(2))*acan(2)
  A(4,3) = -tdif(2)
  A(5,4) = -tdif(1)
  b(1) = tdif(1)*Sdif
  b(2) = 0
  b(3) = asrf*tdir(1)*tdir(2)*Sdir
  b(4) = (1 - tdif(2))*acan(2)*tdir(1)*Sdir
  b(5) = (1 - tdif(1))*acan(1)*Sdif + (1 - tdir(1))*acan(1)*Sdir
  call LUDCMP(5,A,b,x)
  SWout = x(5)
  SWveg(1) = Sdif - x(1) + x(4) - x(5) + (1 - tdir(1))*Sdir
  SWveg(2) = x(1) - x(2) + x(3) - x(4) + tdir(1)*(1 - tdir(2))*Sdir
  SWsub = x(2) + tdir(1)*tdir(2)*Sdir
  SWsrf = (1 - asrf)*SWsub
end if

end subroutine SWRAD

!-----------------------------------------------------------------------
! Thermal properties of snow and soil
!-----------------------------------------------------------------------
subroutine THERMAL(Dsnw,Nsnow,Sice,Sliq,Tsnow,Tsoil,Vsmc,              &
                   csoil,Ds1,gs1,ksnow,ksoil,ks1,Ts1)

use CONSTANTS, only: &
  grav,              &! Acceleration due to gravity (m/s^2)
  hcap_ice,          &! Specific heat capacity of ice (J/K/kg)
  hcap_wat,          &! Specific heat capacity of water (J/K/kg)
  hcon_air,          &! Thermal conductivity of air (W/m/K)
  hcon_ice,          &! Thermal conducivity of ice (W/m/K)
  hcon_wat,          &! Thermal conductivity of water (W/m/K)
  Lf,                &! Latent heat of fusion (J/kg)
  rho_ice,           &! Density of ice (kg/m^3)
  rho_wat,           &! Density of water (kg/m^3)
  Tm                  ! Melting point (K)

use LAYERS, only: &
  Dzsoil,            &! Soil layer thicknesses (m)
  Nsmax,             &! Maximum number of snow layers
  Nsoil               ! Number of soil layers

use PARAMETERS, only: &
  b,                 &! Clapp-Hornberger exponent
  gsat,              &! Surface conductance for saturated soil (m/s)
  hcap_soil,         &! Volumetric heat capacity of dry soil (J/K/m^3)
  hcon_soil,         &! Thermal conductivity of dry soil (W/m/K)
  kfix,              &! Fixed thermal conductivity of snow (W/m/K)
  rhof,              &! Fresh snow density (kg/m^3)
  sathh,             &! Saturated soil water pressure (m)
  Vcrit,             &! Volumetric soil moisture at critical point
  Vsat                ! Volumetric soil moisture at saturation

implicit none

integer, intent(in) :: &
  Nsnow               ! Number of snow layers

real, intent(in) :: &
  Dsnw(Nsmax),       &! Snow layer thicknesses (m)
  Sice(Nsmax),       &! Ice content of snow layers (kg/m^2)
  Sliq(Nsmax),       &! Liquid content of snow layers (kg/m^2)
  Tsnow(Nsmax),      &! Snow layer temperatures (K)
  Tsoil(Nsoil),      &! Soil layer temperatures (K)
  Vsmc(Nsoil)         ! Volumetric soil moisture content in layers

real, intent(out) :: &
  Ds1,               &! Surface layer thickness (m)
  gs1,               &! Surface moisture conductance (m/s)
  ks1,               &! Surface layer thermal conductivity (W/m/K)
  Ts1,               &! Surface layer temperature (K)
  csoil(Nsoil),      &! Areal heat capacity of soil layers (J/K/m^2)
  ksnow(Nsmax),      &! Thermal conductivity of snow layers (W/m/K)
  ksoil(Nsoil)        ! Thermal conductivity of soil layers (W/m/K)

integer :: &
  k                   ! Level counter

real :: &
  dPsidT,            &! d(ice potential)/dT (m/K)
  dthudT,            &! d(unfrozen soil moisture)/dT (1/K)
  hcon_sat,          &! Thermal conductivity of saturated soil (W/m/K)
  Mf,                &! Frozen moisture content of soil layer (kg/m^2)
  Mu,                &! Unfrozen moisture content of soil layer (kg/m^2)
  rhos,              &! Snow density (kg/m^3)
  Smf,               &! Fractional frozen soil moisture content
  Smu,               &! Fractional unfrozen soil moisture content
  sthf,              &! Frozen soil moisture content
  sthu,              &! Unfrozen soil moisure content
  Tc,                &! Soil temperature (C)
  thice,             &! Soil ice saturation at current liquid/ice ratio
  thwat,             &! Soil water saturation at current liquid/ice ratio
  Tmax                ! Maximum temperature for frozen soil moisture (K)

! Thermal conductivity of snow
ksnow = kfix
do k = 1, Nsnow
  rhos = rhof
  if (Dsnw(k) > epsilon(Dsnw)) rhos = (Sice(k) + Sliq(k)) / Dsnw(k)

  ksnow(k) = 2.224*(rhos/rho_wat)**1.885
end do

! Heat capacity and thermal conductivity of soil
dPsidT = - rho_ice*Lf/(rho_wat*grav*Tm)
do k = 1, Nsoil
  csoil(k) = hcap_soil*Dzsoil(k)
  ksoil(k) = hcon_soil
  if (Vsmc(k) > epsilon(Vsmc)) then
    dthudT = 0
    sthu = Vsmc(k)
    sthf = 0
    Tc = Tsoil(k) - Tm
    Tmax = Tm + (sathh/dPsidT)*(Vsat/Vsmc(k))**b
    if (Tsoil(k) < Tmax) then
      dthudT = (-dPsidT*Vsat/(b*sathh)) *  &
               (dPsidT*Tc/sathh)**(-1/b - 1)
      sthu = Vsat*(dPsidT*Tc/sathh)**(-1/b)
      sthu = min(sthu, Vsmc(k))
      sthf = (Vsmc(k) - sthu)*rho_wat/rho_ice
    end if
    Mf = rho_ice*Dzsoil(k)*sthf
    Mu = rho_wat*Dzsoil(k)*sthu
    csoil(k) = hcap_soil*Dzsoil(k) + hcap_ice*Mf + hcap_wat*Mu +       &
               rho_wat*Dzsoil(k)*((hcap_wat - hcap_ice)*Tc + Lf)*dthudT
    Smf = rho_ice*sthf/(rho_wat*Vsat)
    Smu = sthu/Vsat
    thice = 0
    if (Smf > 0) thice = Vsat*Smf/(Smu + Smf) 
    thwat = 0
    if (Smu > 0) thwat = Vsat*Smu/(Smu + Smf)
    hcon_sat = hcon_soil*(hcon_wat**thwat)*(hcon_ice**thice) /         &
              (hcon_air**Vsat)
    ksoil(k) = (hcon_sat - hcon_soil)*(Smf + Smu) + hcon_soil
    if (k == 1) gs1 = gsat*max((Smu*Vsat/Vcrit)**2, 1.)
  end if
end do

! Surface layer
Ds1 = max(Dzsoil(1), Dsnw(1))
Ts1 = Tsoil(1) + (Tsnow(1) - Tsoil(1))*Dsnw(1)/Dzsoil(1)
ks1 = Dzsoil(1)/(2*Dsnw(1)/ksnow(1) + (Dzsoil(1) - 2*Dsnw(1))/ksoil(1))
if (Dsnw(1) > 0.5*Dzsoil(1)) ks1 = ksnow(1)
if (Dsnw(1) > Dzsoil(1)) Ts1 = Tsnow(1)

end subroutine THERMAL

!-----------------------------------------------------------------------
! Solve tridiagonal matrix equation
!-----------------------------------------------------------------------
subroutine TRIDIAG(Nvec,Nmax,a,b,c,r,x)

implicit none

integer, intent(in) :: &
  Nvec,              &! Vector length
  Nmax                ! Maximum vector length

real, intent(in) :: &
  a(Nmax),           &! Below-diagonal matrix elements
  b(Nmax),           &! Diagonal matrix elements
  c(Nmax),           &! Above-diagonal matrix elements
  r(Nmax)             ! Matrix equation rhs

real, intent(out) :: &
  x(Nmax)             ! Solution vector

integer :: n          ! Loop counter 

! Work space   
real :: beta, g(Nvec) 

beta = b(1)
x(1) = r(1) / beta

do n = 2, Nvec
  g(n) = c(n-1) / beta
  beta = b(n) - a(n)*g(n)
  x(n) = (r(n) - a(n)*x(n-1)) / beta
end do

do n = Nvec - 1, 1, -1
  x(n) = x(n) - g(n+1)*x(n+1)
end do

end subroutine TRIDIAG
