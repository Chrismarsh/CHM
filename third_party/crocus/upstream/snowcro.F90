!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!SFX_LIC for details. version 1.
!     ##########################################################################
      SUBROUTINE SNOWCRO(HSNOWRES, TPTIME, OMEB, HIMPLICIT_WIND,      &
                      PPEW_A_COEF, PPEW_B_COEF,                                 &
                      PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,       &
                      PSNOWSWE,PSNOWRHO,PSNOWHEAT,PSNOWALB,                     &
                      PSNOWDIAMOPT,PSNOWSPHERI,PSNOWHIST,PSNOWAGE, PSNOWIMPUR,  &
                      PTSTEP,PPS,PSR,PUNLOAD, PRR,PPSN3L,PRESA_SN,              &
                      PTA,PTG,PSW_RAD,PQA,PVMOD,PWIND_DRIFT,PLW_RAD, PRHOA,     &
                      PUREF,PEXNS,PEXNA,PDIRCOSZW, PSLOPEDIR,                   &
                      PZREF,PZ0,PZ0EFF,PZ0H,PALB,                               &
                      PSOILCOND,PD_G,                                           &
                      PSNOWLIQ,PSNOWTEMP,PSNOWDZ,                               &
                      PTHRUFAL,PGRNDFLUX,PEVAPCOR, PGFLXCOR,                    &
                      PSWNETSNOW,PSWNETSNOWS,PLWNETSNOW,                        &
                      PRNSNOW,PHSNOW,PGFLUXSNOW,                                &
                      PHPSNOW,PLES3L,PLEL3L,PEVAP,PSNDRIFT,PRI,                 &
                      PEMISNOW,PCDSNOW,PUSTAR,PCHSNOW,PSNOWHMASS,PQS,           &
                      PPERMSNOWFRAC,PZENITH,PAZIM,PXLAT,PXLON,PBLOWSNW,         &
                      HSNOWDRIFT,HSNOWFPAPPUS,OSNOWDRIFT_SUBLIM,OSNOW_ABS_ZENITH,&
                      HSNOWMETAMO,HSNOWRAD,OATMORAD,P_DIR_SW, P_SCA_SW,         &
                      PSPEC_ALB, PDIFF_RATIO,PSPEC_TOT,PSNOWFLUX,PIMPWET,PIMPDRY,&
                      HSNOWFALL, HSNOWCOND,HSNOWHOLD,HSNOWCOMP,HSNOWZREF,       &
                      PSNOWMAK, OSNOWCOMPACT_BOOL, OSNOWMAK_BOOL, OSNOWTILLER,  &
                      OSELF_PROD, OSNOWMAK_PROP, PVFRIC_T, PHVEGPOL,            &
                      OSNOWAGING_VAR,PAGINGCOEF, OFOREST, PSFCMELT)
!     ##########################################################################
!
!
!!***** *SNOWCRO*
!!
!!    PURPOSE
!!    -------
!
!     Detailed snowpack scheme Crocus, computationnally based on the
!     3-Layer snow scheme option (Boone and Etchevers 1999)!
!!
!!**  METHOD
!!    ------
!
!     Direct calculation
!
!!    EXTERNAL
!!    --------
!
!     None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCES
!!    ----------
!!
!!    Crocus : Brun et al., 1989 (https://doi.org/10.3189/S0022143000009254)
!!    Crocus : Brun et al., 1992 (https://doi.org/10.3189/S0022143000009552)
!!    Crocus : Vionnet et al., 2012 (https://doi.org/10.5194/gmd-5-773-2012)
!!    New metamorphism : Carmagnola et al., 2014 (https://doi.org/10.5194/tc-8-417-2014)
!!    ESCROC : Lafaysse et al., 2017 (https://doi.org/10.5194/tc-11-1173-2017)
!!
!!
!!    AUTHORS
!!    -------
!!    Initial authors (2011)
!!      A. Boone    * CNRM *
!!      V. Vionnet  * CNRM *
!!      E. Brun     * CNRM *
!!    Main others contributors (global maintenance)
!!      M. Lafaysse * CNRM *  2012-2024
!!      R. Nheili   * CNRM *  2018-2022
!!      M. Fructus  * CNRM *  2018-2024
!!    Other inputs
!!      S. Morin    * CNRM *
!!      M. Dumont   * CNRM  
!!      P. Spandre  * CNRM *
!!      C. Carmagnola * CNRM*
!!      B. Cluzet   * CNRM *
!!      F. Tuzet    * CNRM *
!!      M. Baron    * CNRM *
!!      A. Haddjeri * CNRM *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    7/99
!!      Modified by A.Boone 05/02 (code, not physics)
!!      Modified by A.Boone 11/04 i) maximum density limit imposed (although
!!                                rarely if ever reached), ii) check to
!!                                see if upermost layer completely sublimates
!!                                during a timestep (as snowpack becomes vanishly
!!                                thin), iii) impose maximum grain size limit
!!                                in radiation transmission computation.
!!
!!      Modified by B. Decharme  (03/2009): Consistency with Arpege permanent
!!                                          snow/ice treatment (LGLACIER for alb)
!!      Modified by A. Boone     (04/2010): Implicit coupling and replace Qsat and DQsat
!!                                          by Qsati and DQsati, respectively.
!!      Modified by E. Brun, V. Vionnet, S. Morin (05/2011):
!!                                          Addition of Crocus processes and
!!                                          parametrizations to
!!                                          the SNOW-3L code. This includes the dynamic handling
!!                                          of snow layers and the inclusion of snow metamorphism
!!                                          rules similar to the original Crocus implementation.
!!      Modified by B. Decharme  (09/2012): New wind implicitation
!!
!!      Modified by M. Lafaysse (07/2012) :
!!                                          * Albedo and roughness parametrizations
!!                                            for surface ice over glaciers
!!                                                     MODIF 2012-10-03 : don't modify roughness if implicit coupling
!!                                                                 (test PPEW_A_COEF == 0. )
!!                                          * SNOWCROALB is now called by SNOWCRORAD to remove duplicated code
!!                                          * Parameters for albedo are moved to modd_snow_par
!!                                          * PSNOWAGE is stored as an age
!!                                            (days since snowfall) and not as a date
!!                                            to allow spinup simulations
!!                                          * New rules for optimal discretization of very thick snowpacks
!!                                          * Optional outputs for debugging
!!
!!       Modified by E. Brun and M. Lafaysse (07/2012) :
!!                                          * Implement sublimation in SNOWDRIFT
!!                                          * Flag in namelist to activate SNOWDRIFT and SNOWDRIFT_SUBLIM
!!       Modified by E. Brun and M. Lafaysse (08/2012) :
!!                                          * XUEPSI replaced by 0 in the if statement of case 1.3.3.2 (SNOWCROMETAMO)
!!                                          * If SNOWDRIFT is activated the wind do not modify grain types during snowfall
!!                                            (redundant with snowdrift)
!!       Modified by E. Brun (24/09/2012) :
!!                                          * Correction coupling coefficient for specific humidity in SNOWCROEBUD
!!                                          * PSFCFRZ(:)  = 1.0 for systematic solid/vapor latent fluxes in SNOWCROEBUD
!!       Modified by C. Carmagnola (3/2013):
!!                                          * Dendricity and size replaced by the optical diameter
!!                                          * Test of different evolution laws for the optical diameter
!!
!!       Modified by B. Decharme  (08/2013): Qsat as argument (needed for coupling with atm)
!!                                           add PSNDRIFT
!!
!!       Modified by M. Lafaysse (08/2015): MEB-Crocus coupling
!!       Modified by M. Dumont (11/2015) : atmotartes and spectral outputs
!!       Modified by F. Tuzet (06/2016): Add of a new dimension for impurity: The type of impurity
!!                                       Add of Impurity scavenging during melt
!!       Modified by M. Fructus (09/2018): Snow_sytron can be used with optical diameter
!!       Modified by M. Baron (2021): Fix metamorphism of Carmagnola et al. 2014
!!       Modified by M. Lafaysse (2021): Re-implement ice formation from Quéno et al., 2018
!!       Modified by M. Baron and A. Haddjeri (2024): SnowPappus blowing snow inputs
!!       Modified by M. Lafaysse (09/2024):
!!                                          * Improve MEB-Crocus coupling (snow layering should not change between MEB and Crocus) :
!!                                            call metamorphism, compaction, unloading, freezing rain etc. in that case after energy balance
!!                                            note that this is much more complex and expensive and still partly inaccurate in case of interception
!!                                            but the only way to obtain stable runs of MEB-Crocus at large scale.
!!                                          * Separate vegetation unloading and ice formation due to freezing rain for better stability
!!                                          * Separate in mode_snowcro module several routines for fresh snow / discretization / albedo
!!                                            necessary in MEB preparation
!!                                          * Remove LGLACIER option
!!                                          * Remove numerical inaccuracy in the different computations of INLVLS_OLD
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
!
USE MODD_CSTS, ONLY : XTT, XRHOLW, XLMTT,XLSTT,XLVTT, XCL, XCI, XPI, XRHOLI
USE MODD_SNOW_PAR, ONLY : XZ0ICEZ0SNOW, XRHOTHRESHOLD_ICE, XPERCENTAGEPORE,             &
                          XPERCENTAGEPORE_FRZ, XPERCENTAGEPORE_ICE, XIMPUR_EFOLD,       &
                          XIMPUR_DRY,XIMPUR_WET,XRHO_SNOWMAK, XCROCOEF_FF
USE MODD_SNOW_METAMO
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_PREP_SNOW, ONLY : NIMPUR
USE MODD_CONST_TARTES, ONLY:  XPSNOWG0, XPSNOWY0, XPSNOWW0, XPSNOWB0,NPNBANDS
USE MODD_CONST_ATM, ONLY: JPNBANDS_ATM
!
USE MODE_SNOW3L
USE MODE_SNOWCRO, ONLY : SNOWNLFALL_UPGRID, SNOWNLGRIDFRESH_1D, SNOWCROALB, &
                         SNOWCROUNLOAD, SNOWCROFREEZINGRAIN
USE MODE_TARTES, ONLY : SNOWCRO_TARTES, SURFACE_IMPURITY_REPARTITION
!
USE MODE_THERMOS
!
USE MODE_CRODEBUG
!
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE PARKIND1  ,ONLY : JPRB
#ifdef SFX_OL
USE MODN_IO_OFFLINE,  ONLY : LFORCIMP
#endif
!
!
IMPLICIT NONE
!
!
!*      0.1    declarations of arguments
!
! CHM diagnostic extension: surface ice melted per unit area and second.
REAL, DIMENSION(:), OPTIONAL, INTENT(OUT) :: PSFCMELT
REAL, INTENT(IN)                       :: PTSTEP
!                                      PTSTEP    = time step of the integration
TYPE(DATE_TIME), INTENT(IN)            :: TPTIME      ! current date and time
!
 CHARACTER(LEN=*), INTENT(IN)          :: HSNOWRES
!                                      HSNOWRES  = ISBA-SNOW3L turbulant exchange option
!                                      'DEF' = Default: Louis (ISBA: Noilhan and Mahfouf 1996)
!                                      'RIL' = Limit Richarson number under very stable
!                                              conditions (currently testing)
!                                      'M98'  = Martin et Lejeune 1998 : older computation for turbulent fluxes coefficents in Crocus
LOGICAL, INTENT(IN)                    :: OMEB       ! True = coupled to MEB. This means surface fluxes ae IMPOSED
!                                                    ! as an upper boundary condition to the explicit snow schemes.
!                                                    ! If = False, then energy
!                                                    ! budget and fluxes are computed herein.
!
 CHARACTER(LEN=*), INTENT(IN)          :: HIMPLICIT_WIND   ! wind implicitation option
!                                                          ! 'OLD' = direct
!                                                          ! 'NEW' = Taylor serie, order 1
!
REAL, DIMENSION(:), INTENT(IN)         :: PPS, PTA, PSW_RAD, PQA, PVMOD, PWIND_DRIFT, PLW_RAD, PSR, PUNLOAD, PRR
!                                      PSW_RAD = incoming solar radiation (W/m2)
!                                      PLW_RAD = atmospheric infrared radiation (W/m2)
!                                      PRR     = rain rate [kg/(m2 s)]
!                                      PSR     = snow rate (SWE) [kg/(m2 s)]
!                                      PTA     = atmospheric temperature at level za (K)
!                                      PVMOD   = modulus of the wind parallel to the orography (m/s)
!                                      PWIND_DRIFT   = modulus of the wind under canopy if PVMOD is above canopy, otherwise PWIND_DRIFT = PVMOD (m/s)
!                                      PPS     = surface pressure
!                                      PQA     = atmospheric specific humidity
!                                                at level za
REAL, DIMENSION(:), INTENT(INOUT)         :: PRESA_SN ! Aerodynamic resistance computed externally
!
REAL, DIMENSION(:,:), INTENT(IN)       :: P_DIR_SW, P_SCA_SW ! direct and diffuse spectral irradiance (W/m2/um)
!
REAL, DIMENSION(:,:), INTENT(IN)       :: PIMPWET, PIMPDRY  !Dry and wet deposit coefficient from Forcing File(g/m²/s)
!
REAL, DIMENSION(:), INTENT(IN)       :: PHVEGPOL ! Mean polar vegetation height used to reduce snow compaction and snow drift 
                                                 ! below veg height Only used with options from Royer et al. 2021 
!
REAL, DIMENSION(:), INTENT(IN)         :: PTG, PSOILCOND, PD_G, PPSN3L
!                                      PTG       = Surface soil temperature (effective
!                                                  temperature the of layer lying below snow)
!                                      PSOILCOND = soil thermal conductivity [W/(m K)]
!                                      PD_G      = Assumed first soil layer thickness (m)
!                                                  Used to calculate ground/snow heat flux
!                                      PPSN3L    = snow fraction
!
REAL, DIMENSION(:), INTENT(IN)         :: PZREF, PUREF, PEXNS, PEXNA, PDIRCOSZW, PSLOPEDIR, &
                                          PRHOA, PZ0, PZ0EFF, PALB, PZ0H, PPERMSNOWFRAC, PVFRIC_T
!                                      PZ0EFF        = roughness length for momentum
!                                      PZ0           = grid box average roughness length
!                                      PZ0H          = grid box average roughness length for heat
!                                      PZREF         = reference height of the first
!                                                      atmospheric level
!                                      PUREF         = reference height of the wind
!                                      PRHOA         = air density
!                                      PEXNS         = Exner function at surface
!                                      PEXNA         = Exner function at lowest atmos level
!                                      PDIRCOSZW     = Cosinus of the angle between the
!                                                      normal to the surface and the vertical
!                                      PSLOPEDIR     = Slope direction
!                                      PALB          = soil/vegetation albedo
!                                      PPERMSNOWFRAC = fraction of permanet snow/ice
!                                      PVFRIC_T      = threshold friction velocity for wind induced snow transport (SnowPappus)
!
REAL, DIMENSION(:), INTENT(IN)         :: PPEW_A_COEF, PPEW_B_COEF,                   &
                                          PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF,      &
                                          PPEQ_B_COEF
!                                      PPEW_A_COEF = wind coefficient (m2s/kg)
!                                      PPEW_B_COEF = wind coefficient (m/s)
!                                      PPET_A_COEF = A-air temperature coefficient
!                                      PPET_B_COEF = B-air temperature coefficient
!                                      PPEQ_A_COEF = A-air specific humidity coefficient
!                                      PPEQ_B_COEF = B-air specific humidity coefficient
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PSNOWALB
!                                    PSNOWALB = Prognostic surface snow albedo
!                                               (does not include anything but
!                                               the actual snow cover)
!
REAL, DIMENSION(:,:), INTENT(INOUT)    :: PSNOWHEAT, PSNOWRHO, PSNOWSWE
!                                      PSNOWHEAT = Snow layer(s) heat content (J/m2)
!                                      PSNOWRHO  = Snow layer(s) averaged density (kg/m3)
!                                      PSNOWSWE  = Snow layer(s) Water Equivalent (SWE:kg m-2)
!
REAL, DIMENSION(:,:), INTENT(INOUT)  :: PSNOWDIAMOPT, PSNOWSPHERI, PSNOWHIST
!                                      PSNOWDIAMOPT = Snow layers grain feature 1
!                                      PSNOWSPHERI = Snow layer grain feature 2
!                                      PSNOWHIST  = Snow layer grain historical
!                                                   parameter (only for non
!                                                   dendritic snow)
REAL, DIMENSION(:,:), INTENT(INOUT)    :: PSNOWAGE  ! Snow grain age
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PSNOWIMPUR  ! Snow impurity content (g) (LOCATION,LAYER,NIMPUR)) Impur type :1/BC 2/Dust
!
REAL, DIMENSION(:,:), INTENT(OUT)      :: PSNOWTEMP
!                                      PSNOWTEMP = Snow layer(s) temperature (m)

REAL, DIMENSION(:,:), INTENT(OUT)      :: PSNOWLIQ, PSNOWDZ
!                                      PSNOWLIQ  = Snow layer(s) liquid water content (m)
!                                      PSNOWDZ   = Snow layer(s) thickness (m)
!
REAL, DIMENSION(:), INTENT(OUT)        :: PTHRUFAL,  PEVAPCOR,  PGFLXCOR

!                                      PTHRUFAL  = rate that liquid water leaves snow pack:
!                                                  paritioned into soil infiltration/runoff
!                                                  by ISBA [kg/(m2 s)]
!
!                                      PEVAPCOR  = evaporation/sublimation correction term:
!                                                  extract any evaporation exceeding the
!                                                  actual snow cover (as snow vanishes)
!                                                  and apply it as a surface soil water
!                                                  sink. [kg/(m2 s)]
!                                      PGFLXCOR  = flux correction to underlying soil for vanishing snowpack
!                                                  (to put any energy excess from snow to soil) (W/m2)
!

REAL, DIMENSION(:,:), INTENT(OUT)      :: PSPEC_ALB, PDIFF_RATIO,PSPEC_TOT  !! spectral albedo, diffuse to total irradiance ratio and total incoming spectral irradiance (npoints,nbands)

REAL, DIMENSION(:), INTENT(OUT)        :: PSNOWFLUX
!                                      PSNOWFLUX = heat flux between the surface and sub-surface
!                                                  snow layers (W/m2))

! Diagnostics : on verra plus tard si c'est nécessaire
! REAL, DIMENSION(:), INTENT(OUT)      :: PDELHEATN, PDELHEATN_SFC
!                                      PDELHEATN = total snow heat content change in the surface layer (W m-2)
!                                      PDELHEATN_SFC = total snow heat content change during the timestep (W m-2)

!
REAL, DIMENSION(:), INTENT(INOUT)      :: PGRNDFLUX
!                                      PGRNDFLUX = soil/snow interface heat flux (W/m2)

REAL, DIMENSION(:), INTENT(INOUT)      :: PRNSNOW, PHSNOW, PGFLUXSNOW, PLES3L, PLEL3L, &
                                          PHPSNOW, PCDSNOW, PUSTAR, PEVAP
!                                      PLES3L      = sublimation heat flux from snow (W/m2)
!                                      PLEL3L      = liquid water evaporation flux from snow (always 0 W/m2 in Crocus)
!                                      PHPSNOW     = heat release from rainfall (W/m2)
!                                      PRNSNOW     = net radiative flux from snow (W/m2)
!                                      PHSNOW      = sensible heat flux from snow (W/m2)
!                                      PGFLUXSNOW  = net heat flux from snow (W/m2)
!                                      PCDSNOW     = drag coefficient for momentum over snow
!                                      PUSTAR      = friction velocity over snow (m/s)
!                                      PEVAP       = total evaporative flux (kg/m2/s)

REAL, DIMENSION(:), INTENT(OUT)        :: PSNDRIFT
!                                      PSNDRIFT    = blowing snow sublimation (kg/m2/s)
!
REAL, DIMENSION(:), INTENT(INOUT)      :: PSWNETSNOW, PLWNETSNOW, PSWNETSNOWS
!                                      PSWNETSNOW = net shortwave radiation entering top of snowpack
!                                                  (W m-2) Imposed if MEB=T, diagnosed herein if MEB=F
!                                      PSWNETSNOWS= net shortwave radiation in uppermost layer of snowpack
!                                                  (W m-2) Imposed if MEB=T, diagnosed herein if MEB=F
!                                                   Used for surface energy budget diagnostics
!                                      PLWNETSNOW = net longwave radiation entering top of snowpack
!                                                  (W m-2) Imposed if MEB=T, diagnosed herein if MEB=F
!
REAL, DIMENSION(:), INTENT(INOUT)      :: PCHSNOW,PRI
!                                      PCHSNOW     = drag coefficient for heat over snow
!                                      PRI = Ridcharson number

REAL, DIMENSION(:), INTENT(OUT)        :: PEMISNOW, PSNOWHMASS
!                                      PEMISNOW    = snow surface emissivity
!                                      PSNOWHMASS  = heat content change due to mass
!                                                    changes in snowpack (J/m2): for budget
!                                                    calculations only.
!
REAL, DIMENSION(:), INTENT(OUT)        :: PQS
!                                      PQS = surface humidity
!
REAL, DIMENSION(:), INTENT(IN)         :: PZENITH ! solar zenith angle
REAL, DIMENSION(:), INTENT(IN)         :: PAZIM   ! azimuthal angle      (radian from North, clockwise)
REAL, DIMENSION(:), INTENT(IN)         :: PXLAT,PXLON ! LAT/LON after packing
!
REAL, DIMENSION(:,:), INTENT(INOUT)    :: PBLOWSNW !  Properties of deposited blowing snow (from Sytron or Meso-NH/Crocus)
                                       !    1 : Deposition flux (kg/m2/s)
                                       !    2 : Density of deposited snow (kg/m3)
                                       !    3 : SGRA1 of deposited snow
                                       !    4 : SGRA2 of deposited snow
!
 CHARACTER(4), INTENT(IN)              :: HSNOWDRIFT        ! Snowdrift scheme :
                                       ! Mechanical transformation of snow grain and compaction + effect of wind
                                       ! on falling snow properties
                                       !    'NONE': No snowdrift scheme
                                       !    'DFLT': falling snow falls as purely dendritic
                                       !    'GA01': Gallee et al 2001
                                       !    'VI13': Vionnet et al 2013
                                       !    'PAPP': snowdrift scheme coupled with SnowPappus threshold wind speed
                                       !    'R21F': Royer et al 2021 (Full effects: Increase in Maximum Density and Wind Effect)
                                       !    'R21W': Royer et al 2021 (Increase in Wind_Effect)
                                       !    'R21R': Royer et al 2021 (Increase in Maximum Density)                                       
CHARACTER(4), INTENT(IN)               :: HSNOWFPAPPUS
                                       ! Option to force fresh snow characteristics when snowpappus is used
                                       ! 'GM98' => forces fresh snow characteristics as if HSNOWDRIFT = 'NONE'
                                       ! 'VI13' => """" 'VI13'
                                       ! 'NONE' => no influence ( it is forced when snowpappus is not activated )
!
LOGICAL, INTENT(IN)                    :: OSNOWDRIFT_SUBLIM ! activate sublimation during drift
REAL, DIMENSION (:), INTENT(IN)        :: PSNOWMAK        ! Snowmaking thickness (m)
LOGICAL, INTENT(IN)                    :: OSNOWCOMPACT_BOOL, OSNOWMAK_BOOL, OSNOWTILLER, &
                                          OSELF_PROD, OSNOWMAK_PROP

LOGICAL, INTENT(IN)                    :: OSNOWAGING_VAR ! activate option to read spatially-variable snow aging coefficient for albedo calculation 
REAL, DIMENSION(:), INTENT(IN)         :: PAGINGCOEF ! Snow aging coefficient for albedo param. (only used in OSNOWAGING_VAR = TRUE)
!
!
LOGICAL, INTENT(IN)                    :: OSNOW_ABS_ZENITH ! activate parametrization of solar absorption for polar regions
 CHARACTER(3), INTENT(IN)              :: HSNOWMETAMO, HSNOWRAD, HSNOWFALL, HSNOWCOND, HSNOWHOLD, HSNOWCOMP, HSNOWZREF
LOGICAL, INTENT(IN)                    :: OATMORAD ! activate atmotartes scheme
                                       !-----------------------
                                       ! Metamorphism scheme
                                       ! HSNOWMETAMO=B21 Baron, 2023
                                       ! HSNOWMETAMO=T07 Taillandier et al 2007
                                       ! HSNOWMETAMO=F06 Flanner et al 2006
                                       ! HSNOWMETAMO=S-F Schlef et al 2014
                                       ! HSNOWMETAMO=S-B Schlef et al 2014                                       
                                       !-----------------------
                                       ! Radiative transfer scheme
                                       ! HSNOWRAD=B92 Brun et al 1992
                                       ! HSNOWRAD=T17 (Tuzet et al. 2017) (Libois et al. 2013) TARTES with impurities content scheme
                                       !-----------------------
                                       ! New options for multiphysics version (Cluzet et al 2016)
                                       ! Falling snow scheme
                                       ! HSNOWFALL=V12 Vionnet et al. 2012 from Brun et al. 1989
                                       ! HSNOWFALL=A76 Anderson et al. 1976
                                       ! HSNOWFALL=S02 Lehning el al. 2002
                                       ! HSNOWFALL=P75 Pahaut 1975
                                       ! HSNOWFALL=NZE Constant density 200 kg/m3 (who knows ?)
                                       ! HSNOWFALL=R21 Royer et al. 2021
                                       ! HSNOWFALL=L22 Lackner et al. 2022
                                       ! HSNOWFALL=GW1 
                                       ! HSNOWFALL=GW2                                          
                                       !---------------------
                                       ! Thermal conductivity scheme
                                       ! HSNOWCOND=Y81 default Crocus from Yen et al. 1981
                                       ! HSNOWCOND=I02 ISBA_ES snow conductivity parametrization (Boone et al. 2002)
                                       ! HSNOWCOND=C11 Calonne et al. 2011 snow conductivity parametrization
                                       !-----------------------
                                       ! liquid water content scheme
                                       ! HSNOWHOLD=B92 default Crocus from Brun et al. 1992 or Vionnet et al. 2012
                                       ! HSNOWHOLD=B02 ISBA_ES  parametrization (Boone et al. 2002)
                                       ! HSNOWHOLD=O04 CLM parametrization (Oleson et al 2004)
                                       ! HSNOWHOLD=SPK SNOWPACK aprametrization (Lehning et al 2002)
                                       !-----------------------
                                       ! reference height is constant or variable from the snow surface
                                       ! HSNOWZREF='CST' constant reference height from the snow surface
                                       ! HSNOWZREF='VAR' variable reference height from the snow surface (i.e. constant from the ground)
                                       !-----------------------
!
LOGICAL, INTENT(IN)                    :: OFOREST ! True if Crocus is called in a forested environment. 
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2))        :: ZSNOWSSA_BEFORE, ZSNOWSSA_AFTER,ZSNOWDSSA
REAL, DIMENSION(SIZE(PSNOWRHO,1))                         :: ZIMPUR_NORM !impurities content (g) (npoints,nlayer,ntypes_impurities)
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2))        :: ZSNOWTEMP, ZSCAP, ZSNOWDZO, ZSNOWDZN, ZSCOND, ZRADSINK
!                                                         ZSNOWTEMP  = Snow layer(s) averaged temperature (K)
!                                                         ZSCAP      = Snow layer(s) heat capacity [J/(K m3)]
!                                                         ZSNOWDZN   = Updated snow layer thicknesses (m)
!                                                         ZSCOND     = Snow layer(s) thermal conducivity [W/(m K)]
!                                                         ZRADSINK   = Snow solar Radiation source terms (W/m2)
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2))        :: ZDRYDENSITY ! Density of solid phase kg/m3
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2))        :: ZWHOLDMAX
REAL, DIMENSION(SIZE(PSNOWRHO,1),JPNBANDS_ATM)            :: ZSNOWALB_SP
REAL, DIMENSION(SIZE(PSNOWRHO,1),JPNBANDS_ATM)            :: ZSPEC_DIR, ZSPEC_DIF
!
!For now these values are constant
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2))        :: ZSNOWG0 ! asymmetry parameter of snow grains at nr=1.3 and at non absorbing wavelengths (no unit) (npoints,nlayer)
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2))        :: ZSNOWY0 ! Value of y of snow grains at nr=1.3 (no unit
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2))        :: ZSNOWW0 ! Value of W of snow grains at nr=1.3 (no unit)
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2))        :: ZSNOWB0 ! absorption enhancement parameter of snow grains at nr=1.3 and at non absorbing wavelengths (no unit)
!
!spectral albedo (3 bands for now) :: ready to output if necessary
REAL, DIMENSION(SIZE(PSNOWRHO,1),3) :: ZSPECTRALALBEDO
!
REAL, DIMENSION(SIZE(PTA))          :: ZSNOWBIS
!                                      ZSNOWBIS      = Total snow depth after snowfall
!
REAL, DIMENSION(SIZE(PTA))          :: ZSNOW, ZSFCFRZ, ZTSTERM1, ZTSTERM2, ZCT, ZRA, ZSNOWTEMPO1
!                                      ZSNOW      = Total snow depth (m)
!                                      ZCT        = inverse of the product of snow heat capacity
!                                                   and layer thickness [(m2 K)/J]
!                                      ZRA        = Surface aerodynamic resistance
!                                      ZTSTERM1,ZTSTERM2 = Surface energy budget coefficients
!
!                                      ZSNOWTEMPO1= value of uppermost snow temperature
!                                                   before time integration (K)
!
REAL, DIMENSION(SIZE(PTA))          :: ZRSRA, ZDQSAT, ZQSAT, ZRADXS, ZLIQHEATXS, ZGRNDFLUXI, ZPSN3L
!                                      ZRSRA    = air density over aerodynamic resistance
!                                      ZDQSAT   = derrivative of saturation specific humidity
!                                      ZQSAT    = saturation specific humidity
!                                      ZRADXS   = shortwave radiation absorbed by soil surface
!                                                 (for thin snow sover) (W m-2)
!                                      ZLIQHEATXS = excess snowpack heating for vanishingly thin
!                                                 snow cover: add energy to snow/ground heat
!                                                 flux (W m-2)
!                                      ZGRNDFLUXI= for the case where the ground flux is imposed,
!                                                  this is the actual imposed value.
!                                      ZPSN3L    = snow fraction: different use if MEB "on".
!                                                  In this case, it is only used for Tg update
!                                                  since only this variable has a sub-grid relevance.
!
REAL, DIMENSION(SIZE(PTA))          :: ZUSTAR2_IC, ZTA_IC, ZQA_IC, &
                                       ZPET_A_COEF_T, ZPEQ_A_COEF_T, ZPET_B_COEF_T, ZPEQ_B_COEF_T
!                                      ZUSTAR2_IC    = implicit lowest atmospheric level friction (m2/s2)
!                                      ZTA_IC        = implicit lowest atmospheric level air temperature
!                                      ZQA_IC        = implicit lowest atmospheric level specific humidity
!                                      ZPET_A_COEF_T = transformed A-air temperature coefficient
!                                      ZPET_B_COEF_T = transformed B-air temperature coefficient
!                                      ZPEQ_A_COEF_T = transformed A-air specific humidity coefficient
!                                      ZPEQ_B_COEF_T = transformed B-air specific humidity coefficient
!
REAL, DIMENSION(SIZE(PTA))          :: ZSR ! snow rate with minor correction (include unloading on bare ground)
REAL, DIMENSION(SIZE(PTA))          :: ZUNLOAD ! unloading rate with minor correction (0 on bare ground)
REAL, DIMENSION(SIZE(PTA))          :: ZSNOWRHOF, ZSNOWDZF, ZSNOWDIAMOPTF, ZSNOWSPHERIF, ZSNOWHISTF
REAL, DIMENSION(SIZE(PTA))          :: ZSNOWAGEF
REAL, DIMENSION(SIZE(PTA),NIMPUR)   :: ZSNOWIMPURF

REAL, DIMENSION(SIZE(PTA),NIMPUR)   :: ZDRYCOEF ! Dry deposit coefficient for each kind of impurity (g/m²)
REAL, DIMENSION(SIZE(PTA),NIMPUR)   :: ZWETCOEF ! Wet deposit coefficient for each kind of impurity (g/m²)


! New roughness lengths in case of glaciers without snow.
REAL, DIMENSION(SIZE(PTA))          :: ZZ0_SNOWICE, ZZ0H_SNOWICE, ZZ0EFF_SNOWICE
!
! Reference heights for temperature and wind can be modified depending on snow depth when HSNOWZREF=="VAR".
REAL, DIMENSION(SIZE(PTA))          :: ZZREF,ZUREF
!To control and print eneregy balance
REAL , DIMENSION(SIZE(PTA))         :: ZSUMMASS_INI,ZSUMHEAT_INI,ZSUMMASS_FIN,ZSUMHEAT_FIN
!
REAL, DIMENSION(SIZE(PTA))          :: ZMASSBALANCE, ZENERGYBALANCE, ZEVAPCOR2
!
INTEGER, DIMENSION(SIZE(PTA))       :: INLVLS_USE, INLVLS_USE_OLD ! varying number of effective layers
!
LOGICAL, DIMENSION(SIZE(PTA))       :: GSNOWFALL,GMODIF_MAILLAGE
!                                      GSNOWFALL  = FLAG if snowfall exceed PSNOW/10, used for
!                                                   grid updating.
LOGICAL, DIMENSION(SIZE(PTA))       :: GFRZRAIN ! flag for freezing rain
REAL, DIMENSION(SIZE(PTA))          :: ZSNOWDZ_FRZ ! Depth of freezing layer
REAL, DIMENSION(SIZE(PTA))          :: ZWORK,ZWORK2
!
REAL                                :: ZTSTEPDAYS ! time step in days
!
LOGICAL                             :: GCROINFOPRINT ! print daily informations
LOGICAL                             :: GCRODEBUGPRINT, GCRODEBUGDETAILSPRINT, GCRODEBUGPRINTATM ! print diagnostics for debugging
LOGICAL                             :: GCRODEBUGPRINTBALANCE
LOGICAL                             :: GSUCCESS ! flag to test the success of regridding
!
INTEGER                             :: JJ,JST, JP,JIMP  ! looping indexes
INTEGER                             :: IPRINT  ! gridpoint number to be printed
INTEGER                             :: IDEBUG
!
REAL, DIMENSION(SIZE(PTA),SIZE(PSNOWRHO,2))     :: ZWORK2D
INTEGER :: IMAX_USE ! maximum number of layers over the domain
REAL, DIMENSION(SIZE(PTA)) :: ZICE_BEFORE_GONE, ZSURFACE_LIQ
REAL(KIND=JPRB)                     :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO',0,ZHOOK_HANDLE)
!
IF (LCRODAILYINFO .OR. LCRODEBUG) THEN
  !***************************************PRINT IN**********************************************
  ! Look if we have to print snow profiles for debugging
  GCROINFOPRINT = LCRODAILYINFO .AND. (TPTIME%TIME ==0.0)
  !***************************************PRINT OUT*********************************************
  !***************************************DEBUG IN**********************************************
  GCRODEBUGPRINTBALANCE = ( TPTIME%TDATE%YEAR*10000 + TPTIME%TDATE%MONTH*100 + TPTIME%TDATE%DAY &
                          >= NTIMECRODEBUG ) .AND. &
                          ( TPTIME%TIME/3600. >= NHOURCRODEBUG ) .AND. &
                          ( TPTIME%TDATE%YEAR*10000 + TPTIME%TDATE%MONTH*100 + TPTIME%TDATE%DAY &
                          < NENDCRODEBUG )
  !
  IF (LCRODEBUG) THEN
    GCRODEBUGPRINT        = GCRODEBUGPRINTBALANCE
    GCRODEBUGDETAILSPRINT = LCRODEBUGDETAILS .AND. GCRODEBUGPRINT
    GCRODEBUGPRINTATM     = LCRODEBUGATM .AND. GCRODEBUGPRINT
  ELSE
    GCRODEBUGPRINT        = .FALSE.
    GCRODEBUGDETAILSPRINT = .FALSE.
    GCRODEBUGPRINTATM     = .FALSE.
  END IF

  !
  ! Look if we have to compute and print energy balance control
  GCRODEBUGPRINTBALANCE = LCONTROLBALANCE .AND. GCRODEBUGPRINTBALANCE
  !
  IF ( LCRODEBUG .OR. GCROINFOPRINT .OR. GCRODEBUGPRINTBALANCE ) THEN
    !
    IF ( XLATCRODEBUG >= -90 .AND. XLONCRODEBUG >= -180. ) THEN
      CALL GETPOINT_CRODEBUG(PXLAT,PXLON,IDEBUG)
    ELSE
      IDEBUG = NPOINTCRODEBUG
    END IF
    !
    IF (SIZE(PTA) < IDEBUG) THEN
      GCRODEBUGPRINT        = .FALSE.
      GCRODEBUGDETAILSPRINT = .FALSE.
      GCRODEBUGPRINTATM     = .FALSE.
    END IF
    ! For parallel runs : we just want to do this for the thread where there is this point.
    IF ( XLATCRODEBUG >= -90 .AND. XLONCRODEBUG >= -180. ) THEN
      IF ( ABS( PXLAT(IDEBUG)-XLATCRODEBUG ) + ABS(PXLON(IDEBUG) - XLONCRODEBUG) > 0.001 ) THEN
        GCRODEBUGPRINT        = .FALSE.
        GCRODEBUGDETAILSPRINT = .FALSE.
        GCRODEBUGPRINTATM     = .FALSE.
      END IF
    END IF
    !
  END IF
ELSE
  GCROINFOPRINT = .FALSE.
  GCRODEBUGPRINT = .FALSE.
  GCRODEBUGDETAILSPRINT = .FALSE.
  GCRODEBUGPRINTATM = .FALSE.
  GCRODEBUGPRINTBALANCE = .FALSE.
END IF
!***************************************DEBUG OUT*********************************************
!
IF ( HSNOWRAD=="T17") THEN
  !For now fix constant values
  ZSNOWG0 = XPSNOWG0
  ZSNOWY0 = XPSNOWY0
  ZSNOWW0 = XPSNOWW0
  ZSNOWB0 = XPSNOWB0
  !
  !
END IF
!
ZUSTAR2_IC = 0.0
ZTA_IC     = 0.0
ZQA_IC     = 0.0
!
!
IF (.NOT. (OMEB) )THEN
  PGRNDFLUX  = 0.
  PHSNOW     = 0.
  PRNSNOW    = 0.
  PLES3L     = 0.
  PLEL3L     = 0.
  PHPSNOW    = 0.
  PGFLXCOR   = 0.
END IF

PSNDRIFT   = 0.
PSNOWHMASS = 0.
PEVAPCOR   = 0.
PTHRUFAL   = 0.
ZIMPUR_NORM(:)=0.
!
! pour imprimer des diagnostics sur un des points
IPRINT = 1
!
! - - ---------------------------------------------------
!
!       0.     Initialization
!               --------------
! NOTE that snow layer thickness is used throughout this code: SWE
! is only used to diagnose the thickness at the beginning of this routine
! and it is updated at the end of this routine.
!
! Initialization of the actual number of snow layers, total snow depth
!  and layer thicknesses
!
ZSNOWTEMP(:,:) = 0.
!
INLVLS_USE(:) = 0
DO JST = 1,SIZE(PSNOWSWE(:,:),2)
  DO JJ = 1,SIZE(ZSNOW)
    IF ( PSNOWSWE(JJ,JST)>0. ) THEN
      PSNOWDZ(JJ,JST) = PSNOWSWE(JJ,JST) / PSNOWRHO(JJ,JST)
      INLVLS_USE(JJ) = JST
    ELSE
      PSNOWDZ(JJ,JST) = 0.
    END IF
  END DO  !  end loop snow layers
END DO    ! end loop grid points
!
IMAX_USE = MAXVAL(INLVLS_USE)
! Incrementation of snow layers age
ZTSTEPDAYS = PTSTEP/86400. ! time step in days
! Lafaysse / Cluzet : reimplementation of first Morin/Charrois impuritites content option:
! part of code modified by S. Morin 22/08/2013 on impurities behavior
!
IF ( HSNOWRAD=="T17") THEN
#ifdef SFX_OL
  IF (LFORCIMP) THEN ! Les flux de dépots atmosphériques doivent être en g/m²/s en entrée, format ALADIN
    DO JIMP=1,NIMPUR
      DO JJ = 1,SIZE(ZSNOW)
        ZWETCOEF(JJ,JIMP)=PIMPWET(JJ,JIMP) * PTSTEP         !from g/m²/s to g/m² during the time step
        ZDRYCOEF(JJ,JIMP)=PIMPDRY(JJ,JIMP) * PTSTEP         !from g/m²/s to g/m² during the time step
      END DO
    END DO
  ELSE
    DO JIMP=1,NIMPUR
      DO JJ = 1,SIZE(ZSNOW)
        ZWETCOEF(JJ,JIMP)=XIMPUR_WET(JIMP) * PTSTEP      ! Value defined in the namelist (g m-2 s-1) converted to g m-2 during the time step
        ZDRYCOEF(JJ,JIMP)=XIMPUR_DRY(JIMP)  *PTSTEP    ! Value defined in the namelist (g m-2 s-1) converted to g m-2 during the time step
      END DO
    END DO
  END IF
#else
  DO JIMP=1,NIMPUR
    DO JJ = 1,SIZE(ZSNOW)
      ZWETCOEF(JJ,JIMP)=XIMPUR_WET(JIMP) * PTSTEP      ! Value defined in the namelist (g m-2 s-1) converted to g m-2 during the time step
      ZDRYCOEF(JJ,JIMP)=XIMPUR_DRY(JIMP) * PTSTEP     ! Value defined in the namelist (g m-2 s-1) converted to g m-2 during the time step
    END DO
  END DO
#endif
END IF
!
DO JST = 1,IMAX_USE
  DO JJ = 1,SIZE(PSNOWSWE,1)
    IF (JST <= INLVLS_USE(JJ)) THEN
      PSNOWAGE(JJ,JST)=PSNOWAGE(JJ,JST)+ZTSTEPDAYS  ! this is the classical version where snowage is a real age of snow layers
    END IF
  END DO
END DO

!
!***************************************PRINT IN**********************************************
!
IF (LCRODAILYINFO .OR. LCRODEBUG) THEN
  !Compute total SWE and heat for energy control
  IF ( GCRODEBUGPRINTBALANCE ) THEN
    DO JJ = 1,SIZE(ZSNOW)
      ZSUMMASS_INI(JJ) = SUM(PSNOWSWE (JJ,1:INLVLS_USE(JJ)))
      ZSUMHEAT_INI(JJ) = SUM(PSNOWHEAT(JJ,1:INLVLS_USE(JJ)))
    END DO    ! end loop grid points
  END IF
  !
  !  Print of some simulation characteristics
  IF(GCROINFOPRINT) THEN
    CALL SNOWCROPRINTDATE()
    WRITE(*,FMT="(A12,I3,A12,I4)") 'nlayer:',INLVLS_USE(IDEBUG), ' nbpoints:', SIZE(ZSNOW)
  !   WRITE(*,*) 'PZ0H: ', PZ0H(IDEBUG)
    WRITE(*,*) 'Snow fraction =',PPSN3L(IDEBUG)
  END IF
  !
  !***************************************PRINT OUT*********************************************
  !***************************************DEBUG IN**********************************************
  IF (GCRODEBUGPRINT) THEN
    CALL SNOWCROPRINTDATE()
    CALL SNOWCROPRINTPROFILE("crocus initialization",INLVLS_USE(IDEBUG),LPRINTGRAN,      &
                             PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),   &
                             PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:),&
                             PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:),&
                             HSNOWMETAMO)
  END IF
  !
  IF (GCRODEBUGPRINTATM) THEN
    CALL SNOWCROPRINTATM("forcing data :",PTA(IDEBUG),PQA(IDEBUG),PVMOD(IDEBUG),   &
                         PRR(IDEBUG),PSR(IDEBUG),PSW_RAD(IDEBUG),PLW_RAD(IDEBUG),  &
                         PTG(IDEBUG),PSOILCOND(IDEBUG),PD_G(IDEBUG),PPSN3L(IDEBUG),&
                         PBLOWSNW(IDEBUG,1) )
  END IF
END IF
!***************************************DEBUG OUT********************************************
!
!*       1.     Snow total depth
!               ----------------
!
ZSNOW(:) = 0.
DO JJ = 1,SIZE(ZSNOW)
  ZSNOW(JJ) = SUM(PSNOWDZ(JJ,1:INLVLS_USE(JJ)))
END DO
!
ZSNOWBIS(:) = ZSNOW(:)
!
! First estimate of ZUREF to estimate wind for snowfall et snowdrift routines
IF (HSNOWZREF=="VAR") THEN
  ZUREF(:)=MAX(PUREF(:)-ZSNOW(:),0.2)
ELSE
  ZUREF(:)=PUREF(:)
END IF
!*       2.     Snowfall
!               --------
! Calculate uppermost density and thickness changes due to snowfall,
! and add heat content of falling snow
!
IF (OMEB) THEN
  GFRZRAIN(:) = .FALSE.
ELSE
  DO JJ=1,SIZE(ZSNOW)
  ! VV GFRZRAIN(JJ) = (PTA(JJ) < XTT) .AND. (PRR(JJ)*PTSTEP > XUEPSI) ! Modif for single precision
  GFRZRAIN(JJ) = (PTA(JJ) < XTT) .AND. (PRR(JJ)*PTSTEP > XUEPSI_SMP*917.) !Cannot create ice layer thinner than  XUEPSI_SMP 
! In case of troubles with MEB, following restriction might temporarily help but is not a good option
!.AND. &
!  (ZSNOW(JJ) > 0.02) .AND. (PSR(JJ)*PTSTEP < XUEPSI) .AND. (PUNLOAD(JJ)*PTSTEP < XUEPSI)
  END DO
END IF
!
! Tricky case : if unloading on bare ground, convert unloading in new snow
! because otherwise it is not possible to compute any energy balance (no snow before and after)
! and because regridding routines do not consider this case
! Note that in that case the MEB solving has not seen any snow...
ZSR(:) = PSR(:)
ZUNLOAD(:) = PUNLOAD(:)
!VV IF (OMEB) THEN
IF (OMEB .OR. OFOREST) THEN
  DO JJ=1,SIZE(ZSNOW)
    IF ((INLVLS_USE(JJ) == 0) .AND. (ZUNLOAD(JJ) > 0.) .AND. (ZSR(JJ) < XUEPSI_SMP)) THEN
      ZSR(JJ) = ZSR(JJ)+ZUNLOAD(JJ)
      ZUNLOAD(JJ) = 0.
    END IF
  END DO
END IF
!
IF (GCRODEBUGDETAILSPRINT) THEN

  PRINT*, "natural snowfall: ", ZSR(IDEBUG)>XUEPSI, ZSR(IDEBUG)
  PRINT*, "unloading: ", ZUNLOAD(IDEBUG)> XUEPSI, ZUNLOAD(IDEBUG) 
  PRINT*, "blowing snow: ", PBLOWSNW(IDEBUG,1) > XUEPSI
  PRINT*, "machine-made snow: ", PSNOWMAK(IDEBUG)*XRHO_SNOWMAK/PTSTEP > XUEPSI
  PRINT*, "freezing rain / temperature: ", GFRZRAIN(IDEBUG), PRR(IDEBUG), PTA(IDEBUG)
END IF

! For SVS2 call snow unloading before SNOWFALL
IF(OFOREST) THEN
    CALL SNOWCROUNLOAD(PSNOWDZ, PSNOWRHO, PSNOWHEAT, PSNOWDIAMOPT, PSNOWSPHERI, &
                   PSNOWHIST, PSNOWAGE, PSNOWIMPUR, ZUNLOAD, PTA, PTSTEP, INLVLS_USE)
   ! Update ZSNOW
    DO JJ = 1,SIZE(ZSNOW)
       ZSNOW(JJ) = SUM(PSNOWDZ(JJ,1:INLVLS_USE(JJ)))
    END DO
    ZSNOWBIS(:) = ZSNOW(:)
END IF

!
! save previous number of layers
INLVLS_USE_OLD = INLVLS_USE
ZSNOWDZO(:,:) = PSNOWDZ(:,:)
!
 CALL SNOWNLFALL_UPGRID(PTSTEP,ZSR,PTA,PWIND_DRIFT,ZSNOWBIS,PSNOWRHO,     &
                        PSNOWDZ, PSNOWHEAT,PSNOWHMASS,PSNOWALB,PPERMSNOWFRAC,   &
                        PSNOWDIAMOPT,PSNOWSPHERI,PSNOWHIST,PSNOWAGE,GSNOWFALL,  &
                        ZSNOWDZN, ZSNOWRHOF, ZSNOWDZF,                &
                        ZSNOWDIAMOPTF, ZSNOWSPHERIF, ZSNOWHISTF,                &
                        ZSNOWAGEF,ZWETCOEF, ZSNOWIMPURF,GMODIF_MAILLAGE,INLVLS_USE,       &
                        HSNOWDRIFT,HSNOWFPAPPUS,PZ0EFF,ZUREF,                   &
                        PBLOWSNW, HSNOWFALL,                   &
                        PSNOWMAK, OSNOWMAK_BOOL, OSNOWMAK_PROP,IMAX_USE)
!
!***************************************DEBUG IN**********************************************
IF (GCRODEBUGDETAILSPRINT) THEN
  CALL SNOWCROPRINTPROFILE("after SNOWFALL_UPGRID",INLVLS_USE(IDEBUG),LPRINTGRAN,      &
                           PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),   &
                           PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:),&
                           PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:),&
                           HSNOWMETAMO )
END IF
!***************************************DEBUG OUT**********************************************
!
ZSNOW(:) = ZSNOWBIS(:)
!
!*       3.     Update grid/discretization
!               --------------------------
! Reset grid to conform to model specifications:
!
DO JJ=1,SIZE(ZSNOW)
  !
  IF ( GMODIF_MAILLAGE(JJ) ) THEN
    CALL SNOWNLGRIDFRESH_1D(JJ,ZSNOW(JJ),PSNOWDZ(JJ,:),ZSNOWDZN(JJ,:),PSNOWRHO(JJ,:),    &
                            PSNOWHEAT(JJ,:),PSNOWDIAMOPT(JJ,:),PSNOWSPHERI(JJ,:),           &
                            PSNOWHIST(JJ,:),PSNOWAGE(JJ,:),PSNOWIMPUR(JJ,:,:),GSNOWFALL(JJ),ZSNOWRHOF(JJ),  &
                            ZSNOWDZF(JJ),PSNOWHMASS(JJ),ZSNOWDIAMOPTF(JJ),ZSNOWSPHERIF(JJ), &
                            ZSNOWHISTF(JJ),ZSNOWAGEF(JJ),ZSNOWIMPURF(JJ,:),INLVLS_USE(JJ),INLVLS_USE_OLD(JJ),&
                            GSUCCESS)
! Uncomment this in case of troubles
!    IF (.NOT. GSUCCESS) THEN
!      PRINT*, 'regridding problem in snowcro'
!      PRINT*, 'JJ=', JJ
!      PRINT*, 'GSNOWFALL=',GSNOWFALL(JJ)
!      PRINT*, 'ZSNOWDZF=',ZSNOWDZF(JJ)
!      PRINT*, 'ZSR=',ZSR(JJ)
!      PRINT*, 'PSR=', PSR(JJ)
!      PRINT*, 'ZUNLOAD=',ZUNLOAD(JJ)
!      PRINT*, 'PUNLOAD=', PUNLOAD(JJ)
!      PRINT*, 'INLVLS_USE, INLVLS_USE_OLD=',INLVLS_USE(JJ),INLVLS_USE_OLD(JJ)
!      PRINT*, 'ZSNOWDZO=',ZSNOWDZO(JJ,:)
!      PRINT*, 'ZSNOWDZN=',ZSNOWDZN(JJ,:)
!    END IF
  END IF
  !
END DO
!
!
!***************************************DEBUG IN**********************************************
IF (GCRODEBUGDETAILSPRINT) THEN
  CALL SNOWCROPRINTPROFILE("after SNOWNLGRIDFRESH_1D",INLVLS_USE(IDEBUG),LPRINTGRAN,     &
                           PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                         &
                           PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:),&
                           PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:), &
                           HSNOWMETAMO,HSNOWRAD,PSNOWIMPUR(IDEBUG,:,:))
END IF
!***************************************DEBUG OUT**********************************************
! Add ice layer if freezing rain. In MEB case, this will be done after solving heat diffusion
IF (.NOT. OMEB) THEN
  IF (ANY(GFRZRAIN)) THEN
    CALL SNOWCROFREEZINGRAIN(GFRZRAIN, PSNOWDZ, PSNOWRHO, PSNOWHEAT, PSNOWDIAMOPT, PSNOWSPHERI, &
                             PSNOWHIST, PSNOWAGE, PSNOWIMPUR, PRR, PTA, PTSTEP, INLVLS_USE)
  END IF
END IF

!VV Update total snow depth after ice layer from freezing rain
DO JJ = 1,SIZE(ZSNOW)
  ZSNOW(JJ) = SUM(PSNOWDZ(JJ,1:INLVLS_USE(JJ)))
END DO

!*       4.     Liquid water content and snow temperature
!               -----------------------------------------
!
! First diagnose snow temperatures and liquid
! water portion of the snow from snow heat content:
! update some snow layers parameters after new discretization
!
IMAX_USE = MAXVAL(INLVLS_USE)
! active layers
DO JST = 1, IMAX_USE
  DO JJ = 1,SIZE(ZSNOW)
    IF (JST <= INLVLS_USE(JJ)) THEN
      PSNOWSWE (JJ,JST) = PSNOWDZ(JJ,JST) * PSNOWRHO(JJ,JST)
      !
      ZSCAP    (JJ,JST) = PSNOWRHO(JJ,JST) * XCI
      !
      ZSNOWTEMP(JJ,JST) = XTT + &
                          ( ( PSNOWHEAT(JJ,JST)/PSNOWDZ(JJ,JST) + XLMTT*PSNOWRHO(JJ,JST) )/ZSCAP(JJ,JST) )
      !
      PSNOWLIQ (JJ,JST) = MAX( 0.0, ZSNOWTEMP(JJ,JST)-XTT ) * ZSCAP(JJ,JST) * &
                          PSNOWDZ(JJ,JST) / (XLMTT*XRHOLW)
      !
      !
      ZSNOWTEMP(JJ,JST) = MIN( XTT, ZSNOWTEMP(JJ,JST) )
    END IF
  END DO  !  end loop active snow layers
END DO
  !
  ! unactive layers
DO JJ = 1,SIZE(ZSNOW)
  DO JST = INLVLS_USE(JJ)+1,SIZE(PSNOWSWE,2)
    PSNOWSWE    (JJ,JST) = 0.0
    PSNOWRHO    (JJ,JST) = 999.
    PSNOWDZ     (JJ,JST) = 0.
    PSNOWDIAMOPT(JJ,JST) = 0.
    PSNOWSPHERI (JJ,JST) = 0.
    PSNOWHIST   (JJ,JST) = 0.
    PSNOWAGE    (JJ,JST) = 0.
    PSNOWHEAT   (JJ,JST) = 0.
    ZSNOWTEMP   (JJ,JST) = XTT
    PSNOWLIQ    (JJ,JST) = 0.
  END DO  !  end loop unactive snow layers
  !
END DO    ! end loop grid points
!
!***************************************DEBUG IN**********************************************
IF (GCRODEBUGDETAILSPRINT) THEN
  CALL SNOWCROPRINTPROFILE("after liquid water/temperature diagnostic",                  &
                           INLVLS_USE(IDEBUG),LPRINTGRAN,                                &
                           PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                         &
                           PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:),&
                           PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:), &
                           HSNOWMETAMO )
END IF
!***************************************DEBUG OUT**********************************************
!                -----------------
!        In MEB case, this will be done after solving heat diffusion
!        4.BIS   Snow metamorphism
!
IF (.NOT. OMEB) THEN
  !
  CALL SNOWCROGETSSA(PSNOWDIAMOPT,INLVLS_USE,ZSNOWSSA_BEFORE)
  !
  CALL SNOWCROMETAMO(PSNOWDZ,PSNOWDIAMOPT,PSNOWSPHERI,PSNOWHIST,ZSNOWTEMP,   &
                     PSNOWLIQ,PTSTEP,PSNOWSWE,INLVLS_USE,PSNOWAGE,HSNOWMETAMO)
  !
  CALL SNOWCROGETSSA(PSNOWDIAMOPT,INLVLS_USE,ZSNOWSSA_AFTER)
  !
  ZSNOWDSSA=ZSNOWSSA_AFTER-ZSNOWSSA_BEFORE
  !
  !***************************************DEBUG IN**********************************************
  IF (GCRODEBUGDETAILSPRINT) THEN
    CALL SNOWCROPRINTPROFILE("after SNOWCROMETAMO", INLVLS_USE(IDEBUG),LPRINTGRAN,          &
                              PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                         &
                              PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:),&
                              PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:), &
                              HSNOWMETAMO )
  END IF
  !***************************************DEBUG OUT**********************************************
  !
  !*       5.     Snow Compaction
  !               ---------------
  ! Calculate snow density: compaction/aging: density increases
  !
  CALL SNOWCROCOMPACTN(PTSTEP,PSNOWRHO,PSNOWDZ,ZSNOWTEMP,ZSNOW,                      &
                       PSNOWDIAMOPT,PSNOWSPHERI,PSNOWHIST,PSNOWSWE,PSNOWAGE,         &
                       ZSNOWDSSA,PSNOWLIQ,INLVLS_USE,PDIRCOSZW,                      &
                       HSNOWMETAMO,HSNOWCOMP, OSNOWCOMPACT_BOOL, PHVEGPOL)
  !
  !***************************************DEBUG IN**********************************************
  IF (GCRODEBUGDETAILSPRINT) THEN
    CALL SNOWCROPRINTPROFILE("after SNOWCROCOMPACTN", INLVLS_USE(IDEBUG),LPRINTGRAN,        &
                              PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                         &
                              PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:),&
                              PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:), &
                              HSNOWMETAMO )
  END IF
  !***************************************DEBUG OUT**********************************************
  !
  !*       5.1    Snow Compaction and Metamorphism due to snow drift
  !               ---------------
  IF (HSNOWDRIFT  .NE. 'NONE') THEN
    CALL SNOWDRIFT(PTSTEP, PWIND_DRIFT, PSNOWRHO,PSNOWDZ, ZSNOW, HSNOWMETAMO, HSNOWDRIFT,  &
                   HSNOWCOMP,PSNOWDIAMOPT,PSNOWSPHERI,PSNOWHIST,INLVLS_USE,PTA,PQA,PPS,PRHOA,  &
                   PZ0EFF,ZUREF,OSNOWDRIFT_SUBLIM,PSNDRIFT,PVFRIC_T, PHVEGPOL           )
  END IF
  !***************************************DEBUG IN**********************************************
  IF (GCRODEBUGDETAILSPRINT) THEN
    CALL SNOWCROPRINTPROFILE("after SNOWDRIFT", INLVLS_USE(IDEBUG),LPRINTGRAN,              &
                              PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                         &
                              PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:),&
                              PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:), &
                              HSNOWMETAMO )
  END IF
  !***************************************DEBUG OUT**********************************************
  !
  ! Update snow heat content (J/m2) using dry density instead of total density:
  !
  ! Comment M Lafaysse : Before 2020 a threshold was applied on PSNOWDZ (XSNOWDZMIN=1E-4)
  ! in the computation of ZSCAP. This leads to incorrect energy computation in SNOWCROREFRZ
  ! when refreezing occurs on very thin layers
  !
  DO JST = 1, IMAX_USE
    DO JJ = 1,SIZE(ZSNOW)
      IF (JST <= INLVLS_USE(JJ)) THEN
        ZSCAP(JJ,JST) = ( PSNOWRHO(JJ,JST) - PSNOWLIQ(JJ,JST) * XRHOLW / PSNOWDZ(JJ,JST)) * XCI
        !
        PSNOWHEAT(JJ,JST) = PSNOWDZ(JJ,JST) * &
                            ( ZSCAP(JJ,JST)*(ZSNOWTEMP(JJ,JST)-XTT) - XLMTT*PSNOWRHO(JJ,JST) ) + &
                            XLMTT * XRHOLW * PSNOWLIQ(JJ,JST)
      END IF
    END DO  !  end loop snow layers
  END DO    ! end loop grid points
  !
  !***************************************DEBUG IN**********************************************
  IF (GCRODEBUGDETAILSPRINT) THEN
    CALL SNOWCROPRINTPROFILE("after  update snow heat content", INLVLS_USE(IDEBUG),LPRINTGRAN, &
                              PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                            &
                              PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:),   &
                              PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:),    &
                              HSNOWMETAMO      )
  END IF
  !***************************************DEBUG OUT********************************************
!
END IF
!*       6.     Solar radiation transmission
!               -----------------------------
!
! Heat source (-sink) term due to shortwave
! radiation transmission within the snowpack:
!
!
SELECT CASE (HSNOWRAD)
  CASE ("B92","B93")
    PSNOWIMPUR(:,:,:)=0.
    CALL SNOWCRORAD(TPTIME,OMEB,                          &
                    PSW_RAD,PSNOWALB,PSNOWDZ,PSNOWRHO,    &
                    PALB,PSWNETSNOW,PSWNETSNOWS,          &
                    ZRADSINK,ZRADXS,                      &
                    PSNOWDIAMOPT, PSNOWSPHERI, PSNOWAGE,PPS, &
                    PZENITH, PPERMSNOWFRAC,INLVLS_USE,    &
                    OSNOW_ABS_ZENITH, OSNOWAGING_VAR,     &
                    PAGINGCOEF)
  !
  !
  CASE ("T17")!
    IF (NIMPUR >=1) THEN
      !Calculate the factor to norm the impurity content following parameterization from S. Morin(F.tuzet)
      DO JJ=1, SIZE(ZSNOW)
        ZIMPUR_NORM(JJ)=EXP(-0.5*PSNOWDZ(JJ,1)/XIMPUR_EFOLD) !Initialise the norm
        DO JST=2, IMAX_USE
          IF (JST <= INLVLS_USE(JJ)) THEN
            ZWORK2D (JJ,JST) = EXP(-(SUM(PSNOWDZ(JJ,1:JST-1))+0.5*PSNOWDZ(JJ,JST))/XIMPUR_EFOLD)
            ZIMPUR_NORM(JJ)=ZIMPUR_NORM(JJ)+ ZWORK2D (JJ,JST) !add the contribution of each layer
          END IF
        END DO
      END DO
      ! Increase impurity content following parameterization from S. Morin
      DO JIMP=1,NIMPUR
        DO JJ=1, SIZE(ZSNOW)
          PSNOWIMPUR(JJ,1,JIMP)=PSNOWIMPUR(JJ,1,JIMP)+(ZDRYCOEF(JJ,JIMP)*&
          EXP(-0.5*PSNOWDZ(JJ,1)/XIMPUR_EFOLD))/ZIMPUR_NORM(JJ)
          DO JST=2, IMAX_USE
            IF (JST <= INLVLS_USE(JJ)) THEN
              PSNOWIMPUR(JJ,JST,JIMP)=PSNOWIMPUR(JJ,JST,JIMP)+(ZDRYCOEF(JJ,JIMP)*&
              ZWORK2D (JJ,JST))/ZIMPUR_NORM(JJ)
            END IF
          END DO
        END DO
      END DO
     ! TAke into account wet deposit by rainfall
      DO JIMP=1,NIMPUR
        DO JJ=1, SIZE(ZSNOW)
          IF (PRR(JJ)>XUEPSI .AND. ZSR(JJ)<XUEPSI) THEN
            PSNOWIMPUR(JJ,1,JIMP)=PSNOWIMPUR(JJ,1,JIMP)+ZWETCOEF(JJ,JIMP)
          END IF
        END DO
      END DO
    END IF
    !
    CALL SURFACE_IMPURITY_REPARTITION(PSNOWLIQ,PSNOWSWE,PSNOWIMPUR)

    CALL SNOWCRO_TARTES(PSNOWDIAMOPT, PSNOWSPHERI, PSNOWRHO, PSNOWDZ, ZSNOWG0, ZSNOWY0, ZSNOWW0,   &
                        ZSNOWB0, PSNOWIMPUR, PALB, PSW_RAD, PZENITH, PAZIM, PDIRCOSZW, PSLOPEDIR,  &
                        INLVLS_USE, PSNOWALB, ZRADSINK, ZRADXS, GCRODEBUGDETAILSPRINT, HSNOWMETAMO,&
                        P_DIR_SW, P_SCA_SW, ZSNOWALB_SP, ZSPEC_DIR, ZSPEC_DIF, OATMORAD)
                       
    !    PRINT*,"ZRADSINK",ZRADSINK,"ZRADXS",ZRADXS   , "TIME",tptime
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! activation of spectral outputs
    !
    ! output spectral albedo and diffuse to total irradiance ratio
    IF (SIZE(PSPEC_ALB,2) == JPNBANDS_ATM) THEN
      ! upper condition should be true only if LSPECSNOW = TRUE in namelist,
      ! otherwise PSPEC_ALB and PSPEC_TOT have second dimension fixed to 1
      DO JJ=1,SIZE(PSPEC_ALB,2)
        DO JP=1, SIZE(ZSNOW)
          PSPEC_ALB(JP,JJ)=ZSNOWALB_SP(JP,JJ)   ! Snow spectral albedo
          PSPEC_TOT(JP,:)=ZSPEC_DIR(JP,:) + ZSPEC_DIF(JP,:)  ! Total incident radiation spectrally reparted
          IF ((ZSPEC_DIR(JP,JJ)+ZSPEC_DIF(JP,JJ))>0.) THEN
            PDIFF_RATIO(JP,JJ)=ZSPEC_DIF(JP,JJ)/(ZSPEC_DIR(JP,JJ)+ZSPEC_DIF(JP,JJ))
          ELSE
            PDIFF_RATIO(JP,JJ)=XUNDEF
          END IF
        END DO
      END DO
    ELSE
      PSPEC_ALB(:,:) = XUNDEF
      PSPEC_TOT(:,:) = XUNDEF
      PDIFF_RATIO(:,:) = XUNDEF
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !
  CASE DEFAULT
    CALL ABOR1_SFX("UNKNOWN CSNOWRAD OPTION")
  !
END SELECT
!
!***************************************DEBUG IN**********************************************
IF (GCRODEBUGDETAILSPRINT) THEN
  CALL SNOWCROPRINTPROFILE("after SNOWCRORAD", INLVLS_USE(IDEBUG),LPRINTGRAN,              &
                           PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                           &
                           PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:),  &
                           PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:),   &
                           HSNOWMETAMO,HSNOWRAD,PSNOWIMPUR(IDEBUG,:,:))
END IF
!***************************************DEBUG OUT********************************************
!
!*       7.     Heat transfer and surface energy budget
!               ---------------------------------------
! Snow thermal conductivity:
!
 CALL SNOWCROTHRM(PSNOWRHO,ZSCOND,ZSNOWTEMP,PPS,PSNOWLIQ,HSNOWCOND)
!
!***************************************DEBUG IN**********************************************
IF (GCRODEBUGDETAILSPRINT) THEN
  CALL SNOWCROPRINTPROFILE("after SNOWCROTHRM", INLVLS_USE(IDEBUG),LPRINTGRAN,             &
                           PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                           &
                           PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:),  &
                           PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:),   &
                           HSNOWMETAMO)
END IF
!***************************************DEBUG OUT********************************************
!
! Precipitation heating term:
! Rainfall renders it's heat to the snow when it enters
! the snowpack:
! if freezing rain (PTA<XTT), PHPSNOW=0, no special case necessary
PHPSNOW(:) = PRR(:) * XCL * ( MAX( XTT,PTA(:) ) - XTT )    ! (W/m2)
!
! Surface Energy Budget calculations using ISBA linearized form
! and standard ISBA turbulent transfer formulation
!
IF ( ALL(PPEW_A_COEF==0.) ) THEN
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Modif Matthieu Lafaysse for glaciers
  ! For surface ice, modify roughness lengths
  ! Only if not implicit coupling
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  WHERE( PSNOWRHO(:,1)>XRHOTHRESHOLD_ICE )
    ZZ0_SNOWICE    = PZ0    * XZ0ICEZ0SNOW
    ZZ0H_SNOWICE   = PZ0H   * XZ0ICEZ0SNOW
    ZZ0EFF_SNOWICE = PZ0EFF * XZ0ICEZ0SNOW
  ELSEWHERE
    ZZ0_SNOWICE    = PZ0
    ZZ0H_SNOWICE   = PZ0H
    ZZ0EFF_SNOWICE = PZ0EFF
  END WHERE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE
  ZZ0_SNOWICE    = PZ0
  ZZ0H_SNOWICE   = PZ0H
  ZZ0EFF_SNOWICE = PZ0EFF
END IF

! Recompute this with update snowdepth
IF (HSNOWZREF=="VAR") THEN
  ZZREF(:)=MAX(PZREF(:)-ZSNOW(:),0.2)
  ZUREF(:)=MAX(PUREF(:)-ZSNOW(:),0.2)
ELSE
  ZZREF(:)=PZREF(:)
  ZUREF(:)=PUREF(:)
END IF

IF (OMEB) THEN
  CALL SNOWCROEBUDMEB(PTSTEP,XSNOWDZMIN,                                     &
                      ZSNOWTEMP(:,1),PSNOWDZ(:,1),PSNOWDZ(:,2),              &
                      ZSCOND(:,1),ZSCOND(:,2),ZSCAP(:,1),                    &
                      PRNSNOW,PHSNOW,PLES3L,PLEL3L,ZRADSINK(:,1),PHPSNOW,    &
                      ZCT,ZTSTERM1,ZTSTERM2,PGFLUXSNOW)
ELSE
  CALL SNOWCROEBUD(HSNOWRES, HIMPLICIT_WIND,                                    &
                   PPEW_A_COEF, PPEW_B_COEF,                                    &
                   PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,          &
                   XSNOWDZMIN,                                                  &
                   ZZREF,ZSNOWTEMP(:,1),PSNOWRHO(:,1),PSNOWLIQ(:,1),ZSCAP(:,1), &
                   ZSCOND(:,1),ZSCOND(:,2),                                     &
                   ZUREF,PEXNS,PEXNA,PDIRCOSZW,PVMOD,                           &
                   PLW_RAD,PSW_RAD,PTA,PQA,PPS,PTSTEP,                          &
                   PSNOWDZ(:,1),PSNOWDZ(:,2),PSNOWALB,ZZ0_SNOWICE,              &
                   ZZ0EFF_SNOWICE,ZZ0H_SNOWICE,                                 &
                   ZSFCFRZ,ZRADSINK(:,1),PHPSNOW,                               &
                   ZCT,PEMISNOW,PRHOA,ZTSTERM1,ZTSTERM2,ZRA,PCDSNOW,PCHSNOW,    &
                   ZQSAT, ZDQSAT, ZRSRA, ZUSTAR2_IC, PRI,                       &
                   ZPET_A_COEF_T,ZPEQ_A_COEF_T,ZPET_B_COEF_T,ZPEQ_B_COEF_T ,    &
                   GFRZRAIN,PRR, PRESA_SN)
END IF
!
!***************************************DEBUG IN**********************************************
IF (GCRODEBUGDETAILSPRINT) THEN
  CALL SNOWCROPRINTPROFILE("after SNOWCROEBUD", INLVLS_USE(IDEBUG),LPRINTGRAN,              &
                           PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),        &
                           PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:),     &
                           PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:),     &
                           HSNOWMETAMO)
END IF
!***************************************DEBUG OUT********************************************
!
! Heat transfer: simple diffusion along the thermal gradient
!
ZSNOWTEMPO1(:) = ZSNOWTEMP(:,1) ! save surface snow temperature before update
!
ZGRNDFLUXI(:)  = PGRNDFLUX(:) ! on sauvegarde le flux imposé par MEB
!
CALL SNOWCROSOLVT(OMEB,PTSTEP,XSNOWDZMIN,PSNOWDZ,ZSCOND,ZSCAP,PTG,           &
                  PSOILCOND,PD_G,ZRADSINK,ZCT,ZTSTERM1,ZTSTERM2,             &
                  ZPET_A_COEF_T,ZPEQ_A_COEF_T,ZPET_B_COEF_T,ZPEQ_B_COEF_T,   &
                  ZTA_IC,ZQA_IC,PGRNDFLUX, ZSNOWTEMP ,PSNOWFLUX,             &
                  INLVLS_USE                                                 )
!
!***************************************DEBUG IN**********************************************
IF (GCRODEBUGDETAILSPRINT) THEN
  CALL SNOWCROPRINTPROFILE("after SNOWCROSOLVT", INLVLS_USE(IDEBUG),LPRINTGRAN,            &
                           PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                           &
                           PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:),  &
                           PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:),   &
                           HSNOWMETAMO)
END IF
!***************************************DEBUG OUT********************************************
!
!*       8.     Surface fluxes
!               --------------
!

! Since surface fluxes already computed under MEB option (and already
! recomputed for the case when T>Tf), Only call if MEB not in use:

IF (.NOT. OMEB) THEN
  CALL SNOWCROFLUX(ZSNOWTEMP(:,1),PSNOWDZ(:,1),PEXNS,PEXNA,            &
                   ZUSTAR2_IC,                                         &
                   PTSTEP,PSNOWALB,PSW_RAD,PEMISNOW,PLW_RAD,           &
                   ZTA_IC,ZSFCFRZ,ZQA_IC,PHPSNOW,                      &
                   ZSNOWTEMPO1,PSNOWFLUX,ZCT,ZRADSINK(:,1),            &
                   ZQSAT,ZDQSAT,ZRSRA,                                 &
                   PRNSNOW,PHSNOW,PGFLUXSNOW,PLES3L,PLEL3L,PEVAP,      &
                   PUSTAR, GFRZRAIN,PRR)
END IF
!
! Update snow heat content (J/m2) using dry density with new temperatures to appropriately check
! the full snow cover disappearance with SNOWCROGONE
!
! Introduced by M Lafaysse 28/10/2022 to solve a bug with no detection of vanishing snow
! and moved after snowcroflux 26/01/2023
!
DO JST = 1, IMAX_USE
  DO JJ = 1,SIZE(ZSNOW)
    IF (JST <= INLVLS_USE(JJ)) THEN
      ZSCAP(JJ,JST) = ( PSNOWRHO(JJ,JST) - PSNOWLIQ(JJ,JST) * XRHOLW / PSNOWDZ(JJ,JST)) * XCI
      !
      PSNOWHEAT(JJ,JST) = PSNOWDZ(JJ,JST) * &
                          ( ZSCAP(JJ,JST)*(ZSNOWTEMP(JJ,JST)-XTT) - XLMTT*PSNOWRHO(JJ,JST) ) + &
                          XLMTT * XRHOLW * PSNOWLIQ(JJ,JST)
    END IF
  END DO  !  end loop snow layers
END DO    ! end loop grid points
!
!***************************************DEBUG IN**********************************************
IF (GCRODEBUGDETAILSPRINT) THEN
  CALL SNOWCROPRINTPROFILE("after SNOWCROFLUX", INLVLS_USE(IDEBUG),LPRINTGRAN,            &
                           PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                          &
                           PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:), &
                           PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:),  &
                           HSNOWMETAMO)
END IF
!***************************************DEBUG OUT********************************************
!*       9.     Snow melt
!               ---------
!
! First Test to see if snow pack vanishes during this time step:
!
IF (PRESENT(PSFCMELT)) THEN
  DO JJ=1,SIZE(PTA)
    JST = INLVLS_USE(JJ)
    ZICE_BEFORE_GONE(JJ) = MAX(0., SUM(PSNOWRHO(JJ,1:JST)*PSNOWDZ(JJ,1:JST) - &
        PSNOWLIQ(JJ,1:JST)*XRHOLW)/PTSTEP - MAX(0., PLES3L(JJ))/XLSTT)
  END DO
  PSFCMELT = 0.
END IF
CALL SNOWCROGONE(PTSTEP,PLEL3L,PLES3L,PSNOWRHO,                  &
                 PSNOWHEAT,ZRADSINK,PEVAPCOR,PTHRUFAL,PGRNDFLUX, &
                 PGFLUXSNOW,PSNOWDZ,PSNOWLIQ,ZSNOWTEMP,ZRADXS,   &
                 PRR,INLVLS_USE,GFRZRAIN                         )
!
!***************************************DEBUG IN**********************************************
IF (PRESENT(PSFCMELT)) THEN
  WHERE (INLVLS_USE == 0) PSFCMELT = ZICE_BEFORE_GONE
END IF
IF (GCRODEBUGDETAILSPRINT) THEN
  CALL SNOWCROPRINTPROFILE("after SNOWCROGONE", INLVLS_USE(IDEBUG),LPRINTGRAN,              &
                           PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                            &
                           PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:),   &
                           PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:),    &
                           HSNOWMETAMO)
END IF
!***************************************DEBUG OUT********************************************
!
! Add radiation not absorbed by snow to soil/vegetation interface flux
! (for thin snowpacks):
!
PGRNDFLUX(:) = PGRNDFLUX(:) + ZRADXS(:)
!
! Second Test to see if one or several snow layers vanishe during this time
! step. In such a case, the concerned snow layers are agregated to neighbours
!
CALL SNOWCROLAYER_GONE(PTSTEP,ZSCAP,PSNOWHEAT,ZSNOWTEMP,PSNOWDZ,          &
                       PSNOWRHO,PSNOWLIQ,PSNOWDIAMOPT,PSNOWSPHERI,        &
                       PSNOWHIST,PSNOWAGE,PSNOWIMPUR,PLES3L, INLVLS_USE   )
!
IMAX_USE = MAXVAL(INLVLS_USE)
!
!***************************************DEBUG IN**********************************************
IF (GCRODEBUGDETAILSPRINT) THEN
  CALL SNOWCROPRINTPROFILE("after SNOWCROLAYER_GONE", INLVLS_USE(IDEBUG),LPRINTGRAN,      &
                           PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                          &
                           PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:), &
                           PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:),  &
                           HSNOWMETAMO,HSNOWRAD,PSNOWIMPUR(IDEBUG,:,:))
END IF
!***************************************DEBUG OUT********************************************
!
! For partial melt: transform excess heat content into snow liquid:
!
IF (PRESENT(PSFCMELT)) ZSURFACE_LIQ = PSNOWLIQ(:,1)
CALL SNOWCROMELT(ZSCAP,ZSNOWTEMP,PSNOWDZ,PSNOWRHO,PSNOWLIQ,PSNOWIMPUR,INLVLS_USE)
IF (PRESENT(PSFCMELT)) PSFCMELT = PSFCMELT + MAX(0., PSNOWLIQ(:,1)-ZSURFACE_LIQ)*XRHOLW/PTSTEP
!
!***************************************DEBUG IN**********************************************
IF (GCRODEBUGDETAILSPRINT) THEN
  CALL SNOWCROPRINTPROFILE("after SNOWCROMELT", INLVLS_USE(IDEBUG),LPRINTGRAN,            &
                           PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                          &
                           PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:), &
                           PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:),  &
                           HSNOWMETAMO,HSNOWRAD,PSNOWIMPUR(IDEBUG,:,:))
END IF
!***************************************DEBUG OUT********************************************
!
!*      10.     Snow water flow and refreezing
!               ------------------------------
! Liquid water vertical transfer and possible snowpack runoff
! And refreezing/freezing of meltwater/rainfall (ripening of the snow)
!
CALL SNOWCROREFRZ(PTSTEP,PRR,PSNOWRHO,ZSNOWTEMP,PSNOWDZ,PSNOWLIQ,PSNOWIMPUR,  &
                  PTHRUFAL,ZSCAP,PLEL3L,INLVLS_USE,HSNOWHOLD,IMAX_USE,GFRZRAIN)
!
!***************************************DEBUG IN**********************************************
IF (GCRODEBUGDETAILSPRINT) THEN
  CALL SNOWCROPRINTPROFILE("after SNOWCROREFRZ", INLVLS_USE(IDEBUG),LPRINTGRAN,           &
                           PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                          &
                           PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:), &
                           PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:),  &
                           HSNOWMETAMO,HSNOWRAD,PSNOWIMPUR(IDEBUG,:,:))
END IF
!***************************************DEBUG OUT********************************************
!
! Compute dry density. This quantity is not allowed to change anymore during this time step
! except in the case of regridding (SNOWCROEVAPGONE or SNOWCROFREEZINGRAIN)
DO JST = 1, IMAX_USE
  DO JJ = 1,SIZE(ZSNOW)
    IF (JST <= INLVLS_USE(JJ)) THEN
      ZDRYDENSITY(JJ,JST) = PSNOWRHO(JJ,JST) - PSNOWLIQ(JJ,JST) * XRHOLW / PSNOWDZ(JJ,JST)
    END IF
  END DO
END DO
!
!*      11.     Snow Evaporation/Sublimation mass updates:
!               ------------------------------------------
!
CALL SNOWCROEVAPN(PLES3L,PTSTEP,ZSNOWTEMP(:,1), ZDRYDENSITY(:,1),  &
                  PSNOWDZ(:,1),PSNOWLIQ(:,1:2),PEVAPCOR,PSNOWHMASS )
!
!***************************************DEBUG IN**********************************************
IF (GCRODEBUGDETAILSPRINT) THEN
  CALL SNOWCROPRINTPROFILE("after SNOWCROEVAPN", INLVLS_USE(IDEBUG),LPRINTGRAN,           &
                           PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                          &
                           PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:), &
                           PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:),  &
                           HSNOWMETAMO,HSNOWRAD,PSNOWIMPUR(IDEBUG,:,:))
END IF
!***************************************DEBUG OUT********************************************
!
! If all snow in uppermost layer evaporates/sublimates, re-distribute
! grid (below could be evoked for vanishingly thin snowpacks):
!
CALL SNOWCROEVAPGONE(PSNOWHEAT,PSNOWDZ, ZDRYDENSITY, ZSNOWTEMP,PSNOWLIQ,PSNOWDIAMOPT, &
                     PSNOWSPHERI,PSNOWHIST,PSNOWAGE,INLVLS_USE,HSNOWMETAMO            )
!
!***************************************DEBUG IN**********************************************
IF (GCRODEBUGDETAILSPRINT) THEN
  CALL SNOWCROPRINTPROFILE("after SNOWCROEVAPGONE", INLVLS_USE(IDEBUG),LPRINTGRAN,        &
                           PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                          &
                           PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:), &
                           PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:),  &
                           HSNOWMETAMO,HSNOWRAD,PSNOWIMPUR(IDEBUG,:,:))
END IF
!***************************************DEBUG OUT********************************************
!
!                -----------------
!        4.BIS   Snow metamorphism
!
IF (OMEB) THEN
  ! Adjust liquid water content as unrealistic densities might be possible in the 2 first layers after sublimation
  ! and update enthalpy
  ! Cluzet et al 2016 : different lwc options
  IF ( HSNOWHOLD == 'B92' ) THEN
    DO JST = 1, IMAX_USE
      DO JJ = 1,SIZE(ZSNOW)
        IF (JST <= INLVLS_USE(JJ)) THEN
          ZWHOLDMAX (JJ,JST) = XPERCENTAGEPORE/XRHOLI * (PSNOWDZ(JJ,JST) * &
                               (XRHOLI-ZDRYDENSITY(JJ,JST)))
        END IF
      END DO
    END DO
  ELSEIF ( HSNOWHOLD == 'BFZ' ) THEN
    ZSNOWDZ_FRZ=0.
    DO JST = 1, IMAX_USE
      DO JJ = 1,SIZE(ZSNOW)
        IF (JST <= INLVLS_USE(JJ)) THEN
          ! GFRZRAIN
          IF (GFRZRAIN(JJ)) THEN
            ZSNOWDZ_FRZ(JJ)=PSNOWDZ(JJ,JST)+ZSNOWDZ_FRZ(JJ)
            IF (ZSNOWDZ_FRZ(JJ) <= 0.02) THEN
              ZWHOLDMAX (JJ,JST) = XPERCENTAGEPORE_FRZ/XRHOLI * (PSNOWDZ(JJ,JST) * &
                                   (XRHOLI-ZDRYDENSITY(JJ,JST)))
              CONTINUE
            END IF
          END IF
          ! PSNOWRHO > 700
          IF (PSNOWRHO(JJ,JST)>700.) THEN
            ZWHOLDMAX (JJ,JST) = XPERCENTAGEPORE_ICE/XRHOLI * (PSNOWDZ(JJ,JST) * &
                                (XRHOLI-PSNOWRHO(JJ,JST)) + PSNOWLIQ(JJ,JST)*XRHOLW)
          ELSE
            ZWHOLDMAX (JJ,JST) = XPERCENTAGEPORE/XRHOLI * (PSNOWDZ(JJ,JST) * &
                                 (XRHOLI-ZDRYDENSITY(JJ,JST)))
          END IF
        END IF
      END DO
    END DO
  ELSE IF ( HSNOWHOLD == 'B02' ) THEN
    DO JST = 1, IMAX_USE
      DO JJ = 1,SIZE(ZSNOW)
        IF (JST <= INLVLS_USE(JJ)) THEN
          ZWHOLDMAX(JJ,JST) = SNOW3LHOLD( PSNOWRHO(JJ,JST),PSNOWDZ(JJ,JST))
        END IF
      END DO
    END DO
  ELSE IF ( HSNOWHOLD == 'SPK' ) THEN
    DO JST = 1, IMAX_USE
      DO JJ = 1,SIZE(ZSNOW)
        IF (JST <= INLVLS_USE(JJ)) THEN
          ZWHOLDMAX(JJ,JST) = SNOWSPKHOLD(ZDRYDENSITY(JJ,JST),PSNOWDZ(JJ,JST))
        END IF
      END DO
    END DO
  ELSE IF ( HSNOWHOLD == 'O04' ) THEN
    DO JST = 1, IMAX_USE
      DO JJ = 1,SIZE(ZSNOW)
        IF (JST <= INLVLS_USE(JJ)) THEN
          ZWHOLDMAX(JJ,JST) = SNOWO04HOLD_0D(ZDRYDENSITY(JJ,JST),PSNOWDZ(JJ,JST) )
        END IF
      END DO
    END DO
  END IF
  !
  DO JST = 1, IMAX_USE
    DO JJ = 1,SIZE(ZSNOW)
      IF (JST <= INLVLS_USE(JJ)) THEN
        ZLIQHEATXS(JJ)     = MAX( 0.0, (PSNOWLIQ(JJ,JST) - ZWHOLDMAX(JJ,JST)) * XRHOLW ) * XLMTT/PTSTEP
        PSNOWLIQ  (JJ,JST) = PSNOWLIQ(JJ,JST) - ZLIQHEATXS(JJ)*PTSTEP/(XRHOLW*XLMTT)
        PSNOWLIQ  (JJ,JST) = MAX( 0.0, PSNOWLIQ(JJ,JST) )
        ! Update total density
        PSNOWRHO  (JJ,JST) = ZDRYDENSITY(JJ,JST) + PSNOWLIQ(JJ,JST) * XRHOLW / PSNOWDZ(JJ,JST)    
        PGRNDFLUX (JJ)     = PGRNDFLUX(JJ) + ZLIQHEATXS(JJ)
        PSNOWTEMP (JJ,JST) = ZSNOWTEMP(JJ,JST)
        ! Heat content using total density
        ZSCAP     (JJ,JST) = PSNOWRHO(JJ,JST) * XCI
        PSNOWHEAT (JJ,JST) = PSNOWDZ(JJ,JST) * &
                             ( ZSCAP(JJ,JST)*(PSNOWTEMP(JJ,JST)-XTT) - XLMTT*PSNOWRHO(JJ,JST) ) + &
                             XLMTT * XRHOLW * PSNOWLIQ(JJ,JST)
        !
        PSNOWSWE  (JJ,JST)  = PSNOWDZ(JJ,JST) * PSNOWRHO(JJ,JST)
      END IF
    END DO  !  end loop active snow layers
  END DO
  !
  !***************************************DEBUG IN**********************************************
  IF (GCRODEBUGDETAILSPRINT) THEN
    CALL SNOWCROPRINTPROFILE("after update", INLVLS_USE(IDEBUG),LPRINTGRAN,          &
                       PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                         &
                       PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:),&
                       PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:), &
                       HSNOWMETAMO )
  END IF
  !***************************************DEBUG OUT**********************************************
!                  
  IF (ANY(GFRZRAIN)) THEN
    CALL SNOWCROFREEZINGRAIN(GFRZRAIN, PSNOWDZ, PSNOWRHO, PSNOWHEAT, PSNOWDIAMOPT, PSNOWSPHERI, &
                             PSNOWHIST, PSNOWAGE, PSNOWIMPUR, PRR, PTA, PTSTEP, INLVLS_USE)
    !
    !***************************************DEBUG IN**********************************************
    IF (GCRODEBUGDETAILSPRINT) THEN
      CALL SNOWCROPRINTPROFILE("after SNOWCROFREEZINGRAIN", INLVLS_USE(IDEBUG),LPRINTGRAN,   &
                               PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                         &
                               PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:),&
                               PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:), &
                               HSNOWMETAMO )
    END IF
  !***************************************DEBUG OUT**********************************************
  END IF
!
  CALL SNOWCROUNLOAD(PSNOWDZ, PSNOWRHO, PSNOWHEAT, PSNOWDIAMOPT, PSNOWSPHERI, &
                     PSNOWHIST, PSNOWAGE, PSNOWIMPUR, ZUNLOAD, PTA, PTSTEP, INLVLS_USE)
  !
  ! Update temperature and liquid water and dry density after freezing rain and unloading
  !
  IMAX_USE = MAXVAL(INLVLS_USE)
  ! active layers
  DO JST = 1, IMAX_USE
    DO JJ = 1,SIZE(ZSNOW)
      IF (JST <= INLVLS_USE(JJ)) THEN
        !
        PSNOWSWE (JJ,JST) = PSNOWDZ(JJ,JST) * PSNOWRHO(JJ,JST)
        !
        ZSCAP    (JJ,JST) = PSNOWRHO(JJ,JST) * XCI
        !
        ZSNOWTEMP(JJ,JST) = XTT + &
                            ( ( PSNOWHEAT(JJ,JST)/PSNOWDZ(JJ,JST) + XLMTT*PSNOWRHO(JJ,JST) )/ZSCAP(JJ,JST) )
        !
        PSNOWLIQ (JJ,JST) = MAX( 0.0, ZSNOWTEMP(JJ,JST)-XTT ) * ZSCAP(JJ,JST) * &
                                 PSNOWDZ(JJ,JST) / (XLMTT*XRHOLW)
        !
        ZSNOWTEMP(JJ,JST) = MIN( XTT, ZSNOWTEMP(JJ,JST) )
        !
        ! update dry density 
        ZDRYDENSITY(JJ,JST) = PSNOWRHO(JJ,JST) - PSNOWLIQ(JJ,JST) * XRHOLW / PSNOWDZ(JJ,JST)
      END IF
    END DO  !  end loop active snow layers
  END DO
  !
  !***************************************DEBUG IN**********************************************
  IF (GCRODEBUGDETAILSPRINT) THEN
    CALL SNOWCROPRINTPROFILE("after SNOWCROUNLOAD", INLVLS_USE(IDEBUG),LPRINTGRAN,          &
                              PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                         &
                              PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:),&
                              PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:), &
                              HSNOWMETAMO )
  END IF
  !***************************************DEBUG OUT**********************************************
  !
  CALL SNOWCROGETSSA(PSNOWDIAMOPT,INLVLS_USE,ZSNOWSSA_BEFORE)
  !
  CALL SNOWCROMETAMO(PSNOWDZ,PSNOWDIAMOPT,PSNOWSPHERI,PSNOWHIST,ZSNOWTEMP,     &
                     PSNOWLIQ,PTSTEP,PSNOWSWE,INLVLS_USE,PSNOWAGE,HSNOWMETAMO  )
  !
  CALL SNOWCROGETSSA(PSNOWDIAMOPT,INLVLS_USE,ZSNOWSSA_AFTER)
  !
  ZSNOWDSSA=ZSNOWSSA_AFTER-ZSNOWSSA_BEFORE
  !
  !***************************************DEBUG IN**********************************************
  IF (GCRODEBUGDETAILSPRINT) THEN
    CALL SNOWCROPRINTPROFILE("after SNOWCROMETAMO", INLVLS_USE(IDEBUG),LPRINTGRAN,         &
                             PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                         &
                             PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:),&
                             PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:), &
                             HSNOWMETAMO )
  END IF
  !***************************************DEBUG OUT**********************************************
  !
  !*       5.     Snow Compaction
  !               ---------------
  ! Calculate snow density: compaction/aging: density increases
  !
  CALL SNOWCROCOMPACTN(PTSTEP,PSNOWRHO,PSNOWDZ,ZSNOWTEMP,ZSNOW,               &
                  PSNOWDIAMOPT,PSNOWSPHERI,PSNOWHIST,PSNOWSWE,PSNOWAGE,       &
                  ZSNOWDSSA,PSNOWLIQ,INLVLS_USE,PDIRCOSZW,                    &
                  HSNOWMETAMO,HSNOWCOMP, OSNOWCOMPACT_BOOL, PHVEGPOL)
  !
  !***************************************DEBUG IN**********************************************
  IF (GCRODEBUGDETAILSPRINT) THEN
    CALL SNOWCROPRINTPROFILE("after SNOWCROCOMPACTN", INLVLS_USE(IDEBUG),LPRINTGRAN,        &
                             PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                          &
                             PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:), &
                             PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:),  &
                             HSNOWMETAMO )
  END IF
  !***************************************DEBUG OUT**********************************************
  !
  !*       5.1    Snow Compaction and Metamorphism due to snow drift
  !               ---------------
  !
  IF (HSNOWDRIFT  .NE. 'NONE') THEN
    CALL SNOWDRIFT(PTSTEP, PVMOD, PSNOWRHO,PSNOWDZ, ZSNOW, HSNOWMETAMO, HSNOWDRIFT,     &
             HSNOWCOMP,PSNOWDIAMOPT,PSNOWSPHERI,PSNOWHIST,INLVLS_USE,PTA,PQA,PPS,PRHOA,&
             PZ0EFF,ZUREF,OSNOWDRIFT_SUBLIM,PSNDRIFT,PVFRIC_T,PHVEGPOL)
  END IF
  !***************************************DEBUG IN**********************************************
  IF (GCRODEBUGDETAILSPRINT) THEN
    CALL SNOWCROPRINTPROFILE("after SNOWDRIFT", INLVLS_USE(IDEBUG),LPRINTGRAN,              &
                             PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                          &
                             PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:), &
                             PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:),  &
                             HSNOWMETAMO )
  END IF
  !***************************************DEBUG OUT**********************************************
  !
  ! Update snow heat content (J/m2) using dry density instead of total density:
  !
  ! Comment M Lafaysse : Before 2020 a threshold was applied on PSNOWDZ (XSNOWDZMIN=1E-4)
  ! in the computation of ZSCAP. This leads to incorrect energy computation in SNOWCROREFRZ
  ! when refreezing occurs on very thin layers
  !
  DO JST = 1, IMAX_USE
    DO JJ = 1,SIZE(ZSNOW)
      IF (JST <= INLVLS_USE(JJ)) THEN
        ZSCAP(JJ,JST) = ( PSNOWRHO(JJ,JST) - PSNOWLIQ(JJ,JST) * XRHOLW / PSNOWDZ(JJ,JST)) * XCI
        !
        PSNOWHEAT(JJ,JST) = PSNOWDZ(JJ,JST) * &
                            ( ZSCAP(JJ,JST)*(ZSNOWTEMP(JJ,JST)-XTT) - XLMTT*PSNOWRHO(JJ,JST) ) + &
                            XLMTT * XRHOLW * PSNOWLIQ(JJ,JST)
        ! update dry density 
        ZDRYDENSITY(JJ,JST) = PSNOWRHO(JJ,JST) - PSNOWLIQ(JJ,JST) * XRHOLW / PSNOWDZ(JJ,JST)
      END IF
    END DO  !  end loop snow layers
  END DO    ! end loop grid points
  !
  !***************************************DEBUG IN**********************************************
  IF (GCRODEBUGDETAILSPRINT) THEN
    CALL SNOWCROPRINTPROFILE("after update snow heat content", INLVLS_USE(IDEBUG),LPRINTGRAN,&
                             PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),                           &
                             PSNOWLIQ(IDEBUG,:),PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:),  &
                             PSNOWSPHERI(IDEBUG,:),PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:),   &
                             HSNOWMETAMO )
  END IF
  !***************************************DEBUG OUT********************************************
!
END IF
!
!*      12.     Update surface albedo:
!               ----------------------
! Snow clear sky albedo:
!
IMAX_USE = MAXVAL(INLVLS_USE)
IF ( HSNOWRAD=='B92' .OR. HSNOWRAD=='B93' ) THEN
  CALL SNOWCROALB(PSNOWALB,ZSPECTRALALBEDO,PSNOWDZ(:,1),PSNOWRHO(:,1:2),          &
                  PPERMSNOWFRAC,PSNOWDIAMOPT(:,1),PSNOWSPHERI(:,1),               &
                  PSNOWAGE(:,1),PSNOWDIAMOPT(:,2),PSNOWSPHERI(:,2),PSNOWAGE(:,2), &
                  PPS, PZENITH,OSNOWAGING_VAR, PAGINGCOEF, INLVLS_USE)
  !
  !the albedo is not updated in the case of TARTES scheme
ELSEIF ( HSNOWRAD=='T17') THEN
    DO JJ=1, SIZE(ZSNOW)
      IF (SUM(PSNOWDZ(JJ,1:INLVLS_USE(JJ))) .LT. EPSILON(XUNDEF))THEN
        PSPEC_ALB(JJ,:)=XUNDEF
      END IF
    END DO
END IF
!
!*      13.     Update snow heat content:
!               -------------------------
! Update the heat content (variable stored each time step)
! using current snow temperature and liquid water content:
!
! First, make check to make sure heat content not too large
! (this can result due to signifigant heating of thin snowpacks):
! add any excess heat to ground flux:
!
      ! Cluzet et al 2016 : different lwc options
IF ( HSNOWHOLD == 'B92' ) THEN
  DO JST = 1, IMAX_USE
    DO JJ = 1,SIZE(ZSNOW)
      IF (JST <= INLVLS_USE(JJ)) THEN
        ZWHOLDMAX (JJ,JST) = XPERCENTAGEPORE/XRHOLI * (PSNOWDZ(JJ,JST) * &
        (XRHOLI-ZDRYDENSITY(JJ,JST)))
      END IF
    END DO
  END DO
ELSEIF ( HSNOWHOLD == 'BFZ' ) THEN
  ZSNOWDZ_FRZ=0.
  DO JST = 1, IMAX_USE
    DO JJ = 1,SIZE(ZSNOW)
      IF (JST <= INLVLS_USE(JJ)) THEN
        IF (GFRZRAIN(JJ)) Then
          ZSNOWDZ_FRZ(JJ)=PSNOWDZ(JJ,JST)+ZSNOWDZ_FRZ(JJ)
            IF (ZSNOWDZ_FRZ(JJ) <= 0.02) THEN
               ZWHOLDMAX (JJ,JST) = XPERCENTAGEPORE_FRZ/XRHOLI * (PSNOWDZ(JJ,JST) * &
               (XRHOLI-ZDRYDENSITY(JJ,JST)))
               CONTINUE
            END IF
        END IF
        IF (PSNOWRHO(JJ,JST)>700.) THEN
          ZWHOLDMAX (JJ,JST) = XPERCENTAGEPORE_ICE/XRHOLI * (PSNOWDZ(JJ,JST) * &
          (XRHOLI-PSNOWRHO(JJ,JST)) + PSNOWLIQ(JJ,JST)*XRHOLW)
        ELSE
          ZWHOLDMAX (JJ,JST) = XPERCENTAGEPORE/XRHOLI * (PSNOWDZ(JJ,JST) * &
              (XRHOLI-ZDRYDENSITY(JJ,JST)))
        END IF
      END IF
    END DO
  END DO
ELSE IF ( HSNOWHOLD == 'B02' ) THEN
  DO JST = 1, IMAX_USE
    DO JJ = 1,SIZE(ZSNOW)
      IF (JST <= INLVLS_USE(JJ)) THEN
        ZWHOLDMAX(JJ,JST) = SNOW3LHOLD( PSNOWRHO(JJ,JST),PSNOWDZ(JJ,JST))
      END IF
    END DO
  END DO
ELSE IF ( HSNOWHOLD == 'SPK' ) THEN
  DO JST = 1, IMAX_USE
    DO JJ = 1,SIZE(ZSNOW)
      IF (JST <= INLVLS_USE(JJ)) THEN
        ZWHOLDMAX(JJ,JST) = SNOWSPKHOLD(ZDRYDENSITY(JJ,JST),PSNOWDZ(JJ,JST))
      END IF
    END DO
  END DO
ELSE IF ( HSNOWHOLD == 'O04' ) THEN
  DO JST = 1, IMAX_USE
    DO JJ = 1,SIZE(ZSNOW)
      IF (JST <= INLVLS_USE(JJ)) THEN
        ZWHOLDMAX(JJ,JST) = SNOWO04HOLD_0D(ZDRYDENSITY(JJ,JST),PSNOWDZ(JJ,JST) )
      END IF
    END DO
  END DO
END IF
DO JST = 1, IMAX_USE
  DO JJ = 1,SIZE(ZSNOW)
    IF (JST <= INLVLS_USE(JJ)) THEN
      ZLIQHEATXS(JJ)     = MAX( 0.0, (PSNOWLIQ(JJ,JST) - ZWHOLDMAX(JJ,JST)) * XRHOLW ) * XLMTT/PTSTEP
      PSNOWLIQ  (JJ,JST) = PSNOWLIQ(JJ,JST) - ZLIQHEATXS(JJ)*PTSTEP/(XRHOLW*XLMTT)
      PTHRUFAL (JJ) = PTHRUFAL (JJ) + ZLIQHEATXS(JJ)/XLMTT  ! account for the excess removed above      
      PSNOWLIQ  (JJ,JST) = MAX( 0.0, PSNOWLIQ(JJ,JST) )
      ! Update total density
      PSNOWRHO  (JJ,JST) = ZDRYDENSITY(JJ,JST) + PSNOWLIQ(JJ,JST) * XRHOLW / PSNOWDZ(JJ,JST)    
      PGRNDFLUX (JJ)     = PGRNDFLUX(JJ) + ZLIQHEATXS(JJ)
      PSNOWTEMP (JJ,JST) = ZSNOWTEMP(JJ,JST)
  !   Heat content using total density
      ZSCAP     (JJ,JST) = PSNOWRHO(JJ,JST) * XCI
      PSNOWHEAT (JJ,JST) = PSNOWDZ(JJ,JST) * &
                        ( ZSCAP(JJ,JST)*(PSNOWTEMP(JJ,JST)-XTT) - XLMTT*PSNOWRHO(JJ,JST) ) + &
                        XLMTT * XRHOLW * PSNOWLIQ(JJ,JST)
  !
      PSNOWSWE(JJ,JST)  = PSNOWDZ(JJ,JST) * PSNOWRHO(JJ,JST)
    END IF
  END DO  !  end loop active snow layers
END DO
!
! unactive layers
DO JJ = 1,SIZE(ZSNOW)
  DO JST = INLVLS_USE(JJ)+1,SIZE(PSNOWSWE,2)
    PSNOWSWE(JJ,JST)  = 0.
    PSNOWHEAT(JJ,JST) = 0.
    PSNOWRHO(JJ,JST)  = 999.
    PSNOWTEMP(JJ,JST) = 0.
    PSNOWDZ(JJ,JST)   = 0.
    DO JIMP=1,NIMPUR
        PSNOWIMPUR(JJ,JST,JIMP)=0.
    END DO
  END DO  !  end loop unactive snow layers
!
END DO    ! end loop grid points
!
IF(OMEB)THEN
  !
  ! Final adjustment: If using MEB, then all excess heat (above initial guess ground flux)
  ! is used herein to adjust the ground temperature since for the MEB case, boundary
  ! fluxes are imposed (the energy budget of the ground has already been computed).
  ! This correction is needed since snow-soil flux here is implicit, while guess
  ! used in surface energy soil budget was semi-implicit. This forces the flux
  ! seen by the ground to be *the same* as that leaving the snow (energy conservation).
  ! Also, add any excessing cooling from sublimation as snowpack becomes vanishingly thin.
  ! This is added back to total heat content of the snowpack (and distributed among
  ! all snow layers while conserving total heat content plus correction)
  !
  ZWORK(:)    = 0.
  DO JST=1,IMAX_USE
    DO JJ=1,SIZE(ZSNOW)
      IF (JST <= INLVLS_USE(JJ)) THEN
        ZWORK(JJ)        = ZWORK(JJ) + PSNOWHEAT(JJ,JST)
      END IF
    END DO
  END DO
  !
  DO JJ=1,SIZE(ZSNOW)
    ZWORK2(JJ)              = MIN(0.0, ZWORK(JJ) + PGRNDFLUX(JJ) - ZGRNDFLUXI(JJ))
    PGFLXCOR(JJ)            = MAX(0., ZWORK2(JJ)) ! add any possible (rare!) excess to soil
    IF (ZWORK(JJ) > -1.E-10)THEN
      ZWORK(JJ) = 1.  ! i.e. no modifs to H profile
    ELSE
      ZWORK(JJ) = ZWORK2(JJ)/ZWORK(JJ)
    END IF
  END DO
  !
  ! comment from Matthieu : PGFLXCOR do not have an initial value like in ISBA-ES because there is no problem of excess energy due to melt (dealt with SNOWCROLAYER_GONE)
  !
  DO JST=1,IMAX_USE
    DO JJ=1,SIZE(ZSNOW)
      IF (JST <= INLVLS_USE(JJ)) THEN
        PSNOWHEAT(JJ,JST) = PSNOWHEAT(JJ,JST)*ZWORK(JJ)
      END IF
    END DO
  END DO
  !
  ! Here we impose that the implicit snow T profile found in e_budget_meb
  ! be imposed, while conserving the overall heat content of the snowpack
  ! This results in smoother near surface Ts for certain conditions
  ! and more consistency with sfc fluxes. Total heat content of the snowpack
  ! is unchanged/conserved.
  !
  !
  ! Matthieu : je ne comprends rien à cette partie --> attendre réponse d'Aaron
  !
  !
  !    ZWORK2D(:,:)   = PSNOWHEAT(:,:)
  !    DO JJ=1,INLVLS_USE(JJ) ! this can be 1 to INLVLS...
  !       DO JI=1,INI
  !          ZWORK2D(JI,JJ) = MIN(0., PSNOWDZ(JI,JJ)*( ZSCAP(JI,JJ)*(ZSNOWTEMPO(JI,JJ)-XTT)  &
  !                          - XLMTT*PSNOWRHO(JI,JJ) ) + XLMTT*XRHOLW*PSNOWLIQ(JI,JJ) )
  !       END DO
  !    END DO
  !    ZWORK(:)       = 0.0
  !    ZWORK2(:)      = 0.0
  !    DO JJ=1,INLVLS
  !       DO JI=1,INI
  !          ZWORK(JI)   = ZWORK(JI)  + PSNOWHEAT(JI,JJ)
  !          ZWORK2(JI)  = ZWORK2(JI) + ZWORK2D(JI,JJ)
  !       END DO
  !    END DO
  !    WHERE(ZWORK2(:) > -1.E-10)
  !       ZWORK(:) = 1. ! i.e. no modifs to T profile
  !    ELSEWHERE
  !       ZWORK(:) = ZWORK(:)/ZWORK2(:)
  !    END WHERE
  !    DO JJ=1,1,INLVLS
  !       DO JI=1,INI
  !          PSNOWHEAT(JI,JJ) = ZWORK2D(JI,JJ)*ZWORK(JI) ! possibly change heat profile
  !                                                      ! but conserve total
  !       END DO
  !    END DO
  !
  ! ! these are just updated diagnostics at this point:
  !
  !    PSNOWTEMP(:,:) = XTT + ( ((PSNOWHEAT(:,:)/PSNOWDZ(:,:))                   &
  !                     + XLMTT*PSNOWRHO(:,:))/ZSCAP(:,:) )
  !    PSNOWLIQ(:,:)  = MAX(0.0,PSNOWTEMP(:,:)-XTT)*ZSCAP(:,:)*                  &
  !                     PSNOWDZ(:,:)/(XLMTT*XRHOLW)
  !    PSNOWTEMP(:,:) = MIN(XTT,PSNOWTEMP(:,:))
  !
END IF
!
! print some final diagnostics
! ! ! IF (INLVLS_USE(I)>0) THEN
! ! !  WRITE(*,FMT="(I4,2I4,F4.0,A7,F5.2,A10,F7.1,A11,F6.2,A13,F6.2)")  &
! ! !      TPTIME%TDATE%YEAR,TPTIME%TDATE%MONTH,TPTIME%TDATE%DAY,TPTIME%TIME/3600.,&
! ! !           'HTN= ',SUM(PSNOWDZ(I,1:INLVLS_USE(I))), 'FLUX Sol=', PGRNDFLUX(I),&
! ! !   'Tsurf_sol=',PTG(I)-273.16, 'Tbase_neige=',PSNOWTEMP(I,INLVLS_USE(I))-273.16
! ! ! END IF
!
!***************************************PRINT IN*********************************************
! check suspect low temperature
DO JJ = 1,SIZE(ZSNOW)
!IF(INLVLS_USE(JJ)>0) WRITE(*,*) 'SOL:',JJ,INLVLS_USE(JJ),PGRNDFLUX(JJ),PTG(JJ),&
! PSNOWTEMP(jj,INLVLS_USE(JJ)),PSNOWTEMP(jj,1),PZENITH(JJ)
  DO JST = 1,INLVLS_USE(JJ)
    IF ( PSNOWTEMP(JJ,JST) < 100. ) THEN
      WRITE(*,*) 'pb tempe snow :',PSNOWTEMP(JJ,JST)
      WRITE(*,FMT='("DATE:",2(I2.2,"/"),I4.4,F3.0)')          &
        TPTIME%TDATE%DAY,TPTIME%TDATE%MONTH,TPTIME%TDATE%YEAR,TPTIME%TIME/3600.
      WRITE(*,*) 'point',JJ,"/",SIZE(ZSNOW)
      WRITE(*,*) 'layer',JST
      WRITE(*,*) 'pressure',PPS(JJ)
      WRITE(*,*) 'slope',ACOS(PDIRCOSZW(JJ))*(180./XPI),"deg"
      WRITE(*,*) 'XLAT=',PXLAT(JJ),'XLON=',PXLON(JJ)
      WRITE(*,*) 'solar radiation=',PSW_RAD(JJ)
      WRITE(*,*) 'INLVLS_USE(JJ):',INLVLS_USE(JJ)
      WRITE(*,*) 'DZ',PSNOWDZ(JJ,1:INLVLS_USE(JJ))
      WRITE(*,*) 'RHO',PSNOWRHO(JJ,1:INLVLS_USE(JJ))
      WRITE(*,*) 'TEMP',PSNOWTEMP(JJ,1:INLVLS_USE(JJ))
      WRITE(*,*) 'RADSINK',ZRADSINK(JJ,:)
      WRITE(*,*) 'RADXS',ZRADXS(JJ)
      WRITE(*,*) 'ZENITH',PZENITH(JJ)*(180./XPI)
      CALL ABOR1_SFX('SNOWCRO: erreur tempe snow')
    END IF
  END DO
END DO
!***************************************PRINT OUT*********************************************
!***************************************DEBUG IN*********************************************
!
IF (LCRODAILYINFO .OR. LCRODEBUG) THEN
  IF (GCRODEBUGPRINT) THEN
    CALL SNOWCROPRINTDATE()
    CALL SNOWCROPRINTPROFILE("CROCUS : end of time step",INLVLS_USE(IDEBUG),LPRINTGRAN, &
          PSNOWDZ(IDEBUG,:),PSNOWRHO(IDEBUG,:),PSNOWLIQ(IDEBUG,:),  &
          PSNOWHEAT(IDEBUG,:),PSNOWDIAMOPT(IDEBUG,:),PSNOWSPHERI(IDEBUG,:),                &
          PSNOWHIST(IDEBUG,:),PSNOWAGE(IDEBUG,:), HSNOWMETAMO)
  END IF
  !Control and print energy balance
  IF (GCRODEBUGPRINTBALANCE) THEN
    !
    ZSUMMASS_FIN(IDEBUG) = SUM( PSNOWSWE (IDEBUG,1:INLVLS_USE(IDEBUG)) )
    ZSUMHEAT_FIN(IDEBUG) = SUM( PSNOWHEAT(IDEBUG,1:INLVLS_USE(IDEBUG)) )
    !
    CALL GET_BALANCE(ZSUMMASS_INI(IDEBUG),ZSUMHEAT_INI(IDEBUG),ZSUMMASS_FIN(IDEBUG), &
                     ZSUMHEAT_FIN(IDEBUG),ZSR(IDEBUG),PRR(IDEBUG),PTHRUFAL(IDEBUG),  &
                     PEVAP(IDEBUG),PEVAPCOR(IDEBUG),PGRNDFLUX(IDEBUG),PHSNOW(IDEBUG),&
                     PRNSNOW(IDEBUG),PLEL3L(IDEBUG),PLES3L(IDEBUG),PHPSNOW(IDEBUG),  &
                     PSNOWHMASS(IDEBUG),PSNOWDZ(IDEBUG,1),PTSTEP,                    &
                     ZMASSBALANCE(IDEBUG),ZENERGYBALANCE(IDEBUG),ZEVAPCOR2(IDEBUG)   )
    !
    CALL SNOWCROPRINTBALANCE(ZSUMMASS_INI(IDEBUG),ZSUMHEAT_INI(IDEBUG),ZSUMMASS_FIN(IDEBUG), &
                             ZSUMHEAT_FIN(IDEBUG),ZSR(IDEBUG),PRR(IDEBUG),PTHRUFAL(IDEBUG),  &
                             PEVAP(IDEBUG),PEVAPCOR(IDEBUG),PGRNDFLUX(IDEBUG),PHSNOW(IDEBUG),&
                             PRNSNOW(IDEBUG),PLEL3L(IDEBUG),PLES3L(IDEBUG),PHPSNOW(IDEBUG),  &
                             PSNOWHMASS(IDEBUG),PSNOWDZ(IDEBUG,1),PTSTEP,                    &
                             ZMASSBALANCE(IDEBUG),ZENERGYBALANCE(IDEBUG),ZEVAPCOR2(IDEBUG))
    !
  END IF
  !
  IF (LPSTOPBALANCE) THEN
    !
    ! bilan pour tous points pour test eventuel sur depassement seuil des residus
    DO JJ=1, SIZE(ZSNOW)
      !
      ZSUMMASS_FIN(JJ) = SUM( PSNOWSWE (JJ,1:INLVLS_USE(JJ)) )
      ZSUMHEAT_FIN(JJ) = SUM( PSNOWHEAT(JJ,1:INLVLS_USE(JJ)) )
      !
      CALL GET_BALANCE(ZSUMMASS_INI(JJ),ZSUMHEAT_INI(JJ),ZSUMMASS_FIN(JJ), &
                       ZSUMHEAT_FIN(JJ),ZSR(JJ),PRR(JJ),PTHRUFAL(JJ),      &
                       PEVAP(JJ),PEVAPCOR(JJ),PGRNDFLUX(JJ),PHSNOW(JJ),    &
                       PRNSNOW(JJ),PLEL3L(JJ),PLES3L(JJ),PHPSNOW(JJ),      &
                       PSNOWHMASS(JJ),PSNOWDZ(JJ,1),PTSTEP,                &
                       ZMASSBALANCE(JJ),ZENERGYBALANCE(JJ),ZEVAPCOR2(JJ)   )
      !
    END DO    ! end loop grid points
    !
    CALL SNOWCROSTOPBALANCE(ZMASSBALANCE(:),ZENERGYBALANCE(:))
    !
  END IF
END IF
!***************************************DEBUG OUT********************************************
!
IF (OMEB) THEN
  DO JJ=1, SIZE(ZSNOW)
    IF (INLVLS_USE(JJ) >= 1) THEN
      ZWORK(JJ) = PSNOWTEMP(JJ,1)
    ELSE
      ZWORK(JJ) = XTT ! melting point assumed if snow has disappeared
    END IF
  END DO
  PQS(:) = QSATI(ZWORK(:), PPS(:))
ELSE
  PQS(:) = ZQSAT(:)
END IF
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO',1,ZHOOK_HANDLE)
!
CONTAINS
!
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOWCROCOMPACTN(PTSTEP,PSNOWRHO,PSNOWDZ,                         &
                           PSNOWTEMP,PSNOW,PSNOWDIAMOPT,PSNOWSPHERI,PSNOWHIST, &
                           PSNOWSWE, PSNOWAGE, PSNOWDSSA, PSNOWLIQ,         &
                           INLVLS_USE,PDIRCOSZW,                            &
                           HSNOWMETAMO,HSNOWCOMP, OSNOWCOMPACT_BOOL,        &
                           PHVEGPOL)                    
!
!!    PURPOSE
!!    -------
!     Snow compaction due to overburden and settling.
!     Mass is unchanged: layer thickness is reduced
!     in proportion to density increases. Method
!     directly inherited from Crocus v2.4 and
!     coarsely described in Brun et al., J. Glac 1989 and 1992
!
!     de/e = -sigma/eta * dt
!
!     where e is layer thickness, sigma is the vertical stress, dt is the
!     time step and eta is the snow viscosity
!     * sigma is currently calculated taking into account only the overburden
!     (a term representing "metamorphism stress" in fresh snow may be added
!      in the future)
!     * eta is computed as a function of snowtype, density and temperature
!
!     The local slope is taken into account, through the variable PDIRCOSZW
!     which is directly the cosine of the local slope
!
!
!     HISTORY:
!     Basic structure from ISBA-ES model (Boone and Etchevers, 2001)
!     Implementation of Crocus laws : E. Brun, S. Morin, J.-M. Willemet July 2010.
!     Implementation of slope effect on settling : V. Vionnet, S. Morin May 2011
!
!
!Comments by P.Spandre (30/10/2013) from descriptions in (ref):
!The detailed snowpack scheme Crocus and its implementation in SURFEX v7.2
!by V.Vionnet, 2012

USE MODD_CSTS,     ONLY : XTT, XG
USE MODD_SNOW_PAR, ONLY : XRHOSMAX_ES
USE MODD_SNOW_METAMO
USE MODE_SNOW3L, ONLY : GETGRAINSIZE_B21
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP       ! Time step UNIT : s
REAL, DIMENSION(:), INTENT(IN)      :: PDIRCOSZW    ! cosine of local slope
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWTEMP    ! Snow temperature UNIT : K
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWRHO, PSNOWDZ   ! Density UNIT : kg m-3, Layer thickness UNIT : m
!
REAL, DIMENSION(:), INTENT(OUT)     :: PSNOW        ! Snowheight UNIT : m
!

REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDIAMOPT, PSNOWSPHERI, PSNOWHIST
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWDSSA    !Snowtype variables
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWLIQ     ! Snow liquid water content UNIT ???
INTEGER, DIMENSION(:), INTENT(IN)   :: INLVLS_USE   ! Number of snow layers used
CHARACTER(3), INTENT(IN)            :: HSNOWMETAMO  ! metamorphism scheme
CHARACTER(3), INTENT(IN)            :: HSNOWCOMP    ! compaction option
!
REAL, DIMENSION(:), INTENT(IN)      :: PHVEGPOL        ! Mean polar low vegetation height (m)
!
LOGICAL, INTENT(IN)                 :: OSNOWCOMPACT_BOOL   ! 20160211
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWAGE, PSNOWSWE  ! Age and SWE of snow layers
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWRHO2,    &! work snow density UNIT : kg m-3
                                                      ZVISCOSITY,   &! Snow viscosity UNIT : N m-2 s (= Pa s)
                                                      ZSMASS,       &!overburden mass for a given layer UNIT : kg m-2
                                                      ZSNOW_JST,    &!Snow layer height (Royer et al 2021) 
                                                      ZSMASSCOEFF    ! Coefficient for extra pressure due to grooming (p.s 20150721)
!
REAL,PARAMETER     :: PPK=0.18
REAL,PARAMETER     :: PPB=-6.6E-3
!
REAL,PARAMETER     :: PPA0=0.4
!
INTEGER            :: JJ,JST   ! looping indexes
REAL               :: ZSNOWSIZE ! Size of grain snow in non-dendritic case
!
REAL(KIND=JPRB)    :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 1. Cumulative snow mass (kg/m2):
! --------------------------------
!
! P.S   (ref) p.779 Vertical stress from the weight of the overlying layers for each layer JST (6)
!   ZSMASS(JST)=sigma(i)/gcos@
!   Uppermost layer JST=1 on top (air surface)
!   Loaded with half of its own weight
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROCOMPACTN',0,ZHOOK_HANDLE)
!
ZSMASS(:,1) = 0.0
DO JST = 2,IMAX_USE
  DO JJ = 1,SIZE(PSNOW)
    IF ( JST<=INLVLS_USE(JJ)) THEN
      ZSMASS(JJ,JST) = ZSMASS(JJ,JST-1) + PSNOWDZ(JJ,JST-1) * PSNOWRHO(JJ,JST-1)
    END IF
  END DO
END DO
ZSMASS(:,1) = 0.5 * PSNOWDZ(:,1) * PSNOWRHO(:,1)  ! overburden of half the mass of the uppermost layer applied to itself
!
! 2. Compaction/Settling
! ----------------------

DO JJ = 1,SIZE(PSNOW)
  !
  ! Calculate height of each snow layer (Royer 2021)
  ZSNOW_JST(JJ,1) = PSNOW(JJ) - 0.5*PSNOWDZ(JJ,1)
  DO JST=2,INLVLS_USE(JJ)
     ZSNOW_JST(JJ,JST) = ZSNOW_JST(JJ,JST-1) - 0.5 * (PSNOWDZ(JJ,JST-1) + PSNOWDZ(JJ,JST))
  ENDDO
ENDDO
!
IF ((HSNOWCOMP=="S14")) THEN
  DO JST = 1,IMAX_USE
    DO JJ = 1,SIZE(PSNOW)
      IF ( JST<=INLVLS_USE(JJ)) THEN
        IF ((PSNOWAGE(JJ,JST)<=2.) .AND. (PSNOWRHO(JJ,JST) < XRHOTHRESHOLD_ICE))  THEN
          ZSNOWRHO2(JJ,JST) = PSNOWRHO(JJ,JST) + PSNOWRHO(JJ,JST) * PPB * PSNOWDSSA(JJ,JST) * &
                                             ( XG*PDIRCOSZW(JJ)*ZSMASS(JJ,JST))**PPK
        ELSE
          ZVISCOSITY(JJ,JST) = XVVISC1 * &
                           EXP( XVVISC3*PSNOWRHO(JJ,JST) + XVVISC4*ABS(XTT-PSNOWTEMP(JJ,JST)) ) * &
                           PSNOWRHO(JJ,JST) / XVRO11
        END IF
      END IF
    END DO
  END DO
ELSEIF (HSNOWCOMP=="T11") THEN
  DO JST = 1,IMAX_USE
    DO JJ = 1,SIZE(PSNOW)
      IF ( JST<=INLVLS_USE(JJ)) THEN
        ZVISCOSITY(JJ,JST) = 0.05 * PSNOWRHO(JJ,JST)**(-0.0371*(PSNOWTEMP(JJ,JST)-XTT)+4.4)*(1E-4*EXP(0.018*PSNOWRHO(JJ,JST))+1)
      END IF
    END DO
  END DO
ELSE
  DO JST = 1,IMAX_USE
    DO JJ = 1,SIZE(PSNOW)
      IF ( JST<=INLVLS_USE(JJ)) THEN
           ! Snow viscosity basic equation (depend on temperature and density only):
           ZVISCOSITY(JJ,JST) = XVVISC1 * &
                           EXP( XVVISC3*PSNOWRHO(JJ,JST) + XVVISC4*ABS(XTT-PSNOWTEMP(JJ,JST)) ) * &
                           PSNOWRHO(JJ,JST) / XVRO11
      END IF
    END DO
  END DO
END IF


DO JST = 1,IMAX_USE
  DO JJ = 1,SIZE(PSNOW)
  !
    IF ( JST<=INLVLS_USE(JJ)) THEN
      ! Equations below apply changes to the basic viscosity value, based on snow microstructure properties
      IF ((HSNOWCOMP .EQ. "S14") .AND. PSNOWAGE(JJ,JST)<=2. .AND. (PSNOWRHO(JJ,JST) < XRHOTHRESHOLD_ICE)) CYCLE
      IF ( PSNOWLIQ(JJ,JST)>0.0 ) THEN
        ZVISCOSITY(JJ,JST) = ZVISCOSITY(JJ,JST) / &
                             ( XVVISC5 + XVVISC6*PSNOWLIQ(JJ,JST)/PSNOWDZ(JJ,JST) )
      END IF
      !
      ! 25 Aug 2025 (crocus3.0) Remove this unpublished parameterization slowing the compaction of refrozen crusts
      ! which is not supported by scientific literature and has detrimental impacts on modern diagnostics of mechanical stability
      ! Note that MEPRA indexes are affected by this modification (lower values)
      ! Note also that this parameterization was also removed in the work of Quéno et al. 2018.
      ! Also impacts on snow depths scores.
      !IF( PSNOWLIQ(JJ,JST)/PSNOWDZ(JJ,JST)<=0.01 .AND. NINT(PSNOWHIST(JJ,JST))>=NVHIS2 ) THEN
      !  ZVISCOSITY(JJ,JST) = ZVISCOSITY(JJ,JST) * XVVISC7
      !END IF
      !
      IF ( NINT(PSNOWHIST(JJ,JST))==NVHIS1 ) THEN
        !
        IF ( PSNOWDIAMOPT(JJ,JST)>=XVDIAM6*(4.-PSNOWSPHERI(JJ,JST)) .AND. PSNOWSPHERI(JJ,JST)<XVGRAN6 ) THEN
            ! non dendritic case
            CALL GETGRAINSIZE_B21(PSNOWDIAMOPT(JJ,JST),PSNOWSPHERI(JJ,JST),ZSNOWSIZE)

          ZVISCOSITY(JJ,JST) = ZVISCOSITY(JJ,JST) * &
                               MIN( 4., EXP( MIN( XVDIAM1, &
                                                 ZSNOWSIZE-XVDIAM4 ) / XVDIAM6 ) )
        END IF
        !
      END IF
      ! Increase snow viscosity for snow layer height <= vegetation threshold / M. Barrere
      IF (HSNOWCOMP == 'R21' .OR. HSNOWCOMP == 'R2V' ) THEN
        IF ( PSNOWLIQ(JJ,JST)<=XUEPSI_SMP .AND. PHVEGPOL(JJ) > 0. ) THEN ! only for dry snow layers when shrubs are present (PHVEGPOL>0.)
          IF(ZSNOW_JST(JJ,JST) <= MIN(0.1,PHVEGPOL(JJ) )) THEN
             ZVISCOSITY(JJ,JST) = 100. * ZVISCOSITY(JJ,JST)
          ELSE IF ( ZSNOW_JST(JJ,JST) <= PHVEGPOL(JJ)  ) THEN
             ZVISCOSITY(JJ,JST) = 10. * ZVISCOSITY(JJ,JST)
          ENDIF
        ENDIF
      ENDIF
      !      
      !
      ! Calculate new snow snow density: compaction from weight/over-burden
      ZSNOWRHO2(JJ,JST) = PSNOWRHO(JJ,JST) + PSNOWRHO(JJ,JST) * PTSTEP * &
                                           ( XG*PDIRCOSZW(JJ)*ZSMASS(JJ,JST)/ZVISCOSITY(JJ,JST) )

    END IF
  END DO
END DO
!
DO JST = 1,IMAX_USE
  DO JJ = 1,SIZE(PSNOW)
    IF ( JST<=INLVLS_USE(JJ)) THEN
      !
      !!print*, JJ,JST,  PSNOWAGE(JJ,JST), ZSNOWRHO2(JJ,JST),PSNOWDZ(JJ,JST)
      ! Calculate new grid thickness in response to density change
      PSNOWDZ(JJ,JST) = PSNOWDZ(JJ,JST) * ( PSNOWRHO(JJ,JST)/ZSNOWRHO2(JJ,JST) )
      !
      !  Update density (kg m-3):
      PSNOWRHO(JJ,JST) = ZSNOWRHO2(JJ,JST)
      !
    END IF
    !
  END DO    ! end loop snow layers
  !
END DO      ! end loop grid points
!
!
! 3. Update total snow depth:
! -----------------------------------------------
!

IF (OSNOWCOMPACT_BOOL) THEN
  DO JST = 1,IMAX_USE
    DO JJ = 1,SIZE(PSNOW)
    !
      IF ( JST<=INLVLS_USE(JJ)) THEN
        ZSMASSCOEFF(JJ,JST) = PTSTEP*XG*PDIRCOSZW(JJ)/ZVISCOSITY(JJ,JST)
      END IF
    END DO
  END DO
  CALL SNOWGROOMING(PSNOWDZ,PSNOWSWE,PSNOWAGE,PSNOWRHO,PSNOWDIAMOPT,&
                    PSNOWSPHERI,INLVLS_USE,ZSMASSCOEFF, &
                    OSNOWMAK_BOOL, OSNOWTILLER)
END IF
!
DO JJ = 1,SIZE(PSNOWDZ,1)
  PSNOW(JJ) = SUM( PSNOWDZ(JJ,1:INLVLS_USE(JJ)) )
END DO
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROCOMPACTN',1,ZHOOK_HANDLE)

!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOWCROCOMPACTN


!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOWCROMETAMO(PSNOWDZ,PSNOWDIAMOPT, PSNOWSPHERI,         &
                         PSNOWHIST, PSNOWTEMP, PSNOWLIQ, PTSTEP, &
                         PSNOWSWE,INLVLS_USE, PSNOWAGE, HSNOWMETAMO)
!

!     SNOW METAMORPHISM
!     HISTORICAL VARIABLES
!
!     GRAIN METAMORPHISM ACCORDING TO BRUN ET AL (1992)
!     THE DIFFERENT CASES ARE :
!     1.2 WET SNOW
!     1.3 DRY SNOW
!       1.3.1. LOW      TEMPERATURE GRADIENT
!       1.3.2. MODERATE TEMPERATURE GRADIENT
!       1.3.3. HIGH     TEMPERATURE GRADIENTi
!     THE CASE OF DENTRITIC OR NON DENDRITIC SNOW IS TREATED SEPARATELY
!     THE LIMIT DENTRITIC ==> NON DENDRITIC IS REACHED WHEN SGRAN1>0

!     SNOW SETTLING : VISCOSITY DEPENDS ON THE GRAIN TYPES

!     HISTORICAL VARIABLES (NON DENDRITIC CASE)
!     MSHIST DEFAUT
!        0           CAS NORMAL
!     NVHIS1   1     FACETED CRISTAL
!     NVHIS2   2     LIQUID WATER AND NO FACETED CRISTALS BEFORE
!     NVHIS3   3     LIQUID WATER AND FACETED CRISTALS BEFORE

!     EXTERNES.
!     ---------

!     REFERENCES.
!     -----------

!     AUTEURS.
!     --------
!        ERIC BRUN ET AL. - JOURNAL OF GLACIOLOGY 1989/1992.

!     MODIFICATIONS.
!     --------------
!        08/95: YANNICK DANIELOU - CODAGE A LA NORME DOCTOR.
!        09/96: ERIC MARTIN      - CORRECTION COMMENTAIRES
!        03/06: JM Willemet      - F90 and SI units
!        08/06: JM Willemet      - new formulation for TEL (Mwat/(Mice+Mwat) instead of Mwat/Mice.
!                                  Threshold on the diameter increasing of the wet grains.
!        01/07 : JM Willemet     - CORRECTION DES COUCHES SATUREES SUBISSANT DU TASSEMENT
!                                  CORRECTION ON THE SATURATED LAYERS WHICH ARE SETTLED
!        12/12: CM Carmagnola    - Dendricity and size replaced by the optical diameter
!                                - Test of different evolution laws for the optical diameter
!        08/13: M Lafaysse       - Simplification of historical parameter computation (logicals GNONDENDRITIC, GFACETED, GSPHE_LW)
!        04/21: M Baron          - new option 'B21' for metamorphism, as a proposition of correction of C13
!        04/23: M Lafaysse       - conversion of historical param in integer to remove instabilities
!                                  due to inappropriate comparisons real / integer
USE MODD_SNOW_METAMO
USE MODD_CSTS, ONLY : XTT, XPI, XRHOLW, XRHOLI
USE MODD_SURF_PAR, ONLY : XUNDEF
!
USE MODE_SNOW3L
!
IMPLICIT NONE
!
!     0.1 declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWDZ, PSNOWTEMP, PSNOWLIQ, PSNOWSWE
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDIAMOPT, PSNOWSPHERI, PSNOWHIST
!
REAL, INTENT(IN)                    :: PTSTEP
!
INTEGER, DIMENSION(:), INTENT(IN)   :: INLVLS_USE
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWAGE
!
CHARACTER(3), INTENT(IN)              :: HSNOWMETAMO ! metamorphism scheme
!
!     0.2 declaration of local variables
!
INTEGER, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ISNOWHIST ! historical value stored as integer
REAl :: ZGRADT(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)),ZTELM(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2))
!       ZGRADT = temperature gradient in the layer ( K/m )
REAL :: ZSPHE(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)),ZVAP,ZSNOWSIZE
!       ZSPHE = sphericity ( intermediate variable for calculations )
!       ZVAP = coeff in the rate of evolution of microstructural variables
!       ZSNOWSIZE = snow grain size in the 'Brun92' formalism ( m )
REAL :: ZDANGL, ZSSA, ZSSA0, ZA, ZB, ZC, &
        ZA2, ZB2, ZC2, ZOPTR, ZOPTR0, ZDRDT, ZDSPHESURDT
REAL :: ZVDENT1 (SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)), ZVDENT2(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)),&
         ZCOEF_SPH(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2))
REAL :: ZDENOM1, ZDENOM2, ZFACT1, ZFACT2
REAL :: Z4PI ,ZONETHIRD,ZSIXOFXRHOLI,ZTSTEPHOUR !RAFIFE RENOMMER
INTEGER :: JST,JJ                                !Loop controls
INTEGER :: IDRHO, IDGRAD, IDTEMP           !Indices for values from Flanner 2006
LOGICAL :: GNONDENDRITIC (SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2))
LOGICAL :: GCOND_B21(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2))
LOGICAL :: GCOND_SPH(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2))
!          GCOND_SPH = logical testing if sphericity is allowed to move
LOGICAL :: GCOND_DEND
!          GCOND_DEND = logical testing if snow is dendritic
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!     INITIALISATION
!     --------------
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROMETAMO',0,ZHOOK_HANDLE)
!
! Convert historical values in integer for exact comparisons below
! Strict comparisons between real and integers are extremely dangerous
ISNOWHIST(:,:) = NINT(PSNOWHIST(:,:))
!
ZTELM=0.
ZDANGL=0.
ZCOEF_SPH = 2.
Z4PI=4.*XPI
ZONETHIRD=1./3.
ZTSTEPHOUR=PTSTEP/3600.
ZSIXOFXRHOLI=6./XRHOLI
!
GNONDENDRITIC = .FALSE.
!
!*    1. METAMORPHOSES DANS LES STRATES. / METAMORPHISM
!        -----------------------------------------------
    ! 1.1 INITIALISATION: GRADIENT DE TEMPERATURE / TEMPERATURE GRADIENT
DO JJ = 1,SIZE(PSNOWRHO,1)
  !
  !FOR JST =INLVLS_USE(JJ)
  IF (INLVLS_USE(JJ) == 0) CYCLE
  
  IF (INLVLS_USE(JJ) > 1) THEN
    ZGRADT(JJ,INLVLS_USE(JJ)) = ABS(PSNOWTEMP(JJ,INLVLS_USE(JJ))         &
    - PSNOWTEMP(JJ,INLVLS_USE(JJ)-1))*2. / (PSNOWDZ(JJ,INLVLS_USE(JJ)-1) &
    + PSNOWDZ(JJ,INLVLS_USE(JJ)))
  ENDIF

  !FOR JST =1
  ZGRADT(JJ,1) = ABS(PSNOWTEMP(JJ,1+1) - PSNOWTEMP(JJ,1)  )*2. / (PSNOWDZ(JJ,1) + PSNOWDZ(JJ,1+1))
  !
  !FOR JST ELSE
  DO JST = 2,IMAX_USE
    IF ( JST<=INLVLS_USE(JJ)-1) THEN
      ZGRADT(JJ,JST) = ABS(PSNOWTEMP(JJ,JST+1) - PSNOWTEMP(JJ,JST-1))*2. / &
                (PSNOWDZ(JJ,JST-1) + PSNOWDZ(JJ,JST)*2. + PSNOWDZ(JJ,JST+1))
    END IF
  END DO
END DO
!
!======================Version of metamorphism proposed by M.Baron=============
  !Evolution of sphericity
  DO JST = 1,IMAX_USE
    !
    DO JJ = 1,SIZE(PSNOWRHO,1)
      !
      IF (INLVLS_USE(JJ) == 0) CYCLE
      IF ( JST<=INLVLS_USE(JJ)) THEN
        !
        IF ( PSNOWLIQ(JJ,JST)>XUEPSI ) THEN
          ! 1.2 METAMORPHOSE HUMIDE. / WET SNOW METAMORPHISM
          !
          ! TENEUR EN EAU LIQUIDE / LIQUID WATER CONTENT
          ZTELM(JJ,JST)  = XUPOURC * PSNOWLIQ(JJ,JST) * XRHOLW / PSNOWSWE(JJ,JST)
          !
          ! VITESSES DE DIMINUTION DE LA DENDRICITE / RATE OF THE DENDRICITY DECREASE
          ZVDENT1(JJ,JST) = MAX( XVDENT2 * ZTELM(JJ,JST)**NVDENT1, XVDENT1 * EXP(XVVAP1/XTT) )
          ZVDENT2(JJ,JST) = ZVDENT1(JJ,JST)
          !
          GCOND_B21(JJ,JST) = .TRUE.  !CONDITION POUR LE CALCUL DE SNOWDIAMOPT EN UTILISANT B21
          !
          ZSPHE(JJ,JST) = PSNOWSPHERI(JJ,JST) + ZVDENT2(JJ,JST) * PTSTEP
          ZSPHE(JJ,JST) = MIN(1.,ZSPHE(JJ,JST))
          !
          GCOND_SPH(JJ,JST) = ( ZSPHE(JJ,JST) < 1.-XUEPSI )
          !
        ELSE
          !
          ! CONDITION POUR LE CALCUL DE SNOWDIAMOPT EN UTILISANT B21
          GCOND_B21(JJ,JST) = ( HSNOWMETAMO=='B21' ) .OR. ((HSNOWMETAMO=='S-B') .AND. (PSNOWAGE(JJ,JST)>2.0))
          !
          IF ( ZGRADT(JJ,JST)<XVGRAT1 ) THEN
          ! 1.3.1 METAMORPHOSE SECHE FAIBLE/ DRY LOW GRADIENT (0-5 DEG/M).
          !
            ZVAP = EXP( XVVAP1/PSNOWTEMP(JJ,JST) )
            !
            ! VITESSES DE DIMINUTION DE LA DENDRICITE / RATE OF THE DENDRICITY DECREASE
            ZVDENT1(JJ,JST) = XVDENT1 * ZVAP
            CALL GETGRAINSIZE_B21(PSNOWDIAMOPT(JJ,JST),PSNOWSPHERI(JJ,JST),ZSNOWSIZE)
            IF (ISNOWHIST(JJ,JST) == NVHIS1 .AND. ZSNOWSIZE > XVDIAM2) THEN !(rq : si dendritique, on calcule un snowsize=0 donc FALSE)
            !BIG FACETED GRAINS => SPHERICITY GROWTH IS LIMITED ( source : old code of Crocus, unpublished, probably a "hand
            !adjustment" supposed to prevent the unphysical evolution "DH" => "RG")

              ZVDENT2(JJ,JST)=XVSPHE2*ZVAP*EXP(MIN(0.,XVDIAM3-ZSNOWSIZE)/XVDIAM6)
              ZSPHE(JJ,JST)=PSNOWSPHERI(JJ,JST)+ZVDENT2(JJ,JST)*PTSTEP
              ZSPHE(JJ,JST)=MIN(XVSPHE3,ZSPHE(JJ,JST))!sphericity limited to 0.5
              GCOND_SPH(JJ,JST)=(ZSPHE(JJ,JST)<XVSPHE3-XUEPSI)
              ZCOEF_SPH(JJ,JST) = 2.5
            ELSE
              ZVDENT2(JJ,JST) = XVSPHE2 * ZVAP
              ! CONDITION POUR LE CAS NON DENDRITIQUE SPHERICITE NON LIMITEE

              ZSPHE(JJ,JST) = PSNOWSPHERI(JJ,JST) + ZVDENT2(JJ,JST) * PTSTEP
              ZSPHE(JJ,JST) = MIN(1.,ZSPHE(JJ,JST))
            !
              GCOND_SPH(JJ,JST) = ( ZSPHE(JJ,JST) < 1.-XUEPSI )
            END IF
          !
          ELSE
            ! 1.3.2 METAMORPHOSE SECHE GRADIENT MOYEN / DRY MODERATE (5-15).
            ! 1.3.3 METAMORPHOSE SECHE FORT / DRY HIGH GRADIENT
            !
            ZVAP = EXP( XVVAP1/PSNOWTEMP(JJ,JST) ) * (ZGRADT(JJ,JST))**XVVAP2
            !
            ! VITESSES DE DIMINUTION DE LA DENDRICITE / RATE OF THE DENDRICITY DECREASE
            ZVDENT1(JJ,JST) = XVDENT1 * ZVAP
            ZVDENT2(JJ,JST) = - XVDENT1 * ZVAP
            ! CONDITION POUR LE CAS NON DENDRITIQUE NON COMPLETEMENT ANGULEUX
            ZCOEF_SPH(JJ,JST) = 3.
            !
            ZSPHE(JJ,JST) = PSNOWSPHERI(JJ,JST) + ZVDENT2(JJ,JST) * PTSTEP
            ZSPHE(JJ,JST) = MAX(0.,ZSPHE(JJ,JST))
            !
            GCOND_SPH(JJ,JST) = ( ZSPHE(JJ,JST) > XUEPSI )
          END IF
        END IF
      END IF
    END DO
  END DO
! Evolution of optical diameter for all cases with B21 options and for wet metamorphism with FO6, S-F, T07 options
  DO JST = 1,IMAX_USE
    !
    DO JJ = 1,SIZE(PSNOWRHO,1)
      !
      IF (INLVLS_USE(JJ) == 0) CYCLE
      IF ( JST<=INLVLS_USE(JJ)) THEN
        !
        IF (GCOND_B21(JJ,JST)) THEN
          CALL CHECK_DENDRITIC(PSNOWDIAMOPT(JJ,JST),PSNOWSPHERI(JJ,JST),GCOND_DEND)
          !
        END IF
        IF (GCOND_B21(JJ,JST) .AND. GCOND_DEND) THEN
          !EVOLUTION DU DIAMETRE OPTIQUE DANS LE CAS DENDRITIQUE ( based on work published in Carmagnola et al. 2014 )
          IF ( GCOND_SPH (JJ,JST)) THEN
            !case when sphericity can evolve
            
            ! gcond_sph <=> sph == 1 with liquid water (Eq. 52 second line)
            
            ! This formula corresponds to Eq 47 and Eq 52 for dendritic case because ZVDENT2/ZVDENT1 includes the differences
            ! from the change in the evolution of sphericity as computed above.
            ! liquid case --> ZVDENT2 = ZVDENT1 which comes back to Eq 48.
            
            PSNOWDIAMOPT(JJ,JST) = PSNOWDIAMOPT(JJ,JST) +XVDIAM6 * PTSTEP * &
                                 ( ZVDENT2(JJ,JST)*((PSNOWDIAMOPT(JJ,JST)/XVDIAM6)-1.)/(ZSPHE(JJ,JST)-3.) - &
                                   ZVDENT1(JJ,JST)*(ZSPHE(JJ,JST)-3.) )
          ELSE
            !case when sphericity is stuck on the borders
            ! This is a simplification of Eq. 47 and Eq.52 for the dendritic case because delta S = 0 and S is known
            ! (0, 0.5 or 1 as accounted for by ZCOEF_SPH)
            PSNOWDIAMOPT(JJ,JST) =  PSNOWDIAMOPT(JJ,JST) + XVDIAM6 * PTSTEP * ZVDENT1(JJ,JST) * ZCOEF_SPH(JJ,JST)
          END IF
          !
        ELSEIF (GCOND_B21(JJ,JST) .AND. GCOND_SPH(JJ,JST) ) THEN
          !OPTICAL DIAMETER EVOLUTION IN THE NON-DENDRITIC CASE, WITH SPHERICITY ALLOWED TO EVOLVE
          !( application of the partial derivatives formula with the formulation proposed by M. Baron )
          
          ! Eq 47 (sphericity >0 for dry non-dendritic snow) or Eq 52 depending on ZVDENT2 (sphericity < 1 with liquid water)
          PSNOWDIAMOPT(JJ,JST)=MAX(PSNOWDIAMOPT(JJ,JST)+PTSTEP*ZVDENT2(JJ,JST)*(PSNOWDIAMOPT(JJ,JST)-4*XVDIAM6)&
                  /(1+ZSPHE(JJ,JST)),XVDIAM6*(4.-ZSPHE(JJ,JST)))
        !
        ELSE IF (PSNOWLIQ(JJ,JST) > XUEPSI) THEN
          !CASE WHEN SPHERICITY=1 AND WET SNOW ( Last line of Eq 52)
          !
          ! old incorrect formulation
          ! PSNOWDIAMOPT(JJ,JST) = 2. * ( 3./(Z4PI) * &
          ! ( Z4PI/3. * (PSNOWDIAMOPT(JJ,JST)/2.)**3. + &
          ! ( XVTAIL1 + XVTAIL2 * ZTELM(JJ,JST)**NVDENT1 ) * PTSTEP ) )**(ZONETHIRD)
          !
          !
          PSNOWDIAMOPT(JJ,JST) = PSNOWDIAMOPT(JJ,JST) + (2./XPI) * &
          ( XVTAIL1 + XVTAIL2 * ZTELM(JJ,JST)**NVDENT1 ) / (PSNOWDIAMOPT(JJ,JST) )**2 * PTSTEP
          !
        ELSEIF (GCOND_B21(JJ,JST) .AND. (ZGRADT(JJ,JST) >= XVGRAT2)) THEN
          !CASE WHEN DRY SNOW AND TEMPERATURE GRADIENT >= 15 K/m AND S=0
          ! Cases with S>0 i.e. GCOND_SPH = True are dealt before
          ZDANGL = SNOW3L_MARBOUTY(PSNOWRHO(JJ,JST),PSNOWTEMP(JJ,JST),ZGRADT(JJ,JST))
          PSNOWDIAMOPT(JJ,JST) = PSNOWDIAMOPT(JJ,JST) + 0.5 * ZDANGL * XVFI * PTSTEP
        END IF
        PSNOWSPHERI(JJ,JST)=ZSPHE(JJ,JST)
      END IF
    END DO
  END DO
  !
    !---------------------------------
    !    TAILLANDIER et al. 2007 (T07)
    !
    ! -> Dry snow
    ! -> Evolution of optical diameter
    !---------------------------------
    !
IF ( HSNOWMETAMO=='T07' ) THEN
  !
!   WRITE(*,*) CSNOWMETAMO,': you are using T07 formulation!!'
  !
  DO JST = 1,IMAX_USE
    !
    DO JJ = 1,SIZE(PSNOWRHO,1)
      !
      IF (INLVLS_USE(JJ) == 0) CYCLE
      IF ( JST<=INLVLS_USE(JJ) .AND. PSNOWLIQ(JJ,JST)<=XUEPSI) THEN
          ! Coefficients from Taillander et al. 2007
          ZSSA0 = 6./( XRHOLI*XVDIAM6 ) * 10.
          !
          ZA  =  0.659*ZSSA0 - 27.2 * ( PSNOWTEMP(JJ,JST)-273.15-2.03 )      ! TG conditions
          ZB  = 0.0961*ZSSA0 - 3.44 * ( PSNOWTEMP(JJ,JST)-273.15+1.90 )
          ZC  = -0.341*ZSSA0 - 27.2 * ( PSNOWTEMP(JJ,JST)-273.15-2.03 )
          ZA2 =  0.629*ZSSA0 - 15.0 * ( PSNOWTEMP(JJ,JST)-273.15-11.2 )     ! ET conditions
          ZB2 = 0.0760*ZSSA0 - 1.76 * ( PSNOWTEMP(JJ,JST)-273.15-2.96 )
          ZC2 = -0.371*ZSSA0 - 15.0 * ( PSNOWTEMP(JJ,JST)-273.15-11.2 )
          !
          ! Compute SSA (rate equation with Taylor series)
          ZSSA = 6./( XRHOLI*PSNOWDIAMOPT(JJ,JST) ) * 10.
          !
          ZDENOM1 = (PSNOWAGE(JJ,JST)*24.) + EXP(ZC/ZB)
          ZDENOM2 = (PSNOWAGE(JJ,JST)*24.) + EXP(ZC2/ZB2)
          ZFACT1 = 0.5 + 0.5*TANH( 0.5*(ZGRADT(JJ,JST)-10.) )
          ZFACT2 = 0.5 - 0.5*TANH( 0.5*(ZGRADT(JJ,JST)-10.) )
          ZSSA = ZSSA + (ZTSTEPHOUR) * (                   ZFACT1 * (-ZB/ZDENOM1)    + ZFACT2 * (-ZB2/ZDENOM2)   + &
                                          (ZTSTEPHOUR) * ( ZFACT1 * (ZB/ZDENOM1**2.) + ZFACT2 * (ZB2/ZDENOM2**2.) ) * 1./2. )
          ZSSA = MAX( ZSSA, 8.*10. )
          !
          PSNOWDIAMOPT(JJ,JST) = 6./( XRHOLI*ZSSA ) * 10.
          !
      END IF
    END DO
  END DO
!---------------------------------
!    FLANNER et al. 2006 (F06)
!
! -> Dry snow
! -> Evolution of optical diameter
!---------------------------------
ELSEIF (HSNOWMETAMO=='F06'  .OR. HSNOWMETAMO=='S-F')THEN
  !
  DO JST = 1,IMAX_USE
    DO JJ = 1,SIZE(PSNOWRHO,1)
      IF (INLVLS_USE(JJ) == 0) CYCLE
      IF (  JST<=INLVLS_USE(JJ) .AND. PSNOWLIQ(JJ,JST)<=XUEPSI .AND. &
        ((HSNOWMETAMO=='F06' ) .OR. ( (HSNOWMETAMO=='S-F').AND.((PSNOWAGE(JJ,JST)>2.) &
                                      .OR.(PSNOWRHO(JJ,JST) >= XRHOTHRESHOLD_ICE)) ) ) )THEN
          !
          !  WRITE(*,*) CSNOWMETAMO,': you are using F06 formulation!!'
          !
          ! XDRDT0(dens,gradT,T), XTAU(dens,gradT,T), XKAPPA(dens,gradT,T)
          ! dens: [1-8 <-> 50.-400. kg/m3]
          ! gradT: [1-31 <-> 0.-300. K/m]
          ! T: [1-11 <-> 223.15-273.15 K]
          !
          !  Select indices of density, temperature gradient and temperature
          IDRHO  = MIN( ABS( INT( (PSNOWRHO(JJ,JST)-25.)/50.        ) + 1 ), 8  )
          IDGRAD = MIN( ABS( INT( (ZGRADT(JJ,JST)-5.)/10.+2.                )     ), 31 )
          IDTEMP = MIN( ABS( INT( (PSNOWTEMP(JJ,JST)-225.65 )/5.+2. )     ), 11 )
          IF ( PSNOWTEMP(JJ,JST)<221. ) IDTEMP = 1
          !
          ! Compute SSA
          ZOPTR0 = XVDIAM6/2. * 10.**6.
          ! Matthieu : Add maximum because there can be continuity problems in the
          ! S-F metamorphism option
          ZOPTR  = MAX(PSNOWDIAMOPT(JJ,JST)/2. * 10.**6.,ZOPTR0)
          ZDRDT  = XDRDT0(IDRHO,IDGRAD,IDTEMP) * &
                   EXP(1./XKAPPA(IDRHO,IDGRAD,IDTEMP)* &
                   LOG(XTAU(IDRHO,IDGRAD,IDTEMP) /( ZOPTR - ZOPTR0 + XTAU(IDRHO,IDGRAD,IDTEMP) )))
!          ZDRDT  = XDRDT0(IDRHO,IDGRAD,IDTEMP) * &
!                   ( XTAU(IDRHO,IDGRAD,IDTEMP) / &
!                     ( ZOPTR - ZOPTR0 + XTAU(IDRHO,IDGRAD,IDTEMP) ) )**(1./XKAPPA(IDRHO,IDGRAD,IDTEMP))
          ZOPTR  = ZOPTR + ZDRDT * ZTSTEPHOUR
          ZOPTR  = MIN( ZOPTR, 3./(XRHOLI*2.) * 10.**6.)
          !
          !PSNOWDIAMOPT(JJ,JST) = ZOPTR * 2./10.**6.
          PSNOWDIAMOPT(JJ,JST) = ZOPTR * 0.000002
          !
      END IF
    END DO
  END DO
  !
END IF
IF ( HSNOWMETAMO=='S-F'.OR. HSNOWMETAMO=='S-B')THEN
  !
  DO JST = 1,IMAX_USE
    DO JJ = 1,SIZE(PSNOWRHO,1)
      IF (INLVLS_USE(JJ) == 0) CYCLE
      IF (JST<=INLVLS_USE(JJ) .AND. PSNOWLIQ(JJ,JST)<=XUEPSI &
                              .AND. PSNOWAGE(JJ,JST)<=2. &
                              .AND. PSNOWRHO(JJ,JST) < XRHOTHRESHOLD_ICE) THEN
        ZSSA = ZSIXOFXRHOLI / PSNOWDIAMOPT(JJ,JST)
        ! Equation 5 by Schleef et al 2014 gives -dSSA/dt in m2 kg-1 hour -1:
        ZSSA=ZSSA-(1.1E-6+3.1E-8*(PSNOWTEMP(JJ,JST)-273.15))*(ZSSA**3.1)*ZTSTEPHOUR
        PSNOWDIAMOPT(JJ,JST) = ZSIXOFXRHOLI / ZSSA
      !
      END IF
    END DO
  END DO
  !
END IF
!
!*    2. MISE A JOUR VARIABLES HISTORIQUES (CAS NON DENDRITIQUE).
!        UPDATE OF THE HISTORICAL VARIABLES
!        --------------------------------------------------------
GNONDENDRITIC(:,:) = ( PSNOWDIAMOPT(:,:)>=XVDIAM6*(4.-PSNOWSPHERI(:,:))-XUEPSI )
!
DO JST = 1,IMAX_USE
  !
  DO JJ = 1,SIZE(PSNOWRHO,1)
    !
    IF (INLVLS_USE(JJ) == 0) CYCLE
    IF ( JST<=INLVLS_USE(JJ)) THEN
    !
      IF ( GNONDENDRITIC (JJ,JST)) THEN
        !
        IF (PSNOWSPHERI(JJ,JST)<XVSPHE4 .AND. ISNOWHIST(JJ,JST)==0) THEN
          !
          ISNOWHIST(JJ,JST) = NVHIS1
          !
        ELSEIF ( XVSPHE1-PSNOWSPHERI(JJ,JST)<XVSPHE4 .AND. PSNOWLIQ(JJ,JST)/PSNOWDZ(JJ,JST)>XVTELV1  ) THEN
          !
          IF (ISNOWHIST(JJ,JST)==0)     ISNOWHIST(JJ,JST) = NVHIS2
          IF (ISNOWHIST(JJ,JST)==NVHIS1) ISNOWHIST(JJ,JST) = NVHIS3
          !
        ELSEIF ( PSNOWTEMP(JJ,JST) < XTT ) THEN
          !
          IF(ISNOWHIST(JJ,JST)==NVHIS2) ISNOWHIST(JJ,JST) = NVHIS4
          IF(ISNOWHIST(JJ,JST)==NVHIS3) ISNOWHIST(JJ,JST) = NVHIS5
          !
        END IF
        !
      END IF
    END IF
    !
  END DO
  !
END DO
!
PSNOWHIST(:,:) = ISNOWHIST(:,:) ! for output
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROMETAMO',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWCROMETAMO
!
!####################################################################

!####################################################################
SUBROUTINE SNOWCRORAD(TPTIME, OMEB,                            &
                      PSW_RAD, PSNOWALB, PSNOWDZ,              &
                      PSNOWRHO, PALB, PSWNETSNOW, PSWNETSNOWS, &
                      PRADSINK, PRADXS,                        &
                      PSNOWDIAMOPT, PSNOWSPHERI, PSNOWAGE,PPS,    &
                      PZENITH, PPERMSNOWFRAC,KNLVLS_USE,       &
                      OSNOW_ABS_ZENITH, OSNOWAGING_VAR, PAGINGCOEF)
!
!!    PURPOSE
!!    -------
!     Calculate the transmission of shortwave (solar) radiation
!     through the snowpack (using a form of Beer's Law: exponential
!     decay of radiation with increasing snow depth).
!     Needs a first calculation of the albedo to stay coherent with
!     ISBA-ES ==> make sure to keep SNOWCRORAD coherent with SNOWCROALB
!
USE MODD_SNOW_PAR, ONLY : XVSPEC1, XVSPEC2, XVSPEC3, XVBETA1, &
                          XVBETA2, XVBETA3, XVBETA4, XVBETA5
USE MODD_MEB_PAR,  ONLY : XSW_WGHT_VIS, XSW_WGHT_NIR
USE MODD_TYPE_DATE_SURF, ONLY : DATE_TIME
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(DATE_TIME), INTENT(IN)       :: TPTIME      ! current date and time

LOGICAL,            INTENT(IN)      :: OMEB ! if=T, then uppermost abs is diagnosed
!                                           !       since fluxes known
!
REAL, DIMENSION(:), INTENT(IN)      :: PSW_RAD, PSNOWALB, PALB,PPERMSNOWFRAC
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWRHO, PSNOWDZ
!
LOGICAL, INTENT(IN)                 :: OSNOW_ABS_ZENITH ! parametrization for polar regions (not physic but better results)
!                                                       ! default FALSE
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PSWNETSNOW, PSWNETSNOWS
!
REAL, DIMENSION(:), INTENT(OUT)     :: PRADXS
!
REAL, DIMENSION(:,:), INTENT(OUT)   :: PRADSINK
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWDIAMOPT, PSNOWSPHERI, PSNOWAGE
REAL, DIMENSION(:), INTENT(IN)      :: PPS
INTEGER, DIMENSION(:), INTENT(IN)   :: KNLVLS_USE
REAL, DIMENSION(:), INTENT(IN)      :: PZENITH
REAL, DIMENSION(:), INTENT(IN)      :: PAGINGCOEF
LOGICAL, INTENT(IN)                 :: OSNOWAGING_VAR
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO,1))    :: ZALB_NEW
REAL, DIMENSION(SIZE(PSNOWRHO,1),3)  :: ZALB !albedo 3 bands
!REAL, DIMENSION(SIZE(PSNOWRHO,1),3)  :: ZALB_2B !albedo 2 bands for MEB
REAL, DIMENSION(SIZE(PSNOWRHO,2))    :: ZDIAM
REAL, DIMENSION(SIZE(PSNOWRHO,2),3)  :: ZBETA
REAL, DIMENSION(3) :: ZOPTICALPATH, ZFACT
REAL :: ZPROJLAT
!
INTEGER :: JJ,JST,JB   ! looping indexes
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCRORAD',0,ZHOOK_HANDLE)
!
! 0. Initialization:
! ------------------
!
PRADSINK(:,:) = 0.
!
! 1. Computation of the new albedo (see SNOWCROALB):
! -----------------------------------
!
 CALL SNOWCROALB(ZALB_NEW,ZALB,PSNOWDZ(:,1),PSNOWRHO(:,1:2),                  &
                 PPERMSNOWFRAC,PSNOWDIAMOPT(:,1),PSNOWSPHERI(:,1),               &
                 PSNOWAGE(:,1),PSNOWDIAMOPT(:,2),PSNOWSPHERI(:,2),PSNOWAGE(:,2), &
                 PPS, PZENITH, OSNOWAGING_VAR,PAGINGCOEF,KNLVLS_USE )
!
!!!IF (OMEB) THEN
!!!  ZALB_2B(:,1)=ZALB(:,1)
!!!  ZALB_2B(:,2)=(PSNOWALB(:) - XSW_WGHT_VIS*ZALB_2B(:,1))/XSW_WGHT_NIR
!!!END IF

DO JJ = 1,SIZE(PSW_RAD)
  !
  ZDIAM(:)=PSNOWDIAMOPT(JJ,:)
!!!!!  DO JST = 1,KNLVLS_USE(JJ)
!!!!!    CALL GET_DIAM(PSNOWDIAMOPT(JJ,JST),PSNOWSPHERI(JJ,JST),ZDIAM(JST),HSNOWMETAMO)
!!!!!  END DO    ! end loop snow layers
  !
  ! 2. Extinction of net shortwave radiation
  ! ----------------------------------------
  ! First calculates extinction coefficients fn of grain size and density
  ! then calculates exctinction in the layer and increases optical path length
  !
  ! Coefficient for taking into account the increase of path length of rays
  ! in snow due to zenithal angle
  ZPROJLAT = 1. / MAX( XUEPSI, COS(PZENITH(JJ)) )
  !
  IF (OMEB) THEN
      PRADSINK(JJ,:) = -PSWNETSNOW(JJ)  / ( 1.-ZALB_NEW(JJ) )
  ELSE
      PRADSINK(JJ,:) = -PSW_RAD(JJ) * ( 1.-PSNOWALB(JJ) ) / ( 1.-ZALB_NEW(JJ) )
  END IF
  !
  !   Initialize optical depth
  ZOPTICALPATH(:) = 0.
  !
  DO JST = 1,IMAX_USE
    IF (JST <= KNLVLS_USE(JJ)) THEN
      !
      ZBETA(JST,1) = MAX( XVBETA1 * PSNOWRHO(JJ,JST) / SQRT(ZDIAM(JST)), XVBETA2 )
      ZBETA(JST,2) = MAX( XVBETA3 * PSNOWRHO(JJ,JST) / SQRT(ZDIAM(JST)), XVBETA4 )
      ZBETA(JST,3) = XVBETA5
      !
      ZFACT(:) = 0.
      DO JB = 1,3
        ZOPTICALPATH(JB) = ZOPTICALPATH(JB) + ZBETA(JST,JB) * PSNOWDZ(JJ,JST)
        IF (OSNOW_ABS_ZENITH) THEN
          !This formulation is incorrect but it compensate partly the fact that the albedo formulation does not account for zenithal angle
          ZFACT(JB) = (1.-ZALB(JJ,JB)) * EXP( -ZOPTICALPATH(JB)*ZPROJLAT)
        ELSE
          ZFACT(JB) = (1.-ZALB(JJ,JB)) * EXP( -ZOPTICALPATH(JB) )
        END IF
  !
      END DO
      !
      IF (OMEB) THEN
! ML : This is now possible because snowpack should be identical at MEB and Crocus steps.
! However, take care if snowfall is intercepted... still partly incorrect
        IF (JST==1) THEN
          PRADSINK(JJ,JST) = PSWNETSNOWS(JJ)-PSWNETSNOW(JJ)
        ELSE
          PRADSINK(JJ,JST) = PRADSINK(JJ,JST) * &
                         ( XSW_WGHT_VIS*ZFACT(1) + XSW_WGHT_NIR*ZFACT(2) )
        END IF
      ELSE
      PRADSINK(JJ,JST) = PRADSINK(JJ,JST) * &
                         ( XVSPEC1*ZFACT(1) + XVSPEC2*ZFACT(2) + XVSPEC3*ZFACT(3) )
      END IF
      !
    END IF
  END DO    ! end loop snow layers
  !
  ! For thin snow packs, radiation might reach base of
  ! snowpack and the reflected energy can be absorbed by the bottom of snow layer:
  ! THIS PROCESS IS NOT SIMULATED
  !
  ! 4. Excess radiation not absorbed by snowpack (W/m2)JJ
  ! ----------------------------------------------------
  !
  PRADXS(JJ) = -PRADSINK( JJ,KNLVLS_USE(JJ) )
  !
END DO    !end loop grid points
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCRORAD',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOWCRORAD
!####################################################################
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOWCROEBUD(HSNOWRES, HIMPLICIT_WIND,                                   &
                       PPEW_A_COEF, PPEW_B_COEF,                                   &
                       PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,         &
                       PSNOWDZMIN,                                                 &
                       PZREF,PTS,PSNOWRHO,PSNOWLIQ,PSCAP,PSCOND1,PSCOND2,          &
                       PUREF,PEXNS,PEXNA,PDIRCOSZW,PVMOD,                          &
                       PLW_RAD,PSW_RAD,PTA,PQA,PPS,PTSTEP,                         &
                       PSNOWDZ1,PSNOWDZ2,PALBT,PZ0,PZ0EFF,PZ0H,                    &
                       PSFCFRZ,PRADSINK,PHPSNOW,                                   &
                       PCT,PEMIST,PRHOA,PTSTERM1,PTSTERM2,PRA,PCDSNOW,PCHSNOW,     &
                       PQSAT,PDQSAT,PRSRA,PUSTAR2_IC, PRI,                         &
                       PPET_A_COEF_T,PPEQ_A_COEF_T,PPET_B_COEF_T,PPEQ_B_COEF_T,    &
                       OFRZRAIN,PRR, PRESA_SN)
!
!!    PURPOSE
!!    -------
!     Calculate surface energy budget linearization (terms) and turbulent
!     exchange coefficients/resistance between surface and atmosphere.
!     (Noilhan and Planton 1989; Giordani 1993; Noilhan and Mahfouf 1996)
!
!!    MODIFICATIONS
!!    -------------
!!      Original A. Boone
!!      Modified by E. Brun (24/09/2012) :
!!      * Correction coupling coefficient for specific humidity in SNOWCROEBUD
!!      * PSFCFRZ(:)  = 1.0 for systematic solid/vapor latent fluxes in SNOWCROEBUD
!!      Modified by B. Decharme 09/12  new wind implicitation
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_CSTS,     ONLY : XCPD, XSTEFAN, XLVTT, XLSTT, XCL
USE MODD_SNOW_PAR, ONLY : X_RI_MAX, XEMISSN
!
USE MODE_THERMOS
!
USE MODI_SURFACE_RI
USE MODI_SURFACE_AERO_COND
USE MODI_SURFACE_CD
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                   :: PTSTEP, PSNOWDZMIN
!
 CHARACTER(LEN=*),     INTENT(IN)  :: HSNOWRES ! type of sfc resistance
!                                      DEFAULT=Louis (1979), standard ISBA
!                                      method. Option to limit Ri number
!                                      for very stable conditions
!
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
!                                                     ! 'OLD' = direct
!                                                     ! 'NEW' = Taylor serie, order 1
!
REAL, DIMENSION(:), INTENT(IN)     :: PPEW_A_COEF, PPEW_B_COEF,                   &
                                      PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF,      &
                                      PPEQ_B_COEF
!                                      PPEW_A_COEF = wind coefficient (m2s/kg)
!                                      PPEW_B_COEF = wind coefficient (m/s)
!                                      PPET_A_COEF = A-air temperature coefficient
!                                      PPET_B_COEF = B-air temperature coefficient
!                                      PPEQ_A_COEF = A-air specific humidity coefficient
!                                      PPEQ_B_COEF = B-air specific humidity coefficient
!
REAL, DIMENSION(:), INTENT(IN)     :: PZREF, PTS, PSNOWDZ1, PSNOWDZ2,        &
                                      PRADSINK, PSNOWRHO, PSNOWLIQ, PSCAP,   &
                                      PSCOND1, PSCOND2,                      &
                                      PZ0, PHPSNOW,                          &
                                      PALBT, PZ0EFF, PZ0H
!
REAL, DIMENSION(:), INTENT(IN)     :: PSW_RAD, PLW_RAD, PTA, PQA, PPS, PRHOA
!
REAL, DIMENSION(:), INTENT(IN)     :: PUREF, PEXNS, PEXNA, PDIRCOSZW, PVMOD
!
REAL, DIMENSION(:), INTENT(IN)     :: PRR
REAL, DIMENSION(:), INTENT(INOUT)     :: PRESA_SN
LOGICAL, DIMENSION(:), INTENT(IN)  :: OFRZRAIN
!
REAL, DIMENSION(:), INTENT(OUT)    :: PTSTERM1, PTSTERM2, PEMIST, PRA,         &
                                      PCT, PSFCFRZ, PCDSNOW, PCHSNOW,          &
                                      PQSAT, PDQSAT, PRSRA
!
REAL, DIMENSION(:), INTENT(OUT)    :: PUSTAR2_IC,                      &
                                      PPET_A_COEF_T, PPEQ_A_COEF_T,    &
                                      PPET_B_COEF_T, PPEQ_B_COEF_T
!
REAL, DIMENSION(:), INTENT(OUT)    :: PRI
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PTS))        :: ZAC, ZRI,                        &
                                     ZSCONDA, ZA, ZB, ZC,             &
                                     ZCDN, ZSNOWDZM1, ZSNOWDZM2,      &
                                     ZVMOD, ZUSTAR2, ZTS3, ZLVT,      &
                                     Z_CCOEF
REAL, DIMENSION(SIZE(PTS))        :: ZDENOM, ZNUMER
!
REAL, DIMENSION(SIZE(PTA))          :: ZPHIFRZ ! ZPHIFRZ = instant freezing ratio of supercooled water
INTEGER :: JJ   ! looping indexes
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
! 1. New saturated specific humidity and derrivative:
! ---------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROEBUD',0,ZHOOK_HANDLE)
!
ZRI   (:) = XUNDEF
!
PQSAT (:) = QSATI(PTS(:),PPS(:))
PDQSAT(:) = DQSATI(PTS(:),PPS(:),PQSAT(:))
!
! 2. Surface properties:
! ----------------------
! It might be of interest to use snow-specific roughness
! or a temperature dependence on emissivity:
! but for now, use ISBA defaults.
!
PEMIST(:) = XEMISSN
!
! 2. Computation of resistance and drag coefficient
! -------------------------------------------------
!
 CALL SURFACE_RI(PTS, PQSAT, PEXNS, PEXNA, PTA, PQA,  &
                 PZREF, PUREF, PDIRCOSZW, PVMOD, ZRI  )
!
! Simple adaptation of method by Martin and Lejeune (1998)
! to apply a lower limit to turbulent transfer coef
! by defining a maximum Richarson number for stable
! conditions:
!
IF ( HSNOWRES=='RIL' .or. HSNOWRES=='RI1' .or. HSNOWRES=='RI2' ) THEN
  DO JJ=1,SIZE(ZRI)
    ZRI(JJ) = MIN( X_RI_MAX, ZRI(JJ) )
  END DO
END IF
!
PRI(:) = ZRI(:)
!
! Surface aerodynamic resistance for heat transfers
!
 CALL SURFACE_AERO_COND(ZRI, PZREF, PUREF, PVMOD, PZ0, PZ0H, PRESA_SN, ZAC, PRA, PCHSNOW, HSNOWRES=HSNOWRES)
!
PRSRA(:) = PRHOA(:) / PRA(:)
!
! For atmospheric model coupling:
!
 CALL SURFACE_CD(ZRI, PZREF, PUREF, PZ0EFF, PZ0H, PCDSNOW, ZCDN)
!
!
! Modify flux-form implicit coupling coefficients:
! - wind components:
!
ZNUMER(:) = PCDSNOW(:)*PVMOD(:)
ZDENOM(:) = PRHOA(:) * PCDSNOW(:) * PVMOD(:) * PPEW_A_COEF(:)
IF(HIMPLICIT_WIND=='OLD')THEN
! old implicitation
  ZUSTAR2(:) = ( ZNUMER(:) * PPEW_B_COEF(:) ) / ( 1.0 - ZDENOM(:) )
ELSE
! new implicitation
  ZUSTAR2(:) = ( ZNUMER(:) * ( 2.*PPEW_B_COEF(:) - PVMOD(:) ) ) / ( 1.0 - 2.0*ZDENOM(:) )
END IF
!
ZVMOD(:) = PRHOA(:)*PPEW_A_COEF(:)*ZUSTAR2(:) + PPEW_B_COEF(:)

ZVMOD(:) = MAX( ZVMOD(:),0. )
!
WHERE ( PPEW_A_COEF(:)/= 0. )
  ZUSTAR2(:) = MAX( ( ZVMOD(:) - PPEW_B_COEF(:) ) / (PRHOA(:)*PPEW_A_COEF(:)), 0. )
ENDWHERE
!
! implicit wind friction
ZUSTAR2(:) = MAX( ZUSTAR2(:),0. )
!
PUSTAR2_IC(:) = ZUSTAR2(:)
!
! 3. Calculate linearized surface energy budget components:
! ---------------------------------------------------------
! To prevent numerical difficulties for very thin snow
! layers, limit the grid "thinness": this is important as
! layers become vanishing thin:
!
ZSNOWDZM1(:) = MAX( PSNOWDZ1(:), PSNOWDZMIN )
ZSNOWDZM2(:) = MAX( PSNOWDZ2(:), PSNOWDZMIN )
!
! Surface thermal inertia:
!
PCT(:) = 1.0 / ( PSCAP(:)*ZSNOWDZM1(:) )
!
! Fraction of surface frozen (sublimation) with the remaining
! fraction being liquid (evaporation):
!
PSFCFRZ(:) = 1.0
!
! Thermal conductivity between uppermost and lower snow layers:
!
ZSCONDA(:) = ( ZSNOWDZM1(:)*PSCOND1(:) + ZSNOWDZM2(:)*PSCOND2(:) ) / &
             ( ZSNOWDZM1(:)            + ZSNOWDZM2(:)            )
!
! Transform implicit coupling coefficients:
! Note, surface humidity is 100% over snow.
!
! - specific humidity:
!
Z_CCOEF(:) = 1.0 - PPEQ_A_COEF(:) * PRSRA(:)
!
PPEQ_A_COEF_T(:) = - PPEQ_A_COEF(:) * PRSRA(:) * PDQSAT(:) / Z_CCOEF(:)
!
PPEQ_B_COEF_T(:) = ( PPEQ_B_COEF(:) &
                   - PPEQ_A_COEF(:) * PRSRA(:) * (PQSAT(:) - PDQSAT(:)*PTS(:)) ) / Z_CCOEF(:)
!
! - air temperature:
!   (assumes A and B correspond to potential T):
!
Z_CCOEF(:) = ( 1.0 - PPET_A_COEF(:) * PRSRA(:) ) / PEXNA(:)
!
PPET_A_COEF_T(:) = - PPET_A_COEF(:) * PRSRA(:) / ( PEXNS(:) * Z_CCOEF(:) )
!
PPET_B_COEF_T(:) = PPET_B_COEF(:) / Z_CCOEF(:)
!
!
! Energy budget solution terms:
!
ZTS3(:) = PEMIST(:) * XSTEFAN * PTS(:)**3
ZLVT(:) = (1.-PSFCFRZ(:))*XLVTT + PSFCFRZ(:)*XLSTT
!
ZA(:) = 1./PTSTEP + PCT(:) * &
         ( 4. * ZTS3(:) + PRSRA(:) * ZLVT(:) * ( PDQSAT(:) - PPEQ_A_COEF_T(:) )            &
                        + PRSRA(:) * XCPD * ( (1./PEXNS(:))-(PPET_A_COEF_T(:)/PEXNA(:)) )  &
                        + ( 2*ZSCONDA(:) / ( ZSNOWDZM2(:)+ZSNOWDZM1(:) ) ) )
!
ZB(:) = 1./PTSTEP + PCT(:) * &
         ( 3. * ZTS3(:) + PRSRA(:) * ZLVT(:) * PDQSAT(:) )
!
WHERE(OFRZRAIN(:))
  ZPHIFRZ(:) = ( XTT - PTA(:) ) * XCL / XLMTT
ELSEWHERE
  ZPHIFRZ(:) = 1.
ENDWHERE
ZC(:) = PCT(:) * ( - PRSRA(:) * ZLVT(:) * ( PQSAT(:) - PPEQ_B_COEF_T(:) )  &
                   + PRSRA(:) * XCPD * PPET_B_COEF_T(:) / PEXNA(:)         &
                   + PSW_RAD(:) * (1. - PALBT(:)) + PEMIST(:) * PLW_RAD(:) &
                   + PHPSNOW(:) + PRADSINK(:)                              &
                   + ( 1. - ZPHIFRZ(:) ) * XLMTT * PRR(:))
!
!
! Coefficients needed for implicit solution
! of linearized surface energy budget:
!
PTSTERM2(:) = 2. * ZSCONDA(:) * PCT(:) / ( ZA(:) * (ZSNOWDZM2(:)+ZSNOWDZM1(:) ) )
!
PTSTERM1(:) = ( PTS(:)*ZB(:) + ZC(:) ) / ZA(:)
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROEBUD',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWCROEBUD
!####################################################################
!####################################################################
!####################################################################
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOWCROEBUDMEB(PTSTEP, PSNOWDZMIN,                        &
           PTS, PSNOWDZ1, PSNOWDZ2, PSCOND1, PSCOND2, PSCAP,        &
           PRNSNOW, PHSNOW, PLES3L, PLEL3L, PRADSINK, PHPSNOW,      &
           PCT, PTSTERM1, PTSTERM2, PGFLUXSNOW                      )
!
!!    PURPOSE
!!    -------
!     Calculate surface energy budget with surface fluxes imposed.
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,    INTENT(IN)               :: PTSTEP, PSNOWDZMIN
!
REAL, DIMENSION(:), INTENT(IN)    :: PTS, PSNOWDZ1, PSNOWDZ2, PSCOND1, PSCOND2, PSCAP,  &
                                     PRNSNOW, PHSNOW, PLES3L, PLEL3L, PRADSINK, PHPSNOW
!
REAL, DIMENSION(:), INTENT(OUT)   :: PCT, PTSTERM1, PTSTERM2, PGFLUXSNOW
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PTS))        :: ZSCONDA, ZA, ZB, ZC,             &
                                     ZSNOWDZM1, ZSNOWDZM2
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROEBUDMEB',0,ZHOOK_HANDLE)
!
! Calculate surface energy budget components:
! ---------------------------------------------------------
! To prevent numerical difficulties for very thin snow
! layers, limit the grid "thinness": this is important as
! layers become vanishing thin:
!
ZSNOWDZM1(:)  = MAX(PSNOWDZ1(:), PSNOWDZMIN)
ZSNOWDZM2(:)  = MAX(PSNOWDZ2(:), PSNOWDZMIN)
!
! Surface thermal inertia:
!
PCT(:)        = 1.0/(PSCAP(:)*ZSNOWDZM1(:))
!
! Surface fluxes entering the snowpack (radiative and turbulent):
!
PGFLUXSNOW(:) = PRNSNOW(:) - PHSNOW(:) - PLES3L(:) - PLEL3L(:)
!
! Thermal conductivity between uppermost and lower snow layers:
!
ZSCONDA(:)    = (ZSNOWDZM1(:)+ZSNOWDZM2(:))/                           &
               ((ZSNOWDZM1(:)/PSCOND1(:)) + (ZSNOWDZM2(:)/PSCOND2(:)))
!
!
! Energy budget solution terms (with surface flux imposed):
!
ZB(:)         = 1./PTSTEP ! Equation 14 doc MEB/snow(A. Boone)
!
ZA(:)         = ZB(:) + PCT(:)*(2*ZSCONDA(:)/(ZSNOWDZM2(:)+ZSNOWDZM1(:)))  ! combinaison des équations 15 et 2 doc MEB/snow (A. Boone)
!
ZC(:)         = PCT(:)*( PGFLUXSNOW(:) + PHPSNOW(:) + PRADSINK(:) ) ! Equation 16 doc MEB/snow (A. Boone)
!
! Coefficients needed for implicit solution
! of linearized surface energy budget:
!
PTSTERM2(:)   = 2*ZSCONDA(:)*PCT(:)/(ZA(:)*(ZSNOWDZM2(:)+ZSNOWDZM1(:))) ! Equation 19 doc MEB/snow (A. Boone)
!
PTSTERM1(:)   = (PTS(:)*ZB(:) + ZC(:))/ZA(:) ! Equation 18 doc MEB/snow (A. Boone)
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROEBUDMEB',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWCROEBUDMEB
!
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOWCROSOLVT(OMEB,PTSTEP,PSNOWDZMIN,                &
                        PSNOWDZ,PSCOND,PSCAP,PTG,              &
                        PSOILCOND,PD_G,                        &
                        PRADSINK,PCT,PTERM1,PTERM2,            &
                        PPET_A_COEF_T,PPEQ_A_COEF_T,           &
                        PPET_B_COEF_T,PPEQ_B_COEF_T,           &
                        PTA_IC, PQA_IC,                        &
                        PGBAS,PSNOWTEMP,PSNOWFLUX,             &
                        KNLVLS_USE                             )
!
!!    PURPOSE
!!    -------
!     This subroutine solves the 1-d diffusion of 'ZSNOWTEMP' using a
!     layer averaged set of equations which are time differenced
!     using the backward-difference scheme (implicit).
!     The eqs are solved rapidly by taking advantage of the
!     fact that the matrix is tridiagonal. This is a very
!     general routine and can be used for the 1-d diffusion of any
!     quantity as long as the diffusity is not a function of the
!     quantity being diffused. Aaron Boone 8-98. Soln to the eq:
!
!                 c  dQ    d  k dQ    dS
!                    -- =  --   -- -  --
!                    dt    dx   dx    dx
!
!     where k = k(x) (thermal conductivity), c = c(x) (heat capacity)
!     as an eg. for temperature/heat eq. S is a sink (-source) term.
!     Diffusivity is k/c
!
!!     MODIFICATIONS
!!     -------------
!!      Original A. Boone
!!       05/2011: Brun  Special treatment to tackle the variable number
!!                      of snow layers
!
USE MODD_CSTS, ONLY : XTT
!
USE MODI_TRIDIAG_GROUND_SNOWCRO
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
LOGICAL,            INTENT(IN)      :: OMEB
!
REAL, INTENT(IN)                    :: PTSTEP, PSNOWDZMIN
!
REAL, DIMENSION(:), INTENT(IN)      :: PTG, PSOILCOND, PD_G,        &
                                        PCT, PTERM1, PTERM2

!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWDZ, PSCOND, PSCAP,      &
                                        PRADSINK
!
REAL, DIMENSION(:), INTENT(IN)      :: PPET_A_COEF_T, PPEQ_A_COEF_T, &
                                        PPET_B_COEF_T, PPEQ_B_COEF_T
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWTEMP
!
REAL, DIMENSION(:), INTENT(OUT)     :: PGBAS, PSNOWFLUX, PTA_IC, PQA_IC
!
INTEGER, DIMENSION(:), INTENT(IN)   :: KNLVLS_USE
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWDZ,1),SIZE(PSNOWDZ,2)) :: ZSNOWTEMP, ZDTERM, ZCTERM, &
                                                    ZFRCV, ZAMTRX, ZBMTRX,     &
                                                    ZCMTRX
!
REAL, DIMENSION(SIZE(PSNOWDZ,1),SIZE(PSNOWDZ,2)) :: ZWORK1, ZWORK2, ZDZDIF,    &
                                                    ZSNOWDZM
!
REAL, DIMENSION(SIZE(PSNOWDZ,1),SIZE(PSNOWDZ,2)-1) :: ZSNOWTEMP_M
!
REAL, DIMENSION(SIZE(PTG)) :: ZGBAS, ZSNOWTEMP_DELTA
!
INTEGER :: JJ, JST   ! looping indexes
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROSOLVT',0,ZHOOK_HANDLE)
!
! 0. Initialize:
! ------------------
!
ZSNOWTEMP(:,:) = PSNOWTEMP(:,:)
!
!
! 1. Calculate tri-diagnoal matrix coefficients:
! ----------------------------------------------
! For heat transfer, assume a minimum grid
! thickness (to prevent numerical
! problems for very thin snow cover):
!
DO JJ=1,SIZE(PTG) ! LOOP FOR POINTS
  !
  DO JST = KNLVLS_USE(JJ),1,-1 !PARCOUR LES COUCHES DANS L' ORDRE INVERSE
    !
    ZSNOWDZM(JJ,JST) = MAX( PSNOWDZ(JJ,JST), PSNOWDZMIN ) !APPLICATION DU MINIMUM DE TAILLE DES COUCHES (CREATION DE NEIGE ) M
    !
    ZWORK1(JJ,JST)   = ZSNOWDZM(JJ,JST) * PSCOND(JJ,JST) ! W/K CONDUCTANCE THERMIQUE W/K
    !
    IF ( JST<KNLVLS_USE(JJ) ) THEN ! SI PAS AU CONTACT DU SOL
      !
      ZDZDIF(JJ,JST) = ZSNOWDZM(JJ,JST) + ZSNOWDZM(JJ,JST+1) ! FUSION DE LA COUCHE 1 ET 2 , 2 ET 3 ....
      !
      ZWORK2(JJ,JST) = ZSNOWDZM(JJ,JST+1) * PSCOND(JJ,JST+1) !CONDUCTANCE THERMIQUE DE LA COUCHE SUIVANTEZSNOWHISTZSNOWHIST W/K
      !
    ELSE
      !
      ZDZDIF(JJ,JST) = ZSNOWDZM(JJ,JST) + PD_G(JJ) ! FUSION DE LA DERNIERE COUCHE AVEC LA PREMIER DU SOL
      !
      ZWORK2(JJ,JST) = PD_G(JJ) * PSOILCOND(JJ) ! CONDUCTANCE THERMOIQUE DU SOL
      !
    END IF
    !
    ZDTERM(JJ,JST) = 2.0 * ( ZWORK1(JJ,JST) + ZWORK2(JJ,JST) ) / ZDZDIF(JJ,JST)**2. ![W/KM²]
    !
    ZCTERM(JJ,JST) = PSCAP(JJ,JST) * ZSNOWDZM(JJ,JST) / PTSTEP ! W/KM²
    !
  END DO
  !
END DO
!
! 2. Set up tri-diagonal matrix
! -----------------------------
!
ZAMTRX(:,:) = 0.
ZBMTRX(:,:) = 0.
ZCMTRX(:,:) = 0.
ZFRCV(:,:) = 0.
! Upper BC
!
ZAMTRX(:,1) =  0.0 ! DETAILS ONE DAY IN LAFAYSSE ET AL GMD
ZBMTRX(:,1) =  1. / ( PCT(:)*PTSTEP ) ![W/KM²]
ZCMTRX(:,1) = - PTERM2(:) * ZBMTRX(:,1)
ZFRCV (:,1) =   PTERM1(:) * ZBMTRX(:,1)
!
DO JST = 2,IMAX_USE !NOMBRE MAX DE COUCHE SUR TOUS LES POINTS
  !
  DO JJ = 1, SIZE(PTG) ! POINTS ?
  !
    IF (JST <= KNLVLS_USE(JJ)) THEN
      !
      ! Interior Grid & Lower BC
      ZAMTRX(JJ,JST) = -ZDTERM(JJ,JST-1) !w/Km2
      ZBMTRX(JJ,JST) =  ZCTERM(JJ,JST) + ZDTERM(JJ,JST-1) + ZDTERM(JJ,JST) !w/km 2
      ZFRCV (JJ,JST) =  ZCTERM(JJ,JST)*PSNOWTEMP(JJ,JST) - (PRADSINK(JJ,JST-1)-PRADSINK(JJ,JST)) ! right hand side W/m2
      !
      IF ( JST<KNLVLS_USE(JJ) ) THEN !si pas la derniere couche
        ZCMTRX(JJ,JST) = -ZDTERM(JJ,JST)
      ELSE
        ZCMTRX(JJ,JST) = 0.0
        ZFRCV (JJ,JST) = ZFRCV(JJ,JST) + ZDTERM(JJ,JST)*PTG(JJ)
      END IF
      !
    END IF
  END DO
  !
END DO
!
!
! 4. Compute solution vector
! --------------------------
!
 CALL TRIDIAG_GROUND_SNOWCRO(ZAMTRX,ZBMTRX,ZCMTRX,ZFRCV,ZSNOWTEMP,KNLVLS_USE,IMAX_USE,0)
!
! Heat flux between surface and 2nd snow layers: (W/m2)
!
PSNOWFLUX(:) = ZDTERM(:,1) * ( ZSNOWTEMP(:,1) - ZSNOWTEMP(:,2) )
!
! 5. Snow melt case
! -----------------
! If melting in uppermost layer, assume surface layer
! temperature at freezing point and re-evaluate lower
! snowpack temperatures. This is done as it is most likely
! most signigant melting will occur within a time step in surface layer.
! Surface energy budget (and fluxes) will
! be re-calculated (outside of this routine):
!
! NOTE: if MEB is active, then surface fluxes have been defined outside
! of the snow routine and have been adjusted such that they are evaluated
! at a snow surface temperature no greater than Tf. Thus, the implicit surface temperature
! will likely never greatly exceed Tf (before melt computed and they are adjusted to Tf)
! so we can skip the next block of code when MEB is active.

IF(.NOT.OMEB)THEN
!
  ZCMTRX(:,2) = -ZDTERM(:,2)
  ZFRCV (:,2) =  ZCTERM(:,2)*PSNOWTEMP(:,2) - (PRADSINK(:,1)-PRADSINK(:,2)) + ZDTERM(:,1)*XTT
!
  CALL TRIDIAG_GROUND_SNOWCRO(ZAMTRX(:,2:IMAX_USE),ZBMTRX(:,2:IMAX_USE), &
  ZCMTRX(:,2:IMAX_USE),ZFRCV(:,2:IMAX_USE),ZSNOWTEMP_M,KNLVLS_USE,IMAX_USE,1)
!
! If melting for 2 consecuative time steps, then replace current T-profile
! with one assuming T=Tf in surface layer:
!
  ZSNOWTEMP_DELTA(:) = 0.0
!
  WHERE( ZSNOWTEMP(:,1)>XTT .AND. PSNOWTEMP(:,1)>=XTT )
    PSNOWFLUX(:) = ZDTERM(:,1) * ( XTT-ZSNOWTEMP_M(:,1) )
    ZSNOWTEMP_DELTA(:) = 1.0
  END WHERE
!
  DO JST = 2,IMAX_USE
    DO JJ = 1,SIZE(PTG)
      IF (JST <= KNLVLS_USE(JJ)) THEN
        ZSNOWTEMP(JJ,JST) = ZSNOWTEMP_DELTA(JJ)       * ZSNOWTEMP_M(JJ,JST-1) + &
                          (1.0-ZSNOWTEMP_DELTA(JJ)) * ZSNOWTEMP  (JJ,JST)
      END IF
    END DO
  END DO
!
END IF
!
! 6. Lower boundary flux:
! -----------------------
! NOTE: evaluate this term assuming the snow layer
! can't exceed the freezing point as this adjustment
! is made in melting routine. Then must adjust temperature
! to conserve energy:
!
DO JJ=1, SIZE(PTG)
 ZGBAS(JJ) = ZDTERM(JJ,KNLVLS_USE(JJ)) * ( ZSNOWTEMP(JJ,KNLVLS_USE(JJ))             - PTG(JJ) )
 PGBAS(JJ) = ZDTERM(JJ,KNLVLS_USE(JJ)) * ( MIN( XTT, ZSNOWTEMP(JJ,KNLVLS_USE(JJ)) ) - PTG(JJ) )
 ZSNOWTEMP(JJ,KNLVLS_USE(JJ)) = ZSNOWTEMP(JJ,KNLVLS_USE(JJ)) + &
                                ( ZGBAS(JJ)-PGBAS(JJ) ) / ZCTERM(JJ,KNLVLS_USE(JJ))
END DO
!
! 7. Update temperature profile in time:
! --------------------------------------
!
DO JJ=1, SIZE(PTG)
  PSNOWTEMP(JJ,1:KNLVLS_USE(JJ)) = ZSNOWTEMP(JJ,1:KNLVLS_USE(JJ))
END DO
!
!
! 8. Compute new (implicit) air T and specific humidity
! -----------------------------------------------------
!
IF(.NOT. OMEB) THEN
!
  PTA_IC(:) = PPET_B_COEF_T(:) + PPET_A_COEF_T(:) * PSNOWTEMP(:,1)
!
  PQA_IC(:) = PPEQ_B_COEF_T(:) + PPEQ_A_COEF_T(:) * PSNOWTEMP(:,1)
END IF
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROSOLVT',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWCROSOLVT
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOWCROMELT(PSCAP,PSNOWTEMP,PSNOWDZ,         &
                       PSNOWRHO,PSNOWLIQ,PSNOWIMPUR,KNLVLS_USE     )
!
!!    PURPOSE
!!    -------
!     Calculate snow melt (resulting from surface fluxes, ground fluxes,
!     or internal shortwave radiation absorbtion). It is used to
!     augment liquid water content, maintain temperatures
!     at or below freezing, and possibly reduce the mass
!     or compact the layer(s).
!
USE MODD_CSTS,ONLY : XTT, XLMTT, XRHOLW
!
USE MODE_SNOW3L
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSCAP
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDZ, PSNOWTEMP, PSNOWRHO,   &
                                          PSNOWLIQ
!
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PSNOWIMPUR
!
INTEGER, DIMENSION(:), INTENT(IN)   :: KNLVLS_USE
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZPHASE, ZCMPRSFACT,   &
                                                      ZSNOWLWE,             &
                                                      ZSNOWMELT, ZSNOWTEMP
!
INTEGER :: JJ, JST ! looping indexes
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROMELT',0,ZHOOK_HANDLE)
!
! 0. Initialize:
! ---------------------------
!
 ZPHASE     = 0.0
 ZCMPRSFACT = 0.0
 ZSNOWLWE   = 0.0
 ZSNOWMELT  = 0.0
 ZSNOWTEMP  = 0.0
!
! 1. Determine amount of melt in each layer:
! ------------------------------------------
!
DO JST = 1,IMAX_USE
  DO JJ = 1,SIZE(PSNOWDZ,1)
    !
    IF (JST <= KNLVLS_USE(JJ)) THEN
      !
      ! Total Liquid equivalent water content of snow (m):
      ZSNOWLWE(JJ,JST) = PSNOWRHO(JJ,JST) * PSNOWDZ(JJ,JST) / XRHOLW
      !
      ! Melt snow if excess energy and snow available:
      ! Phase change (J/m2)
!VV      ZPHASE(JJ,JST) = MIN( PSCAP(JJ,JST) * MAX(0.0,PSNOWTEMP(JJ,JST)-XTT) * PSNOWDZ(JJ,JST), &
!VV                           MAX(0.0,ZSNOWLWE(JJ,JST)-PSNOWLIQ(JJ,JST)) * XLMTT * XRHOLW )
!VV   Modification to run Crocus in single precision mode
      ZPHASE(JJ,JST) = MIN( PSCAP(JJ,JST) * MAX(0.0,PSNOWTEMP(JJ,JST)-XTT) * PSNOWDZ(JJ,JST), &
                           MAX(0.0,ZSNOWLWE(JJ,JST)-PSNOWLIQ(JJ,JST)) * XLMTT * XRHOLW-      &
                               XUEPSI_SMP * XLMTT * PSNOWRHO(JJ,JST))                   
      !
      ! Update snow liq water content and temperature if melting:
      ! liquid water available for next layer from melting of snow
      ! which is assumed to be leaving the current layer (m):
      ZSNOWMELT(JJ,JST) = ZPHASE(JJ,JST) / (XLMTT*XRHOLW)
      !
      ! Cool off snow layer temperature due to melt:
      ZSNOWTEMP(JJ,JST) = PSNOWTEMP(JJ,JST) - ZPHASE(JJ,JST) / (PSCAP(JJ,JST)*PSNOWDZ(JJ,JST))
      !
      ! Difference with ISBA_ES: ZMELTXS should never be different of 0.
      ! because of the introduction of the tests in SNOWCROLAYER_GONE
      ! VV PSNOWTEMP(JJ,JST) =  ZSNOWTEMP(JJ,JST)
      !VV Modification to allow model run in single precision
      PSNOWTEMP(JJ,JST) = MIN(ZSNOWTEMP(JJ,JST),XTT)         
      !
      !
      ! Loss of snowpack depth: (m) and liquid equiv (m):
      ! Compression factor for melt loss: this decreases
      ! layer thickness and increases density thereby leaving
      ! total SWE constant.
      !
      ! Difference with ISBA_ES: All melt is considered to decrease the depth
      ! without consideration to the irreducible content
      !
      ZCMPRSFACT(JJ,JST) = ( ZSNOWLWE(JJ,JST) - (PSNOWLIQ(JJ,JST)+ZSNOWMELT(JJ,JST)) ) &
                       / ( ZSNOWLWE(JJ,JST) - PSNOWLIQ(JJ,JST) )
      PSNOWDZ   (JJ,JST) = PSNOWDZ (JJ,JST) * ZCMPRSFACT(JJ,JST)
      PSNOWRHO  (JJ,JST) = ZSNOWLWE(JJ,JST) * XRHOLW / PSNOWDZ(JJ,JST)
      ! scavenging of 20% of melt
      ! 2. Add snow melt to current snow liquid water content:
      ! ------------------------------------------------------
      !
      PSNOWLIQ(JJ,JST) = PSNOWLIQ(JJ,JST) + ZSNOWMELT(JJ,JST)
      !
    END IF
  END DO   ! loop JJ grid points
END DO   ! loop JST active snow layers
!
! The control below should be suppressed after further tests
DO JST = 1,IMAX_USE
  DO JJ = 1,SIZE(PSNOWDZ,1)
     IF (JST <= KNLVLS_USE(JJ) .AND. PSNOWTEMP(JJ,JST)-XTT > XUEPSI) THEN
       WRITE(*,*) 'pb dans MELT PSNOWTEMP(JJ,JST) >XTT:', JJ,JST, PSNOWTEMP(JJ,JST)
       CALL ABOR1_SFX('SNOWCRO: pb dans MELT')
     END IF
  END DO   ! loop JJ grid points
END DO
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROMELT',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWCROMELT
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOWCROREFRZ(PTSTEP,PRR, PSNOWRHO,                           &
                       PSNOWTEMP,PSNOWDZ,PSNOWLIQ, PSNOWIMPUR,PTHRUFAL, &
                        PSCAP, PLEL3L,KNLVLS_USE, HSNOWHOLD,KMAX_USE,   &
                        OFRZRAIN)
!
!!    PURPOSE
!!    -------
!     Calculate any freezing/refreezing of liquid water in the snowpack.
!     Also, calculate liquid water transmission and snow runoff.
!     Refreezing causes densification of a layer.
!
USE MODD_CSTS,     ONLY : XTT, XLMTT, XRHOLW, XCI,XRHOLI
USE MODD_SNOW_PAR, ONLY : XSNOWDMIN, XSCAVEN_COEF
USE MODD_PREP_SNOW, ONLY : NIMPUR
!
USE MODE_SNOW3L
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                      :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)        :: PRR
!
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWDZ, PSNOWTEMP, PSNOWLIQ, PSNOWRHO
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PSNOWIMPUR
!
REAL, DIMENSION(:), INTENT(INOUT)     :: PTHRUFAL
!
 CHARACTER(3), INTENT(IN)              :: HSNOWHOLD
! modifs_EB layers
INTEGER, DIMENSION(:), INTENT(IN)      :: KNLVLS_USE
REAL, DIMENSION(:,:), INTENT(IN)       :: PSCAP
REAL, DIMENSION(:), INTENT(IN)         :: PLEL3L
INTEGER, INTENT(IN)                    :: KMAX_USE
LOGICAL, DIMENSION(:), INTENT(IN)      :: OFRZRAIN
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZPHASE,              &
                                                      ZSNOWLIQ, ZSNOWRHO,  &
                                                      ZWHOLDMAX, ZSNOWDZ,  &
                                                      ZSNOWTEMP,ZSNOWSWE
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),0:SIZE(PSNOWRHO,2),NIMPUR) ::  ZFLOWIMPUR             !Mass of impurity scavenged to layer down(g)
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),0:SIZE(PSNOWRHO,2)) :: ZFLOWLIQ
!

REAL :: ZSOLIDMASS_AFTER_REFR, ZSOLIDMASS_BEFORE_REFR, ZFINALMASS ! kg m-2
!
INTEGER :: JJ, JST , JIMP   ! looping indexes
INTEGER :: INLVLS     ! maximum snow layers number
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROREFRZ',0,ZHOOK_HANDLE)
!
! 0. Initialize:
! --------------
!
INLVLS = SIZE(PSNOWDZ,2)
ZSNOWRHO (:,:) = PSNOWRHO(:,:)
ZSNOWTEMP(:,:) = PSNOWTEMP(:,:)
!
DO JST = 1, IMAX_USE
  DO JJ = 1,SIZE(ZSNOW)
    IF (JST <= INLVLS_USE(JJ)) THEN
      ZDRYDENSITY(JJ,JST) = PSNOWRHO(JJ,JST) -  PSNOWLIQ(JJ,JST) * XRHOLW / PSNOWDZ(JJ,JST)
    END IF
  END DO
END DO
!
IF ( HSNOWHOLD == 'B92' ) THEN
  DO JST=1,KMAX_USE
    DO JJ=1,SIZE(PSNOWDZ,1)
      IF(JST <= KNLVLS_USE(JJ)) THEN
        ZWHOLDMAX(JJ,JST) = XPERCENTAGEPORE/XRHOLI * (PSNOWDZ(JJ,JST) * &
              (XRHOLI-ZDRYDENSITY(JJ,JST)))
      END IF
    END DO
  END DO
ELSEIF ( HSNOWHOLD == 'BFZ' ) THEN
  DO JST=1,KMAX_USE
    DO JJ=1,SIZE(PSNOWDZ,1)
      IF(JST <= KNLVLS_USE(JJ)) THEN
        IF ( PSNOWRHO(JJ,JST)>700. ) THEN
          ZWHOLDMAX(JJ,JST) = XPERCENTAGEPORE_ICE/XRHOLI * (PSNOWDZ(JJ,JST) * &
              (XRHOLI-PSNOWRHO(JJ,JST)) + PSNOWLIQ(JJ,JST)*XRHOLW)
        ELSE
          ZWHOLDMAX(JJ,JST) = XPERCENTAGEPORE/XRHOLI * (PSNOWDZ(JJ,JST) * &
              (XRHOLI-PSNOWRHO(JJ,JST)) + PSNOWLIQ(JJ,JST)*XRHOLW)
        END IF
      END IF
    END DO
  END DO
  IF (ANY(OFRZRAIN)) THEN
    ZSNOWDZ_FRZ(:) = 0
    DO JST=1,KMAX_USE
      DO JJ=1,SIZE(PSNOWDZ,1)
        IF(JST <= KNLVLS_USE(JJ)) THEN
          IF(OFRZRAIN(JJ)) THEN
            ZSNOWDZ_FRZ(JJ)=PSNOWDZ(JJ,JST)+ZSNOWDZ_FRZ(JJ)
            IF (ZSNOWDZ_FRZ(JJ) <= 0.02) THEN
              ZWHOLDMAX(JJ,JST) = XPERCENTAGEPORE_FRZ/XRHOLI * (PSNOWDZ(JJ,JST) * &
                    (XRHOLI-PSNOWRHO(JJ,JST)) + PSNOWLIQ(JJ,JST)*XRHOLW)
            END IF
          END IF
        END IF
      END DO
    END DO
  END IF
ELSEIF ( HSNOWHOLD == 'B02' )THEN
  DO JST=1,KMAX_USE
    DO JJ=1,SIZE(PSNOWDZ,1)
      IF(JST <= KNLVLS_USE(JJ)) THEN
        ZWHOLDMAX(JJ,JST) = SNOW3LHOLD(PSNOWRHO(JJ,JST),PSNOWDZ(JJ,JST))
      END IF
    END DO
  END DO
ELSEIF ( HSNOWHOLD == 'SPK' )THEN
  DO JST=1,KMAX_USE
    DO JJ=1,SIZE(PSNOWDZ,1)
      IF(JST <= KNLVLS_USE(JJ)) THEN
        ZWHOLDMAX(JJ,JST) = SNOWSPKHOLD(ZDRYDENSITY(JJ,JST),PSNOWDZ(JJ,JST))
      END IF
    END DO
  END DO
ELSEIF ( HSNOWHOLD == 'O04' )THEN
  DO JST=1,KMAX_USE
    DO JJ=1,SIZE(PSNOWDZ,1)
      IF(JST <= KNLVLS_USE(JJ)) THEN
        ZWHOLDMAX(JJ,JST) = SNOWO04HOLD(ZDRYDENSITY(JJ,JST),PSNOWDZ(JJ,JST) )
      END IF
    END DO
  END DO
END IF

DO JST=1,KMAX_USE
  DO JJ=1,SIZE(PSNOWDZ,1)
    IF(JST <= KNLVLS_USE(JJ)) THEN
      ZSNOWSWE(JJ,JST) = PSNOWRHO(JJ,JST)*PSNOWDZ(JJ,JST)
    END IF
  END DO
END DO
!
DO JJ = 1,SIZE(PSNOWDZ,1)  ! loop JJ grid points
  !
  ! 1. Increases Liquid Water of top layer from rain
  !    ---------------------------------------------
  !
  !  Rainfall (m) initialises the liquid flow which feeds the top layer
  !  and evaporation/condensation are taken into account
  !
  DO JIMP=1,NIMPUR
    ZFLOWIMPUR(JJ,0,JIMP) = 0.
  END DO
  IF ( KNLVLS_USE(JJ)>0. ) THEN
    IF (OFRZRAIN(JJ)) THEN
      ZFLOWLIQ(JJ,0) = 0.
    ELSE
      ZFLOWLIQ(JJ,0) = PRR(JJ) * PTSTEP / XRHOLW
      ZFLOWLIQ(JJ,0) = MAX(0., ZFLOWLIQ(JJ,0) - PLEL3L(JJ)*PTSTEP/(XLVTT*XRHOLW))
    END IF
  ELSE
    ZFLOWLIQ(JJ,0) = 0
  END IF
END DO
!
DO JST=1,IMAX_USE ! loop JST active snow layers
  DO JJ = 1,SIZE(PSNOWDZ,1)  ! loop JJ grid points
    IF (JST<=KNLVLS_USE(JJ))THEN
      !
      !
      ! 2. Increases Liquid Water from the upper layers flow (or rain for top layer)
      !    -----------------------------
      PSNOWLIQ(JJ,JST) = PSNOWLIQ(JJ,JST) + ZFLOWLIQ(JJ,JST-1)
      ZSNOWSWE(JJ,JST) = ZSNOWSWE(JJ,JST) + XRHOLW*ZFLOWLIQ(JJ,JST-1) ! Actualise the swe to take into account the new water arrival
      !Used to scavenge impurities
      ! 3. Freezes liquid water in any cold layers
      !    ---------------------------------------
      !
      ! Calculate the maximum possible refreezing
      ZPHASE(JJ,JST) = MIN( PSCAP(JJ,JST)* MAX(0.0, XTT - ZSNOWTEMP(JJ,JST)) * PSNOWDZ(JJ,JST), &
                            PSNOWLIQ(JJ,JST) * XLMTT * XRHOLW )
      !
      ! Reduce liquid content if freezing occurs:
      ZSNOWLIQ(JJ,JST) = PSNOWLIQ(JJ,JST) - ZPHASE(JJ,JST)/(XLMTT*XRHOLW)
      !
      ! Warm layer and reduce liquid if freezing occurs:
      ZSNOWDZ(JJ,JST) = MAX(XSNOWDMIN/INLVLS, PSNOWDZ(JJ,JST))
      !
      !
      ! Difference with ISBA-ES: a possible cooling of current refreezing water
      !                          is taken into account to calculate temperature change
      ZSOLIDMASS_BEFORE_REFR =  ( ZSNOWRHO(JJ,JST) * ZSNOWDZ(JJ,JST) - &
                                ( PSNOWLIQ(JJ,JST) - ZFLOWLIQ(JJ,JST-1) ) * XRHOLW )
      ZSOLIDMASS_AFTER_REFR =  ( ZSNOWRHO(JJ,JST) * ZSNOWDZ(JJ,JST) - &
                               ( ZSNOWLIQ(JJ,JST) - ZFLOWLIQ(JJ,JST-1) ) * XRHOLW )
      ! Note that it can be easily shown that ZSOLIDMASS_AFTER_REFR = ZSOLIDMASS_BEFORE_REFR + ( PSNOWLIQ - ZSNOWLIQ )  * XRHOLW
      ! (PSNOWLIQ - ZSNOWLIQ) * XRHOLW represents the refrozen mass (kg/m2)
      !
      ! Energy conservation applied on solid mass after refreezing:
      ! ( ZSNOWTEMP - PSNOWTEMP ) * ZSOLIDMASS_AFTER_REFR * XCI + XLMTT * ( ZSOLIDMASS_AFTER_REFR - ZSOLIDMASS_BEFORE_REFR ) + ( XTT - ZSNOWTEMP ) * ( ZSOLIDMASS_AFTER_REFR - ZSOLIDMASS_BEFORE_REFR ) * XCI = 0
      ! First term corresponds to heating of mass ZSOLIDMASS_AFTER_REFR
      ! Second term corresponds to refreezing of mass ( ZSOLIDMASS_AFTER_REFR - ZSOLIDMASS_BEFORE_REFR ) with XLMTT * ( ZSOLIDMASS_AFTER_REFR - ZSOLIDMASS_BEFORE_REFR ) = ZPHASE
      ! Third term corresponds to the cooling of mass ( ZSOLIDMASS_AFTER_REFR - ZSOLIDMASS_BEFORE_REFR ) for the equilibrium of the solid fraction
      ! After developing and simplifying this equation, we obtain the temperature evolution equation:
      PSNOWTEMP(JJ,JST) = XTT + ( ZSNOWTEMP(JJ,JST)-XTT )*ZSOLIDMASS_BEFORE_REFR/ZSOLIDMASS_AFTER_REFR + &
                          ZPHASE(JJ,JST)/( XCI*ZSOLIDMASS_AFTER_REFR )
      !
      ! 4. Calculate flow from the excess of holding capacity
      !    --------------------------------------------------------------
      !
      ! Any water in excess of the maximum holding space for liquid water
      ! amount is drained into next layer down.
      ZFLOWLIQ(JJ,JST) = MAX( 0., ZSNOWLIQ(JJ,JST)-ZWHOLDMAX(JJ,JST) )
      !
      !Compute the mass of impurity scavenged to layer down. Scavenging coefficient is the percentage of impurity scavenged
      ! the ratio flowliq/swe is used to determine the proportion of water leaving the layer
      !and the same proportion of impurity is supposed to leave the layer as well.
      DO JIMP=1,NIMPUR
        ZFLOWIMPUR(JJ,JST,JIMP) = MAX( 0., XSCAVEN_COEF(JIMP)*PSNOWIMPUR(JJ,JST,JIMP)*&
        ((ZFLOWLIQ(JJ,JST)*XRHOLW)/ZSNOWSWE(JJ,JST)))
      END DO
      ZSNOWLIQ(JJ,JST) = ZSNOWLIQ(JJ,JST) - ZFLOWLIQ(JJ,JST)
      !
      ! 5. Density is adjusted to conserve the mass
      !    --------------------------------------------------------------
      ZFINALMASS =  ( ZSNOWRHO(JJ,JST) * PSNOWDZ(JJ,JST) - ( ZFLOWLIQ(JJ,JST) - ZFLOWLIQ(JJ,JST-1) ) * XRHOLW )
      !
      ZSNOWRHO(JJ,JST) = ZFINALMASS / ZSNOWDZ(JJ,JST)
      ! keeps snow denisty below ice density
      IF ( ZSNOWRHO(JJ,JST)>XRHOLI ) THEN
        PSNOWDZ (JJ,JST) = PSNOWDZ(JJ,JST) * ZSNOWRHO(JJ,JST) / XRHOLI
        ZSNOWRHO(JJ,JST) = XRHOLI
      END IF
      !
      ! 6. Update thickness and density and any freezing:
      !    ----------------------------------------------
      PSNOWRHO(JJ,JST) = ZSNOWRHO(JJ,JST)
      PSNOWLIQ(JJ,JST) = ZSNOWLIQ(JJ,JST)
      !
   END IF
  END DO
END DO   ! loop JST active snow layers
  !
  ! Any remaining throughflow after freezing is available to
  ! the soil for infiltration or surface runoff (m).
  ! I.E. This is the amount of water leaving the snowpack:
  ! Rate water leaves the snowpack [kg/(m2 s)]:
  !
DO JJ = 1,SIZE(PSNOWDZ,1)  ! loop JJ grid points
  PTHRUFAL(JJ)  = PTHRUFAL(JJ) + ZFLOWLIQ(JJ,KNLVLS_USE(JJ)) * XRHOLW / PTSTEP
  !
END DO ! loop JJ grid points
!
! Impurity scavenging
DO JIMP=1,NIMPUR
  DO JJ = 1,SIZE(PSNOWDZ,1)
    DO JST=1,KNLVLS_USE(JJ)
        PSNOWIMPUR(JJ,JST,JIMP)=PSNOWIMPUR(JJ,JST,JIMP) + ZFLOWIMPUR(JJ,JST-1,JIMP)
        PSNOWIMPUR(JJ,JST,JIMP)=PSNOWIMPUR(JJ,JST,JIMP)-ZFLOWIMPUR(JJ,JST,JIMP)
    END DO
  END DO
END DO
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROREFRZ',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWCROREFRZ
!####################################################################
SUBROUTINE GET_RHO(PRHO_IN,PDZ,PSNOWLIQ,PFLOWLIQ,PRHO_OUT)
!
USE MODD_CSTS,     ONLY : XRHOLW
!
IMPLICIT NONE
!
REAL, INTENT(IN)  :: PRHO_IN, PDZ, PSNOWLIQ,PFLOWLIQ
REAL, INTENT(OUT) :: PRHO_OUT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:GET_RHO',0,ZHOOK_HANDLE)
!
PRHO_OUT =  ( PRHO_IN * PDZ - ( PSNOWLIQ - PFLOWLIQ ) * XRHOLW )
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:GET_RHO',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_RHO
!####################################################################
!####################################################################
SUBROUTINE SNOWCROFLUX(PSNOWTEMP,PSNOWDZ,PEXNS,PEXNA,          &
                       PUSTAR2_IC,                             &
                       PTSTEP,PALBT,PSW_RAD,PEMIST,            &
                       PLW_RAD,PTA,PSFCFRZ,PQA,PHPSNOW,        &
                       PSNOWTEMPO1,PSNOWFLUX,PCT,PRADSINK,     &
                       PQSAT,PDQSAT,PRSRA,                     &
                       PRN,PH,PGFLUX,PLES3L,PLEL3L,PEVAP,      &
                       PUSTAR, OFRZRAIN, PRR                   )
!
!!    PURPOSE
!!    -------
!     Calculate the surface fluxes (atmospheric/surface).
!     (Noilhan and Planton 1989; Noilhan and Mahfouf 1996)
!
USE MODD_CSTS,ONLY : XTT
!
USE MODE_THERMOS
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)      :: PSNOWDZ, PSNOWTEMPO1, PSNOWFLUX, PCT, &
                                        PRADSINK, PEXNS, PEXNA
!
REAL, DIMENSION(:), INTENT(IN)      :: PALBT, PSW_RAD, PEMIST, PLW_RAD,      &
                                        PTA, PSFCFRZ, PQA,                    &
                                        PHPSNOW, PQSAT, PDQSAT, PRSRA,        &
                                        PUSTAR2_IC, PRR
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PSNOWTEMP
!
REAL, DIMENSION(:), INTENT(OUT)     :: PRN, PH, PGFLUX, PLES3L, PLEL3L,      &
                                        PEVAP, PUSTAR
LOGICAL, DIMENSION(:), INTENT(IN)      :: OFRZRAIN
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWDZ))      :: ZEVAPC, ZSNOWTEMP
REAL :: ZSMSNOW, ZGFLUX
!
INTEGER :: JJ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROFLUX',0,ZHOOK_HANDLE)
!
! 0. Initialize:
! --------------
!
! 1. Flux calculations when melt not occuring at surface (W/m2):
! --------------------------------------------------------------
!
DO JJ = 1,SIZE(PALBT)
  !
  CALL GET_FLUX(PALBT(JJ),PEMIST(JJ),PSW_RAD(JJ),PLW_RAD(JJ),   &
                PEXNS(JJ),PEXNA(JJ),PTA(JJ),PQA(JJ),PRSRA(JJ), &
                PQSAT(JJ),PDQSAT(JJ),PSFCFRZ(JJ),PHPSNOW(JJ),   &
                PSNOWTEMP(JJ),PSNOWTEMPO1(JJ),                  &
                PRN(JJ),PH(JJ),ZEVAPC(JJ),                      &
                PLES3L(JJ),PLEL3L(JJ),ZGFLUX,PRR(JJ),OFRZRAIN(JJ))
  !
  IF ( PSNOWTEMP(JJ)>XTT ) THEN
    !
    IF ( PSNOWTEMPO1(JJ)<XTT ) THEN
      !
      ! 2. Initial melt adjustment
      ! --------------------------
      ! If energy avalabile to melt snow, then recalculate fluxes
      ! at the freezing point and add residual heat to layer
      ! average heat.
      !
      ! A) If temperature change is > 0 and passes freezing point this timestep,
      !    then recalculate fluxes at freezing point and excess energy
      !    will be used outside of this routine to change snow heat content:
      !
      ! WRITE (*,*) 'attention test LFLUX traitement XTT supprime!'
      !
      CALL GET_FLUX(PALBT(JJ),PEMIST(JJ),PSW_RAD(JJ),PLW_RAD(JJ),   &
                    PEXNS(JJ),PEXNA(JJ), PTA(JJ),PQA(JJ),PRSRA(JJ), &
                    PQSAT(JJ),PDQSAT(JJ),PSFCFRZ(JJ),PHPSNOW(JJ),   &
                    XTT,PSNOWTEMPO1(JJ),                            &
                    PRN(JJ),PH(JJ),ZEVAPC(JJ),                      &
                    PLES3L(JJ),PLEL3L(JJ),PGFLUX(JJ),PRR(JJ),OFRZRAIN(JJ))
      !
      ZSMSNOW = ZGFLUX - PGFLUX(JJ)
      !
      ! This will be used to change heat content of snow:
      ZSNOWTEMP(JJ) = PSNOWTEMP(JJ) - ZSMSNOW * PTSTEP * PCT(JJ)
      !
    ELSE
      !
      ! 3. Ongoing melt adjustment: explicit solution
      ! ---------------------------------------------
      !    If temperature change is 0 and at freezing point this timestep,
      !    then recalculate fluxes and surface temperature *explicitly*
      !    as this is *exact* for snow at freezing point (Brun, Martin)
      !
      CALL GET_FLUX(PALBT(JJ),PEMIST(JJ),PSW_RAD(JJ),PLW_RAD(JJ),   &
                    PEXNS(JJ),PEXNA(JJ), PTA(JJ),PQA(JJ),PRSRA(JJ), &
                    PQSAT(JJ),PDQSAT(JJ),PSFCFRZ(JJ),PHPSNOW(JJ),   &
                    XTT,XTT,                                        &
                    PRN(JJ),PH(JJ),ZEVAPC(JJ),                      &
                    PLES3L(JJ),PLEL3L(JJ),PGFLUX(JJ),PRR(JJ),OFRZRAIN(JJ))
      !
      ZSNOWTEMP(JJ) = XTT + PTSTEP * PCT(JJ) * ( PGFLUX(JJ) + PRADSINK(JJ) - PSNOWFLUX(JJ) )
      !
    END IF
    !
  ELSE
    !
    ZSNOWTEMP(JJ) = PSNOWTEMP(JJ)
    !
    PGFLUX(JJ) = ZGFLUX
    !
  END IF
  !
END DO
!
! 4. Update surface temperature:
! ------------------------------
!
PSNOWTEMP(:) = ZSNOWTEMP(:)
!
! 5. Final evaporative flux (kg/m2/s)
!
PEVAP(:) = ZEVAPC(:)
!WRITE(*,*) 'Flux PGFLUX:',PGFLUX(:)
!WRITE(*,*) 'Flux PRN:',PRN(:)
!WRITE(*,*) 'Flux PH:',PH(:)
!!WRITE(*,*) 'Flux surface:',PHPSNOW(:)
!!WRITE(*,*) 'Flux surface2:',PEXNS(:)
!!WRITE(*,*) 'Flux surface2:',PEXNA(:)
!WRITE(*,*) 'Flux PUSTAR2_IC:',PUSTAR2_IC(:)
!WRITE(*,*) 'Flux PTSTEP,PALBT:',PTSTEP,PALBT(:)
!!WRITE(*,*) 'Flux surface2:',PSW_RAD(:)
!!WRITE(*,*) 'Flux surface2:',PEMIST(:)
!!WRITE(*,*) 'Flux surface3:',PLW_RAD(:)
!!WRITE(*,*) 'Flux surface3:',PTA(:)
!!WRITE(*,*) 'Flux surface3:',PSFCFRZ(:)
!!WRITE(*,*) 'Flux surface3:',PQA(:)
!WRITE(*,*) 'Flux PSNOWTEMPO1:',PSNOWTEMPO1(:)
!WRITE(*,*) 'Flux PSNOWFLUX:',PSNOWFLUX(:)
!WRITE(*,*) 'Flux PCT:',PCT(:)
!!WRITE(*,*) 'Flux surface3:',PRADSINK(:)
!WRITE(*,*) 'Flux PQSAT:',PQSAT(:)
!WRITE(*,*) 'Flux PDQSAT:',PDQSAT(:)
!WRITE(*,*) 'Flux PRSRA:',PRSRA(:)
!WRITE(*,*) 'Flux PLES3L:',PLES3L(:)
!!WRITE(*,*) 'Flux surface4:',PLEL3L(:)
!WRITE(*,*) 'Flux PEVAP:',PEVAP(:)
!!WRITE(*,*) 'Flux surface4:',PUSTAR(:)

!
! 5. Friction velocity
! --------------------
!
PUSTAR(:) = SQRT(PUSTAR2_IC(:))
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROFLUX',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWCROFLUX
!####################################################################
SUBROUTINE GET_FLUX(PALBT,PEMIST,PSW_RAD,PLW_RAD,PEXNS,PEXNA,   &
                    PTA,PQA,PRSRA,PQSAT,PDQSAT,PSFCFRZ,PHPSNOW, &
                    PSNOWTEMP,PSNOWTEMPO1,                      &
                    PRN,PH,PEVAPC,PLES3L,PLEL3L,PGFLUX,PRR,OFRZRAIN)
!
USE MODD_CSTS,ONLY : XSTEFAN, XCPD, XLSTT, XLVTT, XLMTT
!
IMPLICIT NONE
!
REAL, INTENT(IN) :: PALBT, PEMIST
REAL, INTENT(IN) :: PSW_RAD, PLW_RAD
REAL, INTENT(IN) :: PEXNS, PEXNA
REAL, INTENT(IN) :: PTA, PQA, PRSRA, PQSAT, PDQSAT, PSFCFRZ, PHPSNOW
REAL, INTENT(IN) :: PSNOWTEMP,PSNOWTEMPO1
REAL, INTENT(IN) :: PRR
LOGICAL, INTENT(IN) :: OFRZRAIN
REAL, INTENT(OUT):: PRN, PH, PEVAPC, PLES3L, PLEL3L, PGFLUX
!
REAL :: ZLE, ZDELTAT, ZLWUPSNOW, ZSNOWTO3
REAL :: ZPHIFRZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:GET_FLUX',0,ZHOOK_HANDLE)
!
ZSNOWTO3  = PSNOWTEMPO1**3.  ! to save some CPU time, store this
!
ZDELTAT   = PSNOWTEMP - PSNOWTEMPO1   ! surface T time change:
!
ZLWUPSNOW = PEMIST * XSTEFAN * ZSNOWTO3 * ( PSNOWTEMPO1 + 4.*ZDELTAT )
!
PRN       = ( 1.-PALBT )*PSW_RAD + PEMIST*PLW_RAD - ZLWUPSNOW
!
PH        = PRSRA * XCPD * ( PSNOWTEMP/PEXNS - PTA/PEXNA )
!
PEVAPC    = PRSRA * ( (PQSAT - PQA) + PDQSAT*ZDELTAT )
!
PLES3L    = PSFCFRZ      * XLSTT * PEVAPC
!
PLEL3L    = (1.-PSFCFRZ) * XLVTT * PEVAPC
!
ZLE       = PLES3L + PLEL3L
!
IF (OFRZRAIN) THEN
  ZPHIFRZ = (XTT-PTA)*XCL/XLMTT
  PGFLUX    = PRN - PH - ZLE + PHPSNOW + (1.-ZPHIFRZ)*XLMTT*PRR
ELSE
  PGFLUX    = PRN - PH - ZLE + PHPSNOW
END IF
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:GET_FLUX',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_FLUX
!
!####################################################################
!####################################################################
SUBROUTINE SNOWCROEVAPN(PLES3L,PTSTEP,PSNOWTEMP,PDRYDENSITY,&
                       PSNOWDZ,PSNOWLIQ,PEVAPCOR,PSNOWHMASS        )
!
!!    PURPOSE
!!    -------
!     Remove mass from uppermost snow layer in response to
!     evaporation (liquid) and sublimation.
!
!!     MODIFICATIONS
!!     -------------
!!      Original A. Boone
!!      05/2011: E. Brun  Takes only into account sublimation and solid
!!                         condensation. Evaporation and liquid condensation
!!                         are taken into account in SNOWCROREFRZ
!
USE MODD_CSTS,     ONLY : XLSTT, XLMTT, XCI, XTT, XRHOLW
USE MODD_SNOW_PAR, ONLY : XRHOSMIN_ES
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)      :: PSNOWTEMP
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWLIQ
!
REAL, DIMENSION(:), INTENT(IN)      :: PLES3L   ! (W/m2)
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PDRYDENSITY ! Dry density of the first layer
REAL, DIMENSION(:), INTENT(INOUT)   :: PSNOWDZ, PSNOWHMASS, PEVAPCOR
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PLES3L))       :: ZSNOWEVAPS,          &
                                       ZSNOWDZ, ZEVAPCOR
!                                      ZEVAPCOR = for vanishingy thin snow cover,
!                                                 allow any excess evaporation
!                                                 to be extracted from the soil
!                                                 to maintain an accurate water
!                                                 balance [kg/(m2 s)]
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROEVAPN',0,ZHOOK_HANDLE)
!
! 0. Initialize:
! --------------
!
ZEVAPCOR  (:) = 0.0
ZSNOWEVAPS(:) = 0.0
ZSNOWDZ   (:) = 0.0
!
DO JJ = 1,SIZE(PSNOWTEMP)

  IF(PSNOWDZ(JJ)>0.0) THEN
   !
   ! 1. Sublimation/condensation of snow ice
   ! ----------------------------------------
   ! Reduce layer thickness and total snow depth
   ! if sublimation: add to correction term if potential
   ! sublimation exceeds available snow cover.
   !
   ! Mass change occurs at constant dry density
   ! Corresponding depth change below
     IF(PLES3L(JJ)>0) THEN ! Sublimation occurs
        ZSNOWEVAPS(JJ) = PLES3L(JJ) * PTSTEP / ( XLSTT*PDRYDENSITY(JJ) )
     ELSE ! Solid condensation occurs
        ! New mass is added at a minimum dry density of 50 kg/m3. 
        ZSNOWEVAPS(JJ) = PLES3L(JJ) * PTSTEP / ( XLSTT*MAX(XRHOSMIN_ES,PDRYDENSITY(JJ)) )
     ENDIF

     ZSNOWDZ(JJ)    = PSNOWDZ(JJ) - ZSNOWEVAPS(JJ)

     IF(PLES3L(JJ)>0) THEN ! Sublimation occurs 
        ! Ensure thickness of upper snow layer remains positive or zero      
        PSNOWDZ(JJ)    = MAX( 0.0, ZSNOWDZ(JJ) )
        ZEVAPCOR(JJ)   = ZEVAPCOR(JJ) + MAX(0.0,-ZSNOWDZ(JJ)) * PDRYDENSITY(JJ) / PTSTEP
     ELSE  ! Solid condensation occurs
        ! Update dry density       
        PDRYDENSITY(JJ) = (PDRYDENSITY(JJ)*PSNOWDZ(JJ) + MAX(XRHOSMIN_ES,PDRYDENSITY(JJ))* ABS(ZSNOWEVAPS(JJ))) / &
                             ZSNOWDZ(JJ) 
        PSNOWDZ(JJ)    = MAX( 0.0, ZSNOWDZ(JJ) )  
     ENDIF
        
     !
     !
     ! Total heat content change due to snowfall and sublimation (added here)JJ
     ! (for budget calculations)JJ
     !
     PSNOWHMASS(JJ) = PSNOWHMASS(JJ) &
                  - PLES3L(JJ) * (PTSTEP/XLSTT) * ( XCI * (PSNOWTEMP(JJ)-XTT) - XLMTT )

  ENDIF
ENDDO
!
! exceeding liquid water transferred to the next layer in case of total sublimation
WHERE (PSNOWDZ == 0.)
  PSNOWLIQ(:,2) = PSNOWLIQ(:,2) + PSNOWLIQ(:,1)
END WHERE
!
! 3. Update evaporation correction term:
! --------------------------------------
!
PEVAPCOR(:) = PEVAPCOR(:) + ZEVAPCOR(:)
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROEVAPN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOWCROEVAPN
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOWCROGONE(PTSTEP,PLEL3L,PLES3L,PSNOWRHO,                     &
                       PSNOWHEAT,PRADSINK_2D,PEVAPCOR,PTHRUFAL,PGRNDFLUX, &
                       PGFLUXSNOW,PSNOWDZ,PSNOWLIQ,PSNOWTEMP,PRADXS,      &
                       PRR,KNLVLS_USE, OFRZRAIN                           )
!
!!    PURPOSE
!!    -------
!     Account for the case when the last trace of snow melts
!     during a time step: ensure mass and heat balance of
!     snow AND underlying surface.
!     Original A. Boone
!     05/2011: E. Brun  Takes into account sublimation and PGRNDFLUX
!                       Adds rain and evaporation/liquid condensation
!                       in PTHRUFAL
!
USE MODD_CSTS,ONLY : XTT, XLSTT, XLVTT
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)      :: PLEL3L, PLES3L, PGFLUXSNOW,PRR
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PRADSINK_2D
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWRHO, PSNOWHEAT
!
LOGICAL, DIMENSION(:), INTENT(IN) :: OFRZRAIN
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PGRNDFLUX, PRADXS
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDZ, PSNOWLIQ, PSNOWTEMP
!
REAL, DIMENSION(:), INTENT(OUT)     :: PTHRUFAL   ! melt water [kg/(m2 s)]
!
REAL, DIMENSION(:), INTENT(OUT)     :: PEVAPCOR   ! [kg/(m2 s)]
!                                      PEVAPCOR = for vanishingy thin snow cover,
!                                                 allow any excess evaporation
!                                                 to be extracted from the soil
!                                                 to maintain an accurate water
!                                                 balance.
!
INTEGER, DIMENSION(:), INTENT(INOUT)      :: KNLVLS_USE
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PLES3L))       :: ZRADSINK
REAL, DIMENSION(SIZE(PLES3L))       :: ZSNOWHEATC
REAL, DIMENSION(SIZE(PLES3L))       :: ZSUBLIM ! mass of snow that sublimates kg m-2
REAL, DIMENSION(SIZE(PLES3L))       :: ZSNOWHEATSUBLIM ! enthalpy of snow that sublimates W m-2
INTEGER, DIMENSION(SIZE(PLES3L))    :: ISNOWGONE_DELTA
REAL :: ZSNOWSWE, ZSNOWSWEC
!
LOGICAL, DIMENSION(SIZE(PLES3L))    :: GALLSUBLIM, GALLMELT
INTEGER  :: JJ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROGONE',0,ZHOOK_HANDLE)
!
! 0. Initialize:
! --------------
!
PEVAPCOR(:) = 0.0
PTHRUFAL(:) = 0.0
!
DO JJ = 1,SIZE(ZRADSINK)
  ZRADSINK  (JJ) = PRADSINK_2D(JJ,INLVLS_USE(JJ))
  ZSNOWHEATC(JJ) = SUM(PSNOWHEAT(JJ,1:INLVLS_USE(JJ))) !total heat content (J m-2)
END DO
!
ISNOWGONE_DELTA(:) = 1
!
! 1. Test to see if snow vanishes:
! ---------------------------------------
! If so, set thicknesses (and therefore mass and heat) and liquid content
! to zero, and adjust fluxes of water, evaporation and heat into underlying
! surface.
!
! takes into account the heat content corresponding to the occasional
! sublimation  and then PGRNDFLUX
!
! Mass that should be sublimated if available
ZSUBLIM(:) = MAX( 0., PLES3L(:)*PTSTEP ) / XLSTT
!
! Enthalpy of sublimated snow (cumulate enthalpy of layers from the surface until
! sublimation is reached)
DO JJ = 1,SIZE(ZRADSINK)
  ZSNOWSWEC = 0.
  ZSNOWHEATSUBLIM(JJ) = 0.
  DO JST =1, KNLVLS_USE(JJ)
    ZSNOWSWE = PSNOWDZ(JJ,JST)*PSNOWRHO(JJ,JST)-PSNOWLIQ(JJ,JST)*XRHOLW
    ZSNOWSWEC = ZSNOWSWEC + ZSNOWSWE
    IF (ZSNOWSWEC < ZSUBLIM(JJ)) THEN
      ! Layer fully sublimated
      ZSNOWHEATSUBLIM(JJ) = ZSNOWHEATSUBLIM(JJ) + PSNOWHEAT(JJ,JST)
    ELSE
      ! Layer partially sublimated
      ZSNOWHEATSUBLIM(JJ) = ZSNOWHEATSUBLIM(JJ) &
      + PSNOWHEAT(JJ,JST)*(ZSNOWSWE+ZSUBLIM(JJ)-ZSNOWSWEC)/ZSNOWSWE
      EXIT
    END IF
  END DO
END DO
! Enthalpy of not sublimated snow
ZSNOWHEATC(:) = ZSNOWHEATC(:) - ZSNOWHEATSUBLIM
! Logical true if all snow is sublimated
GALLSUBLIM(:) = ZSNOWHEATC(:) >= 0.
! Logical true if all non-sublimated snow melts
GALLMELT(:)   = PGFLUXSNOW(:)+ZRADSINK(:)-PGRNDFLUX(:) >= (-ZSNOWHEATC(:)/PTSTEP)
!
! Where snow disappears, transfer energy and mass to soil and remove all snow layers
WHERE (GALLMELT .OR. GALLSUBLIM)
  PGRNDFLUX(:)       = PGFLUXSNOW(:) + (ZSNOWHEATC(:)/PTSTEP)
  PEVAPCOR (:)       = PLES3L(:)/XLSTT
  PRADXS   (:)       = 0.0
  ISNOWGONE_DELTA(:) = 0          ! FLAG...if=0 then snow vanishes, else=1
END WHERE
!
! 2. Final update of snow state and computation of corresponding flow
!    Only if snow vanishes
! -----------------------------
!
!
DO JJ=1, SIZE(ZRADSINK)
  !
  IF(ISNOWGONE_DELTA(JJ) == 0 ) THEN
    PTHRUFAL(JJ) = PTHRUFAL(JJ) + &
                   SUM( PSNOWRHO(JJ,1:INLVLS_USE(JJ))*PSNOWDZ(JJ,1:INLVLS_USE(JJ)) ) / PTSTEP
    ! takes into account rain and condensation/evaporation
    ! if freezing rain, rainfall not added to water flow)
    IF ( OFRZRAIN(JJ) ) THEN
      PTHRUFAL(JJ) = PTHRUFAL(JJ) - PLEL3L(JJ)/XLVTT
    ELSE
      PTHRUFAL(JJ) = PTHRUFAL(JJ) + PRR(JJ) - PLEL3L(JJ)/XLVTT
    END IF
    PSNOWTEMP(JJ,:) = XTT
    PSNOWDZ  (JJ,:) = 0.
    PSNOWLIQ (JJ,:) = 0.
    INLVLS_USE(JJ) = 0
  END IF
  !
END DO
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROGONE',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWCROGONE
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOWCROEVAPGONE(PSNOWHEAT,PSNOWDZ,PDRYDENSITY,PSNOWTEMP,PSNOWLIQ,      &
                           PSNOWDIAMOPT,PSNOWSPHERI,PSNOWHIST,PSNOWAGE,           &
                           KNLVLS_USE, HSNOWMETAMO)
!
!!    PURPOSE
!!    -------
!
!     If all snow in uppermost layer evaporates/sublimates, re-distribute
!     grid (below assumes very thin snowpacks so layer-thicknesses are
!     constant).
!     Original A. Boone
!     05/2011: E. Brun  Takes into account previous changes in the energy
!                       content
!     04/2023: M. Lafaysse Bug fix: microstructure was not updated
!
!
USE MODD_CSTS,     ONLY : XTT, XRHOLW, XLMTT, XCI
USE MODD_SNOW_PAR, ONLY : XRHOSMIN_ES, XSNOWDMIN, XRHOSMAX_ES
USE MODE_SNOW3L
USE MODD_SNOW_METAMO
USE MODD_TYPE_DATE_SURF
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PDRYDENSITY! snow dry density profile                (kg/m3)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWDZ    ! snow layer thickness profile        (m)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWHEAT  ! snow heat content/enthalpy          (J/m2)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWDIAMOPT ! snow grain parameter 1              (-)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWSPHERI ! snow grain parameter 2              (-)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWHIST  ! snow grain historical variable      (-)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWAGE  ! Snow grain age
!
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWTEMP  ! snow temperature profile            (K)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWLIQ   ! snow liquid water profile           (m)
!
INTEGER, DIMENSION(:), INTENT(INOUT)      :: KNLVLS_USE
CHARACTER(3), INTENT(IN)              :: HSNOWMETAMO ! metamorphism scheme
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWDZ,1),SIZE(PSNOWDZ,2))   :: ZSNOWRHO   ! snow total density profile                (kg/m3)
REAL, DIMENSION(SIZE(PSNOWDZ,1)) :: ZSNOWHEAT_1D ! total heat content                (J/m2)
REAL, DIMENSION(SIZE(PSNOWDZ,1)) :: ZSNOWRHO_1D  ! total snowpack average density    (kg/m3)
REAL, DIMENSION(SIZE(PSNOWDZ,1)) :: ZSNOW        ! total snow depth                  (m)
REAL, DIMENSION(SIZE(PSNOWDZ,1)) :: ZSCAP        ! Snow layer heat capacity          (J/K/m3)
INTEGER, DIMENSION(SIZE(PSNOWDZ,1)) :: INDENT       ! Number of dendritic layers        (-)
INTEGER, DIMENSION(SIZE(PSNOWDZ,1)) :: INVIEU       ! Number of non dendritic layers    (-)
REAL, DIMENSION(SIZE(PSNOWDZ,1)) :: ZSNOWAGE_1D  ! total snowpack average
!age (days)
REAL, DIMENSION(SIZE(PSNOWDZ,1),SIZE(PSNOWDZ,2)) :: ZSNOWDIAMOPTN,           &
                                                      ZSNOWSPHERIN,ZSNOWHISTN !
INTEGER :: JJ, JST          ! loop control
!
LOGICAL :: GVANISHSOMEWHERE ! Logical to avoid computations when possible
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROEVAPGONE',0,ZHOOK_HANDLE)
!
! Initialize:
!
ZSNOWHEAT_1D(:) = 0.
ZSNOW(:)        = 0.
ZSNOWRHO_1D(:)  = 0.
INDENT(:)       = 0
INVIEU(:)       = 0
ZSNOWAGE_1D(:)  = 0.
ZSCAP(:)        = 0.
!
! First, determine where uppermost snow layer has completely
! evaporated/sublimated (as it becomes thin):
GVANISHSOMEWHERE = .FALSE.
DO JJ = 1,SIZE(PSNOWRHO,1)
   !
   IF ( PSNOWDZ(JJ,1)==0.0 ) THEN
    IF ( KNLVLS_USE(JJ)==1 ) THEN
      KNLVLS_USE(JJ) = 0
    ELSE            
     !
     GVANISHSOMEWHERE = .TRUE.
     !
     DO JST = 2,IMAX_USE
       IF (JST <= KNLVLS_USE(JJ)) THEN
         ! Update total density 
         ZSNOWRHO(JJ,JST)  = PDRYDENSITY(JJ,JST) +  PSNOWLIQ(JJ,JST) * XRHOLW / PSNOWDZ(JJ,JST)
         !
         ZSNOWHEAT_1D(JJ) = ZSNOWHEAT_1D(JJ) + PSNOWDZ(JJ,JST) * &
                            ( ZSNOWRHO(JJ,JST)*XCI * (ZSNOWTEMP(JJ,JST)-XTT) &
                              - XLMTT * PSNOWRHO(JJ,JST) ) &
                            + XLMTT * XRHOLW * PSNOWLIQ(JJ,JST)
         ZSNOW       (JJ) = ZSNOW       (JJ) + PSNOWDZ(JJ,JST)
         ZSNOWRHO_1D (JJ) = ZSNOWRHO_1D (JJ) + PSNOWDZ(JJ,JST) * ZSNOWRHO(JJ,JST)
         ZSNOWAGE_1D (JJ) = ZSNOWAGE_1D (JJ) + PSNOWDZ(JJ,JST) * ZSNOWRHO(JJ,JST) * PSNOWAGE(JJ,JST)
         !
         IF (  PSNOWDIAMOPT(JJ,JST)<XVDIAM6*(4.-PSNOWSPHERI(JJ,JST))-XUEPSI  ) THEN   ! Dendritic snow
           INDENT(JJ) = INDENT(JJ) + 1
         ELSE                                    ! Non dendritic snow
           INVIEU(JJ) = INVIEU(JJ) + 1
         END IF
         !
       END IF
     END DO
    ENDIF
    !
   END IF
   !
END DO
!
ZSNOWRHO_1D(:) = ZSNOWRHO_1D (:) / MAX( XSNOWDMIN, ZSNOW(:) )
ZSNOWAGE_1D(:) = ZSNOWAGE_1D (:) / MAX( XSNOWDMIN, ZSNOW(:) * ZSNOWRHO_1D(:) )
ZSNOWRHO_1D(:) = MAX( XRHOSMIN_ES, MIN( XRHOSMAX_ES, ZSNOWRHO_1D(:) ) )
!
! Where uppermost snow layer has vanished, redistribute vertical
! snow mass and heat profiles (and associated quantities):
!
IF (GVANISHSOMEWHERE) THEN ! Following is only useful in case of full sublimation of surface layer
 CALL SNOW3LAVGRAIN(PSNOWDIAMOPT,PSNOWSPHERI,PSNOWHIST,                 &
                    ZSNOWDIAMOPTN,ZSNOWSPHERIN,ZSNOWHISTN,INDENT,INVIEU,KNLVLS_USE)
!
 DO JJ=1,SIZE(PSNOWRHO,1)
  !
  IF( ZSNOW(JJ)/=0.0 ) THEN
    !
    PSNOWDZ  (JJ,1:KNLVLS_USE(JJ)) = ZSNOW(JJ) / KNLVLS_USE(JJ)
    PSNOWHEAT(JJ,1:KNLVLS_USE(JJ)) = ZSNOWHEAT_1D(JJ) / KNLVLS_USE(JJ)
    ZSNOWRHO (JJ,1:KNLVLS_USE(JJ)) = ZSNOWRHO_1D(JJ)
    !
    ! Bug fix: following instructions were missing between 2011 and 2023
    PSNOWDIAMOPT (JJ,1:KNLVLS_USE(JJ)) = ZSNOWDIAMOPTN(JJ,1:KNLVLS_USE(JJ))
    PSNOWSPHERI (JJ,1:KNLVLS_USE(JJ)) = ZSNOWSPHERIN(JJ,1:KNLVLS_USE(JJ))
    PSNOWHIST (JJ,1:KNLVLS_USE(JJ)) = ZSNOWHISTN(JJ,1:KNLVLS_USE(JJ))    
    !
    ZSCAP(JJ) = ZSNOWRHO_1D(JJ) * XCI
    !
    DO JST = 1,IMAX_USE
      IF (JST <= KNLVLS_USE(JJ)) THEN
        !
        PSNOWTEMP(JJ,JST) = XTT + ( ( (PSNOWHEAT(JJ,JST)/PSNOWDZ(JJ,JST))   &
                                      + XLMTT*PSNOWRHO(JJ,JST) ) / ZSCAP(JJ) )
        PSNOWTEMP(JJ,JST) = MIN( XTT, PSNOWTEMP(JJ,JST) )
        !
        PSNOWLIQ (JJ,JST) = MAX( 0.0, PSNOWTEMP(JJ,JST)-XTT ) * ZSCAP(JJ) * &
                                      PSNOWDZ(JJ,JST) / (XLMTT*XRHOLW)
        !
        PDRYDENSITY(JJ,JST) = ZSNOWRHO(JJ,JST) -  PSNOWLIQ(JJ,JST) * XRHOLW / PSNOWDZ(JJ,JST)
      END IF
    END DO
    !
  END IF
  !
 END DO
END IF
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROEVAPGONE',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWCROEVAPGONE
!
!####################################################################
!###################################################################
SUBROUTINE SNOWDRIFT(PTSTEP,PVMOD,PSNOWRHO,PSNOWDZ,PSNOW,HSNOWMETAMO,HSNOWDRIFT,        &
                     HSNOWCOMP,PSNOWDIAMOPT,PSNOWSPHERI,PSNOWHIST,KNLVLS_USE, &
                     PTA,PQA,PPS,PRHOA,PZ0EFF,PUREF,             &
                     OSNOWDRIFT_SUBLIM,PSNDRIFT, PVFRIC_T,PHVEGPOL  )
!
!!    PURPOSE
!!    -------
!     Snow compaction  and metamorphism due to drift
!     Mass is unchanged: layer thickness is reduced
!     in proportion to density increases. Method inspired from
!     Brun et al. (1997) and Guyomarch
!
!     - computes a mobility index of each snow layer from its grains, density
!                 and history
!     - computes a drift index of each layer from its mobility and wind speed
!     - computes a transport index with an exponential decay taking into
!                 account its depth and the mobility of upper layers
!     - increases density and changes grains in case of transport
!
!     HISTORY:
!     Basic parameterization from Crocus/ARPEGE Coupling (1997)
!     Implementation in V5
!     Insertion in V6 of grains type evolution in case of dendritic snow (V.
!     Vionnet)
!     07/2012 (for V7.3): E. Brun, M. Lafaysse : optional sublimation of drifted snow
!     2012-09-20 : bug correction : ZFF was not computed if LSNOWDRIFT_SUBLIM=FALSE.
!
!     2014-02-05 V. Vionnet: systematic use of 5m wind speed to compute drift index
!     2014-06-03 M. Lafaysse: threshold on PZ0EFF
!     2021       M. Baron: new option B21
!     2023-03-23 M. Lafaysse: microstructure evolution laws rewritten according to
!                the new GMD paper (including bug fix in the non dendritic-case)
!     2023-11    M. Baron : new option HSNOWDRIFT=PAPP intended to be used with snowpappus

USE MODD_CSTS,ONLY : XTT, XKARMAN
USE MODE_THERMOS

USE MODD_SNOW_PAR, ONLY : XVTIME, XVROMAX, XVROMIN,XVROMAX_R21, XVMOB1, XZ0SN,  &
                          XVMOB2, XVMOB3, XVMOB4, XVDRIFT1, XVDRIFT2, XVDRIFT3, &
                          XVSIZEMIN, XCOEF_EFFECT,XCOEF_EFFECT_R21, XQS_REF,    &
                          XCROCOEF_FF
USE MODE_SNOW3L, ONLY : GETGRAINSIZE_B21
USE MODD_SNOW_METAMO, ONLY : NVHIS2
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)      :: PTA, PQA, PPS, PRHOA
!
REAL, DIMENSION(:), INTENT(IN)      :: PVMOD
!
INTEGER, DIMENSION(:), INTENT(IN)   :: KNLVLS_USE
!
REAL, DIMENSION(:),INTENT(IN)       :: PZ0EFF,PUREF
!
CHARACTER(3),INTENT(IN)             :: HSNOWMETAMO, HSNOWCOMP
!
CHARACTER(4),INTENT(IN)             :: HSNOWDRIFT
!
LOGICAL,INTENT(IN)                  :: OSNOWDRIFT_SUBLIM
!
REAL, DIMENSION(:), INTENT(IN)      :: PHVEGPOL
!
REAL, DIMENSION(:), INTENT(IN)      :: PVFRIC_T 
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWRHO, PSNOWDZ,PSNOWDIAMOPT, &
                                       PSNOWSPHERI,PSNOWHIST
REAL, DIMENSION(:), INTENT(OUT)     :: PSNOW
REAL, DIMENSION(:), INTENT(OUT)     :: PSNDRIFT !blowing snow sublimation (kg/m2/s)
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWRHO2, &
                                                      ZSNOW_JST  ! Snow layer height (m) Royer et al 2021
REAL, DIMENSION(SIZE(PSNOWRHO,1)                 ) :: ZSNOWDZ1
!
REAL, DIMENSION(SIZE(PSNOWRHO,1))   :: ZQSATI, ZFF ! QS wrt ice, gust speed
!
REAL     :: ZZ0EFF
!
REAL     :: ZPROFEQU, ZRMOB, ZRDRIFT, ZRT, ZDRO, ZDGR1, ZDGR2
REAL     :: ZVT ! 5m wind speed threshold for surface
!transport
REAL     :: ZQS_EFFECT ! effect of QS on snow
REAL     :: ZWIND_EFFECT ! effect of wind on snow
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2))     :: ZDRIFT_EFFECT ! effect of QS and wind on snow
! transformation
REAL     :: ZQS !Blowing snow sublimation (kg/m2/s)
REAL     :: ZRHI, ZFACT
REAL     :: ZSNOWSIZE ! Size of grain snow in non-dendritic case
REAL     :: ZDEND ! Dendricity
REAL     :: ZDDIAM ! Variation of optical diameter
REAL     :: ZDSPHER ! Variation of sphericity
REAL     :: ZSPHERP1 ! sphericity + 1 to avoid multiple computations
!
INTEGER  :: JJ,JST   ! looping indexes
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Reference height for wind speed used to dertermine the occurrence of blowing snow
REAL, PARAMETER :: PPHREF_WIND=5.
REAL, PARAMETER :: PPHREF_MIN=PPHREF_WIND/2.
LOGICAL         :: GDENDRITIC
!                  GDENDRITIC = to test if the layer is dendritic
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWDRIFT',0,ZHOOK_HANDLE)
!
! 0. Initialization:
! ------------------
!
PSNDRIFT(:) = 0.0
ZSNOWDZ1(:) = PSNOWDZ(:,1)
!
DO JST = 1,IMAX_USE
  DO JJ = 1,SIZE(PSNOW)
    IF (JST <= KNLVLS_USE(JJ)) THEN
      ZSNOWRHO2(JJ,JST) = PSNOWRHO(JJ,JST)
    END IF
  END DO
END DO
!
IF ( OSNOWDRIFT_SUBLIM ) THEN
  ZQSATI(:) = QSATI( PTA(:),PPS(:) )
END IF
!
! 1. Computation of drift and induced settling and metamorphism
! ------------------
!
DO JJ=1, SIZE(PSNOW)
  !
  IF (HSNOWCOMP == 'R21' .OR. HSNOWCOMP == 'R2D') THEN
    ! Calculate the height of each snow layer (Royer et al 2021)
    ZSNOW_JST(JJ,1) = PSNOW(JJ) - 0.5*PSNOWDZ(JJ,1)
    DO JST=2, KNLVLS_USE(JJ)
       ZSNOW_JST(JJ,JST) = ZSNOW_JST(JJ,JST-1) - 0.5 * (PSNOWDZ(JJ,JST-1) + PSNOWDZ(JJ,JST))
    ENDDO
  ENDIF
  !  
  ! gust speed at 5m above the snowpack
  ! Computed from PVMOD at PUREF (m) assuming a log profile in the SBL
  ! and a roughness length equal to PZ0EFF
  ZZ0EFF=MIN(PZ0EFF(JJ),PUREF(JJ)*0.5,PPHREF_MIN)
  ZFF(JJ) = XCROCOEF_FF*PVMOD(JJ)*LOG(PPHREF_WIND/ZZ0EFF)/LOG(PUREF(JJ)/ZZ0EFF)
  ! XCROCOEF_FF = 1.0 by default from M.Baron, V.Vionnet and M. Lafaysse
  ! It can change for arctic simulations. In ISBA_ES, XCOEF_FF is 1.25
  !
  ! initialization decay coeff
  ZPROFEQU = 0.
  !
  DO JST = 1,IMAX_USE
    IF (JST <= KNLVLS_USE(JJ)) THEN
      ! calcul de l'indice de mobilité
      IF (HSNOWDRIFT == 'PAPP') THEN
        !conversion de la vitesse de friction seuil de SnowPappus en indice de mobilité
        !vitesse de friction => vitesse du vent à 5 m (eq. du vent log) => indice de mobilité 
        !(application de l'eq. 11 de Vionnet et al. 2012 en considérant SI=0 <=> U=Ut)
        ZRMOB = 2.868 *(5./XZ0SN)**(-XVDRIFT2*PVFRIC_T(JJ)/XKARMAN)-1
      ELSE
        ZFACT = 1.25 - 1.25 * ( MAX( PSNOWRHO(JJ,JST), XVROMIN ) - XVROMIN )/1000./XVMOB1
        !
        IF ( PSNOWDIAMOPT(JJ,JST)<XVDIAM6*(4.-PSNOWSPHERI(JJ,JST))-XUEPSI ) THEN
          ! dendritic case
          ZRMOB = 0.34 * ( 0.5 + 0.75 * &
                 ( PSNOWDIAMOPT(JJ,JST)/XVDIAM6-4.+PSNOWSPHERI(JJ,JST) )/( PSNOWSPHERI(JJ,JST)-3. ) &
                               - 0.5 * PSNOWSPHERI(JJ,JST) ) + &
                  0.66 * ZFACT
        ELSE
          ! non dendritic case
          IF ((HSNOWMETAMO == 'B21') .OR. (HSNOWMETAMO == 'F06').OR. &
             (HSNOWMETAMO == 'S-F').OR.(HSNOWMETAMO == 'T07'))  THEN
             CALL GETGRAINSIZE_B21(PSNOWDIAMOPT(JJ,JST),PSNOWSPHERI(JJ,JST),ZSNOWSIZE)
          ELSEIF (PSNOWDIAMOPT(JJ,JST) - 0.0004*(1+PSNOWSPHERI(JJ,JST)) >= 0) THEN
              ZSNOWSIZE  = 2.*PSNOWDIAMOPT(JJ,JST)/(1+  PSNOWSPHERI(JJ,JST))
          ELSEIF (ABS(PSNOWSPHERI(JJ,JST))< XUEPSI) THEN
              ZSNOWSIZE = 0.0008
          ELSE
              ZSNOWSIZE = (1./PSNOWSPHERI(JJ,JST))*(PSNOWDIAMOPT(JJ,JST)-0.0004*(1-PSNOWSPHERI(JJ,JST)))
          END IF
  
          ZRMOB = 0.34 * ( XVMOB2 - XVMOB2 * PSNOWSPHERI(JJ,JST) &
                                  - XVMOB3 * ZSNOWSIZE ) + &
                  0.66 * ZFACT
        END IF
      END IF
      !
      ! correction in case of former wet snow
      IF ( NINT(PSNOWHIST(JJ,JST)) >= NVHIS2 ) ZRMOB = MIN(ZRMOB, XVMOB4)
      !
      ! computation of drift index supposing no overburden snow
      ZRDRIFT = ZRMOB - ( XVDRIFT1 * EXP( -XVDRIFT2*ZFF(JJ) ) - 1.)
      !
      ! Effect of polar low vegetation that prevent the effect of snowdrift
      IF (HSNOWCOMP == 'R21' .OR. HSNOWCOMP == 'R2D') THEN 
        IF (ZSNOW_JST(JJ,JST) <= PHVEGPOL(JJ)) THEN
             ZRDRIFT = 0.
        ENDIF
      ENDIF
      ! modif_EB exit loop if there is no drift
      IF ( ZRDRIFT<=0. ) EXIT
      !
      ! update the decay coeff by half the current layer
      ZPROFEQU = ZPROFEQU + 0.5 * PSNOWDZ(JJ,JST) * ( XVDRIFT3 - ZRDRIFT )
      ! computation of the drift index inclunding the decay by overburden snow
      ZRT = MAX( 0., ZRDRIFT * EXP( -ZPROFEQU * 10.) )
      !
      IF ( OSNOWDRIFT_SUBLIM .AND. JST==1 ) THEN
        !Specific case for blowing snow sublimation
        ! computation of wind speed threshold QSATI and RH withe respect to ice
        ZVT  = -LOG( (ZRMOB+1.)/XVDRIFT1 ) / XVDRIFT2
        ZRHI = PQA(JJ) / ZQSATI(JJ)
        ! computation of sublimation rate according to Gordon's PhD
        ZQS = 0.0018 * (XTT/PTA(JJ))**4. * ZVT * PRHOA(JJ) * ZQSATI(JJ) * &
              (1.-ZRHI) * (ZFF(JJ)/ZVT)**3.6
        ! WRITE(*,*) 'surface Vt vent*coef  ZRDRIFT ZRMOB :',ZVT,&
        ! ZFF(JJ),ZRDRIFT,ZRMOB
        ! WRITE(*,*) 'V>Vt ZQS   :',ZQS
        ! surface depth decrease in case of blowing snow sublimation
        ! WRITE(*,*) 'V>Vt DSWE DZ Z:',- MAX(0.,ZQS)*PTSTEP/COEF_FF,
        ! - MAX(0.,ZQS)*PTSTEP/COEF_FF/PSNOWRHO(JJ,JST),PSNOWDZ(JJ,JST)
        ! 2 lignes ci-dessous a valider pour avoir sublim drift
        !VV PSNOWDZ(JJ,JST) = MAX( 0.5*PSNOWDZ(JJ,JST), &
        !VV                       PSNOWDZ(JJ,JST) - MAX(0.,ZQS) * PTSTEP/PSNOWRHO(JJ,JST) )
        PSNOWDZ(JJ,JST) = MAX( 0.5*PSNOWDZ(JJ,JST), &
                               PSNOWDZ(JJ,JST) - MAX(0.,ZQS) * PTSTEP/(XCROCOEF_FF*PSNOWRHO(JJ,JST)) )
        PSNDRIFT(JJ) = (ZSNOWDZ1(JJ)-PSNOWDZ(JJ,JST))*PSNOWRHO(JJ,JST)/PTSTEP
      ELSE
        ZQS = 0.
      END IF
      !
      ZQS_EFFECT    = MIN( 3., MAX( 0.,ZQS )/XQS_REF ) * ZRT
      IF (HSNOWDRIFT=='R21W' .OR. HSNOWDRIFT=='R21F') THEN
        ZWIND_EFFECT = XCOEF_EFFECT_R21 * ZRT
      ELSE
        ZWIND_EFFECT  = XCOEF_EFFECT * ZRT
      ENDIF      
      !ZDRIFT_EFFECT(JJ,JST) = ( ZQS_EFFECT + ZWIND_EFFECT ) * PTSTEP / XVTIME
      ZDRIFT_EFFECT(JJ,JST) = ( ZQS_EFFECT + ZWIND_EFFECT ) * PTSTEP / (XVTIME*XCROCOEF_FF)
      ! WRITE(*,*) 'ZQS_EFFECT,ZWIND_EFFECT,ZDRIFT_EFFECT:',ZQS_EFFECT,ZWIND_EFFECT,ZDRIFT_EFFECT
      !
      ! settling by wind transport only in case of not too dense snow
      IF (HSNOWDRIFT=='R21R' .OR. HSNOWDRIFT=='R21F') THEN 
        IF( PSNOWRHO(JJ,JST) <= XVROMAX_R21 ) THEN 
          ZDRO = ZDRIFT_EFFECT(JJ,JST) * ( XVROMAX_R21 - PSNOWRHO(JJ,JST) )
          PSNOWRHO(JJ,JST) = MIN( XVROMAX_R21 , PSNOWRHO(JJ,JST) + ZDRO )
          PSNOWDZ (JJ,JST) = PSNOWDZ(JJ,JST) * ZSNOWRHO2(JJ,JST) / PSNOWRHO(JJ,JST)
        ENDIF
      ELSE      
        IF( PSNOWRHO(JJ,JST) < XVROMAX ) THEN
          ZDRO = ZDRIFT_EFFECT(JJ,JST) * ( XVROMAX - PSNOWRHO(JJ,JST) )
          PSNOWRHO(JJ,JST) = MIN( XVROMAX , PSNOWRHO(JJ,JST) + ZDRO )
          PSNOWDZ (JJ,JST) = PSNOWDZ(JJ,JST) * ZSNOWRHO2(JJ,JST) / PSNOWRHO(JJ,JST)
        END IF
      END IF


        !metamorphism with new proposition from M.Baron
        !
        ! Equation 69
        ! ds/dt = (1-s)/tau
        ZDSPHER = 1-PSNOWSPHERI(JJ,JST)
        !
        CALL CHECK_DENDRITIC( PSNOWDIAMOPT(JJ,JST),PSNOWSPHERI(JJ,JST),GDENDRITIC)
        !
        IF (GDENDRITIC) THEN
          ! Dendricity
          ZDEND = ( PSNOWDIAMOPT(JJ,JST)/XVDIAM6-4.+PSNOWSPHERI(JJ,JST) )/ &
                                                    (PSNOWSPHERI(JJ,JST)-3.)
          ! Equation 68
          ! ddend/dt = -dend/(2tau) ; ddiam/ddend = 1E-4(s-3)
          ! ds/dt = (1-s)/tau ; ddiam/ds = 1E-4 (dend - 1)
          ! This gives Eq. E4 which can be simplified in Eq. 68 implemented below
          ZDDIAM = ZDRIFT_EFFECT(JJ,JST) * XVDIAM6 * &
                   ( (2.5 - 1.5 * PSNOWSPHERI(JJ,JST)) * ZDEND - 1. + PSNOWSPHERI(JJ,JST) )
          !
          PSNOWDIAMOPT(JJ,JST) = MAX(PSNOWDIAMOPT(JJ,JST) + ZDDIAM, XVDIAM6)
          ! can never decrease below 0.1 mm
        ELSE
          !non-dendritic case ( new proposition from M.Baron, adapted from Vionnet 2012 )
          ! Equation 68
          ! dg/dt = -5E-4/tau ; ddiam/dg = (S-1)/2
          ! ds/dt = (1-s)/tau ; ddiam/ds = (g-4E-4)/2 = (d-4E-4)/(1+s)
          ZSPHERP1 = PSNOWSPHERI(JJ,JST) + 1
          !
          ZDDIAM=ZDRIFT_EFFECT(JJ,JST)* (-5.*XVDIAM6 * ZSPHERP1/2. + &
                  (PSNOWDIAMOPT(JJ,JST)-4*XVDIAM6)*ZDSPHER/ZSPHERP1)
          !
          PSNOWDIAMOPT(JJ,JST)= MAX(PSNOWDIAMOPT(JJ,JST)+ZDDIAM, XVDIAM6*(4-PSNOWSPHERI(JJ,JST))) !evolution of optical diameter ( proposed by M. Baron )
          ! Can not decrease below the dendritic threshold
        END IF
        !
        PSNOWSPHERI(JJ,JST) = MIN(1.,PSNOWSPHERI(JJ,JST)+ZDRIFT_EFFECT(JJ,JST)*ZDSPHER)
        ! update the decay coeff by half the current layer
        ZPROFEQU = ZPROFEQU + 0.5 * PSNOWDZ(JJ,JST) * ( XVDRIFT3 - ZRDRIFT )
        !
      !
    END IF
    !
  END DO  ! snow layers loop
  !
END DO    ! grid points loop
!
!
! 2. Update total snow depth:
! -----------------------------------------------
!
! Compaction of total snowpack depth
!
DO JJ = 1,SIZE(PSNOWDZ,1)
  PSNOW(JJ) = SUM( PSNOWDZ(JJ,1:KNLVLS_USE(JJ)) )
END DO
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWDRIFT',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWDRIFT
!####################################################################
!####################################################################
SUBROUTINE SNOWCROLAYER_GONE(PTSTEP,PSCAP,PSNOWHEAT,PSNOWTEMP,PSNOWDZ,        &
                             PSNOWRHO,PSNOWLIQ,PSNOWDIAMOPT,PSNOWSPHERI,      &
                             PSNOWHIST,PSNOWAGE,PSNOWIMPUR,PLES3L,KNLVLS_USE  )
!
!
!!    PURPOSE
!     Account for the case when one or several snow layers melt
!     during a time step:
!     in that case, merge these layers with the underlying layer
!     except for the bottom layer which is merged to the abovelying layer
!     energy and mass are conserved
!     a new merged layer keeps the grain, histo and age properties of the
!     non-melted layer
!
USE MODD_CSTS,ONLY : XTT, XLMTT, XRHOLW, XCI
!
USE MODE_SNOW3L
!
USE MODD_PREP_SNOW, ONLY : NIMPUR
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                     :: PTSTEP
!
REAL, DIMENSION(:,:), INTENT(INOUT)  :: PSCAP
!
REAL, DIMENSION(:,:), INTENT(INOUT)  :: PSNOWDZ, PSNOWTEMP, PSNOWRHO, PSNOWLIQ, PSNOWHEAT
REAL, DIMENSION(:,:), INTENT(INOUT)  :: PSNOWDIAMOPT,PSNOWSPHERI,PSNOWHIST,PSNOWAGE
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PSNOWIMPUR
!
INTEGER, DIMENSION(:), INTENT(INOUT) :: KNLVLS_USE !
!
REAL, DIMENSION(:), INTENT(IN) :: PLES3L
!
!*      0.2    declarations of local variables
!
REAL :: ZHEAT, ZMASS, ZDZ, ZLIQ, ZSNOWLWE
REAL, DIMENSION(NIMPUR) :: ZSNOWIMPUR
!
INTEGER :: JJ,JST,JST_1, JST_2, JST_MAX, IDIFF_LAYER ! loop counter
INTEGER :: ID_1, ID_2
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROLAYER_GONE',0,ZHOOK_HANDLE)
!
DO JJ=1,SIZE(PSNOWRHO,1)  ! loop on gridpoints
  !
  JST_MAX = KNLVLS_USE(JJ)
  !
  IDIFF_LAYER = 0  ! used as shift counter of previously melted layers
  !
  DO JST_1 = JST_MAX,1-JST_MAX,-1 ! loop on 2 x layers in case of multi melt
    !
    JST = JST_1 + IDIFF_LAYER
    !
    ! Merge is possible only in case of 2 active layers or more
    IF ( JST>=1 .AND. KNLVLS_USE(JJ)>1 ) THEN
      !
      ! Total Liquid equivalent water content of snow (m):
      ZSNOWLWE = PSNOWRHO(JJ,JST) * PSNOWDZ(JJ,JST) / XRHOLW
      !
      ! Consideration of sublimation if any
      IF ( JST==1 ) ZSNOWLWE = ZSNOWLWE - MAX( 0., PLES3L(JJ)*PTSTEP/(XLSTT*XRHOLW) )
      !
      ! Test if avalaible energy exceeds total latent heat
      !
      ! Modif ML : aggregate if the resulting layer would be of depth < 1.E-8 m after melting
      ! (for conservation because otherwiseat the next time step these layers
      ! would have been ignored in SNOWNLGRIDFRESH_1D)
      !
      IF ( PSCAP(JJ,JST) * MAX( 0.0, PSNOWTEMP(JJ,JST)-XTT ) * PSNOWDZ(JJ,JST) >=  &
           ( ( ZSNOWLWE-PSNOWLIQ(JJ,JST) ) * XLMTT * XRHOLW ) - &
!VV           XUEPSI * XLMTT * PSNOWRHO(JJ,JST) ) THEN
           1.2*XUEPSI_SMP * XLMTT * PSNOWRHO(JJ,JST) ) THEN
        !
        IF ( JST==KNLVLS_USE(JJ) ) THEN
          ID_1 = JST-1
          ID_2 = JST
        ELSE
          ID_1 = JST
          ID_2 = JST + 1
        END IF
        !
        ! Case of a total melt of the bottom layer: merge with above layer
        !        which keeps its grain, histo and age properties
        ZHEAT = 0.
        ZMASS = 0.
        ZDZ   = 0.
        ZLIQ  = 0.
        DO JIMP=1,NIMPUR
          ZSNOWIMPUR(JIMP)=0.
        END DO
        DO JST_2 = ID_1,ID_2
          ZHEAT = ZHEAT + &
                  PSNOWDZ(JJ,JST_2) * &
                  ( PSCAP(JJ,JST_2)*( PSNOWTEMP(JJ,JST_2)-XTT ) - XLMTT*PSNOWRHO(JJ,JST_2) ) + &
                  XLMTT * XRHOLW * PSNOWLIQ(JJ,JST_2)
          ZMASS = ZMASS + PSNOWDZ(JJ,JST_2) * PSNOWRHO(JJ,JST_2)
          ZDZ   = ZDZ   + PSNOWDZ(JJ,JST_2)
          ZLIQ  = ZLIQ  + PSNOWLIQ(JJ,JST_2)
          ! Compute the total amount of impurity present in the melting layer + the over/underlaying layer
          DO JIMP=1,NIMPUR
            ZSNOWIMPUR(JIMP)= ZSNOWIMPUR(JIMP)+PSNOWIMPUR(JJ,JST_2,JIMP)
          END DO
        END DO
        !
        PSNOWDZ  (JJ,ID_1) = ZDZ
        PSNOWRHO (JJ,ID_1) = ZMASS / ZDZ
        PSNOWLIQ (JJ,ID_1) = ZLIQ
        !
        ! Temperature of the merged layer is deduced from the heat content
        !     Comment M Lafaysse : Before 2020 a threshold was applied on PSNOWDZ (XSNOWDZMIN=1E-4)
        !     in the computation of PSCAP. This leads to incorrect energy computation in SNOWCROREFRZ
        !     when refreezing occurs on very thin layers
        PSCAP    (JJ,ID_1) = ( PSNOWRHO(JJ,ID_1) - &
                               PSNOWLIQ(JJ,ID_1) * XRHOLW / &
                               PSNOWDZ(JJ,ID_1) ) * XCI
        PSNOWHEAT(JJ,ID_1) = ZHEAT
        PSNOWTEMP(JJ,ID_1) = XTT + &
          ( ( ( ( ZHEAT - XLMTT*XRHOLW*PSNOWLIQ(JJ,ID_1) ) / PSNOWDZ(JJ,ID_1) ) + &
              XLMTT*PSNOWRHO(JJ,ID_1) ) &
            / PSCAP(JJ,ID_1) )
        ! The section on impurity management is volountary put in the loop JST/=KNLVLS_USE(JJ) as explained bellow.
        !
        IF( JST/=KNLVLS_USE(JJ) ) THEN
          !
          PSNOWDIAMOPT(JJ,JST) = PSNOWDIAMOPT(JJ,JST+1)
          PSNOWSPHERI (JJ,JST) = PSNOWSPHERI (JJ,JST+1)
          PSNOWHIST   (JJ,JST) = PSNOWHIST   (JJ,JST+1)
          PSNOWAGE    (JJ,JST) = PSNOWAGE    (JJ,JST+1)
          ! The impurity content of the underlaying layer is equal to the sum of its current content+ the content of the melting layer.
          ! We put this instruction in the loop JST/=KNLVLS_USE(JJ) because in the case of the bottom layer, we don't want a transfer
          ! of the impurity content to the overlaying layer (not physical). In that particular case the impurity content is discarded by Crocus.
          DO JIMP=1,NIMPUR
            PSNOWIMPUR(JJ,JST,JIMP)= ZSNOWIMPUR(JIMP)
          END DO
          !
          ! Shift the above layers
          DO JST_2 = JST+1,KNLVLS_USE(JJ)-1
            PSNOWHEAT   (JJ,JST_2) = PSNOWHEAT   (JJ,JST_2+1)
            PSNOWTEMP   (JJ,JST_2) = PSNOWTEMP   (JJ,JST_2+1)
            PSCAP       (JJ,JST_2) = PSCAP       (JJ,JST_2+1)
            PSNOWDZ     (JJ,JST_2) = PSNOWDZ     (JJ,JST_2+1)
            PSNOWRHO    (JJ,JST_2) = PSNOWRHO    (JJ,JST_2+1)
            PSNOWLIQ    (JJ,JST_2) = PSNOWLIQ    (JJ,JST_2+1)
            PSNOWDIAMOPT(JJ,JST_2) = PSNOWDIAMOPT(JJ,JST_2+1)
            PSNOWSPHERI (JJ,JST_2) = PSNOWSPHERI (JJ,JST_2+1)
            PSNOWHIST   (JJ,JST_2) = PSNOWHIST   (JJ,JST_2+1)
            PSNOWAGE    (JJ,JST_2) = PSNOWAGE    (JJ,JST_2+1)

            DO JIMP=1,NIMPUR
              PSNOWIMPUR (JJ, JST_2,JIMP) =PSNOWIMPUR(JJ,JST_2+1,JIMP)
            END DO
            
          END DO !  loop JST_2
          !
          ! Update the shift counter IDIFF_LAYER
          IDIFF_LAYER = IDIFF_LAYER + 1
          !
        END IF ! end test of bottom layer
        !
        ! Decrease the number of active snow layers
        KNLVLS_USE(JJ) = KNLVLS_USE(JJ) - 1
        !
      END IF ! end test on availibility of energy
      !
    END IF ! end test on the number of remaining active layers
    !
  END DO ! end loop on the snow layers
  !
END DO ! end loop gridpoints
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROLAYER_GONE',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWCROLAYER_GONE
!####################################################################
!###################################################################
!####################################################################
!###################################################################
SUBROUTINE SNOWCROPRINTPROFILE(HINFO,KLAYERS,OPRINTGRAN,PSNOWDZ,PSNOWRHO, &
                               PSNOWLIQ,PSNOWHEAT,PSNOWDIAMOPT,   &
                               PSNOWSPHERI,PSNOWHIST,PSNOWAGE,HSNOWMETAMO, &
                               HSNOWRAD,PSNOWIMPUR)
!
! Matthieu Lafaysse 08/06/2012
! This routine prints the snow profile of a given point for debugging
!
!to compute SSA
USE MODD_CSTS, ONLY : XRHOLI, XTT, XLMTT, XCI
!
IMPLICIT NONE
!
 CHARACTER(*), INTENT(IN) :: HINFO
LOGICAL,       INTENT(IN) :: OPRINTGRAN
INTEGER,       INTENT(IN) :: KLAYERS
REAL, DIMENSION(:), INTENT(IN) :: PSNOWDZ,PSNOWRHO,PSNOWLIQ, &
                                  PSNOWHEAT,PSNOWDIAMOPT,PSNOWSPHERI,     &
                                  PSNOWHIST,PSNOWAGE
CHARACTER(3), INTENT(IN)       :: HSNOWMETAMO
CHARACTER(3), INTENT(IN),OPTIONAL       :: HSNOWRAD
REAL, DIMENSION(:,:), INTENT(IN),OPTIONAL :: PSNOWIMPUR
!
REAL, DIMENSION(KLAYERS) :: ZSNOWSSA
REAL, DIMENSION(KLAYERS) :: ZSCAP, ZSNOWTEMP, ZSNOWLIQ
LOGICAL::GPRINTIMPUR
CHARACTER(LEN=12), PARAMETER:: CCT = "------------"
!
INTEGER :: JST
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROPRINTPROFILE',0,ZHOOK_HANDLE)
!
WRITE(*,*)
WRITE(*,*)TRIM(HINFO)
!
IF (PRESENT(PSNOWIMPUR).AND. NIMPUR > 0) THEN
  GPRINTIMPUR= (HSNOWRAD=='T17')
ELSE
  GPRINTIMPUR=.FALSE.
END IF

! Compute snow temperature from enthalpy
ZSCAP     (1:KLAYERS) = PSNOWRHO(1:KLAYERS) * XCI
IF (ALL(PSNOWDZ(1:KLAYERS)>0.) .AND. ALL(ZSCAP(1:KLAYERS)>0.)) THEN
  ZSNOWTEMP (1:KLAYERS) = XTT + &
                          ( ( PSNOWHEAT(1:KLAYERS)/PSNOWDZ(1:KLAYERS) + XLMTT*PSNOWRHO(1:KLAYERS) )/ZSCAP(1:KLAYERS) )
  ZSNOWLIQ (1:KLAYERS) = MAX( 0.0, ZSNOWTEMP(1:KLAYERS)-XTT ) * ZSCAP(1:KLAYERS) * &
                          PSNOWDZ(1:KLAYERS) / (XLMTT*XRHOLW)
  ZSNOWTEMP(1:KLAYERS) = MIN( XTT, ZSNOWTEMP(1:KLAYERS) )
ELSE
  PRINT*, "WARNING: UNABLE TO COMPUTE ZSNOWTEMP"
  ZSNOWTEMP (1:KLAYERS) = -999.
END IF

IF (OPRINTGRAN) THEN
  !
  ! Compute SSA from SNOWDIAMOPT and SNOWSPHERI
  !
  IF (KLAYERS>0) THEN
    WHERE(PSNOWDIAMOPT(1:KLAYERS)>0)
      ZSNOWSSA(:) = 6. / (XRHOLI*PSNOWDIAMOPT(1:KLAYERS))
    ELSEWHERE
      ZSNOWSSA(:) = -999
    ENDWHERE
  END IF
  !
  !
  WRITE(*,'(10(A12,"|"))')CCT,CCT,CCT,CCT,CCT,CCT,CCT,CCT,CCT,CCT
  IF (GPRINTIMPUR) THEN
    WRITE(*,'(10(A12,"|"))')"PSNOWDZ","PSNOWRHO","PSNOWLIQ","ZSNOWLIQ","ZSNOWTEMP","PSNOWHEAT",&
        "PSNOWDIAMOPT","PSNOWSPHERI","PSNOWSSA","PSNOWIMPUR"
  ELSE
    WRITE(*,'(10(A12,"|"))')"PSNOWDZ","PSNOWRHO","PSNOWLIQ","ZSNOWLIQ","ZSNOWTEMP","PSNOWHEAT",&
        "PSNOWDIAMOPT","PSNOWSPHERI","PSNOWSSA","PSNOWAGE"
  END IF
  WRITE(*,'(10(A12,"|"))')CCT,CCT,CCT,CCT,CCT,CCT,CCT,CCT,CCT,CCT
  DO JST = 1,KLAYERS
    IF (GPRINTIMPUR) THEN
      WRITE(*,'(10(ES12.3,"|")," L",I2.2)') PSNOWDZ(JST),PSNOWRHO(JST),    &
                                          PSNOWLIQ(JST),ZSNOWLIQ(JST),ZSNOWTEMP(JST),PSNOWHEAT(JST),PSNOWDIAMOPT(JST), &
                                          PSNOWSPHERI(JST),ZSNOWSSA(JST),PSNOWIMPUR(JST,1),JST
    ELSE
      WRITE(*,'(10(ES12.3,"|")," L",I2.2)') PSNOWDZ(JST),PSNOWRHO(JST),    &
                                          PSNOWLIQ(JST),ZSNOWLIQ(JST),ZSNOWTEMP(JST),PSNOWHEAT(JST),PSNOWDIAMOPT(JST), &
                                          PSNOWSPHERI(JST),ZSNOWSSA(JST),PSNOWAGE(JST),JST
    END IF
  END DO
  WRITE(*,'(10(A12,"|"))')"-------------","-------------","-------------",&
        "-------------","-------------","-------------","-------------",&
        "-------------"
  !
ELSE
  !
  WRITE(*,'(6(A12,"|"))')CCT,CCT,CCT,CCT,CCT,CCT
  WRITE(*,'(6(A12,"|"))')"PSNOWDZ","PSNOWRHO","PSNOWLIQ","ZSNOWLIQ","ZSNOWTEMP","PSNOWHEAT"
  WRITE(*,'(6(A12,"|"))')"------------","------------","------------",&
        "------------"
  DO JST = 1,KLAYERS
    WRITE(*,'(6(ES12.3,"|")," L",I2.2)') PSNOWDZ(JST),PSNOWRHO(JST),&
                                         PSNOWLIQ(JST),ZSNOWLIQ(JST),ZSNOWTEMP(JST),PSNOWHEAT(JST),JST
  END DO
  WRITE(*,'(6(A12,"|"))')CCT,CCT,CCT,CCT,CCT,CCT
  !
END IF
!
WRITE(*,*)
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROPRINTPROFILE',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWCROPRINTPROFILE
!####################################################################
!###################################################################
SUBROUTINE SNOWCROPRINTATM(CINFO,PTA,PQA,PVMOD,PRR,PSR,PSW_RAD,PLW_RAD, &
                           PTG, PSOILCOND,PD_G,PPSN3L,PBLOWSNW          )

! Matthieu Lafaysse 08/06/2012
! This routine prints the atmospheric forcing of a given point for debugging
! and ground data

IMPLICIT NONE

 CHARACTER(*), INTENT(IN) :: CINFO
REAL,          INTENT(IN) :: PTA,PQA,PVMOD,PRR,PSR,PSW_RAD,PLW_RAD
REAL,          INTENT(IN) :: PTG, PSOILCOND, PD_G, PPSN3L,PBLOWSNW
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROPRINTATM',0,ZHOOK_HANDLE)
!
 CALL SNOWCROPRINTDATE()
!
WRITE(*,*)
WRITE(*,*)TRIM(CINFO)
WRITE(*,'(4(A12,"|"))')"------------","------------","------------",&
"------------"
WRITE(*,'(4(A12,"|"))')"PTA","PQA","PRR","PSR"
WRITE(*,'(4(A12,"|"))')"------------","------------","------------",&
"------------"
WRITE(*,'(4(ES12.3,"|")," meteo1")')PTA,PQA,PRR,PSR
WRITE(*,'(4(A12,"|"))')"------------","------------","------------",&
"------------"
WRITE(*,'(4(A12,"|"))')"------------","------------","------------",&
"------------"
WRITE(*,'(4(A12,"|"))')"PSW_RAD","PLW_RAD","PVMOD","PBLOWSNW"
WRITE(*,'(4(A12,"|"))')"------------","------------","------------",&
"------------"
WRITE(*,'(4(ES12.3,"|")," meteo2")')PSW_RAD,PLW_RAD,PVMOD,PBLOWSNW
WRITE(*,'(4(A12,"|"))')"------------","------------","------------",&
"------------"
WRITE(*,*)
WRITE(*,*)"Ground :"
WRITE(*,'(4(A12,"|"))')"------------","------------","------------",&
"------------"
WRITE(*,'(4(A12,"|"))')"PTG","PSOILCOND","PD_G","PPSN3L"
WRITE(*,'(4(A12,"|"))')"------------","------------","------------",&
"------------"
WRITE(*,'(4(ES12.3,"|")," soil")')PTG,PSOILCOND,PD_G,PPSN3L
WRITE(*,'(4(A12,"|"))')"------------","------------","------------",&
"------------"
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROPRINTATM',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWCROPRINTATM
!
!####################################################################
SUBROUTINE SNOWCROSTOPBALANCE(PMASSBALANCE,PENERGYBALANCE)
!
USE MODE_CRODEBUG, ONLY : XWARNING_MASSBALANCE, XWARNING_ENERGYBALANCE
!
USE MODI_ABOR1_SFX
!
! stop if energy and mass balances are not closed
!
IMPLICIT NONE
!
REAL , DIMENSION(:), INTENT(IN) :: PMASSBALANCE, PENERGYBALANCE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROSTOPBALANCE',0,ZHOOK_HANDLE)
!
IF ( ANY( PMASSBALANCE   > XWARNING_MASSBALANCE   ) ) &
        CALL ABOR1_SFX("SNOWCRO: WARNING MASS BALANCE !")
IF ( ANY( PENERGYBALANCE > XWARNING_ENERGYBALANCE ) ) &
        CALL ABOR1_SFX("SNOWCRO: WARNING ENERGY BALANCE !")
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROSTOPBALANCE',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWCROSTOPBALANCE
!
!###################################################################
SUBROUTINE SNOWCROPRINTBALANCE(PSUMMASS_INI,PSUMHEAT_INI,PSUMMASS_FIN,PSUMHEAT_FIN, &
                               PSR,PRR,PTHRUFAL,PEVAP,PEVAPCOR,PGRNDFLUX,PHSNOW,    &
                               PRNSNOW,PLEL3L,PLES3L,PHPSNOW,PSNOWHMASS,PSNOWDZ,    &
                               PTSTEP,PMASSBALANCE,PENERGYBALANCE,PEVAPCOR2         )
!
! Matthieu Lafaysse / Eric Brun 03/10/2012
! Print energy and mass balances.
!
IMPLICIT NONE
!
REAL, INTENT(IN) :: PSUMMASS_INI,PSUMHEAT_INI,PSUMMASS_FIN,PSUMHEAT_FIN
REAL, INTENT(IN) :: PSR,PRR,PTHRUFAL,PEVAP,PEVAPCOR
REAL, INTENT(IN) :: PGRNDFLUX,PHSNOW,PRNSNOW,PLEL3L,PLES3L,PHPSNOW,PSNOWHMASS
REAL, INTENT(IN) :: PSNOWDZ !first layer
REAL, INTENT(IN) :: PTSTEP !time step
REAL, INTENT(IN) :: PMASSBALANCE, PENERGYBALANCE, PEVAPCOR2
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROPRINTBALANCE',0,ZHOOK_HANDLE)
!
WRITE(*,*) ' '
WRITE(*,FMT='(A1,67("+"),A1)') "+","+"
!
 CALL SNOWCROPRINTDATE()
!
WRITE(*,*) ' '
!
! print des residus de bilan et des differents termes pour le point
WRITE (*,FMT="(A25,1x,E17.10)") 'final mass (kg/m2) =' , PSUMMASS_FIN
WRITE (*,FMT="(A25,1x,E17.10)") 'final energy (J/m2) =', ZSUMHEAT_FIN
WRITE(*,*) ' '
!
WRITE(*,FMT="(A25,1x,E17.10)") 'mass balance (kg/m2) =', PMASSBALANCE
!
WRITE(*,*) ' '
WRITE(*,FMT="(A35)") 'mass balance contribution (kg/m2) '
WRITE(*,FMT="(A51,1x,E17.10)") 'delta mass:', (PSUMMASS_FIN-PSUMMASS_INI)
WRITE(*,FMT="(A51,1x,E17.10)") 'hoar or condensation (>0 towards snow):', -PEVAP * PTSTEP
WRITE(*,FMT="(A51,1x,E17.10)") 'rain:',    PRR * PTSTEP
WRITE(*,FMT="(A51,1x,E17.10)") 'snow:',    PSR * PTSTEP
WRITE(*,FMT="(A51,1x,E17.10)") 'run-off:', PTHRUFAL * PTSTEP
WRITE(*,FMT="(A51,1x,E17.10)") 'evapcor:', PEVAPCOR * PTSTEP
!
WRITE(*,FMT='(A1,55("-"),A1)')"+","+"
WRITE(*,*) ' '
!
WRITE(*,FMT="(A25,4(1x,E17.10))") 'energy balance (W/m2)=',PENERGYBALANCE
!
WRITE(*,*) ' '
WRITE(*,FMT="(A55)") 'energy balance contribution (W/m2) >0 towards snow :'
WRITE(*,FMT="(A51,1x,E17.10)") 'delta heat:', (ZSUMHEAT_FIN-ZSUMHEAT_INI)/PTSTEP
WRITE(*,FMT="(A51,1x,E17.10)") 'radiation (LW + SW):', PRNSNOW
WRITE(*,FMT="(A51,1x,E17.10)") 'sensible flux :',      -PHSNOW
WRITE(*,FMT="(A51,1x,E17.10)") 'ground heat flux :',   -PGRNDFLUX
WRITE(*,FMT="(A51,1x,E17.10)") 'liquid latent flux:',  -PLEL3L
WRITE(*,FMT="(A51,1x,E17.10)") 'solid latent flux:',   -PLES3L
WRITE(*,FMT="(A51,1x,E17.10)") 'rain sensible heat:',  PHPSNOW
WRITE(*,FMT="(A51,1x,E17.10)") 'snowfall/hoar heat (sensible + melt heat):', PSNOWHMASS/PTSTEP
WRITE(*,FMT="(A51,1x,E17.10)") 'evapcor:', PEVAPCOR2
WRITE(*,FMT='(A1,67("+"),A1)')"+","+"
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROPRINTBALANCE',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWCROPRINTBALANCE
!
!####################################################################
SUBROUTINE GET_BALANCE(PSUMMASS_INI,PSUMHEAT_INI,PSUMMASS_FIN,PSUMHEAT_FIN, &
                       PSR,PRR,PTHRUFAL,PEVAP,PEVAPCOR,PGRNDFLUX,PHSNOW,    &
                       PRNSNOW,PLEL3L,PLES3L,PHPSNOW,PSNOWHMASS,PSNOWDZ,    &
                       PTSTEP,PMASSBALANCE,PENERGYBALANCE,PEVAPCOR2         )
!
IMPLICIT NONE
!
REAL, INTENT(IN) :: PSUMMASS_INI,PSUMHEAT_INI,PSUMMASS_FIN,PSUMHEAT_FIN
REAL, INTENT(IN) :: PSR,PRR,PTHRUFAL,PEVAP,PEVAPCOR
REAL, INTENT(IN) :: PGRNDFLUX,PHSNOW,PRNSNOW,PLEL3L,PLES3L,PHPSNOW,PSNOWHMASS
REAL, INTENT(IN) :: PSNOWDZ !first layer
REAL, INTENT(IN) :: PTSTEP !time step
!
REAL, INTENT(OUT)  :: PMASSBALANCE, PENERGYBALANCE, PEVAPCOR2
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:GET_BALANCE',0,ZHOOK_HANDLE)
!
PMASSBALANCE = PSUMMASS_FIN - PSUMMASS_INI - &
               ( PSR + PRR - PTHRUFAL - PEVAP + PEVAPCOR ) * PTSTEP
!
PEVAPCOR2 = PEVAPCOR * PSNOWDZ / MAX( XUEPSI,PSNOWDZ ) *  &
           ( ABS(PLEL3L) * XLVTT / MAX( XUEPSI,ABS(PLEL3L) ) + &
             ABS(PLES3L) * XLSTT / MAX( XUEPSI,ABS(PLES3L) ) )
!
PENERGYBALANCE = ( PSUMHEAT_FIN-PSUMHEAT_INI ) / PTSTEP - &
                 ( -PGRNDFLUX - PHSNOW + PRNSNOW - PLEL3L - PLES3L + PHPSNOW ) - &
                 PSNOWHMASS / PTSTEP - PEVAPCOR2
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:GET_BALANCE',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_BALANCE
!
!###################################################################
SUBROUTINE SNOWCROPRINTDATE()
!
IMPLICIT NONE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROPRINTDATE',0,ZHOOK_HANDLE)
!
WRITE(*,FMT='(I4.4,2("-",I2.2)," Hour=",F5.2)') &
  TPTIME%TDATE%YEAR, TPTIME%TDATE%MONTH, TPTIME%TDATE%DAY, TPTIME%TIME/3600.
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROPRINTDATE',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWCROPRINTDATE
!####################################################################
!###################################################################
SUBROUTINE SNOWCROGETSSA(PSNOWDIAMOPT,KLAYERS,PSNOWSSA)
!to compute SSA
USE MODD_CSTS, ONLY : XRHOLI
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:), INTENT(IN) :: PSNOWDIAMOPT
INTEGER,DIMENSION(:), INTENT(IN) :: KLAYERS
!
REAL, DIMENSION(:,:),INTENT(OUT) :: PSNOWSSA

!
INTEGER :: JJ,JST
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROGETSSA',0,ZHOOK_HANDLE)
!
!
PSNOWSSA = -999
DO JST = 1,IMAX_USE
  DO JJ=1,SIZE(PSNOWDIAMOPT,1)
    IF (JST<=KLAYERS(JJ)) THEN
      IF(PSNOWDIAMOPT(JJ,JST)>0)  &
        PSNOWSSA(JJ,JST) = 6. / (XRHOLI*PSNOWDIAMOPT(JJ,JST))
    END IF
  END DO
END DO
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROGETSSA',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWCROGETSSA
!####################################################################
!###################################################################
SUBROUTINE SNOWGROOMING(PSNOWDZ,PSNOWSWE,PSNOWAGE,PSNOWRHO, PSNOWDIAMOPT,&
                        PSNOWSPHERI, INLVLS_USE, PSMASSCOEFF,            &
                        OSNOWMAK_BOOL, OSNOWTILLER)
!!
!!    PURPOSE
!!    -------
!!    Apply grooming on snow density profile
!!    Initial implementation by P.spandre 2013/12/04
!!    Bug fixes from M. Lafaysse / C. Carmagnola 2023/03/03
!!
USE MODD_SNOW_METAMO
USE MODD_SNOW_PAR, ONLY : XSM_END, XFREQ_GRO
!
IMPLICIT NONE
!
!      0.1    declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWSWE, PSMASSCOEFF
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDZ, PSNOWAGE, PSNOWDIAMOPT, &
                                       PSNOWSPHERI, PSNOWRHO
!
INTEGER, DIMENSION(:), INTENT(IN)   :: INLVLS_USE
!
LOGICAL, INTENT(IN)                 :: OSNOWMAK_BOOL, OSNOWTILLER
!
!      0.2    declarations of local variables
!
INTEGER                             :: JJ,JST,II,SMT    ! looping indexes
!
REAL                                :: ZSNOWAGEB, ZSNOWDIAMOPTB, ZSNOWSPHERIB, &
                                       ZSNOWRHOB, ZSNOWRHOC, ZSNOWDZC
REAL, PARAMETER                     :: VSWE = 35.0
!                                      threshold for SWE max impacted layers
!
LOGICAL, DIMENSION(SIZE(PTA))       :: LTIMECOMPACT
LOGICAL                             :: PMONTH
LOGICAL                             :: PDAY
REAL, DIMENSION(SIZE(PSNOWDZ,1),SIZE(PSNOWDZ,2)) :: ZSNOWCOMPACT, ZSNOWRHOG, ZSNOWDZG
!                                                   extra pressure due to grooming kg m-2
!                                                   density recalculation variable
!                                                   layer thickness variable
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWGROOMING',0,ZHOOK_HANDLE)
!
II = 1
SMT = 1   ! SMT = 1 <=> GROOMING ONLY
          ! SMT = 3 <=> GROOMING + SNOWMAKING
PMONTH = .TRUE. !
!
DO JJ=1, SIZE(PSNOWDZ,1)
  !
  ! A. Timing conditions
  ! A.1. Month condition : grooming possible from november to april
  !
  ! Check if GROOMING ONLY or GROOM.+SM
  IF (OSNOWMAK_BOOL) THEN
    SMT = 3
  END IF
  !
  ! Check the closing date for that track
  IF (TPTIME%TDATE%MONTH == XSM_END(SMT) .and. TPTIME%TDATE%DAY > XSM_END(SMT+1)) THEN
    PMONTH = .FALSE.
  END IF
  !
  IF (TPTIME%TDATE%MONTH > XSM_END(SMT)) THEN
    PMONTH = .FALSE.
  END IF
  !
  ! Grooming in NOV and DEC always TRUE
  IF (TPTIME%TDATE%MONTH > 10.) THEN
    PMONTH = .TRUE.
  END IF
  !
  ! A.2. Daily condition : grooming possible from 6.15pm to 8.45pm
  IF (TPTIME%TIME > 64800. .and. TPTIME%TIME < 75600.) THEN
    PDAY = .TRUE.
  ELSE
    PDAY = .FALSE.
  END IF
  !
  ! A.3. Option : Grooming one day out of two (no grooming on even days)
  ! ---> Removed to model daily grooming, ps 2015/05/03
  IF (MOD(TPTIME%TDATE%DAY, XFREQ_GRO) /= 0.) THEN
    PDAY = .FALSE.
  END IF
  !
  ! A.3b Option : Grooming every day when snow is produced on this track
  IF (OSNOWMAK_BOOL) THEN
    PDAY = .TRUE.
  END IF
  !
  ! A.4. Night snowfall => Grooming possible in the morning as well (from 6.15am to 8.45am)
  IF (PSNOWAGE(JJ,1) < 0.5 .and. TPTIME%TIME > 21600. .and. TPTIME%TIME < 32400.) THEN
    PDAY = .TRUE.
  END IF
  !
  ! A.5. Day snowfall => Grooming possible even if not an even day => Grooming possible from 6.15pm to 8.45pm
  IF (PSNOWAGE(JJ,1) < 0.5 .and. TPTIME%TIME > 64800. .and. TPTIME%TIME < 75600.) THEN
    PDAY = .TRUE.
  END IF
  !
  ! A.6. Boolean from timing conditions
  IF (PDAY .and. PMONTH) THEN
    LTIMECOMPACT(JJ) = .TRUE.
  ELSE
    LTIMECOMPACT(JJ) = .FALSE.
  END IF
  !
  ! A.7. Conditions on a minimum threshold for grooming : 20kg/m2 i.e. 20cm of fresh snow (100kg/m3)
  IF (SUM(PSNOWSWE(JJ,1:INLVLS_USE(JJ))) < 20.) THEN
    LTIMECOMPACT(JJ) = .FALSE.
  END IF
  !
  IF (LTIMECOMPACT(JJ)) THEN
    DO JST=1,INLVLS_USE(JJ)
      ! B. Overburden weight by the grooming machine (static load)
      ! 500 on the first 50 kg/m2 from the surface
      ! Linear decrease from 500 to 0 between 50 and 150 kg/m2 (Fig. 4 Spandre et al. 2016)
      ! 0 deeper than 150 kg/m2
      ZSNOWCOMPACT(JJ,JST) = MIN(500.,MAX((150.-SUM(PSNOWSWE(JJ,1:JST)))*5., 0.0))
      !
      ! C. Recalculation of density of layers after overburden weight was applied
      ZSNOWRHOC = PSNOWRHO(JJ,JST) + PSNOWRHO(JJ,JST)*PSMASSCOEFF(JJ,JST)*ZSNOWCOMPACT(JJ,JST)
      ZSNOWDZC = PSNOWDZ(JJ,JST)*PSNOWRHO(JJ,JST)/ZSNOWRHOC
      !
      PSNOWRHO(JJ,JST) = ZSNOWRHOC
      PSNOWDZ(JJ,JST) = ZSNOWDZC
    END DO
  !
  END IF
  !
  !
  !------------------ TILLING OPTION -----------------------------
  !
  ! Tilling occurs between 8pm and 8.45pm i.e. 4 time steps
  IF (LTIMECOMPACT(JJ) .and. OSNOWTILLER .and. TPTIME%TIME < 75600. .and. TPTIME%TIME >= 72000.) THEN
    DO JST=1,INLVLS_USE(JJ)
      ! Depth of tiller penetration based on SWE
      IF (SUM(PSNOWSWE(JJ,1:JST)) > 35.0) THEN
        II = JST
        EXIT
      END IF
    END DO
    !
    ZSNOWAGEB = SUM(PSNOWAGE(JJ,1:II)*PSNOWSWE(JJ,1:II))/SUM(PSNOWSWE(JJ,1:II))
    ZSNOWDIAMOPTB = SUM(PSNOWDIAMOPT(JJ,1:II)*PSNOWSWE(JJ,1:II))/SUM(PSNOWSWE(JJ,1:II))
    ZSNOWSPHERIB = SUM(PSNOWSPHERI(JJ,1:II)*PSNOWSWE(JJ,1:II))/SUM(PSNOWSWE(JJ,1:II))
    !
    ! Optical Diameter
    ZSNOWDIAMOPTB = MAX(ZSNOWDIAMOPTB, (ZSNOWDIAMOPTB*4+0.00026)/5)
    ! Sphericity
    ZSNOWSPHERIB = MAX(ZSNOWSPHERIB, (ZSNOWSPHERIB*4+90./100)/5)
    !
    DO JST=1,II
      PSNOWAGE(JJ,JST) = ZSNOWAGEB
      PSNOWDIAMOPT(JJ,JST) = ZSNOWDIAMOPTB
      PSNOWSPHERI(JJ,JST) = ZSNOWSPHERIB
    END DO
    !
    ! DENSITY RE-CALCULATION
    ! Average density
    ZSNOWRHOB = SUM(PSNOWRHO(JJ,1:II)*PSNOWSWE(JJ,1:II))/SUM(PSNOWSWE(JJ,1:II))
    ! Recalculated average density
    ZSNOWRHOB = MAX(ZSNOWRHOB, (ZSNOWRHOB*4.+450.)/5)
    !
    DO JST=1,II
      ZSNOWRHOG(JJ,JST) = ZSNOWRHOB
      ZSNOWDZG(JJ,JST) = PSNOWDZ(JJ,JST)*PSNOWRHO(JJ,JST)/ZSNOWRHOG(JJ,JST)
      !
      PSNOWRHO(JJ,JST) = ZSNOWRHOG(JJ,JST)
      PSNOWDZ(JJ,JST)=ZSNOWDZG(JJ,JST)
    !
    END DO
  END IF
  !
  !------------------ END OF TILLING OPTION ----------------------
  !
END DO
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWGROOMING',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWGROOMING
!
END SUBROUTINE SNOWCRO
