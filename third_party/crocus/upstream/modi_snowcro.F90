!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!auto_modi:spll_snowcro.D
MODULE MODI_SNOWCRO
INTERFACE

      SUBROUTINE SNOWCRO(HSNOWRES, TPTIME, OMEB, HIMPLICIT_WIND,      &
                      PPEW_A_COEF, PPEW_B_COEF,                                 &
                      PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,       &
                      PSNOWSWE,PSNOWRHO,PSNOWHEAT,PSNOWALB,                     &
                      PSNOWDIAMOPT,PSNOWSPHERI,PSNOWHIST,PSNOWAGE, PSNOWIMPUR,  &
                      PTSTEP,PPS,PSR,PUNLOAD, PRR,PPSN3L,PRESA_SV,              &
                      PTA,PTG,PSW_RAD,PQA,PVMOD,PWIND_DRIFT, PLW_RAD, PRHOA,    &
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

USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
IMPLICIT NONE
! CHM diagnostic extension: surface ice melted per unit area and second.
REAL, DIMENSION(:), OPTIONAL, INTENT(OUT) :: PSFCMELT
REAL, INTENT(IN)                       :: PTSTEP
TYPE(DATE_TIME), INTENT(IN)            :: TPTIME      ! current date and time
 CHARACTER(LEN=*), INTENT(IN)          :: HSNOWRES
LOGICAL, INTENT(IN)                    :: OMEB       ! True = coupled to MEB. This means surface fluxes ae IMPOSED
 CHARACTER(LEN=*), INTENT(IN)          :: HIMPLICIT_WIND   ! wind implicitation option
REAL, DIMENSION(:), INTENT(IN)         :: PPS, PTA, PSW_RAD, PQA, PVMOD, PWIND_DRIFT, PLW_RAD, PSR, PRR, PUNLOAD
REAL, DIMENSION(:,:), INTENT(IN)       :: P_DIR_SW, P_SCA_SW ! direct and diffuse spectral irradiance (W/m2/um)
REAL, DIMENSION(:,:), INTENT(IN)       :: PIMPWET, PIMPDRY  !Dry and wet deposit coefficient from Forcing File(g/m²/s)
REAL, DIMENSION(:), INTENT(IN)         :: PTG, PSOILCOND, PD_G, PPSN3L
REAL, DIMENSION(:), INTENT(IN)         :: PZREF, PUREF, PEXNS, PEXNA, PDIRCOSZW,PSLOPEDIR, PRHOA, PZ0, PZ0EFF, &
                                       PALB, PZ0H, PPERMSNOWFRAC,PVFRIC_T
REAL, DIMENSION(:), INTENT(IN)         :: PPEW_A_COEF, PPEW_B_COEF,                   &
                                        PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF,      &
                                        PPEQ_B_COEF
REAL, DIMENSION(:), INTENT(INOUT)      :: PSNOWALB
REAL, DIMENSION(:,:), INTENT(INOUT)    :: PSNOWHEAT, PSNOWRHO, PSNOWSWE
REAL, DIMENSION(:,:), INTENT(INOUT)  :: PSNOWDIAMOPT, PSNOWSPHERI, PSNOWHIST
REAL, DIMENSION(:,:), INTENT(INOUT)    :: PSNOWAGE  ! Snow grain age
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PSNOWIMPUR  ! Snow impurity content (g) (LOCATION,LAYER,NIMPUR)) Impur type :1/BC 2/Dust
REAL, DIMENSION(:,:), INTENT(OUT)    :: PSNOWTEMP
REAL, DIMENSION(:,:), INTENT(OUT)      :: PSNOWLIQ, PSNOWDZ
REAL, DIMENSION(:), INTENT(OUT)        :: PTHRUFAL,  PEVAPCOR,  PGFLXCOR
REAL, DIMENSION(:,:), INTENT(OUT)      :: PSPEC_ALB, PDIFF_RATIO,PSPEC_TOT  !! spectral albedo, diffuse to total irradiance ratio and total incoming spectral irradiance (npoints,nbands)
REAL, DIMENSION(:), INTENT(OUT)        ::   PSNOWFLUX
REAL, DIMENSION(:), INTENT(INOUT)      :: PGRNDFLUX
REAL, DIMENSION(:), INTENT(INOUT)      :: PRNSNOW, PHSNOW, PGFLUXSNOW, PLES3L, PLEL3L, &
                                          PHPSNOW, PCDSNOW, PUSTAR, PEVAP
REAL, DIMENSION(:), INTENT(OUT)        :: PSNDRIFT
REAL, DIMENSION(:), INTENT(INOUT)      :: PSWNETSNOW, PLWNETSNOW, PSWNETSNOWS
REAL, DIMENSION(:), INTENT(INOUT)      :: PCHSNOW,PRI
REAL, DIMENSION(:), INTENT(OUT)        :: PEMISNOW, PSNOWHMASS
REAL, DIMENSION(:), INTENT(OUT)        ::  PQS
REAL, DIMENSION(:), INTENT(IN)         :: PZENITH ! solar zenith angle
REAL, DIMENSION(:), INTENT(IN)         :: PAZIM   ! azimuthal angle      (radian from North, clockwise)
REAL, DIMENSION(:), INTENT(IN)         :: PXLAT,PXLON ! LAT/LON after packing
REAL, DIMENSION(:,:), INTENT(INOUT)    :: PBLOWSNW !  Properties of deposited blowing snow (from Sytron or Meso-NH/Crocus)
 CHARACTER(4), INTENT(IN)              :: HSNOWDRIFT        ! Snowdrift scheme :
LOGICAL, INTENT(IN)                    :: OSNOWDRIFT_SUBLIM ! activate sublimation during drift
REAL, DIMENSION (:), INTENT(IN)        ::  PSNOWMAK        ! Snowmaking thickness (m)
LOGICAL, INTENT(IN)                    :: OSNOWCOMPACT_BOOL, OSNOWMAK_BOOL, OSNOWTILLER, &
                                          OSELF_PROD, OSNOWMAK_PROP
LOGICAL, INTENT(IN)                    :: OSNOW_ABS_ZENITH ! activate parametrization of solar absorption for polar regions
 CHARACTER(3), INTENT(IN)              :: HSNOWMETAMO, HSNOWRAD, HSNOWFALL, HSNOWCOND, HSNOWHOLD, HSNOWCOMP, HSNOWZREF
CHARACTER(4), INTENT(IN)               :: HSNOWFPAPPUS
LOGICAL, INTENT(IN)                    :: OATMORAD ! activate atmotartes scheme
REAL, DIMENSION(:), INTENT(IN)         :: PRESA_SV ! 
REAL, DIMENSION(:), INTENT(IN)         :: PHVEGPOL ! Mean polar vegetation height
REAL, DIMENSION(:), INTENT(IN)         :: PAGINGCOEF ! Snow aging coefficient
LOGICAL, INTENT(IN)                    :: OSNOWAGING_VAR
LOGICAL, INTENT(IN)                    :: OFOREST
END SUBROUTINE SNOWCRO
END INTERFACE
END MODULE MODI_SNOWCRO
