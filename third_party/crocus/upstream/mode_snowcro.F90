!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!SFX_LIC for details. version 1.
!     ###################
      MODULE MODE_SNOWCRO
!     ###################
!
!!****  *MODE_SNOWCRO * - contains Crocus routines that must also be called in MEB
!!
!!
!!
!!    AUTHOR
!!    ------
!!      Original        31/07/24
!!      M Lafaysse
!!      from previous code of snowcro.F90 by V. Vionnet and E. Brun
!!
!!    MODIFICATIONS
!!    -------------
!!
!----------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
!-------------------------------------------------------------------------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
CONTAINS
!
!####################################################################
!
SUBROUTINE SNOWNLFALL_UPGRID(PTSTEP,PSR,PTA,PVMOD,           &
                             PSNOW,PSNOWRHO,PSNOWDZ,PSNOWHEAT,PSNOWHMASS,                 &
                             PSNOWALB,PPERMSNOWFRAC,PSNOWDIAMOPT,PSNOWSPHERI, PSNOWHIST,  &
                             PSNOWAGE,OSNOWFALL,PSNOWDZN,PSNOWRHOF,PSNOWDZF,     &
                             PSNOWDIAMOPTF,PSNOWSPHERIF,PSNOWHISTF,PSNOWAGEF,             &
                             PWETCOEF, PSNOWIMPURF,                                       &
                             OMODIF_GRID,KNLVLS_USE,HSNOWDRIFT,HSNOWFPAPPUS,PZ0EFF,PUREF, &
                             PBLOWSNW, HSNOWFALL,                         &
                             PSNOWMAK, OSNOWMAK_BOOL, OSNOWMAK_PROP, KMAX_USE)
!
!!    PURPOSE
!!    -------
! Adds new snowfall and updates the vertical grid in order to keep an
! optimal discertisation
!
!!    AUTHOR
!!    ------
!!      E. Brun           * Meteo-France *
!!
!
!!
!!     MODIFICATIONS
!!    ------
!!
!!     2014-02-05 V. Vionnet: wind speed in the parameterization for new snow
!!                            density and characteristic of grains of new snow
!!                            are taken at a reference height
!!     2014-06-03 M. Lafaysse : threshold on PZ0EFF
!!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_CSTS,     ONLY : XLMTT, XTT, XCI
!
USE MODD_SNOW_METAMO, ONLY : XNDEN1, XNDEN2, XNDEN3, &
                             XNSPH1, XNSPH2, XNSPH3, XNSPH4, NVHIS2,&
                             XUEPSI, XVDIAM3, XVDIAM6, XUEPSI_SMP
USE MODD_PREP_SNOW, ONLY : NIMPUR
!
!
USE MODD_SNOW_PAR, ONLY : XRHOSMIN_ES, XSNOWDMIN, XANSMAX, XSNOWCRITD,             &
                          XDZMIN_TOP, XDZMIN_TOP_BIS, XDZMIN_BOT, XSPLIT_COEF,     &
                          XAGREG_COEF_1, XAGREG_COEF_2, XDZ1, XDZ2, XDZ3, XDZ3_BIS,&
                          XDZ4, XDZ5, XDZ_BASE, XDZ_INTERNAL, XSCALE_CM,           &
                          XDZMAX_INTERNAL, XDZMIN_TOP_EXTREM, XSNOWFALL_THRESHOLD, &
                          XRATIO_NEWLAYER, XDEPTH_THRESHOLD1, XDEPTH_THRESHOLD2,   &
                          XDEPTH_SURFACE, XDIFF_1, XDIFF_MAX, XSCALE_DIFF,         &
                          XSNOWFALL_A_SN, XSNOWFALL_B_SN, XSNOWFALL_C_SN,          &
                          XSNOWFALL_A_SN_P75, XSNOWFALL_B_SN_P75, XSNOWFALL_C_SN_P75,&
                          XSNOWFALL_A_SN_R21, XSNOWFALL_B_SN_R21, XSNOWFALL_C_SN_R21,&
                          XSNOWFALL_A_SN_L22, XSNOWFALL_B_SN_L22, XSNOWFALL_C_SN_L22,&
                          XSNOWFALL_A_SN_GW1, XSNOWFALL_B_SN_GW1, XSNOWFALL_C_SN_GW1,&
                          XSNOWFALL_A_SN_GW2, XSNOWFALL_B_SN_GW2, XSNOWFALL_C_SN_GW2,&                          
                          XRHOS_A76_1, XRHOS_A76_2, XRHOS_A76_3, XRHOS_S02_1,      &
                          XRHOS_S02_2, XRHOS_S02_3, XRHOS_S02_4, XRHOS_S02_5,      &
                          XRHOS_S02_6, XRHO_SNOWMAK, XRHOTHRESHOLD_ICE
!
USE MODD_BLOWSNW_SURF
!
USE MODE_SNOW3L
!
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                     :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)       :: PSR, PTA, PVMOD, PPERMSNOWFRAC
!
REAL, DIMENSION(:),INTENT(IN)       :: PZ0EFF,PUREF
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PSNOW, PSNOWALB
!
REAL, DIMENSION(:,:), INTENT(IN)     :: PSNOWRHO, PSNOWDZ, PSNOWHEAT
!
REAL, DIMENSION(:), INTENT(OUT)      :: PSNOWHMASS
!
REAL, DIMENSION(:,:), INTENT(IN)     :: PSNOWDIAMOPT, PSNOWSPHERI, PSNOWHIST, PSNOWAGE
!
REAL, DIMENSION(:,:), INTENT(INOUT)  :: PBLOWSNW
!
LOGICAL, DIMENSION(:), INTENT(OUT) :: OSNOWFALL
!
! Fresh snow characteristics
REAL, DIMENSION(:), INTENT(OUT)      :: PSNOWRHOF, PSNOWDZF
REAL, DIMENSION(:), INTENT(OUT)      :: PSNOWDIAMOPTF, PSNOWSPHERIF, PSNOWHISTF
REAL, DIMENSION(:), INTENT(OUT)      :: PSNOWAGEF
REAL, DIMENSION(:,:), INTENT(IN)       :: PWETCOEF
REAL, DIMENSION(:,:), INTENT(OUT)      :: PSNOWIMPURF
! New vertical grid
REAL, DIMENSION(:,:), INTENT(OUT)    :: PSNOWDZN
!
LOGICAL, DIMENSION(:), INTENT(OUT)   :: OMODIF_GRID
!
!
INTEGER, DIMENSION(:), INTENT(INOUT) :: KNLVLS_USE
!
CHARACTER(4), INTENT(IN)            :: HSNOWDRIFT        ! Snowdrift scheme :
                                      ! Mechanical transformation of snow grain and compaction + effect of wind
                                      ! on falling snow properties
                                      !    'NONE': No snowdrift scheme
                                      !    'DFLT': falling snow falls as purely dendritic
                                      !    'GA01': Gallee et al 2001
                                      !    'VI13': Vionnet et al 2013
                                      !    'PAPP': snowdrift scheme coupled with snowpappus + failing snow properties as 'NONE'
                                      !    'R21F': Royer et al 2021 (Increase in Maximum Density and Wind Effect)
                                      !    'R21W': Royer et al 2021 (Increase in Wind_Effect)
                                      !    'R21R': Royer et al 2021 (Increase in Maximum Density)                                      
CHARACTER(4), INTENT(IN)            :: HSNOWFPAPPUS
                                      ! Option to force fresh snow characteristics when snowpappus is used
                                      ! 'GM98' => forces fresh snow characteristics as if HSNOWDRIFT = 'NONE'
                                      ! 'VI13' => """" 'VI13'
                                      ! 'NONE' => no influence ( it is forced when snowpappus is not activated )
!
CHARACTER(3), INTENT(IN)              :: HSNOWFALL   ! snowfall density scheme Cluzet et al 2016
!
! Snowmaking option by p.spandre 20160211
REAL, DIMENSION (:), INTENT(IN)      :: PSNOWMAK
LOGICAL, INTENT(IN)                  :: OSNOWMAK_BOOL, OSNOWMAK_PROP ! if MM Snow production
INTEGER, INTENT(INOUT) :: KMAX_USE
!
!*      0.2    declarations of local variables
!
!
LOGICAL, DIMENSION(SIZE(PTA))       :: GAGREG_SURF
!
REAL, DIMENSION(SIZE(PTA))          :: ZSNOWFALL, ZSNOWTEMP, ZSCAP
! Snowmaking option by p.spandre 28/01/2014
REAL, DIMENSION(SIZE(PTA))          :: ZPSR_SNOWMAK
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZDZOPT
!
REAL :: ZZ0EFF
!
REAL :: ZSNOW_UPPER, ZSNOW_UPPER2 ! snow depth treatednormally (<= XDEPTH_SURFACE)
REAL :: ZCOEF_DEPTH !coefficient for repartition of deep snow above 3 meters
REAL :: ZTHICKNESS_INTERMEDIATE, ZTHICKNESS2
REAL :: ZPENALTY, ZDIFTYPE_INF, ZDIFTYPE_SUP, ZCRITSIZE_INF
REAL :: ZSNOW2L, ZCOEF
!
INTEGER :: INB_DEEP_LAYER, INB_UPPER_LAYER !separation between deep and upper layers
                                           ! if snow depth below XDEPTH_SURFACE then INB_DEEP_LAYER=0
INTEGER :: INB_MIN_LAYERS    ! why this test ?
INTEGER :: INB_INTERMEDIATE  ! number of intermediate layers (constant optimal gridding)
INTEGER :: IEND_INTERMEDIATE ! layer indice for bottom of intermediate layers
INTEGER :: JSTDEEP, JSTEND
INTEGER :: JST_1, JJ_A_AGREG_SUP, JJ_A_AGREG_INF
INTEGER :: INLVLSMIN, INLVLSMAX, JJ, JST, JIMP
!
! Coefficient to adjust wind speed at the height used in the parameterization
! for:
!         - density of new snow
!         - sphericity and dendricity of new snow
! Default values : 10 m for new snow (Pahaut, 1976) and 5 m for characteristics
! of snow grains (Guyomarc'h et Merindol, 1998)
REAL, PARAMETER                    :: PPHREF_WIND_RHO   = 10.
REAL, PARAMETER                    :: PPHREF_WIND_GRAIN = 5.
REAL, PARAMETER                    :: PPHREF_WIND_MIN = MIN(PPHREF_WIND_RHO,PPHREF_WIND_GRAIN)*0.5
REAL, DIMENSION(SIZE(PTA))         :: ZWIND_RHO
REAL, DIMENSION(SIZE(PTA))         :: ZWIND_GRAIN
!
REAL(KIND=JPRB)                  :: ZHOOK_HANDLE
!
!*      1.0   Initialization and snowage calculation for the present date
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWNLFALL_UPGRID',0,ZHOOK_HANDLE)
!
INLVLSMAX = SIZE (PSNOWRHO(:,:),2)
INLVLSMIN = 3
!
ZSNOWTEMP(:) = XTT
ZSNOWFALL(:) = 0.0 !Matthieu Lafaysse 21/09/2012
!
OSNOWFALL  (:) =.FALSE.
GAGREG_SURF(:) =.FALSE.
!
PSNOWHMASS (:) = 0.0
PSNOWRHOF  (:) = 0.0
PSNOWDZF   (:) = 0.0
PSNOWDIAMOPTF(:) = 0.0
PSNOWSPHERIF(:) = 0.0
PSNOWHISTF (:) = 0.0
DO JIMP=1,NIMPUR
  PSNOWIMPURF(:,JIMP) =0.0
ENDDO
PSNOWDZN (:,:) = PSNOWDZ(:,:)
ZPSR_SNOWMAK(:)= 0.0    !Pierre Spandre 28/01/2014
!
OMODIF_GRID(:) = .FALSE.
!
!************************************************************************************
!*      1.1   Calculation of the optimal vertical grid size ZDZOPT
!             as a function of maximum number of layers and of current
!             snow depth (modified 05/06/2012 by Matthieu Lafaysse)
!
! KNLVLS_USE(JJ) > INB_MIN_LAYERS =>
! KNLVLS_USE(JJ) > 2 + INLVLSMAX/3 =>
! ( KNLVLS_USE(JJ) + INLVLSMAX ) / 6 > (2 + INLVLSMAX/3 + INLVLSMAX) / 6 =>
! INB_DEEP_LAYER > (2 + 4*INLVLSMAX/3 ) / 6 >= 1
INB_MIN_LAYERS = 2 + INLVLSMAX/3
!
DO JJ = 1,SIZE(PSNOW(:))
  !
  IF ( PSNOW(JJ)>XDEPTH_THRESHOLD2 .AND. KNLVLS_USE(JJ)>INB_MIN_LAYERS ) THEN
    ! for very thick snowpack with enough snow layers
    ! special treatment
    ! we put the highest thickness in the lowest layers
    ! about 1/3 of layers for all snow except XDEPTH_SURFACE=3 first meters
    !
    !number of "deep layers"
    INB_DEEP_LAYER  = ( KNLVLS_USE(JJ) + INLVLSMAX ) / 6
    !
    !number of "upper layers"
    INB_UPPER_LAYER = KNLVLS_USE(JJ) - INB_DEEP_LAYER
    !
    !thickness of "upper layers"
    ZSNOW_UPPER = XDEPTH_SURFACE
    !
    !Arithmetic serie : 1+2+3+...+INB_DEEP_LAYER=INB_DEEP_LAYER*(INB_DEEP_LAYER+1)/2
    ZCOEF_DEPTH = ( PSNOW(JJ) - XDEPTH_SURFACE ) * 2. / ( (INB_DEEP_LAYER+1) * INB_DEEP_LAYER )
    !
    ! deep layers optimal thickness :
    ! increasing thickness with depth
    DO JSTDEEP = 1,INB_DEEP_LAYER
      JST = INB_UPPER_LAYER + JSTDEEP
      ZDZOPT(JJ,JST) = ZCOEF_DEPTH * JSTDEEP
      !This sum is equal to PSNOW(JJ)-XDEPTH_SURFACE
    ENDDO
    !
  ELSE
    !
    INB_UPPER_LAYER = KNLVLS_USE(JJ)
    !
    ZSNOW_UPPER = PSNOW(JJ)
    !
  END IF
  !
  !on force le ZDZOPT des 3 premières couches à ZSNOW_UPPER/3 maximum, chacune.
  ! => si on n'a qu'une couche, ZDZOPT(1) = ZSNOW_UPPER/3
  ! quel que soit INB_UPPER_LAYER
  !
  ZSNOW_UPPER2 = ZSNOW_UPPER / MAX( INLVLSMIN, INB_UPPER_LAYER )
  !
  ZDZOPT(JJ,1) = MIN( XDZ1, ZSNOW_UPPER2 )
  IF ( KNLVLS_USE(JJ)>=2 ) ZDZOPT(JJ,2) = MIN( XDZ2, ZSNOW_UPPER2 )
  IF ( KNLVLS_USE(JJ)>=3 ) ZDZOPT(JJ,3) = MIN( XDZ3, ZSNOW_UPPER2 )
  !
  IF ( INB_UPPER_LAYER>0 ) THEN
    !
    ZSNOW_UPPER2 = ZSNOW_UPPER / INB_UPPER_LAYER
    !
    ! dans ce cas, à partir de la 3ème couche, on prend la fraction du nombre de
    ! couches supérieures total, pour les couches jusqu'à 5
    !
    !ML : replace > by >= on 12-12-20 because the last layer was not initialised in case of thick snowpacks
    IF ( INB_UPPER_LAYER>=3 ) ZDZOPT(JJ,3) = MIN( XDZ3_BIS, ZSNOW_UPPER2 )
    IF ( INB_UPPER_LAYER>=4 ) ZDZOPT(JJ,4) = MIN( XDZ4    , ZSNOW_UPPER2 )
    IF ( INB_UPPER_LAYER>=5 ) ZDZOPT(JJ,5) = MIN( XDZ5    , ZSNOW_UPPER2 )
    !
    IF ( INB_UPPER_LAYER==KNLVLS_USE(JJ) ) THEN
      ! si on n'a pas de couches profondes
      !
      ! dans ce cas, on reprend ZSNOW_UPPER/3 maximum pour la dernière couche
      !
      ! last layer of' upper layers' : normal case : thin layer
      ZDZOPT(JJ,INB_UPPER_LAYER) = MIN( XDZ_BASE, ZSNOW_UPPER/MAX(INLVLSMIN,INB_UPPER_LAYER) )
      !
      ! ZTHICKNESS_INTERMEDIATE contient ce qu'il reste d'épaisseur disponible
      ! dans les couches supérieures
      !remaining snow for remaining layers
      ZTHICKNESS_INTERMEDIATE = ZSNOW_UPPER - SUM ( ZDZOPT ( JJ , 1: MIN ( 5, INB_UPPER_LAYER - 1 ))) &
                                            - ZDZOPT ( JJ,INB_UPPER_LAYER )

      IF ( ZSNOW_UPPER<=XDEPTH_THRESHOLD1 .OR. INB_UPPER_LAYER<8 ) THEN
        INB_INTERMEDIATE  = INB_UPPER_LAYER - 6
        IEND_INTERMEDIATE = INB_UPPER_LAYER - 1
      ELSE
        ! si INB_UPPER_LAYER>=8, les avant et avant-dernière couches ne sont pas
        ! considérées commes intermédiaires
        INB_INTERMEDIATE  = INB_UPPER_LAYER - 8
        IEND_INTERMEDIATE = INB_UPPER_LAYER - 3
        ! dans ce cas, on garde un peu d'épaisseur pour les deux couches restantes
        IF ( INB_INTERMEDIATE>0 ) THEN
          ZTHICKNESS_INTERMEDIATE = ZTHICKNESS_INTERMEDIATE * INB_INTERMEDIATE / FLOAT(INB_INTERMEDIATE+1)
        END IF
      END IF
      !
    ELSE
      ! si on a des couches profondes, les couches intermédiaires sont celles
      ! qui restent quand on a enlevé les 5 premières des couches supérieures
      !
      ! case with very thick snowpacks :
      ! the last layer of upper layers is not an exception
      ZTHICKNESS_INTERMEDIATE = ZSNOW_UPPER - SUM(ZDZOPT(JJ,1:5))
      INB_INTERMEDIATE  = INB_UPPER_LAYER - 5
      IEND_INTERMEDIATE = INB_UPPER_LAYER
      !
    END IF
    !
    ! For thick snowpack : add maximum value of optimal thickness to avoid too
    ! large differencies between layers
    IF ( INB_INTERMEDIATE>0 ) THEN
      !
      ZTHICKNESS2 = MAX( XDZ_INTERNAL, ZTHICKNESS_INTERMEDIATE/INB_INTERMEDIATE )
      !
      JSTEND = MIN( IEND_INTERMEDIATE,10 )
      DO JST = 6,JSTEND
        ZDZOPT(JJ,JST) = MIN( XDZMAX_INTERNAL(JST-5), ZTHICKNESS2 )
      END DO
      !
      IF ( IEND_INTERMEDIATE>10 ) THEN
        DO JST = 11,IEND_INTERMEDIATE
          ZDZOPT(JJ,JST) = ZTHICKNESS2
        END DO
      END IF
      !
    END IF
    !
    IF ( ZSNOW_UPPER>=XDEPTH_THRESHOLD1 .AND. INB_UPPER_LAYER>=8 ) THEN
      !Linear interpolation of optimal thickness between layers N-3 and N :
      ZDZOPT(JJ,INB_UPPER_LAYER-2) = 0.34*ZDZOPT(JJ,INB_UPPER_LAYER) + &
                                     0.66*ZDZOPT(JJ,INB_UPPER_LAYER-3)
      ZDZOPT(JJ,INB_UPPER_LAYER-1) = 0.66*ZDZOPT(JJ,INB_UPPER_LAYER) + &
                                     0.34*ZDZOPT(JJ,INB_UPPER_LAYER-3)
    ENDIF
    !
  END IF
  !
END DO
!
!************************************************************************************
!
!*      2.0   Fresh snow characteristics
!
!
!
! Heat content of newly fallen snow (J/m2):
! NOTE for now we assume the snowfall has
! the temperature of the snow surface upon reaching the snow.
! This is done as opposed to using the air temperature since
! this flux is quite small and has little to no impact
! on the time scales of interest. If we use the above assumption
! then, then the snowfall advective heat flux is zero.
!!
!
DO JJ = 1,SIZE(PSNOW(:))
!Snowmaking option 2014/01/28 : calculation of snowmaking rate, =0 if no snowmaking (PSNOWMAK(jj)=0), =XPSR_SNOWMAK otherwise
  IF (OSNOWMAK_BOOL) ZPSR_SNOWMAK(JJ) = PSNOWMAK(JJ)*XRHO_SNOWMAK/PTSTEP
!
!VV  IF (PSR(JJ)>XUEPSI .OR. PBLOWSNW(JJ,1) > XUEPSI  .OR. ZPSR_SNOWMAK(JJ) > XUEPSI ) THEN
   IF (PSR(JJ)>XUEPSI_SMP*200./PTSTEP .OR. PBLOWSNW(JJ,1) > XUEPSI  .OR. ZPSR_SNOWMAK(JJ) > XUEPSI) THEN          
    !
    ! newly fallen snow characteristics:
    !Case of new snowfall on a previously snow-free surface
    IF ( KNLVLS_USE(JJ)>0 ) THEN
      ZSCAP    (JJ) = XCI*PSNOWRHO(JJ,1)
      ZSNOWTEMP(JJ) = XTT + ( PSNOWHEAT(JJ,1) + XLMTT*PSNOWRHO(JJ,1)*PSNOWDZ(JJ,1) ) / &
                            ( ZSCAP(JJ) * MAX( XSNOWDMIN/INLVLSMAX, PSNOWDZ(JJ,1) ) )
    ELSE  ! case with bare ground
      ZSNOWTEMP(JJ) = PTA(JJ)
    ENDIF
    ZSNOWTEMP(JJ) = MIN( XTT, ZSNOWTEMP(JJ) )
    !
    ! ok dans ttes versions
    ! Debut modifs par VV
    !
    !
    ! Wind speeds at reference heights for new snow density and charactristics of
    ! grains of new snow
    ! Computed from PVMOD at PUREF (m) assuming a log profile in the SBL
    ! and a roughness length equal to PZ0EFF
    !PZ0EFF
    ZZ0EFF=MIN(PZ0EFF(JJ),PUREF(JJ)*0.5,PPHREF_WIND_MIN)

    ZWIND_RHO(JJ)   = PVMOD(JJ)*LOG(PPHREF_WIND_RHO/ZZ0EFF)/          &
                               LOG(PUREF(JJ)/ZZ0EFF)
    ZWIND_GRAIN(JJ) = PVMOD(JJ)*LOG(PPHREF_WIND_GRAIN/ZZ0EFF)/        &
                               LOG(PUREF(JJ)/ZZ0EFF)
    
    PSNOWHMASS(JJ) = (PSR(JJ)+PBLOWSNW(JJ,1)+ZPSR_SNOWMAK(JJ))*&
                    (XCI*(ZSNOWTEMP(JJ)-XTT)-XLMTT)*PTSTEP  !20160211
!VV    IF (PSR(JJ)>XUEPSI) THEN 
    IF (PSR(JJ)>XUEPSI_SMP*200./PTSTEP) THEN !VV Modification for single precision            
      !
      !! Cluzet et al 2016
      !! implementation of different parametrical options for fresh snow density.
      !! Be careful to the time-validity of the options A76(<2h) and S02(<1h) as well as the range of densities. Refer to Lehning et al. 2002 SNOWPACKIII, Anderson 76 and Pahaut 1975
      IF ( HSNOWFALL == 'V12' ) THEN ! Crocus original law
          PSNOWRHOF (JJ) = MAX( XRHOSMIN_ES, XSNOWFALL_A_SN + &
                                         XSNOWFALL_B_SN * ( PTA(JJ)-XTT ) + &
                                         XSNOWFALL_C_SN * SQRT(ZWIND_RHO(JJ) ) )
      ELSEIF( HSNOWFALL == 'P75') THEN ! Pahaut original law quoted by Brun 1989 but with different X_SNOWFALL_BSN
          PSNOWRHOF (JJ) = MAX( XRHOSMIN_ES, XSNOWFALL_A_SN_P75 + &
                                         XSNOWFALL_B_SN_P75 * ( PTA(JJ)-XTT ) + &
                                         XSNOWFALL_C_SN_P75 * SQRT(ZWIND_RHO(JJ) ) )
      ELSEIF( HSNOWFALL == 'R21') THEN ! Royer et al. 2021 (Doubled wind speed)
          PSNOWRHOF (JJ) = MAX( XRHOSMIN_ES, XSNOWFALL_A_SN_R21 + &
                                         XSNOWFALL_B_SN_R21 * ( PTA(JJ)-XTT ) + &
                                         XSNOWFALL_C_SN_R21 * SQRT(ZWIND_RHO(JJ) ) )
      ELSEIF( HSNOWFALL == 'L22') THEN ! Lackner et al. 2022 (Doubled density, Increased wind speed by 5)
          PSNOWRHOF (JJ) = MAX( XRHOSMIN_ES, XSNOWFALL_A_SN_L22 + &
                                         XSNOWFALL_B_SN_L22 * ( PTA(JJ)-XTT ) + &
                                         XSNOWFALL_C_SN_L22 * SQRT(ZWIND_RHO(JJ) ) )
      ELSEIF( HSNOWFALL == 'GW1') THEN ! GW1 (XSNOWFALL_C_SN * 1.5) 
          PSNOWRHOF (JJ) = MAX( XRHOSMIN_ES, XSNOWFALL_A_SN_GW1 + &
                                         XSNOWFALL_B_SN_GW1 * ( PTA(JJ)-XTT ) + &
                                         XSNOWFALL_C_SN_GW1 * SQRT(ZWIND_RHO(JJ) ) )
      ELSEIF( HSNOWFALL == 'GW2') THEN ! GW2 (XSNOWFALL_C_SN * 1) 
          PSNOWRHOF (JJ) = MAX( XRHOSMIN_ES, XSNOWFALL_A_SN_GW2 + &
                                         XSNOWFALL_B_SN_GW2 * ( PTA(JJ)-XTT ) + & 
                                         XSNOWFALL_C_SN_GW2 * SQRT(ZWIND_RHO(JJ) ) )
      ELSEIF ( HSNOWFALL == 'S02') THEN ! SNOWPACK 2014 law  min wind speed = 2m/s
          IF (PTA(JJ) > 259.15) THEN
              PSNOWRHOF (JJ)=EXP(( XRHOS_S02_1 + XRHOS_S02_2 * (PTA(JJ)-XTT) +&
                                   XRHOS_S02_3 + XRHOS_S02_4 * ASIN( SQRT(XRHOS_S02_5)) +&
                                   XRHOS_S02_6 * LOG10( MAX( SQRT(ZWIND_RHO(JJ) ), 2. )))&
                                   *LOG(10.))
          ELSE
              PSNOWRHOF (JJ)=EXP(( XRHOS_S02_1 + XRHOS_S02_2 * (PTA(JJ)-XTT) +&
                                   XRHOS_S02_4 * ASIN( SQRT(XRHOS_S02_5)) +&
                                   XRHOS_S02_6 * LOG10( MAX( SQRT(ZWIND_RHO(JJ) ), 2. )))&
                                   *LOG(10.))
          ENDIF
      ELSEIF ( HSNOWFALL == 'A76') THEN ! Anderson 76 law
        !VVV IF(PTA(JJ) - XTT + XRHOS_A76_3 < 0.) THEN
        IF(PTA(JJ) - XTT + XRHOS_A76_3 < XUEPSI_SMP) THEN !VV Modification for single precision                
          PSNOWRHOF (JJ) = XRHOS_A76_1
        ELSE
          PSNOWRHOF (JJ) = XRHOS_A76_1 + MAX( EXP(1.5*LOG(XRHOS_A76_2*( PTA(JJ) - XTT + XRHOS_A76_3 ))),0.  )
          !Nota Cluzet floating-point exception à cette ligne.... essai de réécriture de la puissance sous forme logarithmique. en fait, un log(<0) apparaît lorsqu'il neige pour TA<-15°C...
        ENDIF
      ELSEIF ( HSNOWFALL == 'NZE' ) THEN
          PSNOWRHOF (JJ) = 200.
      END IF
    ENDIF
    !
    !
    !  Density of accumulated snow (falling+blowing snow) : weighted average of
    !  PSNOWRHOF and density of accumulated snow
    !
    IF( PBLOWSNW(JJ,1) > XUEPSI) THEN
      PSNOWRHOF(JJ) = (PSNOWRHOF(JJ)*PSR(JJ) + PBLOWSNW(JJ,2) * PBLOWSNW(JJ,1))/ &
                       (PSR(JJ)+PBLOWSNW(JJ,1))
    ENDIF
!
!!modifs par VV
!
    IF (OSNOWMAK_PROP .and. ZPSR_SNOWMAK(JJ)>XUEPSI) THEN
      PSNOWRHOF(JJ) = ((PSR(JJ)+PBLOWSNW(JJ,1))*PSNOWRHOF(JJ)+ ZPSR_SNOWMAK(JJ)*XRHO_SNOWMAK)/ &  ! Additionnal boolean to use modified properties of machine made snow (MMS) or not p.spandre 2014/07/15
      (PSR(JJ)+PBLOWSNW(JJ,1)+ZPSR_SNOWMAK(JJ))     ! NB : ZPSR_SNOWMAK = XPSR_SNOWMAK si prod de neige. =0 sinon.
    ENDIF
    ZSNOWFALL(JJ) = (PSR(JJ)+PBLOWSNW(JJ,1)+ZPSR_SNOWMAK(JJ)) * PTSTEP / PSNOWRHOF(JJ)  ! snowfall thickness (m)
!
!End of Snowmaking option
!! 20160211
    PSNOW     (JJ) = PSNOW(JJ) + ZSNOWFALL(JJ)
    PSNOWDZF  (JJ) = ZSNOWFALL(JJ)
    !
    IF (  HSNOWDRIFT=='DFLT' ) THEN
      PSNOWDIAMOPTF(JJ) = XVDIAM6
      PSNOWSPHERIF(JJ) = 0.5
    ELSE IF ( HSNOWDRIFT=='GA01' ) THEN
      ! 2nd Option : deposited grains have a thresold wind speed equal to
      ! the current 5m wind speed (cf Gallee et al, 2001)
      ZCOEF = MIN(MAX(2.868*EXP(-0.085*PVMOD(JJ))-1.,0.),1.)
      PSNOWSPHERIF(JJ) = 1. - 0.495 * ZCOEF
      PSNOWDIAMOPTF(JJ) = XVDIAM6 * &
                      ( ZCOEF + ( 1.- ZCOEF ) * &
                                ( 3.*PSNOWSPHERIF(JJ) + 4.*(1.-PSNOWSPHERIF(JJ)) ) )
      !
    ELSE IF ( HSNOWDRIFT=='VI13' .OR. HSNOWFPAPPUS=='VI13' .OR. (HSNOWDRIFT== 'R21F') .OR. (HSNOWDRIFT=='R21W') .OR. (HSNOWDRIFT=='R21R') ) THEN
      ! 3rd Option : parameterization of Vionnet et al (2013) that allows
      ! simulatneous snow transport and snowfall for wind speed higher than 6 m/s
      !PSNOWSPHERIF(JJ) = MIN(MAX(0.14/4.*(ZWIND_GRAIN(JJ)-2.)+0.5,0.5),0.9)
      !ZCOEF =  MAX(MIN(-0.07*(ZWIND_GRAIN(JJ)-2.)+1.,1.),0.2)
      ! Expression (B.3) and (B.4), Vionnet et al (2013):
      PSNOWSPHERIF(JJ) = MIN(MAX(0.035*ZWIND_GRAIN(JJ)+0.43,0.5),0.9)
      ZCOEF =  MIN(MAX(1.14-0.07*ZWIND_GRAIN(JJ),0.2),1.)
      PSNOWDIAMOPTF(JJ) = XVDIAM6 * &
                      ( ZCOEF + ( 1.- ZCOEF ) * &
                                ( 3.*PSNOWSPHERIF(JJ) + 4.*(1.-PSNOWSPHERIF(JJ)) ) )
    ELSE IF ( HSNOWDRIFT=='NONE' .OR. HSNOWDRIFT=='PAPP' .OR. HSNOWFPAPPUS=='GM98' ) THEN
      PSNOWSPHERIF(JJ) = MIN( MAX( 0.0795*ZWIND_GRAIN(JJ)+0.38, 0.5), 0.9 )
      ZCOEF = MAX( MIN( -0.173*ZWIND_GRAIN(JJ)+1.29, 0.2 ), 1.)
      PSNOWDIAMOPTF(JJ) = XVDIAM6 * &
                      ( ZCOEF + ( 1.- ZCOEF ) * &
                                ( 3.*PSNOWSPHERIF(JJ) + 4.*(1.-PSNOWSPHERIF(JJ)) ) )
    END IF
    !
    ! Additionnal boolean to use modified properties of machine made snow or not p.spandre 2014/07/15
    ! Weighted mean of snowfall + blowing snow + snow_making snow
    IF (OSNOWMAK_PROP) THEN
      ! SPECIFICATION OPT DIAM
      PSNOWDIAMOPTF(JJ)=(PSNOWDIAMOPTF(JJ)*PSR(JJ) + PBLOWSNW(JJ,3)*PBLOWSNW(JJ,1) + XVDIAM3*ZPSR_SNOWMAK(JJ))/ &
      (PSR(JJ)+PBLOWSNW(JJ,1)+ZPSR_SNOWMAK(JJ))
      ! SPECIFICATION SPHERICITY
      PSNOWSPHERIF(JJ)=(PSNOWSPHERIF(JJ)*PSR(JJ) + PBLOWSNW(JJ,4)*PBLOWSNW(JJ,1) + 0.9*ZPSR_SNOWMAK(JJ))/ &
      (PSR(JJ)+PBLOWSNW(JJ,1)+ZPSR_SNOWMAK(JJ))
    ELSE IF ( PBLOWSNW(JJ,1) > XUEPSI) THEN
      PSNOWDIAMOPTF(JJ) = (PSNOWDIAMOPTF(JJ) * PSR(JJ) + PBLOWSNW(JJ,3) * PBLOWSNW(JJ,1)) / (PSR(JJ) + PBLOWSNW(JJ,1))
      PSNOWSPHERIF(JJ) = (PSNOWSPHERIF(JJ) * PSR(JJ) + PBLOWSNW(JJ,4) * PBLOWSNW(JJ,1)) / (PSR(JJ) + PBLOWSNW(JJ,1))
    END IF
    !
    PSNOWHISTF (JJ) = 0.0
    IF (NIMPUR>1 .AND. PSR(JJ)>XUEPSI) THEN
      DO JIMP=1,NIMPUR
        PSNOWIMPURF(JJ,JIMP)=PWETCOEF(JJ,JIMP)
      ENDDO
    ENDIF

    PSNOWAGEF  (JJ) = 0.0
    
    OSNOWFALL  (JJ) = .TRUE.
    OMODIF_GRID(JJ) = .TRUE.
    !
  ENDIF
  !
ENDDO
!
! VV WHERE( OSNOWFALL(:) .AND. ABS(PSNOW(:)-ZSNOWFALL(:))< XUEPSI )
! VV Modification for single precision
WHERE( OSNOWFALL(:) .AND. ABS(PSNOW(:)-ZSNOWFALL(:))< XUEPSI_SMP )
  PSNOWALB(:) = XANSMAX
END WHERE
!
! Computation of the new grid size
! It starts with successive exclusive cases
! Each case is described inside the corresponding condition
!
! cases with fresh snow
!
DO JJ=1,SIZE(PSNOW(:)) ! grid point loop
  !
  IF( .NOT.OSNOWFALL(JJ) .AND. PSNOW(JJ)>=XSNOWCRITD .AND. KNLVLS_USE(JJ)>=INLVLSMIN ) THEN
    !
    ! no fresh snow + deep enough snowpack + enough snow layers ==> no change
    !
  ELSEIF( PSNOW(JJ)<XSNOWCRITD .OR. KNLVLS_USE(JJ)<INLVLSMIN .OR. PSNOW(JJ)==ZSNOWFALL(JJ) ) THEN
    !
    ! too shallow snowpack or too few layers or only fresh snow
    ! ==> uniform grid and identical snow layers / number depends on snow depth
    OMODIF_GRID(JJ) = .TRUE.
    KNLVLS_USE (JJ) = MAX( INLVLSMIN, MIN( INLVLSMAX, INT(PSNOW(JJ)*XSCALE_CM) ) )
    KMAX_USE = MAX(KMAX_USE,KNLVLS_USE (JJ))
    PSNOWDZN(JJ,1:KNLVLS_USE(JJ)) = PSNOW(JJ) / KNLVLS_USE(JJ)
    PSNOWDZN(JJ,KNLVLS_USE(JJ) + 1:KMAX_USE) = 0.
    !
  ELSE
    !
    ! fresh snow over snow covered ground + enough snow layers
    OMODIF_GRID(JJ) = .TRUE.
    IF (PSNOWDZ(JJ,1)<   ZDZOPT(JJ,1)) THEN ! test for optimization because following is expensive
      ZDIFTYPE_SUP = SNOW3LDIFTYP( PSNOWDIAMOPT(JJ,1),PSNOWDIAMOPTF(JJ), &
                                 PSNOWSPHERI(JJ,1),PSNOWSPHERIF(JJ), &
                                 PSNOWHIST(JJ,1),PSNOWHISTF(JJ),&
                                 PSNOWRHO(JJ,1),PSNOWRHOF(JJ),&
                                 PSNOWAGE(JJ,1),PSNOWAGEF(JJ))
    ELSE
      ZDIFTYPE_SUP = 0 ! do not care about the value because will never be used
    ENDIF
    IF ( ( ZDIFTYPE_SUP<XDIFF_1        .AND. PSNOWDZ(JJ,1)<   ZDZOPT(JJ,1) ) .OR. &
         (  (PSR(JJ)+PBLOWSNW(JJ,1)) <XSNOWFALL_THRESHOLD .AND. PSNOWDZ(JJ,1)<2.*ZDZOPT(JJ,1) ) .OR. &
         ((PSNOWDZ(JJ,1)<XDZMIN_TOP_EXTREM) .AND. (PSNOWRHO(JJ,1)<XRHOTHRESHOLD_ICE))) THEN
      !
      ! Fresh snow is similar to a shallow surface layer (< ZDZOPT)
      ! or snowfall is very low and the surface layer not too deep (< 2*ZDZOPT) [NEW CONDITION 11/2012]
      ! or the surface layer is extremely thin (< XDZMIN_TOP_EXTREM) [NEW CONDITION 11/2012] and surface is not ice [NEW CONDITION 12/2017]
      ! The two new conditions are necessary for forcings with very low precipitation
      ! (e.g. ERA interim reanalyses, or climate models)
      ! ==> fresh snow is agregated to the surface layer
      !
      PSNOWDZN(JJ,1) = PSNOWDZ(JJ,1) + PSNOWDZF(JJ)
      DO JST = KMAX_USE,2,-1
        IF (JST <=KNLVLS_USE(JJ)) PSNOWDZN(JJ,JST) = PSNOWDZ(JJ,JST)
      ENDDO
      !
    ELSEIF ( KNLVLS_USE(JJ)<INLVLSMAX ) THEN
      !
      ! fresh snow is too different from the surface or the surface is too deep
      ! and there is room for extra layers ==> we create a new layer
      KNLVLS_USE(JJ)=KNLVLS_USE(JJ)+1
      KMAX_USE = MAX(KMAX_USE,KNLVLS_USE (JJ))
      !
      IF ( PSNOWDZF(JJ)>XRATIO_NEWLAYER*PSNOWDZ(JJ,2) .OR. &
           PSNOWRHO(JJ,1)>=XRHOTHRESHOLD_ICE .OR.  PSNOWRHOF(JJ)>=XRHOTHRESHOLD_ICE ) THEN
        !
        ! Snowfall is sufficient to create a new layer not lower than 1/10 of the second layer
        ! or snowfall directly over ice
        ! or freezing rain
        PSNOWDZN(JJ,1) = PSNOWDZF(JJ)
        DO JST =  KMAX_USE,2,-1
          IF (JST <=KNLVLS_USE(JJ)) PSNOWDZN(JJ, JST) = PSNOWDZ(JJ,JST-1)
        ENDDO
        !
      ELSE
        ! The ratio would be lower than 1/10 : [NEW : 11/2012]
        ! aggregate a part of the old layer with fresh snow to limit the ratio to 1/10.
        ZSNOW2L = PSNOWDZF(JJ) + PSNOWDZ(JJ,1)
        PSNOWDZN(JJ,1) = XRATIO_NEWLAYER      * ZSNOW2L
        PSNOWDZN(JJ,2) = (1.-XRATIO_NEWLAYER) * ZSNOW2L
        DO JST = KMAX_USE,3,-1
          IF (JST <=KNLVLS_USE(JJ)) PSNOWDZN(JJ,JST) = PSNOWDZ(JJ,JST-1)
        ENDDO
        !
      ENDIF
      !
    ELSE
      !
      ! fresh snow is too different from the surface or the surface is too deep
      ! and there is no room for extra layers
      ! ==> we agregate internal most similar snowlayers and create a new surface layer
      JJ_A_AGREG_SUP = 1
      JJ_A_AGREG_INF = 2
      !
      ZDIFTYPE_INF  = SNOW3LDIFTYP( PSNOWDIAMOPT(JJ,1),PSNOWDIAMOPTF(JJ), &
                                          PSNOWSPHERI(JJ,1),PSNOWSPHERIF(JJ), &
                                          PSNOWHIST(JJ,1),PSNOWHISTF(JJ),&
                                          PSNOWRHO(JJ,1), PSNOWRHOF(JJ),&
                                          PSNOWAGE(JJ,1),PSNOWAGEF(JJ))
            !
      ZCRITSIZE_INF = XSCALE_DIFF * ( PSNOWDZ(JJ,1)  /ZDZOPT(JJ,1)    + &
                                          PSNOWDZ(JJ,2)/ZDZOPT(JJ,2) )
      ZPENALTY = ZDIFTYPE_INF + ZCRITSIZE_INF
      DO JST = 1,KMAX_USE-1
        IF (JST <=KNLVLS_USE(JJ)-1) THEN
          !
          ZCRITSIZE_INF = XSCALE_DIFF * ( PSNOWDZ(JJ,JST)  /ZDZOPT(JJ,JST)    + &
                                          PSNOWDZ(JJ,JST+1)/ZDZOPT(JJ,JST+1) )
          !
          ZDIFTYPE_INF  = SNOW3LDIFTYP( PSNOWDIAMOPT(JJ,JST+1),PSNOWDIAMOPT(JJ,JST), &
                                        PSNOWSPHERI(JJ,JST+1),PSNOWSPHERI(JJ,JST), &
                                        PSNOWHIST(JJ,JST+1),PSNOWHIST(JJ,JST),&
                                        PSNOWRHO(JJ,JST+1), PSNOWRHO(JJ,JST),&
                                        PSNOWAGE(JJ,JST+1),PSNOWAGE(JJ,JST))
          !
          IF ( ZDIFTYPE_INF+ZCRITSIZE_INF<ZPENALTY ) THEN
            ZPENALTY = ZDIFTYPE_INF + ZCRITSIZE_INF
            JJ_A_AGREG_SUP = JST
            JJ_A_AGREG_INF = JST + 1
          ENDIF
          !
        ENDIF
      ENDDO
      ! agregation of the similar layers and shift of upper layers
      PSNOWDZN(JJ,JJ_A_AGREG_INF) = PSNOWDZ(JJ,JJ_A_AGREG_INF) + PSNOWDZ(JJ,JJ_A_AGREG_SUP)
      DO JST = JJ_A_AGREG_SUP,2,-1
        PSNOWDZN(JJ,JST) = PSNOWDZ(JJ,JST-1)
      ENDDO
      PSNOWDZN(JJ,1) = PSNOWDZF(JJ)
      !
      ! Limit the ratio between the new layer and the one beneath (ratio 1/10)
      ! [NEW : 11/2012]
      IF( PSNOWDZN(JJ,1)<XRATIO_NEWLAYER*PSNOWDZN(JJ,2) ) THEN
        ZSNOW2L = PSNOWDZN(JJ,1) + PSNOWDZN(JJ,2)
        PSNOWDZN(JJ,1) = XRATIO_NEWLAYER      * ZSNOW2L
        PSNOWDZN(JJ,2) = (1.-XRATIO_NEWLAYER) * ZSNOW2L
      ENDIF
      !
    ENDIF
    !
  ENDIF ! end of the case with fresh snow
  !
ENDDO ! end loop grid points
!
! cases with no fresh snow and no previous grid resize
!
IF ( INLVLSMIN==INLVLSMAX ) THEN ! specific case with INLSVSMIN = INLVLSMAX  (INLVLS)
  !
  ! check if surface layer depth is too small
  ! in such a case looks for an other layer to be split
  DO JJ = 1,SIZE(PSNOW(:)) ! loop grid points
    !
    IF (OMODIF_GRID(JJ)) CYCLE
    IF ( PSNOWDZ(JJ,1)<XDZMIN_TOP ) THEN
      CALL GET_SNOWDZN_DEB(INLVLSMAX,PSNOWDZ(JJ,:),ZDZOPT(JJ,:),PSNOWDZN(JJ,:))
      GAGREG_SURF(JJ) = .TRUE.
      OMODIF_GRID(JJ) = .TRUE.
    !ENDIF
    ! check if bottom layer depth is too small
    ! in such a case agregation with upper layer and
    ! looks for an other layer to be splitted
    ELSEIF(PSNOWDZ(JJ,INLVLSMAX)<XDZMIN_TOP ) THEN
      OMODIF_GRID(JJ) = .TRUE.
      CALL GET_SNOWDZN_END(INLVLSMAX,PSNOWDZ(JJ,:),ZDZOPT(JJ,:),PSNOWDZN(JJ,:))
    ENDIF
    !
  ENDDO ! end grid points loop
  !
ENDIF  ! end specific case INLSVSMIN = INLVLSMAX
!
! case without new snowfall and INVLSMAX > INLVLSMIN
!
DO JJ=1,SIZE(PSNOW(:))
  !
  IF (OMODIF_GRID(JJ)) CYCLE
  IF( .NOT.OSNOWFALL(JJ) .AND. PSNOW(JJ)>XSNOWCRITD ) THEN
    ! check if surface layer depth is too small
    ! in such a case agregation with layer beneath unless it is a single snow layer on an ice layer
    ! in case of reaching INLVLSMIN, looks for an other layer to be splitted
    IF(PSNOWDZ(JJ,1)<XDZMIN_TOP_BIS ) THEN ! case shallow surface layer
      !
      OMODIF_GRID(JJ) = .TRUE.
      !
      IF( KNLVLS_USE(JJ)>INLVLSMIN ) THEN ! case minimum not reached
        IF ((.NOT.((PSNOWRHO(JJ,1)<XRHOTHRESHOLD_ICE).AND.(PSNOWRHO(JJ,2)>=XRHOTHRESHOLD_ICE))) .AND. &
           (.NOT.((PSNOWRHO(JJ,1)>=XRHOTHRESHOLD_ICE).AND.(PSNOWRHO(JJ,2)<XRHOTHRESHOLD_ICE)))) THEN
          ! if it is not a single snow layer on an ice layer and not a single ice layer on snow layer
          KNLVLS_USE(JJ) = KNLVLS_USE(JJ) - 1
          PSNOWDZN(JJ,1) = PSNOWDZ(JJ,1) + PSNOWDZ(JJ,2)
          DO JST = 2, KMAX_USE
            IF (JST <=KNLVLS_USE(JJ)) &
            PSNOWDZN(JJ,JST) = PSNOWDZ(JJ,JST+1)
          ENDDO
          PSNOWDZN(JJ,KNLVLS_USE(JJ)+1) = 0.
        ENDIF
      ELSE ! case minimum reached
        CALL GET_SNOWDZN_DEB(KNLVLS_USE(JJ),PSNOWDZ(JJ,:),ZDZOPT(JJ,:),PSNOWDZN(JJ,:))
      ENDIF ! end case minimum reached end case shallow surface layer
      !
      GAGREG_SURF(JJ) = .TRUE.
    !
    ! check if bottom layer depth is too small
    ! in such a case agregation with above layer
    ! in case of reaching INLVLSMIN, looks for an other layer to be splitted
    ! case shallow bottom layer
    ELSEIF(PSNOWDZ(JJ,KNLVLS_USE(JJ))<XDZMIN_TOP .AND. .NOT.GAGREG_SURF(JJ) ) THEN
      !
      OMODIF_GRID(JJ) = .TRUE.
      !
      IF ( KNLVLS_USE(JJ)>INLVLSMIN ) THEN ! case minimum not reached
        KNLVLS_USE(JJ) = KNLVLS_USE(JJ) - 1
        PSNOWDZN(JJ,KNLVLS_USE(JJ)) = PSNOWDZ(JJ,KNLVLS_USE(JJ)) + PSNOWDZ(JJ,KNLVLS_USE(JJ)+1)
        PSNOWDZN(JJ,KNLVLS_USE(JJ)+1) = 0.
      ELSE  ! case minimum reached
        CALL GET_SNOWDZN_END(KNLVLS_USE(JJ),PSNOWDZ(JJ,:),ZDZOPT(JJ,:),PSNOWDZN(JJ,:))
      ENDIF ! end case minimum reached end case shallow surface layer
      !
    ENDIF
    !
  ENDIF
ENDDO ! end grid points loop
!
! case whithout new snow fall and without a previous grid resize
! looks for a thick layer to be splitted according to its depth and to
! the optimal grid size
DO JJ = 1,SIZE(PSNOW(:))
  !
  IF (OMODIF_GRID(JJ)) CYCLE
  IF (KNLVLS_USE(JJ)<INLVLSMAX-3 )THEN
    !
    DO JST = 1,INLVLSMAX-4
      !
      IF ( JST<=KNLVLS_USE(JJ) ) THEN !.AND. .NOT.OMODIF_GRID(JJ) ) THEN
        !
        IF( PSNOWDZ(JJ,JST) > &
            ( XSPLIT_COEF - FLOAT( INLVLSMAX-KNLVLS_USE(JJ) )/MAX( 1., FLOAT( INLVLSMAX-INLVLSMIN ) ) ) &
              * ZDZOPT(JJ,JST) ) THEN
          !
          DO JST_1 = KMAX_USE+1,JST+2,-1
            IF ( JST_1<=KNLVLS_USE(JJ)+1) THEN
              PSNOWDZN(JJ,JST_1) = PSNOWDZ(JJ,JST_1-1)
              ZDZOPT  (JJ,JST_1) = ZDZOPT (JJ,JST_1-1)
            ENDIF
          ENDDO
          !
          ! generale case : old layer divided in two equal layers
          IF ( JST/=1 .OR. PSNOWDZ(JJ,JST)<3.*ZDZOPT(JJ,1) ) THEN
            PSNOWDZN(JJ,JST+1) = 0.5*PSNOWDZ(JJ,JST)
            PSNOWDZN(JJ,JST)   = PSNOWDZN(JJ,JST+1)
          ELSE
            ! if thick surface layer : force the surface layer to this value to avoid successive resizing
            ! [NEW : 11/2012]
            PSNOWDZN(JJ,1) = 1.5 * ZDZOPT(JJ,1)
            PSNOWDZN(JJ,2) = PSNOWDZ(JJ,JST) - PSNOWDZN(JJ,1)
          ENDIF
          !
          KNLVLS_USE (JJ) = KNLVLS_USE(JJ) + 1
          KMAX_USE = MAX(KMAX_USE,KNLVLS_USE (JJ))
          OMODIF_GRID(JJ) = .TRUE.
          EXIT
          !
        ENDIF
        !
      ENDIF
      !
    ENDDO
    !
  ENDIF
  !
ENDDO
!
! case whithout new snow fall and without a previous grid resize
! looks for a deep layer to be agregated to the layer beneath if similar
! according to its depth and to the optimal grid size
!
!NB : allow these changes for 5 layers and more [NEW] (before : 6 layers)
!
DO JJ = 1,SIZE(PSNOW(:))
  !
  IF (OMODIF_GRID(JJ) )  CYCLE
    !
    IF (KNLVLS_USE(JJ)>INLVLSMIN+1) THEN
      DO JST = 2,KMAX_USE
        !
        IF ( JST<=KNLVLS_USE(JJ)-1 ) THEN !.AND. .NOT.OMODIF_GRID(JJ) ) THEN
          !
          !IF ((PSNOWRHO(JJ,JST+1)>=XRHOTHRESHOLD_ICE).AND.(PSNOWRHO(JJ,JST)<XRHOTHRESHOLD_ICE)) CYCLE ! never mix snow and ice
          !IF ((PSNOWRHO(JJ,JST+1)< XRHOTHRESHOLD_ICE).AND.(PSNOWRHO(JJ,JST)>= XRHOTHRESHOLD_ICE))CYCLE
          !
          IF (PSNOWDZ(JJ,JST) + PSNOWDZ(JJ,JST+1) < &
                    XAGREG_COEF_2 * MAX( ZDZOPT(JJ,JST),ZDZOPT(JJ,JST+1) ) )THEN
            ZDIFTYPE_INF = SNOW3LDIFTYP( PSNOWDIAMOPT(JJ,JST+1),PSNOWDIAMOPT(JJ, JST), &
                                         PSNOWSPHERI(JJ,JST+1),PSNOWSPHERI(JJ, JST), &
                                         PSNOWHIST(JJ,JST+1),PSNOWHIST(JJ, JST),&
                                         PSNOWRHO(JJ,JST+1), PSNOWRHO(JJ, JST),&
                                         PSNOWAGE(JJ,JST+1),PSNOWAGE(JJ,JST))
            ZDIFTYPE_INF = MAX( XDIFF_1, MIN( XDIFF_MAX, ZDIFTYPE_INF ) )
            !
            IF( PSNOWDZ(JJ,JST) < ZDZOPT(JJ,JST) * XAGREG_COEF_1 / ZDIFTYPE_INF) THEN
              !
              PSNOWDZN(JJ,JST) = PSNOWDZ(JJ,JST) + PSNOWDZ(JJ,JST+1)
              ZDZOPT  (JJ,JST) = ZDZOPT(JJ,JST+1)
              DO JST_1 = JST+1,KMAX_USE-1
                IF (JST_1<= KNLVLS_USE(JJ)-1) THEN
                  PSNOWDZN(JJ,JST_1) = PSNOWDZ(JJ,JST_1+1)
                  ZDZOPT  (JJ,JST_1) = ZDZOPT (JJ,JST_1+1)
                ENDIF
              ENDDO
              KNLVLS_USE(JJ) = KNLVLS_USE(JJ)-1
              PSNOWDZN(JJ, KNLVLS_USE(JJ) + 1) = 0.
              OMODIF_GRID(JJ)=.TRUE.
              EXIT
              !
            ENDIF
            !
          ENDIF
          !
        ENDIF
        !
      ENDDO
      !
    ENDIF
  !
ENDDO
!
! [NEW : 11/2012]
! In case of very low snow fall checks if a new internal snow layer is too shallow
! even if a the grid has already been resized in this time step
! starts from bottom to INLVS_USE-3 until old and new grid differ
DO JJ = 1,SIZE(PSNOW(:))
  !
  IF ( .NOT.OSNOWFALL(JJ) .OR. KNLVLS_USE(JJ)<INLVLSMIN+3 ) CYCLE ! go to next point
  !
  !VV IF( ABS( PSNOWDZN(JJ,KNLVLS_USE(JJ)) - PSNOWDZ(JJ,KNLVLS_USE(JJ)) ) > XUEPSI ) CYCLE ! go to next point
  IF( ABS( PSNOWDZN(JJ,KNLVLS_USE(JJ)) - PSNOWDZ(JJ,KNLVLS_USE(JJ)) ) > XUEPSI_SMP ) CYCLE ! go to next point  
  !
  ! bottom layer
  IF( PSNOWDZN(JJ,KNLVLS_USE(JJ))<XDZMIN_TOP ) THEN ! case shallow bottom layer
    !
    KNLVLS_USE(JJ) = KNLVLS_USE(JJ)-1
    PSNOWDZN(JJ,KNLVLS_USE(JJ)) = PSNOWDZN(JJ,KNLVLS_USE(JJ)) + PSNOWDZN(JJ,KNLVLS_USE(JJ)+1)
    PSNOWDZN(JJ,KNLVLS_USE(JJ)+1) = 0.
    !
  ELSE
    !
    ! internal layer
    DO JST =KMAX_USE-1,4,-1
      IF (JST <=KNLVLS_USE(JJ)-1) THEN
        !
!VV        IF ( ABS( PSNOWDZN(JJ,JST) - PSNOWDZ(JJ,JST) ) > XUEPSI ) EXIT ! old/new grid differ ==> go to next grid point
        IF ( ABS( PSNOWDZN(JJ,JST) - PSNOWDZ(JJ,JST) ) > XUEPSI_SMP ) EXIT ! old/new grid differ ==> go to next grid point        
        !
        IF ( PSNOWDZN(JJ,JST)> 0.001 ) CYCLE
        !
        IF (PSNOWRHO(JJ,JST)>=XRHOTHRESHOLD_ICE .AND. (PSNOWRHO(JJ,JST-1)<XRHOTHRESHOLD_ICE)) CYCLE ! never mix snow and ice
        !
        ! If an internal layer is too shallow, it is merged with the upper layer
        PSNOWDZN(JJ,JST-1) = PSNOWDZN(JJ,JST) + PSNOWDZN(JJ,JST-1)
        KNLVLS_USE(JJ)   = KNLVLS_USE(JJ) - 1
        !
        ! shifts the lower layers
        DO JST_1 = JST,KMAX_USE!KNLVLS_USE(JJ)
          IF (JST_1 <=KNLVLS_USE(JJ)) THEN
            PSNOWDZN(JJ,JST_1) = PSNOWDZN(JJ,JST_1+1)
            ZDZOPT  (JJ,JST_1) = ZDZOPT (JJ,JST_1+1)
          ENDIF
        ENDDO
        PSNOWDZN(JJ,KNLVLS_USE(JJ)+1) = 0.
        !
        EXIT ! goto to next grid point
        !
      ENDIF
    ENDDO ! end loop internal layers
    !
  ENDIF
  !
ENDDO ! end grid loops for checking shallow layers
!
!final check of the consistensy of the new grid size
!
! Following could be removed to reduce computing time if no debugging mode
DO JJ = 1,SIZE(PSNOW(:))
  !
!VV  IF ( ABS( SUM( PSNOWDZN(JJ,1:KNLVLS_USE(JJ)) ) - PSNOW(JJ) ) > XUEPSI ) THEN
  IF ( ABS( SUM( PSNOWDZN(JJ,1:KNLVLS_USE(JJ)) ) - PSNOW(JJ) ) > XUEPSI_SMP * MAX(1.0, PSNOW(JJ))) THEN          
    !
    WRITE(*,*) 'error in grid resizing', JJ, KNLVLS_USE(JJ), SUM( PSNOWDZN(JJ,1:KNLVLS_USE(JJ)) ),  &
                                         PSNOW(JJ), SUM( PSNOWDZN(JJ,1:KNLVLS_USE(JJ)) )-PSNOW(JJ), &
                                         ZSNOWFALL(JJ)
    WRITE( *,*) 'JJ , PSNOWDZ(JJ):',JJ ,  PSNOWDZ(JJ,:)
    WRITE( *,*) 'JJ , PSNOWDZN(JJ):',JJ , PSNOWDZN(JJ,:)
    !
    CALL ABOR1_SFX("SNOWCRO: error in grid resizing")
    !
  ENDIF
  !
ENDDO
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWNLFALL_UPGRID',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWNLFALL_UPGRID

!###############################################################################
SUBROUTINE GET_SNOWDZN_DEB(KNLVLS,PSNOWDZ,PDZOPT,PSNOWDZN)
!
USE MODD_SNOW_PAR, ONLY : XDZMIN_TOP, XDZMIN_BOT
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KNLVLS
REAL, DIMENSION(:), INTENT(IN)  :: PSNOWDZ, PDZOPT
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWDZN
!
REAL :: ZPENALTY, ZCRITSIZE
INTEGER :: JJ_A_DEDOUB, JST
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:GET_SNOWDZN_DEB',0,ZHOOK_HANDLE)
!
ZPENALTY = PSNOWDZ(2) / PDZOPT(2)
IF( PSNOWDZ(2)<XDZMIN_TOP ) ZPENALTY = 0.
JJ_A_DEDOUB = 2
!
DO JST = 3,KNLVLS
  ZCRITSIZE = PSNOWDZ(JST) / PDZOPT(JST)
  IF ( JST==KNLVLS .AND. PSNOWDZ(JST)<XDZMIN_BOT ) ZCRITSIZE = 0.
  IF ( ZCRITSIZE>ZPENALTY ) THEN
    ZPENALTY    = ZCRITSIZE
    JJ_A_DEDOUB = JST
  ENDIF
ENDDO
!
IF ( JJ_A_DEDOUB==2 ) THEN ! case splitted layer == 2
  PSNOWDZN(1) = 0.5 * ( PSNOWDZ(1) + PSNOWDZ(2) )
  PSNOWDZN(2) = PSNOWDZN(1)
ELSE ! case splitted layer =/ 2
  PSNOWDZN(1) = PSNOWDZ(1) + PSNOWDZ(2)
  DO JST = 2,JJ_A_DEDOUB-2
    PSNOWDZN(JST) = PSNOWDZ(JST+1)
  ENDDO
  PSNOWDZN(JJ_A_DEDOUB-1) = 0.5 * PSNOWDZ(JJ_A_DEDOUB)
  PSNOWDZN(JJ_A_DEDOUB)   = PSNOWDZN(JJ_A_DEDOUB-1)
ENDIF ! end case splitted layer =/ 2
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:GET_SNOWDZN_DEB',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_SNOWDZN_DEB
!
!###############################################################################
SUBROUTINE GET_SNOWDZN_END(KNLVLS,PSNOWDZ,PDZOPT,PSNOWDZN)
!
USE MODD_SNOW_PAR, ONLY : XDZMIN_TOP, XDZMIN_BOT
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KNLVLS
REAL, DIMENSION(:), INTENT(IN)  :: PSNOWDZ, PDZOPT
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWDZN
!
REAL :: ZPENALTY, ZCRITSIZE
INTEGER :: JJ_A_DEDOUB, JST
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:GET_SNOWDZN_END',0,ZHOOK_HANDLE)
!
ZPENALTY = PSNOWDZ(KNLVLS-2) / PDZOPT(KNLVLS-2)
JJ_A_DEDOUB = KNLVLS - 2
!
DO JST = MAX(1,KNLVLS-3),1,-1
  ZCRITSIZE = PSNOWDZ(JST) / PDZOPT(JST)
  IF ( JST==1 .AND. PSNOWDZ(JST)<XDZMIN_BOT ) ZCRITSIZE = 0.
  IF ( ZCRITSIZE>ZPENALTY ) THEN
    ZPENALTY    = ZCRITSIZE
    JJ_A_DEDOUB = JST
  ENDIF
ENDDO
!
IF ( JJ_A_DEDOUB==KNLVLS-1 ) THEN ! case splitted layer == 2
  PSNOWDZN(KNLVLS)   = 0.5 * (PSNOWDZ(KNLVLS-1)+PSNOWDZ(KNLVLS))
  PSNOWDZN(KNLVLS-1) = PSNOWDZN(KNLVLS)
ELSE ! case splitted layer =/ 2
  PSNOWDZN(KNLVLS) = PSNOWDZ(KNLVLS-1) + PSNOWDZ(KNLVLS)
  DO JST = KNLVLS-1,JJ_A_DEDOUB+2,-1
    PSNOWDZN(JST) = PSNOWDZ(JST-1)
  ENDDO
  PSNOWDZN(JJ_A_DEDOUB+1) = 0.5 * PSNOWDZ(JJ_A_DEDOUB)
  PSNOWDZN(JJ_A_DEDOUB  ) = PSNOWDZN(JJ_A_DEDOUB+1)
ENDIF ! end case splitted layer =/ 2
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:GET_SNOWDZN_END',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_SNOWDZN_END
!
!###############################################################################
!################################################################################
!################################################################################
!
SUBROUTINE SNOWNLGRIDFRESH_1D (KJ,PSNOW,PSNOWDZ,PSNOWDZN,                  &
                               PSNOWRHO,PSNOWHEAT,PSNOWDIAMOPT,PSNOWSPHERI,   &
                               PSNOWHIST,PSNOWAGE,PSNOWIMPUR, OSNOWFALL,   &
                               PSNOWRHOF, PSNOWDZF,PSNOWHEATF,PSNOWDIAMOPTF, &
                               PSNOWSPHERIF, PSNOWHISTF,PSNOWAGEF, PSNOWIMPURF, &
                               KNLVLS_USE, KNLVLS_USE_OLD, OSUCCESS)
!
!!    PURPOSE
!!    -------
!     Snow mass,heat and characteristics redistibution in case of
!     grid resizing. Total mass and heat content of the overall snowpack
!     unchanged/conserved within this routine.
!     Grain size and type of mixed layers is deduced from the conservation
!     of the average optical size
!
!!    AUTHOR
!!    ------
!!      E. Brun           * Meteo-France *
!!
!
USE MODD_SNOW_PAR, ONLY : XD1,XD2,XD3,XX,XVALB5,XVALB6
USE MODD_PREP_SNOW, ONLY : NIMPUR
USE MODE_SNOW3L, ONLY : GET_MASS_HEAT
USE MODD_SNOW_METAMO, ONLY : XUEPSI
!
IMPLICIT NONE
!
!
!*      0.1    declarations of arguments
!
INTEGER, INTENT(IN)               :: KJ
REAL, INTENT(IN)                  :: PSNOW
!
REAL, DIMENSION(:), INTENT(INOUT) :: PSNOWHEAT, PSNOWRHO, PSNOWDZ,     &
                                     PSNOWDZN, PSNOWDIAMOPT, PSNOWSPHERI, &
                                     PSNOWHIST
REAL, DIMENSION(:), INTENT(INOUT) :: PSNOWAGE
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWIMPUR

REAL,  INTENT(IN)                 :: PSNOWRHOF, PSNOWDZF,PSNOWHEATF,   &
                                     PSNOWDIAMOPTF,PSNOWSPHERIF, PSNOWHISTF
!
REAL, INTENT(IN)                  :: PSNOWAGEF
REAL,DIMENSION(:),INTENT(IN)         ::PSNOWIMPURF
!
INTEGER, INTENT(IN)               :: KNLVLS_USE, KNLVLS_USE_OLD ! new and old number of active snow layers
!
LOGICAL, INTENT(IN)               :: OSNOWFALL
!
LOGICAL, INTENT(OUT)              :: OSUCCESS
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO,1)+1) :: ZSNOWRHOO,ZSNOWDIAMOPTO,ZSNOWSPHERIO, &
                                       ZSNOWHEATO,ZSNOWHISTO,ZSNOWDZO,    &
                                       ZSNOWZTOP_OLD,ZSNOWZBOT_OLD
REAL,DIMENSION(SIZE(PSNOWRHO,1)+1)  :: ZSNOWAGEO
!
REAL, DIMENSION(SIZE(PSNOWRHO,1)) :: ZSNOWRHON,ZSNOWDIAMOPTN,ZSNOWSPHERIN,   &
                                     ZSNOWHEATN,ZSNOWHISTN,               &
                                     ZSNOWZTOP_NEW,ZSNOWZBOT_NEW
REAL,DIMENSION(SIZE(PSNOWRHO,1)) ::ZSNOWAGEN
!
REAL,DIMENSION(SIZE(PSNOWRHO,1)+1,NIMPUR) :: ZSNOWIMPURO

REAL,DIMENSION(SIZE(PSNOWRHO,1),NIMPUR)   :: ZSNOWIMPURN
!
REAL :: ZMASTOTN, ZMASTOTO, ZSNOWHEAN, ZSNOWHEAO
REAL :: ZPSNOW_OLD, ZPSNOW_NEW
!
INTEGER :: INLVLS_OLD, INLVLS_NEW
INTEGER :: JST, JIMP
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWNLGRIDFRESH_1D',0,ZHOOK_HANDLE)
!
! 0. Initialization:
! ------------------
!
OSUCCESS = .TRUE.
!
! starts by checking the consistency between both vertical grid sizes
INLVLS_NEW = KNLVLS_USE
INLVLS_OLD = -1
!
ZPSNOW_NEW = 0.
ZPSNOW_OLD = 0.
!
! Compute the old number of layers
!
! New total snowdepth
DO JST = 1,INLVLS_NEW
  ZPSNOW_NEW = ZPSNOW_NEW + PSNOWDZN(JST)
ENDDO
!
INLVLS_OLD = KNLVLS_USE_OLD
!
IF ( OSNOWFALL ) INLVLS_OLD = INLVLS_OLD + 1
!
ZPSNOW_OLD = PSNOW
ZPSNOW_NEW = ZPSNOW_OLD
!
! initialization of variables describing the initial snowpack + new snowfall
!
IF ( OSNOWFALL ) THEN
  ! Layers 2:JST of the newsnowpack take properties of layers 1:JST-1 of the old snowpack
  DO JST = 2,INLVLS_OLD
    ZSNOWDZO     (JST)   = PSNOWDZ     (JST-1)
    ZSNOWRHOO    (JST)   = PSNOWRHO    (JST-1)
    ZSNOWHEATO   (JST)   = PSNOWHEAT   (JST-1)
    ZSNOWDIAMOPTO(JST)   = PSNOWDIAMOPT(JST-1)
    ZSNOWSPHERIO (JST)   = PSNOWSPHERI (JST-1)
    ZSNOWHISTO   (JST)   = PSNOWHIST   (JST-1)
    ZSNOWAGEO    (JST)   = PSNOWAGE    (JST-1)
  ENDDO
  ! The new layer takes properties of fresh snow
  ZSNOWDZO     (1) = PSNOWDZF
  ZSNOWRHOO    (1) = PSNOWRHOF
  ZSNOWHEATO   (1) = PSNOWHEATF
  ZSNOWDIAMOPTO(1) = PSNOWDIAMOPTF
  ZSNOWSPHERIO (1) = PSNOWSPHERIF
  ZSNOWHISTO   (1) = PSNOWHISTF
  ZSNOWAGEO    (1) = PSNOWAGEF
  
  DO JIMP=1,NIMPUR
     DO JST = 2,INLVLS_OLD
      ZSNOWIMPURO(JST,JIMP)=PSNOWIMPUR(JST-1,JIMP)
     ENDDO
      ZSNOWIMPURO (1,JIMP)=PSNOWIMPURF(JIMP)
  ENDDO
ELSE
  ! first init without any change of the properties
  DO JST = 1,INLVLS_OLD
    ZSNOWDZO     (JST) = PSNOWDZ     (JST)
    ZSNOWRHOO    (JST) = PSNOWRHO    (JST)
    ZSNOWHEATO   (JST) = PSNOWHEAT   (JST)
    ZSNOWDIAMOPTO(JST) = PSNOWDIAMOPT(JST)
    ZSNOWSPHERIO (JST) = PSNOWSPHERI  (JST)
    ZSNOWHISTO   (JST) = PSNOWHIST   (JST)
    ZSNOWAGEO    (JST) = PSNOWAGE    (JST)
  ENDDO
  !IMPURITIES
  DO JIMP=1,NIMPUR
     DO JST = 1,INLVLS_OLD
       ZSNOWIMPURO (JST,JIMP)=PSNOWIMPUR(JST,JIMP)
     ENDDO
  ENDDO
ENDIF
!
! 1. Calculate vertical grid limits (m):
! --------------------------------------
!
ZSNOWZTOP_OLD(1) = ZPSNOW_OLD
ZSNOWZTOP_NEW(1) = ZPSNOW_NEW
!
DO JST = 1,INLVLS_OLD
  IF ( JST>1 ) ZSNOWZTOP_OLD(JST) = ZSNOWZBOT_OLD(JST-1)
  ZSNOWZBOT_OLD(JST) = ZSNOWZTOP_OLD(JST) - ZSNOWDZO(JST)
ENDDO
!
DO JST = 1,INLVLS_NEW
  IF ( JST>1 ) ZSNOWZTOP_NEW(JST) = ZSNOWZBOT_NEW(JST-1)
  ZSNOWZBOT_NEW(JST) = ZSNOWZTOP_NEW(JST) - PSNOWDZN(JST)
ENDDO
!

! Check consistency
IF ( ABS(ZSNOWZBOT_OLD(INLVLS_OLD)) > 0.00001 ) THEN
   WRITE (*,*) 'Error bottom OLD'
   OSUCCESS = .FALSE.
END IF
!
ZSNOWZBOT_OLD(INLVLS_OLD) = 0.
!
! Check consistency
IF ( ABS(ZSNOWZBOT_NEW(INLVLS_NEW)) > 0.00001 ) THEN
 WRITE (*,*) 'Error bottom NEW'
 OSUCCESS = .FALSE.
END IF
!
ZSNOWZBOT_NEW(INLVLS_NEW) = 0.
!
! 3. Calculate mass, heat, charcateristics mixing due to vertical grid resizing:
! --------------------------------------------------------------------
!
! loop over the new snow layers
! Summ or avergage of the constituting quantities of the old snow layers
! which are totally or partially inserted in the new snow layer
 CALL GET_MASS_HEAT(KJ,INLVLS_NEW,INLVLS_OLD,                                &
                    ZSNOWZTOP_OLD,ZSNOWZTOP_NEW,ZSNOWZBOT_OLD,ZSNOWZBOT_NEW, &
                    ZSNOWRHOO,ZSNOWDZO,ZSNOWDIAMOPTO,ZSNOWSPHERIO,ZSNOWHISTO,   &
                    ZSNOWAGEO,ZSNOWIMPURO,ZSNOWHEATO,                        &
                    ZSNOWRHON,PSNOWDZN,ZSNOWDIAMOPTN,ZSNOWSPHERIN,ZSNOWHISTN,   &
                    ZSNOWAGEN,ZSNOWIMPURN,ZSNOWHEATN             )
!
! check of consistency between new and old snowpacks

ZSNOWHEAN  = SUM(ZSNOWHEATN(1:INLVLS_NEW))
ZMASTOTN   = SUM(ZSNOWRHON(1:INLVLS_NEW)*PSNOWDZN(1:INLVLS_NEW))
ZPSNOW_NEW = SUM(PSNOWDZN(1:INLVLS_NEW))
!
ZSNOWHEAO  = SUM(ZSNOWHEATO(1:INLVLS_OLD))
ZMASTOTO   = SUM(ZSNOWRHOO(1:INLVLS_OLD)*ZSNOWDZO(1:INLVLS_OLD))
ZPSNOW_OLD = SUM(ZSNOWDZO(1:INLVLS_OLD))

IF ( INLVLS_OLD >= 1 .AND. (ABS(ZSNOWHEAN-ZSNOWHEAO)/ABS(ZSNOWHEAO))>1E-5) THEN
   OSUCCESS = .FALSE.
  WRITE(*,*) 'Warning diff heat', (ABS(ZSNOWHEAN-ZSNOWHEAO)/ABS(ZSNOWHEAO))
  WRITE(*,*) 'new value ',ZSNOWHEAN,' old value ',ZSNOWHEAO
ENDIF
IF ( INLVLS_OLD >= 1 .AND. (ABS(ZMASTOTN-ZMASTOTO)/ABS(ZMASTOTO))>1E-5) THEN
  OSUCCESS = .FALSE.
  WRITE(*,*) 'Warning diff mass', (ABS(ZMASTOTN-ZMASTOTO)/ABS(ZMASTOTO))
  WRITE(*,*) 'new value ',ZMASTOTN,' old value ',ZMASTOTO
  !
    WRITE(*,*)'INLVLS_OLD=',INLVLS_OLD  
    WRITE(*,*)'INLVLS_NEW=',INLVLS_NEW
    WRITE(*,*)'PSNOWDZ=',PSNOWDZ
    WRITE(*,*)'PSNOWDZF=',PSNOWDZF
    WRITE(*,*)'PSNOWDZN=',PSNOWDZN

ENDIF
IF ( INLVLS_OLD >= 1 .AND. (ABS(ZPSNOW_NEW-ZPSNOW_OLD)/ABS(ZPSNOW_OLD))>1E-5) THEN
   OSUCCESS = .FALSE.
  WRITE(*,*) 'Warning diff depth', (ABS(ZPSNOW_NEW-ZPSNOW_OLD)/ABS(ZPSNOW_OLD))
  WRITE(*,*) 'new value ',ZPSNOW_NEW,' old value ',ZPSNOW_OLD
ENDIF
!
! 5. Update mass (density and thickness) and heat:
! ------------------------------------------------
!
PSNOWDZ     (:) = PSNOWDZN     (:)
!
PSNOWRHO    (:) = ZSNOWRHON    (:)
PSNOWHEAT   (:) = ZSNOWHEATN   (:)
PSNOWDIAMOPT(:) = ZSNOWDIAMOPTN(:)
PSNOWSPHERI (:) = ZSNOWSPHERIN (:)
PSNOWHIST   (:) = ZSNOWHISTN   (:)
!
PSNOWAGE  (:) =  ZSNOWAGEN (:)
DO JIMP=1,NIMPUR
  PSNOWIMPUR(:,JIMP)=ZSNOWIMPURN(:,JIMP)
ENDDO
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWNLGRIDFRESH_1D',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWNLGRIDFRESH_1D
!####################################################################
SUBROUTINE SNOWCROUNLOAD(PSNOWDZ, PSNOWRHO, PSNOWHEAT, PSNOWDIAMOPT, PSNOWSPHERI, &
                         PSNOWHIST, PSNOWAGE, PSNOWIMPUR, PUNLOAD, PTA, PTSTEP, KNLVLS_USE)
!
! Matthieu Lafaysse
! 4th september 2024
!
! For snow under forest, unloading is added to the snowpack AFTER the thermal
! resolution to prevent discreapancies between MEB and Crocus resolution
!
USE MODD_CSTS,     ONLY : XLMTT, XTT, XCI
USE MODD_SNOW_PAR, ONLY : XSNOWDMIN, XSNOWFALL_THRESHOLD, XDZMIN_TOP_EXTREM, &
                          XRHOTHRESHOLD_ICE, XDZ1, XDIFF_1
USE MODD_SNOW_METAMO, ONLY : XUEPSI, XUEPSI_SMP
USE MODE_SNOW3L, ONLY : SNOW3LDIFTYP
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDZ, PSNOWRHO, PSNOWHEAT, PSNOWDIAMOPT, PSNOWSPHERI,&
                                       PSNOWHIST, PSNOWAGE ! snow prognostic variables
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PSNOWIMPUR
REAL, DIMENSION(:), INTENT(IN) :: PUNLOAD ! unloading rate
REAL, DIMENSION(:), INTENT(IN) :: PTA ! air temperature
REAL, INTENT(IN) :: PTSTEP ! time step
INTEGER, DIMENSION(:), INTENT(INOUT) :: KNLVLS_USE
!
INTEGER :: INLVLSMAX
!
! Properties of unloaded snow
REAL, PARAMETER :: PPRHO_UNLOAD = 200. ! density
REAL, PARAMETER :: PPDIAMOPT_UNLOAD = 2.E-4 ! optical diameter, values from Bouchard et al. (2024) for snow older than 3.1 days
REAL, PARAMETER :: PPSNOWSPHER_UNLOAD = 0.95 ! sphericity
!
REAL, DIMENSION(SIZE(PSNOWDZ, 1), SIZE(PSNOWDZ, 2)) :: ZSNOWDZO, ZSNOWDZN ! new snow profile after unloading
REAL, DIMENSION(SIZE(PSNOWDZ, 1)) :: ZSNOW, ZSNOWTEMP, ZSNOWHEATF, ZDZUNLOAD, ZSCAP ! temperature of new snow
REAL, DIMENSION(SIZE(PSNOWIMPUR, 3)) :: ZSNOWIMPURF
INTEGER, DIMENSION(SIZE(PSNOWDZ, 1)) :: INLVLS_USE_OLD
REAL :: ZDIFTYPE_SUP
LOGICAL :: GSUCCESS
!
INTEGER :: JJ
!
INLVLSMAX = SIZE(PSNOWDZ, 2)
!
ZSNOWDZN(:,:) = 0.
ZSNOWDZO(:,:) = PSNOWDZ
!
ZDZUNLOAD = PUNLOAD * PTSTEP / PPRHO_UNLOAD
!
ZSNOWIMPURF(:) = 0
!
INLVLS_USE_OLD = KNLVLS_USE
!
DO JJ=1, SIZE(PSNOWDZ, 1)
!VV  IF (ZDZUNLOAD(JJ) < XUEPSI) THEN
  IF (ZDZUNLOAD(JJ) < XUEPSI_SMP) THEN
    CYCLE
  ELSE
    ZSNOW(JJ) = SUM(PSNOWDZ(JJ, 1:KNLVLS_USE(JJ))) + ZDZUNLOAD(JJ)
    !Case of unloading on a previously snow surface
    IF ( KNLVLS_USE(JJ)>0 ) THEN
      ZSCAP    (JJ) = XCI*PSNOWRHO(JJ,1)
      ZSNOWTEMP(JJ) = XTT + ( PSNOWHEAT(JJ,1) + XLMTT*PSNOWRHO(JJ,1)*PSNOWDZ(JJ,1) ) / &
                            ( ZSCAP(JJ) * MAX( XSNOWDMIN/INLVLSMAX, PSNOWDZ(JJ,1) ) )
    ELSE  ! case with bare ground
      ZSNOWTEMP(JJ) = PTA(JJ)
    ENDIF
    ZSNOWTEMP(JJ) = MIN( XTT, ZSNOWTEMP(JJ) )
    ZSNOWHEATF(JJ) = PUNLOAD(JJ)*(XCI*(ZSNOWTEMP(JJ)-XTT)-XLMTT)*PTSTEP
    ! Measure difference between unloaded snow and surface layer
    IF (PSNOWDZ(JJ,1) < XDZ1) THEN ! test for optimization
      ZDIFTYPE_SUP = SNOW3LDIFTYP( PSNOWDIAMOPT(JJ,1), PPDIAMOPT_UNLOAD,  &
                                   PSNOWSPHERI(JJ,1), PPSNOWSPHER_UNLOAD, &
                                   PSNOWHIST(JJ,1), 0., PSNOWRHO(JJ,1),   &
                                   PPRHO_UNLOAD, PSNOWAGE(JJ,1), 0.)
    ELSE
      ZDIFTYPE_SUP = 0.
    ENDIF
    !
    IF ( KNLVLS_USE(JJ) > 0 .AND. ( KNLVLS_USE(JJ) == INLVLSMAX .OR. &
        ( ZDIFTYPE_SUP < XDIFF_1 .AND. PSNOWDZ(JJ,1) < XDZ1 ) .OR. &
        ( ZDZUNLOAD(JJ) < XSNOWFALL_THRESHOLD .AND. PSNOWDZ(JJ,1) < 2.*XDZ1 ) .OR. &
        ( PSNOWDZ(JJ,1) < XDZMIN_TOP_EXTREM .AND. PSNOWRHO(JJ,1) < XRHOTHRESHOLD_ICE))) THEN
      ! IF no bare ground and
      ! too many snow layers or
      ! unloaded snow similar to thin surface layer or
      ! very low rate of unloading and not too deep surface layer or
      ! extremely thin surface layer which is not ice
      ! THEN
      ! aggregation with previous layer
      ZSNOWDZN(JJ, 1) = ZDZUNLOAD(JJ) + PSNOWDZ(JJ, 1)
      ZSNOWDZN(JJ, 2:KNLVLS_USE (JJ)) = PSNOWDZ(JJ, 2:KNLVLS_USE (JJ))
      !
    ELSE
      ! other cases: new layer
      ZSNOWDZN(JJ, 1) = ZDZUNLOAD(JJ)
      ZSNOWDZN(JJ, 2:KNLVLS_USE (JJ) + 1) = PSNOWDZ(JJ, 1:KNLVLS_USE (JJ))
      KNLVLS_USE (JJ) = KNLVLS_USE (JJ) + 1
    ENDIF
    ! add unloaded snow
    CALL SNOWNLGRIDFRESH_1D(JJ,ZSNOW(JJ),PSNOWDZ(JJ,:),ZSNOWDZN(JJ,:), &
                            PSNOWRHO(JJ,:),PSNOWHEAT(JJ,:),PSNOWDIAMOPT(JJ,:),PSNOWSPHERI(JJ,:), &
                            PSNOWHIST(JJ,:),PSNOWAGE(JJ,:),PSNOWIMPUR(JJ,:,:), .TRUE.,   &
                            PPRHO_UNLOAD, ZDZUNLOAD(JJ), ZSNOWHEATF(JJ),PPDIAMOPT_UNLOAD, &
                            PPSNOWSPHER_UNLOAD, 0.,0., ZSNOWIMPURF, &
                            KNLVLS_USE(JJ), INLVLS_USE_OLD(JJ), GSUCCESS)
    !IF (.NOT. GSUCCESS) THEN
    !  PRINT*, 'regridding problem in unloading'
    !  PRINT*, 'ZDZUNLOAD=', ZDZUNLOAD(JJ)
    !  PRINT*, 'JJ=', JJ
    !  PRINT*, 'ZSNOWDZO=',ZSNOWDZO(JJ,:)
    !  PRINT*, 'ZSNOWDZN=',ZSNOWDZN(JJ,:)
    !END IF
  END IF
END DO
!
!
!
END SUBROUTINE SNOWCROUNLOAD
!####################################################################
SUBROUTINE SNOWCROFREEZINGRAIN(OFRZRAIN, PSNOWDZ, PSNOWRHO, PSNOWHEAT, PSNOWDIAMOPT, PSNOWSPHERI, &
                         PSNOWHIST, PSNOWAGE, PSNOWIMPUR, PRR, PTA, PTSTEP, KNLVLS_USE)
!
! Matthieu Lafaysse
! 5th september 2024
!
! Freezing rain is added to the snowpack AFTER the thermal
! resolution to prevent discreapancies between MEB and Crocus resolution
!
USE MODD_CSTS,        ONLY : XLMTT, XTT, XCI, XRHOLI
USE MODD_SNOW_PAR,    ONLY : XSNOWDMIN, XSNOWFALL_THRESHOLD, XDZMIN_TOP_EXTREM, XDZ1, XDIFF_1
USE MODD_SNOW_METAMO, ONLY : XUEPSI, NVHIS2
USE MODE_SNOW3L,      ONLY : SNOW3LDIFTYP
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDZ, PSNOWRHO, PSNOWHEAT, PSNOWDIAMOPT, PSNOWSPHERI,&
                                       PSNOWHIST, PSNOWAGE ! snow prognostic variables
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PSNOWIMPUR
REAL, DIMENSION(:), INTENT(IN) :: PRR ! freezing precipitation rate
REAL, DIMENSION(:), INTENT(IN) :: PTA ! air temperature
REAL, INTENT(IN) :: PTSTEP ! time step
INTEGER, DIMENSION(:), INTENT(INOUT) :: KNLVLS_USE
!
INTEGER :: INLVLSMAX
!
! Properties of unloaded snow
REAL, PARAMETER :: PPRHO_FRZ = 917. ! density
REAL, PARAMETER :: PPDIAMOPT_FRZ = 2.E-3 ! optical diameter
REAL, PARAMETER :: PPSNOWSPHER_FRZ = 1.0 ! sphericity
!
REAL, DIMENSION(SIZE(PSNOWDZ, 1), SIZE(PSNOWDZ, 2)) :: ZSNOWDZO, ZSNOWDZN ! new snow profile after unloading
REAL, DIMENSION(SIZE(PSNOWDZ, 1)) :: ZSNOW, ZSNOWTEMP, ZSNOWHEATF, ZDZICE ! temperature of new snow
REAL, DIMENSION(SIZE(PSNOWIMPUR, 3)) :: ZSNOWIMPURF
INTEGER, DIMENSION(SIZE(PSNOWDZ, 1)) :: INLVLS_USE_OLD
LOGICAL :: GSUCCESS
REAL :: ZDIFTYPE_SUP
!
LOGICAL, DIMENSION(:), INTENT(IN) :: OFRZRAIN
INTEGER :: JJ
!
INLVLSMAX = SIZE(PSNOWDZ, 2)
!
ZSNOWDZN(:,:) = 0.
ZSNOWDZO(:,:) = PSNOWDZ
!
ZSNOWIMPURF(:) = 0
!
INLVLS_USE_OLD = KNLVLS_USE
!
DO JJ=1, SIZE(PSNOWDZ, 1)
  IF (.NOT. OFRZRAIN(JJ)) THEN
    CYCLE
  ELSE
    ZDZICE(JJ) = PRR(JJ) * PTSTEP / XRHOLI
    ZSNOW(JJ) = SUM(PSNOWDZ(JJ, 1:KNLVLS_USE(JJ))) + ZDZICE(JJ)
    !Case of unloading on a previously snow surface
    ZSNOWTEMP(JJ) = XTT
    ZSNOWHEATF(JJ) = PRR(JJ) * ( XCI * ( ZSNOWTEMP(JJ)-XTT ) - XLMTT ) * PTSTEP
!
! Measure difference between unloaded snow and surface layer
    IF(PSNOWDZ(JJ,1)<   XDZ1) THEN ! test for optimization   
      ZDIFTYPE_SUP = SNOW3LDIFTYP( PSNOWDIAMOPT(JJ,1),PPDIAMOPT_FRZ, &
                                 PSNOWSPHERI(JJ,1),PPSNOWSPHER_FRZ, &
                                 PSNOWHIST(JJ,1),REAL(NVHIS2),&
                                 PSNOWRHO(JJ,1),XRHOLI,&
                                 PSNOWAGE(JJ,1),0.)
    ELSE
      ZDIFTYPE_SUP = 0.
    ENDIF
!
    IF ((KNLVLS_USE(JJ) > 0) .AND. &
        ((KNLVLS_USE(JJ) == INLVLSMAX) .OR. &
        ( ZDIFTYPE_SUP<XDIFF_1        .AND. PSNOWDZ(JJ,1)<   XDZ1 ) .OR. &
        ((PSNOWDZ(JJ,1)<XDZMIN_TOP_EXTREM) ))) THEN
        ! too many snow layers
        ! or freezing rain similar to thin surface layer
        ! or extremely thin surface layer
        ! but no bare ground
        !
        ! aggregation with previous layer
        ZSNOWDZN(JJ, 1) = ZDZICE(JJ) + PSNOWDZ(JJ, 1)
        ZSNOWDZN(JJ, 2:KNLVLS_USE (JJ)) = PSNOWDZ(JJ, 2:KNLVLS_USE (JJ))

    ELSE
      ! other cases: new layer
      ZSNOWDZN(JJ, 1) = ZDZICE(JJ)
      ZSNOWDZN(JJ, 2:KNLVLS_USE (JJ) + 1) = PSNOWDZ(JJ, 1:KNLVLS_USE (JJ))
      KNLVLS_USE (JJ) = KNLVLS_USE (JJ) + 1
    ENDIF
    ! add unloaded snow
    CALL SNOWNLGRIDFRESH_1D(JJ,ZSNOW(JJ),PSNOWDZ(JJ,:),ZSNOWDZN(JJ,:), &
                            PSNOWRHO(JJ,:),PSNOWHEAT(JJ,:),PSNOWDIAMOPT(JJ,:),PSNOWSPHERI(JJ,:), &
                            PSNOWHIST(JJ,:),PSNOWAGE(JJ,:),PSNOWIMPUR(JJ,:,:), .TRUE.,   &
                            XRHOLI, ZDZICE(JJ), ZSNOWHEATF(JJ),PPDIAMOPT_FRZ, &
                            PPSNOWSPHER_FRZ, REAL(NVHIS2),0., ZSNOWIMPURF, &
                            KNLVLS_USE(JJ), INLVLS_USE_OLD(JJ), GSUCCESS)
!    IF (.NOT. GSUCCESS) THEN
!      PRINT*, 'regridding problem in snowcrofreezingrain'
!      PRINT*, 'JJ=', JJ
!      PRINT*, 'ZDZICE=', ZDZICE(JJ)
!      PRINT*, 'ZSNOWDZO=',ZSNOWDZO(JJ,:)
!      PRINT*, 'ZSNOWDZN=',ZSNOWDZN(JJ,:)
!    END IF
  END IF
END DO
!
END SUBROUTINE SNOWCROFREEZINGRAIN
!####################################################################
SUBROUTINE SNOWCROALB(PALBEDOSC,PSPECTRALALBEDO,PSNOWDZ,          &
                      PSNOWRHO,PPERMSNOWFRAC,                     &
                      PSNOWDIAMOPT_TOP,PSNOWSPHERI_TOP,PSNOWAGE_TOP, &
                      PSNOWDIAMOPT_BOT,PSNOWSPHERI_BOT,PSNOWAGE_BOT, &
                      PPS, PZENITH, OSNOWAGING_VAR,                  &
                      PAGINGCOEF, KNLVLS_USE        )                       
!
!!    PURPOSE
!!    -------
!     Calculate the snow surface albedo. Use the method of original
!     Crocus which considers a specified spectral distribution of solar
!     solar radiation (to be replaced by an input forcing when available)
!     In addition to original crocus, the top 2 surface snow layers are
!     considered in the calculation, using an arbitrary weighting, in order
!     to avoid time discontinuities due to layers agregation
!     Ageing depends on the presence of permanent snow cover
!
USE MODD_SNOW_PAR, ONLY : XANSMIN,                             &
                          XVRPRE1,XVRPRE2,XVAGING_NOGLACIER,   &
                          XVSPEC1,XVSPEC2,                     &
                          XVSPEC3, XVW1,XVW2,XVD1,XVD2,        &
                          XRHOTHRESHOLD_ICE
!
USE MODE_SNOW3L
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN)    :: PSNOWDZ,PPERMSNOWFRAC
!
REAL,DIMENSION(:,:), INTENT(IN)   :: PSNOWRHO ! For now only the 2 first layers are required
!
REAL, DIMENSION(:), INTENT(INOUT) :: PALBEDOSC
!
REAL, DIMENSION(:,:), INTENT(OUT) :: PSPECTRALALBEDO   ! Albedo in the different spectral bands
!
REAL, DIMENSION(:), INTENT(IN)    :: PSNOWDIAMOPT_TOP,PSNOWSPHERI_TOP,PSNOWAGE_TOP, &
                                     PSNOWDIAMOPT_BOT,PSNOWSPHERI_BOT,PSNOWAGE_BOT, PPS
INTEGER, DIMENSION(:), INTENT(IN) :: KNLVLS_USE
!
REAL, DIMENSION(:), INTENT(IN)    :: PZENITH ! solar zenith angle for future use
!
LOGICAL, INTENT(IN)               :: OSNOWAGING_VAR    ! True = Spatially-variable snow aging coefficnent is used
                                                       ! False = Value  provided in namelist are used (XVAGING_NOGLACIER and XVAGING_GLACIER )
REAL, DIMENSION(:), INTENT(IN)    :: PAGINGCOEF ! Spatially variable snow aging coefficient 
!
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(3,SIZE(PSNOWRHO,1)) :: ZALB_TOP, ZALB_BOT
!
REAL, DIMENSION(SIZE(PSNOWRHO,1))   :: ZMIN, ZMAX
REAL, DIMENSION(SIZE(PSNOWRHO,1))   :: ZFAC_TOP, ZFAC_BOT
!
REAL, DIMENSION(SIZE(PALBEDOSC))  :: ZVAGE1
!
INTEGER         :: JJ   ! looping indexes
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROALB',0,ZHOOK_HANDLE)
!
! 0. Initialize:
! ------------------
!
IF(OSNOWAGING_VAR) THEN
  ZVAGE1(:)  = PAGINGCOEF(:) ! Snow darkening parameter varies in space   
                             ! according to climatology of LAP-deposition on snow
                             ! (Gaillard et al., 2025)  
ELSE
  ZVAGE1(:)  = XVAGING_NOGLACIER! Snow darkening parameter fixed in space 
ENDIF 

DO JJ=1, SIZE(PALBEDOSC)
  !case with snow on the ground
  IF (PSNOWRHO(JJ,1)<XRHOTHRESHOLD_ICE ) THEN
    !
    CALL GET_ALB_SNOW(PPS(JJ),ZVAGE1(JJ),PSNOWDIAMOPT_TOP(JJ),&
                      PSNOWAGE_TOP(JJ),ZALB_TOP(:,JJ),        &
                      OSNOWAGING_VAR)
  ELSE
    !
    CALL GET_ALB_ICE(ZALB_TOP(:,JJ))
  ENDIF
  !
  IF (PSNOWRHO(JJ,2)<XRHOTHRESHOLD_ICE .AND. KNLVLS_USE(JJ) .NE. 1) THEN
    !
    CALL GET_ALB_SNOW(PPS(JJ),ZVAGE1(JJ),PSNOWDIAMOPT_BOT(JJ),&
               MIN(365.,PSNOWAGE_BOT(JJ)),ZALB_BOT(:,JJ),     &
                OSNOWAGING_VAR)
  ELSE
    !
    CALL GET_ALB_ICE(ZALB_BOT(:,JJ))
  ENDIF
    !

  IF (KNLVLS_USE(JJ)==1)  ZALB_BOT(:,JJ) = ZALB_TOP(:,JJ)

  !
ENDDO
DO JJ=1, SIZE(PALBEDOSC)
  ! computation of spectral albedo over 3 bands taking into account the respective
  ! depths of top layers
  ZMIN(JJ) = MIN( 1., PSNOWDZ(JJ)/XVD1 )
  ZMAX(JJ) = MAX( 0., (PSNOWDZ(JJ)-XVD1)/XVD2 )
  ZFAC_TOP(JJ) = XVW1 * ZMIN(JJ) + XVW2 * MIN( 1., ZMAX(JJ) )
  ZFAC_BOT(JJ) = XVW1 * ( 1. - ZMIN(JJ) ) + XVW2 * ( 1. - MIN( 1., ZMAX(JJ) ) )
  PSPECTRALALBEDO(JJ,1) = ZFAC_TOP(JJ) * ZALB_TOP(1,JJ) + ZFAC_BOT(JJ) * ZALB_BOT(1,JJ)
  PSPECTRALALBEDO(JJ,2) = ZFAC_TOP(JJ) * ZALB_TOP(2,JJ) + ZFAC_BOT(JJ) * ZALB_BOT(2,JJ)
  PSPECTRALALBEDO(JJ,3) = ZFAC_TOP(JJ) * ZALB_TOP(3,JJ) + ZFAC_BOT(JJ) * ZALB_BOT(3,JJ)
  !
  ! arbitrarily specified spectral distribution
  ! to be changed when solar radiation distribution is an input variable
  PALBEDOSC(JJ) = XVSPEC1 * PSPECTRALALBEDO(JJ,1) + &
                  XVSPEC2 * PSPECTRALALBEDO(JJ,2) + &
                  XVSPEC3 * PSPECTRALALBEDO(JJ,3)
  !
ENDDO ! end loop grid points
!
DO JJ=1, SIZE(PALBEDOSC)
  IF ( KNLVLS_USE(JJ)==0 )  &
    ! case with no snow on the ground
    PALBEDOSC(JJ) = XANSMIN
ENDDO

!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:SNOWCROALB',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOWCROALB
!####################################################################
SUBROUTINE GET_ALB_SNOW(PPS_IN,PVAGE1,PSNOWDIAMOPT,PSNOWAGE,PALB,OSNOWAGING_VAR)
!
USE MODD_SNOW_PAR, ONLY : XVALB2, XVALB3, XVALB4, XVALB5, &
                          XVALB6, XVALB7, XVALB8, XVALB9, &
                          XVALB10, XVALB11, XVDIOP1,      &
                          XVRPRE1, XVRPRE2, XVPRES1
!
IMPLICIT NONE
!
REAL, INTENT(IN) :: PPS_IN
REAL, INTENT(IN) :: PVAGE1
REAL, INTENT(IN) :: PSNOWDIAMOPT, PSNOWAGE
REAL, DIMENSION(3), INTENT(OUT) :: PALB
LOGICAL,INTENT(IN)  :: OSNOWAGING_VAR
!
REAL :: ZDIAM, ZDIAM_SQRT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:GET_ALB_SNOW',0,ZHOOK_HANDLE)
!

! Normal case (snow)
ZDIAM=PSNOWDIAMOPT
ZDIAM_SQRT = SQRT(ZDIAM)
PALB(1) = MIN( XVALB2 - XVALB3*ZDIAM_SQRT, XVALB4 )
PALB(2) = MAX( 0.3, XVALB5 - XVALB6*ZDIAM_SQRT )
ZDIAM   = MIN( ZDIAM, XVDIOP1 )
ZDIAM_SQRT = SQRT(ZDIAM)
PALB(3) = MAX( 0., XVALB7*ZDIAM - XVALB8*ZDIAM_SQRT + XVALB9 )
! AGE CORRECTION ONLY FOR VISIBLE BAND
!VV PALB(1) = MAX( XVALB11, PALB(1) - MIN( MAX(PPS_IN/XVPRES1,XVRPRE1), XVRPRE2 ) * &
!VV                 XVALB10 * PSNOWAGE / PVAGE1 )
IF(OSNOWAGING_VAR) THEN 
   ! Elevation effects are not taken into account 
   ! See Gaillard et al. (2025, TC)
   PALB(1) = MAX( XVALB11, PALB(1) - XVALB10 * PSNOWAGE / PVAGE1 )
ELSE
   !PALB(1) = MAX( XVALB11, PALB(1) - MIN( MAX(PPS_IN/XVPRES1,XVRPRE1), XVRPRE2 ) * &
   !        XVALB10 * PSNOWAGE / PVAGE1 )
   !VV Remove temporary the elevation effect for a fair comparaison with Manon approach
   PALB(1) = MAX( XVALB11, PALB(1) - XVALB10 * PSNOWAGE / PVAGE1 )
ENDIF         
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:GET_ALB_SNOW',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_ALB_SNOW
!####################################################################
SUBROUTINE GET_ALB_ICE(PALB)
!
USE MODD_SNOW_PAR, ONLY : XALBICE1, XALBICE2, XALBICE3
!
IMPLICIT NONE
!
REAL, DIMENSION(3), INTENT(OUT) :: PALB
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:GET_ALB_ICE',0,ZHOOK_HANDLE)
!
! Prescribed spectral albedo for surface ice
  PALB(1) = XALBICE1
  PALB(2) = XALBICE2
  PALB(3) = XALBICE3
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO:GET_ALB_ICE',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_ALB_ICE
!####################################################################
END MODULE MODE_SNOWCRO

