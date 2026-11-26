!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-----------------------------------------------------------------
      MODULE MODD_BLOWSNW_SURF
!
!     Contains parameters for blowing snow simulation 
!      
        IMPLICIT NONE
!
!     Snow density : cf Liston et Sturm (1998)          
      REAL, PARAMETER   :: XRHO_DEP=280.0 
!
!     Deposited snow grains are round and small : s=1,r=0.4mm  
      REAL,PARAMETER    :: XSDEP_GRA1=99.0
      REAL,PARAMETER    :: XSDEP_GRA2=0.0004 
!      
!     Minimal mean radius (um) 
      REAL,PARAMETER    :: XINIRADIUS_SNW = 5.e-6
!     Minimum allowed number concentration (#/m3)
      REAL,PARAMETER    :: XN0MIN_SNW    =  1     

!     Parameters used in KC02 parameterization for settling velocity 
      REAL,PARAMETER    :: XC0     = 0.29
      REAL,PARAMETER    :: XDELTA0 = 9.06
 
!      Parameters used in Mitchell (96) parameterization for settling velocity  
      REAL,PARAMETER    :: XAM1     = 0.04394
      REAL,PARAMETER    :: XAM2     = 0.06049
      REAL,PARAMETER    :: XAM3     = 0.2072
      REAL,PARAMETER    :: XBM1     = 0.970
      REAL,PARAMETER    :: XBM2     = 0.831
      REAL,PARAMETER    :: XBM3     = 0.638
      REAL,PARAMETER    :: XBESTL_1 = 10.0
      REAL,PARAMETER    :: XBESTL_2 = 585.

!     Grain size distribution is a two parameter gamma distribution 
      REAL              :: XEMIRADIUS_SNW
      REAL              :: XEMIALPHA_SNW 
!      
      LOGICAL          :: LSNOW_SALT  ! Flag to fix snow concentration in the saltation layer      
      ! Snow concentration in the saltation layer (kg_{snow}/m3_{air}) if LSNOW_SALT = T
      REAL             :: XCONC_SALT 
      !
      LOGICAL          :: LSNOW_PRECIP  ! Flag to impose uniform and constant precipitation rate
      ! Fixed snow precipitation rate SWE (kg/ (m2 s)) if LSNOW_PRECIP = T
      REAL             :: XSNOW_PRECIP 
      LOGICAL          :: LBLOWSNW_ADV       ! Flag to account for advection effects on 
                                             ! total mass and number in Canopy variables and 
                                             ! to compute divergence of
                                             ! saltation flux
      LOGICAL          :: LBLOWSNW_CANOSUBL  ! Flag to activate sublimation of Canopy variables

      LOGICAL          :: LBLOWSNW_CANODIAG ! Flag to get additional diagnostic at Canopy levels : mean radius 
                                ! and number and mass variables in _/m3

      LOGICAL          :: LSNOW_WIND  ! Flag to activate effects of snow particle on wind profile
      ! Increase in surface rougnhess due to saltation following Dover (1993) z0_s = z0_ns*(u*/uth*)²
      REAL             :: XSNOW_ROUGHNESS         ! = 0 not activated; =1 activated   
      ! Buoyancy effect in presence of suspended particles of blowing snow.
      REAL             :: XSNOW_BUOYANCY          ! = 0 not activated; =1 activated  
      CHARACTER(LEN=4)       :: CSNOW_SALT ! Paramaterization to compute particle flux in saltation
!             'POME': Pomeroy and Gray, 1990
!             'SORE': Sorensen (1991) : used at SLF before Doorshot's model
!             'MANN': Concentration in the saltation layer is computed according to Mann
!                     parameterization for particle number density just above the ground (10 cm) 
      CHARACTER(LEN=4)    :: CSNOW_SEDIM ! Paramaterization to compute settling velocity
!              'CARR': follow Carrier's drag coefficient
!              'TABC': based on tabuleted values of Carrier's formulation for alpha=3
!              'MITC': based on Mitchell's formulation for settling velocity of ice spheres
!              'NONE': sedimentation is desactivated (for test only!)  
!
      REAL            :: XRSNOW_SBL ! Ratio between diffusion coefficient for scalar
                          ! variables and blowing snow variables
                          ! RSNOW_SBL = KSCA/KSNOW = 4.
                          ! See Vionnet (PhD, 2012, In French) and Vionnet et al (TC, 2014)
                          ! for a complete dicsussion
! 
     END MODULE MODD_BLOWSNW_SURF 

