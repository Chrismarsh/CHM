!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!   ######################################################################
    SUBROUTINE SURFACE_AERO_COND(PRI, PZREF, PUREF, PVMOD, PZ0,    &
                                     PZ0H, PRESA_SN, PAC, PRA, PCH    ,HSNOWRES       )
!   ######################################################################
!
!!****  *SURFACE_AERO_COND*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the drag coefficients for heat and momentum near the ground
!         
!     
!!**  METHOD
!!    ------
!
!
!
!    1 and 2 : computation of relative humidity near the ground
!
!    3 : richardson number
!
!    4 : the aerodynamical resistance for heat transfers is deduced
!
!    5 : the drag coefficient for momentum ZCD is computed
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    MODD_CST
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!      V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/01/98 
!!                  02/04/01 (P Jabouille) limitation of Z0 with 0.5 PUREF
!!                  05/2016 M. Lafaysse - B. Cluzet : implement Martin and Lejeune 1998 formulation for multiphysics
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XKARMAN
USE MODI_WIND_THRESHOLD
!
USE MODE_THERMOS
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
REAL, DIMENSION(:), INTENT(IN)    :: PRI      ! Richardson number
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD    ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PZREF    ! reference height of the first
                                              ! atmospheric level
REAL, DIMENSION(:), INTENT(IN)    :: PUREF    ! reference height of the wind
                                              ! NOTE this is different from ZZREF
                                              ! ONLY in stand-alone/forced mode,
                                              ! NOT when coupled to a model (MesoNH)
REAL, DIMENSION(:), INTENT(IN)    :: PZ0      ! roughness length for momentum
REAL, DIMENSION(:), INTENT(IN)    :: PZ0H     ! roughness length for heat
REAL, DIMENSION(:), INTENT(INOUT)    :: PRESA_SN ! aerodynamic surface resistance from external scheme (0 if not used) 

!
REAL, DIMENSION(:), INTENT(OUT)   :: PAC      ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PRA      ! aerodynamical resistance
REAL, DIMENSION(:), INTENT(OUT)   :: PCH      ! drag coefficient for heat
!
CHARACTER(LEN=3), INTENT(IN), OPTIONAL :: HSNOWRES !surface exchange coefficient option 
!
!
!*      0.2    declarations of local variables
!
!
REAL, DIMENSION(SIZE(PRI)) :: ZZ0, ZZ0H, ZMU, ZFH, ZCHSTAR, &
                              ZPH, ZCDN, ZSTA, ZDI, ZWORK1, &
                              ZWORK2, ZWORK3, ZVMOD,        &
                              ZMARTIN, ZCDN_M98
!
CHARACTER(LEN=3) :: YSNOWRES
!
INTEGER         :: INI, JI
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
! Functions:
!
REAL :: ZX, CHSTAR, PH
!
CHSTAR(ZX) = 3.2165 + 4.3431*ZX + 0.5360*ZX*ZX - 0.0781*ZX*ZX*ZX
PH    (ZX) = 0.5802 - 0.1571*ZX + 0.0327*ZX*ZX - 0.0026*ZX*ZX*ZX
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SURFACE_AERO_COND',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!*       1.     Surface aerodynamic resistance for heat transfers
!               -------------------------------------------------
!
INI = SIZE(PRI)
!
IF(PRESENT(HSNOWRES))THEN
   YSNOWRES=HSNOWRES
ELSE
   YSNOWRES='DEF'
ENDIF
!
ZVMOD(:) = WIND_THRESHOLD(PVMOD(:),PUREF(:))
!
DO JI=1,INI
!
  ZZ0 (JI) = MIN(PZ0(JI),PUREF(JI)*0.5)
  ZZ0H(JI) = MIN(ZZ0(JI),PZ0H(JI))
  ZZ0H(JI) = MIN(ZZ0H(JI),PZREF(JI)*0.5)
!
  ZWORK1(JI)=LOG( PUREF(JI)/ZZ0(JI) )
  ZWORK2(JI)=PZREF(JI)/ZZ0H(JI)
  ZWORK3(JI)=ZVMOD(JI)*ZVMOD(JI)
  ZMU   (JI) = MAX(LOG(ZZ0(JI)/ZZ0H(JI)),0.0)
  ZFH   (JI) = ZWORK1(JI) / LOG(ZWORK2(JI))
!
  ZCHSTAR(JI) = CHSTAR(ZMU(JI))
  ZPH    (JI) = PH(ZMU(JI))
!
ENDDO
!
IF(YSNOWRES=='RIL' .OR. YSNOWRES=='DEF' .OR. YSNOWRES=='RI1' .OR. YSNOWRES=='RI2') THEN
!
  DO JI=1,INI
!
     ZCDN(JI) = (XKARMAN/ZWORK1(JI))**2
     ZSTA(JI) = PRI(JI)*ZWORK3(JI)
!
     IF (PRI(JI)< 0.0) THEN
        ZDI(JI) = 1.0 / ( ZVMOD(JI)+ZCHSTAR(JI)*ZCDN(JI)*15.*ZWORK2(JI)**ZPH(JI)*ZFH(JI)*SQRT(-ZSTA(JI)) ) 
        PAC(JI) = ZCDN(JI)*(ZVMOD(JI)-15.0*ZSTA(JI)*ZDI(JI))*  ZFH(JI)
     ELSE
        ZDI(JI) = SQRT(ZWORK3(JI) + 5.0* ZSTA(JI))
        PAC(JI) = ZCDN(JI)*ZVMOD(JI)/(1.0+15.0*ZSTA(JI)*ZDI(JI)/ZWORK3(JI)/ZVMOD(JI))*ZFH(JI)    
     ENDIF
!
     IF (PRESA_SN(JI) < 0) THEN  ! Use formulation from SURFEX
        PRA(JI) = 1. / PAC(JI) 
     ELSE ! Use resistance calculated outside SURFEX
        PRA(JI) = PRESA_SN(JI)
        PAC(JI) = 1. / PRA(JI) 
      ENDIF
     PCH(JI) = 1. / (PRA(JI) * ZVMOD(JI))

     ! output RESA_SN
     PRESA_SN(JI) = PRA(JI)
!
  ENDDO
!
ELSEIF (YSNOWRES=='M98')THEN
!
!   Martin and Lejeune 1998 ; Cluzet et al 2016
!
  DO JI=1,INI
!
     ZCDN_M98(JI)= XKARMAN*XKARMAN/(LOG(PUREF(JI)/ZZ0(JI))*LOG(PZREF(JI)/ZZ0(JI)))
!
     IF (PRI(JI)<0.0) THEN
        IF (ZCDN_M98(JI)==0.) THEN
           ZMARTIN(JI) = 1.0
        ELSE
           ZMARTIN(JI) = 1. + (7./  ( 0.83*(ZCDN_M98(JI))**(-0.62) )  ) &
                       * LOG(1.-0.83*(ZCDN_M98(JI))**(-0.62)*PRI(JI))
        ENDIF          
     ELSE
        IF (PRI(JI)<0.2) THEN
           ZMARTIN(JI) = MAX(0.75,( 1. - 5.*PRI(JI))**2)!
           !Nota B. Cluzet : le min servait à empêcher le CH de remonter pour 
           !des Ri >0.4 c'est une erreur car cela seuille le Ch dès Ri =0 au lieu de le faire dès Ri=0.026
        ELSE
           ZMARTIN(JI) = 0.75
        ENDIF
     ENDIF
!    
     IF (PRESA_SN(JI) < 0) THEN ! Use formulation from SURFEX
       PCH(JI) = ZMARTIN(JI) * ZCDN_M98(JI)
       PRA(JI)=1./(PCH(JI)*ZVMOD(JI))  ! Nota B. Cluzet : checked in noilhan and mahfouf seems ok
     ELSE ! In the forest, resistance calculated in drag_svs2
       PRA(JI) = PRESA_SN(JI)
       PCH(JI)=1./(PRA(JI)*ZVMOD(JI))  ! Nota B. Cluzet : checked in noilhan and mahfouf seems ok
     ENDIF

     PAC(JI) = 1./(PRA(JI))

     ! output RESA_SN
     PRESA_SN(JI) = PRA(JI)
!
  ENDDO
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('SURFACE_AERO_COND',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SURFACE_AERO_COND
