!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!auto_modi:spll_surface_aero_cond.D
MODULE MODI_SURFACE_AERO_COND
INTERFACE
    SUBROUTINE SURFACE_AERO_COND(PRI, PZREF, PUREF, PVMOD, PZ0,&
                                     PZ0H, PRSURF, PAC, PRA, PCH,HSNOWRES)
IMPLICIT NONE
REAL, DIMENSION(:), INTENT(IN)    :: PRI      ! Richardson number
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD    ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PZREF    ! reference height of the first
REAL, DIMENSION(:), INTENT(IN)    :: PUREF    ! reference height of the wind
REAL, DIMENSION(:), INTENT(IN)    :: PZ0      ! roughness length for momentum
REAL, DIMENSION(:), INTENT(IN)    :: PZ0H     ! roughness length for heat
REAL, DIMENSION(:), INTENT(OUT)   :: PAC      ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PRA      ! aerodynamical resistance
REAL, DIMENSION(:), INTENT(OUT)   :: PCH      ! drag coefficient for heat
CHARACTER(LEN=3), INTENT(IN)  ::HSNOWRES !surface exchange coefficient option 
REAL, DIMENSION(:)                :: PRSURF     ! aerodynamic surface resistance if Ta and HR above canopy


END SUBROUTINE SURFACE_AERO_COND
END INTERFACE
END MODULE MODI_SURFACE_AERO_COND
