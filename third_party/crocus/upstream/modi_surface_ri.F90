!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!auto_modi:spll_surface_ri.D
MODULE MODI_SURFACE_RI
INTERFACE
    SUBROUTINE SURFACE_RI(PTG, PQS, PEXNS, PEXNA, PTA, PQA,   &
                               PZREF, PUREF, PDIRCOSZW, PVMOD, PRI )  
IMPLICIT NONE
REAL, DIMENSION(:), INTENT(IN)    :: PTG      ! surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PQS      ! surface specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS    ! surface exner function
REAL, DIMENSION(:), INTENT(IN)    :: PTA      ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA      ! specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PEXNA    ! exner function
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD    ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PZREF    ! reference height of the first
REAL, DIMENSION(:), INTENT(IN)    :: PUREF    ! reference height of the wind
REAL, DIMENSION(:), INTENT(IN)    :: PDIRCOSZW! Cosine of the angle between
REAL, DIMENSION(:), INTENT(OUT)   :: PRI      ! Richardson number
END SUBROUTINE SURFACE_RI
END INTERFACE
END MODULE MODI_SURFACE_RI
