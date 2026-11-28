!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!auto_modi:spll_wind_threshold.D
MODULE MODI_WIND_THRESHOLD
INTERFACE
    FUNCTION WIND_THRESHOLD(PWIND,PUREF) RESULT(PWIND_NEW)
IMPLICIT NONE
REAL, DIMENSION(:), INTENT(IN)   :: PWIND      ! wind
REAL, DIMENSION(:), INTENT(IN)   :: PUREF      ! forcing level
REAL, DIMENSION(SIZE(PWIND))     :: PWIND_NEW  ! modified wind
END FUNCTION WIND_THRESHOLD
END INTERFACE
END MODULE MODI_WIND_THRESHOLD
