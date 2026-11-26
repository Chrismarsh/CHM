!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!auto_modi:spll_abor1_sfx.D
MODULE MODI_ABOR1_SFX
INTERFACE
      SUBROUTINE ABOR1_SFX(YTEXT)
IMPLICIT NONE
 CHARACTER(LEN=*),  INTENT(IN)  :: YTEXT

END SUBROUTINE ABOR1_SFX
END INTERFACE

END MODULE MODI_ABOR1_SFX


!     #############################################################
      SUBROUTINE ABOR1_SFX(YTEXT)
!     #############################################################
!
!!****  *ABOR1_SFX* - abor1 subroutine
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      P. Le Moigne   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/2008 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=*),  INTENT(IN)  :: YTEXT
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
WRITE(*,*)YTEXT

STOP
!
END SUBROUTINE ABOR1_SFX
