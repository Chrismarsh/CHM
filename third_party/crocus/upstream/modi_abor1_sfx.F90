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
USE ISO_C_BINDING, ONLY: C_INT
USE ISO_FORTRAN_ENV, ONLY: OUTPUT_UNIT, ERROR_UNIT
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=*),  INTENT(IN)  :: YTEXT
!
INTERFACE
  SUBROUTINE C_EXIT_IMMEDIATELY(STATUS) BIND(C, NAME="_Exit")
    IMPORT C_INT
    INTEGER(C_INT), VALUE :: STATUS
  END SUBROUTINE
END INTERFACE
!
! In CHM a process exit handler can enter collective HDF5 cleanup while the
! other MPI ranks are still in the timestep. Let the launcher see the failure.
FLUSH(OUTPUT_UNIT)
WRITE(ERROR_UNIT,'(A)') TRIM(YTEXT)
FLUSH(ERROR_UNIT)
CALL C_EXIT_IMMEDIATELY(1_C_INT)
!
END SUBROUTINE ABOR1_SFX
