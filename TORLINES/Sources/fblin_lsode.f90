!-----------------------------------------------------------------------
!     Function:      fblin_lsode
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/16/2012
!     Description:   This subroutine calculates the RHS of the ODE for 
!                    field line following (LSODE).
!-----------------------------------------------------------------------
      SUBROUTINE fblin_lsode(neq,phi,q,qdot)
!-----------------------------------------------------------------------
!     Libraries None
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Input Variables
!          phi        phi angle
!          q          (q(1),q(2)) = (R,Z)
!          qdot       dq/dt
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER          :: neq
      DOUBLE PRECISION :: phi, q(2), qdot(2)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      CALL fblin(phi,q,qdot)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fblin_lsode
