      SUBROUTINE chisq_jstar (sigma, bmax_b, bmin_b,
     1          bmod_b, num, nopt, mskip, sqrt_nsurf)
      USE stel_kinds
      USE chisq_mod
      USE boozer_params, ONLY: nu2_b, nv_boz
      USE vparams, ONLY: twopi
      USE optim, ONLY: bigno
      USE optim_params, ONLY: NumJstar
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nopt, mskip
      INTEGER :: num
      REAL(rprec), INTENT(in) :: bmax_b(*), bmin_b(*), bmod_b(*),
     1    sigma(1,*), sqrt_nsurf
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p5 = 0.5_dp, zero = 0, one = 1
      REAL(rprec), PARAMETER :: max_variation = 2.e-2_dp
      REAL(rprec), PARAMETER :: epl = 0.05_dp, epu = 0.05_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iep, nJstar, n
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: trapJS
      REAL(rprec) :: bmin_global, bmax_global, epsmu, sj, avg_Jstar
      LOGICAL, DIMENSION(:), ALLOCATABLE  :: ljstar
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL j_star
C-----------------------------------------------
      IF (ALL(ABS(sigma(1,1:NumJstar)) .ge. bigno)) RETURN

!        CODE GENERATED BY D. SPONG, 6/97
!        COMPUTE J_star AT NUMJSTAR VALUES OF ep/mu RANGING FROM SLIGHTLY ABOVE
!        THE TRAPPED-PASSING BOUNDARY TO SLIGHTLY BELOW THE
!        DEEPLY TRAPPED-FORBIDDEN BOUNDARY.  THE PARAMETERS epl AND epu
!        DETERMINE DISTANCE TO THESE BOUNDARIES.

      IF (nopt .gt. 0) THEN
         ALLOCATE (ljstar(nu2_b), trapJS(nu2_b))

         bmin_global = MINVAL(bmin_b(:nu2_b))
         bmax_global = MAXVAL(bmax_b(:nu2_b))

         DO iep = 1,NumJstar

            IF (ABS(sigma(1,iep)) .ge. bigno) CYCLE

            sj = REAL(iep - 1,rprec)/(NumJstar - 1)
            epsmu = bmin_global*(one + epl) + sj*(
     1         bmax_global*(one - epu) - bmin_global*(one + epl))
            CALL j_star (bmod_b, bmin_b, bmax_b, epsmu, trapJs,
     1          nv_boz, nu2_b)

            ljstar = trapJs(:nu2_b) .gt. zero
            nJstar = count(ljstar)
            avg_Jstar = SUM(trapJs, mask=ljstar)
            IF (nJstar .gt. 0) avg_Jstar = avg_Jstar/nJstar
            WHERE (.not.ljstar) trapJs = avg_Jstar
!
!          Target all non-zero Jstars to (d Jstar/du) = 0
!
           DO n = 1, nu2_b, mskip
              num = num+1
              index_array(num) = ivar_jstar
              wegt(num) = max_variation * avg_Jstar
     1                   * sqrt_nsurf * sigma(1,iep)
              chisq_target(num) = avg_Jstar
              chisq_match(num) = trapJs(n)
            END DO
         END DO

         DEALLOCATE (ljstar, trapJS)

      ELSE
         DO iep = 1,NumJstar
            IF (ABS(sigma(1,iep)) .ge. bigno) CYCLE
            DO n = 1, nu2_b, mskip
               num = num + 1
               IF (nopt .eq. -2) chisq_descript(num) = 
     1                           descript(ivar_jstar)
            END DO
         END DO
      END IF

      END SUBROUTINE chisq_jstar