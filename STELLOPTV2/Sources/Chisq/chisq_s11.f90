!-----------------------------------------------------------------------
!     Subroutine:    chisq_s11
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          08/02/2019
!     Description:   Calculates the S11 susceptance matrix target.
!                    For details see
!                       Strand and Houlberg, Phys. Plasmas 2001, Vol 8,
!                              No. 6, p. 2782
!-----------------------------------------------------------------------
      SUBROUTINE chisq_s11(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_utils
      USE EZspline_obj
      USE EZspline
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      REAL(rprec), INTENT(in)    ::  target(nprof)
      REAL(rprec), INTENT(in)    ::  sigma(nprof)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!        lreset_s    Gets set to true if using R,PHI,Z Specification
!        ik          Dummy index
!        ier         Error Flag
!        iota_val    Holds profile evaulation
!-----------------------------------------------------------------------
      LOGICAL ::  lreset_s = .true.
      INTEGER ::  ik, ier
      REAL(rprec) :: iota_val, s11, s12, s21, s22
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0 ) RETURN
      ik = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'S11 ',ik,4
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  S11  S'
      IF (niter >= 0) THEN
         ! GET s if necessary
         DO ik = 1, nprof
            IF (sigma(ik) >= bigno) CYCLE
            IF (s_s11(ik) <= 1.0 .and. s_s11(ik) >= 0.0) THEN
               ier = 0
               CALL get_equil_sus(s_vaciota(ik),s11,s12,s21,s22,ier)
            ELSE
               s11 = 0.0
            END IF
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)
            sigmas(mtargets)  = sigma(ik)
            vals(mtargets)    = s11
            IF (iflag == 1) WRITE(iunit_out,'(4ES22.12E3)') target(ik),sigma(ik),s11,s_s11(ik)
         END DO
      ELSE
         DO ik = 1, nprof
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets) = jtarget_s11
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_s11
