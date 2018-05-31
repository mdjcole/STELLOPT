!-----------------------------------------------------------------------
!     Subroutine:    stellopt_prof_to_vmec(iflag)
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/20/2012
!     Description:   This subroutine takes various profile arrays and
!                    constructs the am_aux and ac_aux arrays.  It begins
!                    by redefining the s vector over 25 points
!                    equidistant in toroidal flux.  Then it splines
!                    each quantity to the new grid before recomposing
!                    the quantities onto the new s mesh.
!                    It has also been modified to setup the AH and AT
!                    arrays in order to properly define toroidal
!                    rotation.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_prof_to_vmec(file_str,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stellopt_vars
      USE stellopt_targets
      USE equil_utils
      USE vmec_input, ONLY: am_aux_s, am_aux_f, ac_aux_s, ac_aux_f, &
                            ah_aux_s, ah_aux_f, at_aux_s, at_aux_f, &
                            pmass_type, pcurr_type, ph_type, pt_type
      USE safe_open_mod
      IMPLICIT NONE
      
!-----------------------------------------------------------------------
!     Input Variables
!        iflag   Error flag
!----------------------------------------------------------------------
!      CHARACTER(256),INTENT(in) :: file_str   !Changed due to bounds check issue in init_stellopt
      CHARACTER(*),INTENT(in) :: file_str
      INTEGER, INTENT(inout)   :: iflag
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      LOGICAL ::  lnew_am, lnew_ac, lnew_at, lnew_ah, lno_file
      INTEGER ::  dex_ne, dex_te, dex_ti, dex_zeff, &
                  dex_beamj, dex_bootj, dex_am, dex_ac, &
                  dex_ah, dex_at
      INTEGER ::  ier,ik, iunit, dex
      REAL(rprec), PARAMETER :: ec  = 1.60217653D-19
      REAL(rprec) :: emis_xics_temp
      REAL(rprec), ALLOCATABLE :: ne_temp(:), zeff_temp(:), te_temp(:), ti_temp(:),&
                                  th_temp(:), beamj_temp(:), bootj_temp(:), dump_temp(:)
      
      INTEGER, PARAMETER :: new_prof=50
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      lno_file = .false.
      IF (iflag==-327) lno_file = .true.
      iflag = 0
      ! Setup the internal STELLOPT arrays
      IF (EZspline_allocated(phi_spl)) CALL EZspline_free(phi_spl,iflag)
      IF (EZspline_allocated(ne_spl)) CALL EZspline_free(ne_spl,iflag)
      IF (EZspline_allocated(te_spl)) CALL EZspline_free(te_spl,iflag)
      IF (EZspline_allocated(ti_spl)) CALL EZspline_free(ti_spl,iflag)
      IF (EZspline_allocated(th_spl)) CALL EZspline_free(th_spl,iflag)
      IF (EZspline_allocated(nustar_spl)) CALL EZspline_free(nustar_spl,iflag)
      IF (EZspline_allocated(emis_xics_spl)) CALL EZspline_free(emis_xics_spl,iflag)
      IF (EZspline_allocated(zeff_spl)) CALL EZspline_free(zeff_spl,iflag)
      IF (EZspline_allocated(omega_spl)) CALL EZspline_free(omega_spl,iflag)
      dex = MINLOC(phi_aux_s(2:),DIM=1)
      IF (dex > 4) THEN
         CALL EZspline_init(phi_spl,dex,bcs0,iflag)
         phi_spl%x1 = phi_aux_s
         phi_spl%isHermite = 1
         CALL EZspline_setup(phi_spl,phi_aux_f,ier)
      END IF
      dex = MINLOC(ne_aux_s(2:),DIM=1)
      IF (dex > 4) THEN
         CALL EZspline_init(ne_spl,dex,bcs0,iflag)
         ne_spl%x1 = ne_aux_s(1:dex)
         ne_spl%isHermite = 1
         CALL EZspline_setup(ne_spl,ne_aux_f,ier)
      END IF
      dex = MINLOC(te_aux_s(2:),DIM=1)
      IF (dex > 4) THEN
         CALL EZspline_init(te_spl,dex,bcs0,iflag)
         te_spl%x1 = te_aux_s(1:dex)
         te_spl%isHermite = 1
         CALL EZspline_setup(te_spl,te_aux_f,ier)
      END IF
      dex = MINLOC(ti_aux_s(2:),DIM=1)
      IF (dex > 4) THEN
         CALL EZspline_init(ti_spl,dex,bcs0,iflag)
         ti_spl%x1 = ti_aux_s(1:dex)
         ti_spl%isHermite = 1
         CALL EZspline_setup(ti_spl,ti_aux_f,ier)
      END IF
      dex = MINLOC(th_aux_s(2:),DIM=1)
      IF (dex > 4) THEN
         CALL EZspline_init(th_spl,dex,bcs0,iflag)
         th_spl%x1 = th_aux_s(1:dex)
         th_spl%isHermite = 1
         CALL EZspline_setup(th_spl,phi_aux_f,ier)
      END IF
      dex = MINLOC(nustar_s(2:),DIM=1)
      IF (dex > 4) THEN
         CALL EZspline_init(nustar_spl,dex,bcs0,iflag)
         nustar_spl%x1 = nustar_s(1:dex)
         nustar_spl%isHermite = 1
         CALL EZspline_setup(nustar_spl,nustar_f,ier)
      END IF
      dex = MINLOC(zeff_aux_s(2:),DIM=1)
      IF (dex > 4) THEN
         CALL EZspline_init(zeff_spl,dex,bcs0,iflag)
         zeff_spl%x1 = zeff_aux_s(1:dex)
         zeff_spl%isHermite = 1
         CALL EZspline_setup(zeff_spl,zeff_aux_f,ier)
      END IF
      dex = MINLOC(ah_aux_s(2:),DIM=1)
      IF (dex > 4) THEN
         CALL EZspline_init(omega_spl,dex,bcs0,iflag)
         omega_spl%x1 = ah_aux_s(1:dex)
         omega_spl%isHermite = 1
         CALL EZspline_setup(omega_spl,ah_aux_f,ier)
      END IF
      dex = MINLOC(emis_xics_s(2:),DIM=1)
      IF (dex > 4) THEN
         CALL EZspline_init(emis_xics_spl,dex,bcs0,iflag)
         emis_xics_spl%x1 = emis_xics_s(1:dex)
         emis_xics_spl%isHermite = 1
         CALL EZspline_setup(emis_xics_spl,emis_xics_f,ier)
      END IF
      
      
      ! First get indexes of component arrays
      dex_ne    = MINLOC(ne_aux_s(2:),DIM=1)
      dex_te    = MINLOC(te_aux_s(2:),DIM=1)
      dex_ti    = MINLOC(ti_aux_s(2:),DIM=1)
      dex_zeff  = MINLOC(zeff_aux_s(2:),DIM=1)
      dex_beamj = MINLOC(beamj_aux_s(2:),DIM=1)
      dex_bootj = MINLOC(bootj_aux_s(2:),DIM=1)
      dex_ah    = MINLOC(ah_aux_s(2:),DIM=1)
      dex_at    = MINLOC(at_aux_s(2:),DIM=1)
      ! Now check if we need new am and ac arrays
      lnew_am = .false.; lnew_ac = .false.
      lnew_ah = .false.; lnew_at = .false.
      IF (ANY(ne_opt /= 0.0_rprec)) lnew_am = .true.
      IF (ANY(te_opt /= 0.0_rprec)) lnew_am = .true.
      IF (ANY(ti_opt /= 0.0_rprec)) lnew_am = .true.
      IF (ANY(th_opt /= 0.0_rprec)) lnew_am = .true.
      IF (ANY(bootj_aux_f /= 0.0_rprec)) lnew_ac = .true.
      IF (ANY(beamj_aux_f /= 0.0_rprec)) lnew_ac = .true.
      IF (dex_te > 3) lnew_am = .true.
      IF (dex_ti > 3) lnew_am = .true.
      !IF (dex_beamj > 1) lnew_ac = .true.
      !IF (dex_bootj > 1) lnew_ac = .true.
      IF (TRIM(equil_type)=='flow' .or. TRIM(equil_type)=='satire') THEN
         lnew_ah = .true.
         lnew_at = .true.
      END IF
      ! Now create the new arrays
      IF (lnew_am) THEN
         iunit = 66
         IF (.not. lno_file) CALL safe_open(iunit,iflag,TRIM('tprof.'//TRIM(file_str)),'unknown','formatted')
         !WRITE(iunit,*) am_aux_s
         am_aux_s(:) = -1
         DO ik = 1, new_prof
            am_aux_s(ik) = REAL(ik-1)/REAL(new_prof-1)
         END DO
         am_aux_s(1) = 0.0_rprec
         am_aux_s(new_prof) = 1.0_rprec
         dex_am = MINLOC(am_aux_s(2:),DIM=1)
         !WRITE(iunit,*) am_aux_s
         ! NE
         ALLOCATE(ne_temp(1:dex_am))
         ne_temp = 1.0_rprec
         DO ik = 1, dex_am
            ier = 0
            CALL get_equil_ne(am_aux_s(ik),TRIM(ne_type),ne_temp(ik),ier)
            IF (ier /= 0) ne_temp(ik) = 1.0D18
         END DO
         ier = 0
         ! TE
         ALLOCATE(te_temp(1:dex_am))
         te_temp = 0.0_rprec
         DO ik = 1, dex_am
            ier = 0
            CALL get_equil_te(am_aux_s(ik),TRIM(te_type),te_temp(ik),ier)
         END DO
         ! TI
         ALLOCATE(ti_temp(1:dex_am))
         ti_temp = 0.0_rprec
         DO ik = 1, dex_am
            ier = 0
            CALL get_equil_ti(am_aux_s(ik),TRIM(ti_type),ti_temp(ik),ier)
         END DO
         ! ZEFF
         ALLOCATE(zeff_temp(1:dex_am))
         zeff_temp = 1.0_rprec
         DO ik = 1, dex_am
            ier = 0
            CALL get_equil_zeff(am_aux_s(ik),TRIM(zeff_type),zeff_temp(ik),ier)
         END DO
         pmass_type = 'akima_spline'
         ! Create p
         DO ik = 1, dex_am
            am_aux_f(ik) = ne_temp(ik) * ec * (te_temp(ik) + ti_temp(ik)/zeff_temp(ik))
            IF (am_aux_f(ik) < 0.0_rprec) am_aux_f(ik) = 0.0_rprec
         END DO
         !WRITE(iunit,*) am_aux_s
         !WRITE(iunit,*) am_aux_f
         ! Now make table
         IF (.not. lno_file) THEN
            WRITE(iunit,*) 's     ne     te     ti     zeff     p'
            DO ik = 1, dex_am
               WRITE(iunit,'(7es12.4)') am_aux_s(ik),ne_temp(ik),te_temp(ik),ti_temp(ik),zeff_temp(ik),am_aux_f(ik)
            END DO
            CLOSE(iunit)
         END IF
      END IF
      !  SINCE ROTATION does not depend on anything I don't think we need to do
      !  this.  -SAL 12/2/13
      !IF (lnew_ah) THEN
      !   iunit = 66
      !   CALL safe_open(iunit,iflag,TRIM('rotprof.'//TRIM(proc_string)),'unknown','formatted')
      !   ah_aux_s(:) = -1
      !   DO ik = 1, new_prof
      !      ah_aux_s(ik) = REAL(ik-1)/REAL(new_prof-1)
      !   END DO
      !   ah_aux_s(1) = 0.0
      !   ah_aux_s(new_prof) = 1.0
      !   dex_ah = MINLOC(ah_aux_s(2:),DIM=1)
      !   ph_type = 'akima_spline'
      !   DO ik = 1, dex_ah
      !      ier = 0
      !      CALL get_equil_omega(ah_aux_s(ik),ah_aux_f(ik),ier)
      !   END DO
      !   ! Now make table
      !   WRITE(iunit,*) 's     rot'
      !   DO ik = 1, dex_am
      !      WRITE(iunit,'(7es12.4)') ah_aux_s(ik),ah_aux_f(ik)
      !   END DO
      !   CLOSE(iunit)
      !END IF
      IF (lnew_at) THEN
         !iunit = 66
         !CALL safe_open(iunit,iflag,TRIM('rotprof.'//TRIM(proc_string)),'unknown','formatted')
         at_aux_s(:) = -1
         DO ik = 1, new_prof
            at_aux_s(ik) = REAL(ik-1)/REAL(new_prof-1)
         END DO
         at_aux_s(1) = 0.0_rprec
         at_aux_s(new_prof) = 1.0_rprec
         dex_at = MINLOC(at_aux_s(2:),DIM=1)
         pt_type = 'akima_spline'
         DO ik = 1, dex_at
            ier = 0
            CALL get_equil_ti(at_aux_s(ik),TRIM(ti_type),at_aux_f(ik),ier)
         END DO
         ! Now make table
         !WRITE(iunit,*) 's     rot'
         !DO ik = 1, dex_am
         !   WRITE(iunit,'(7es12.4)') ah_aux_s(ik),ah_aux_f(ik)
         !END DO
         !CLOSE(iunit)
      END IF
         
      IF (lnew_ac) THEN
         iunit = 66
         IF (.not. lno_file) CALL safe_open(iunit,iflag,TRIM('jprof.'//TRIM(file_str)),'unknown','formatted')
         !WRITE(iunit,*) ac_aux_s
         ac_aux_s(:) = -1
         DO ik = 1, new_prof
            ac_aux_s(ik) = REAL(ik-1)/REAL(new_prof-1)
         END DO
         ac_aux_s(1) = 0.0_rprec
         ac_aux_s(new_prof) = 1.0_rprec
         dex_ac = MINLOC(ac_aux_s(2:),DIM=1)
         !WRITE(iunit,*) ac_aux_s
         ! BEAMJ
         ALLOCATE(beamj_temp(1:dex_ac))
         beamj_temp = 0.0_rprec
         DO ik = 1, dex_ac
            ier = 0
            CALL get_equil_beamj(ac_aux_s(ik),beamj_temp(ik),ier)
         END DO
         ! BEAMJ
         ALLOCATE(bootj_temp(1:dex_ac))
         bootj_temp = 0.0_rprec
         DO ik = 1, dex_ac
            ier = 0
            CALL get_equil_bootj(ac_aux_s(ik),bootj_temp(ik),ier)
         END DO
         pcurr_type = 'akima_spline_ip'
         ! Create j
         DO ik = 1, dex_ac
            ac_aux_f(ik) = beamj_temp(ik) + bootj_temp(ik)
         END DO
         !WRITE(iunit,*) ac_aux_s
         !WRITE(iunit,*) ac_aux_f
         IF (.not. lno_file) THEN
            WRITE(iunit,*) 's    j_beam     j_bs     j_tot'
            DO ik = 1, dex_ac
               WRITE(iunit,'(4es12.4)') ac_aux_s(ik),beamj_temp(ik),bootj_temp(ik),ac_aux_f(ik)
            END DO
            CLOSE(iunit)
         END IF
      END IF
      
      ! Free Arrays
      IF (ALLOCATED(ne_temp)) DEALLOCATE(ne_temp)
      IF (ALLOCATED(zeff_temp)) DEALLOCATE(zeff_temp)
      IF (ALLOCATED(te_temp)) DEALLOCATE(te_temp)
      IF (ALLOCATED(ti_temp)) DEALLOCATE(ti_temp)
      IF (ALLOCATED(bootj_temp)) DEALLOCATE(bootj_temp)
      IF (ALLOCATED(beamj_temp)) DEALLOCATE(beamj_temp)

      ! Output a diagnostic table
      ! Now make table
      IF (lnew_am) THEN
         iunit = 68
         IF (.not. lno_file) THEN
            CALL safe_open(iunit,iflag,TRIM('dprof.'//TRIM(file_str)),'unknown','formatted')
            WRITE(iunit,*) 's     emis_xics'
            DO ik = 1, dex_am
               CALL get_equil_emis_xics(am_aux_s(ik),TRIM(emis_xics_type),emis_xics_temp,ier)
               WRITE(iunit,'(2es12.4)') am_aux_s(ik),emis_xics_temp
            END DO
            CLOSE(iunit)
         END IF
      END IF
      
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_prof_to_vmec