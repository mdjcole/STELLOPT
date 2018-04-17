!-----------------------------------------------------------------------
!     Subroutine:    stellopt_regcoil_chi2_b
!     Authors:       J.C.Schmitt (Auburn/PPPL) jcschmitt@auburn.edu
!     Date:          2017-2018
!     Description:   This subroutine calls the coil regularization code
!                    REGCOIL in 'target sqrt(<K^2>)' mode
!                    
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_regcoil_chi2_b(lscreen, iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars
      USE equil_utils
      use vparams, only: my_mpol => mpol_rcws, my_ntor => ntor_rcws

!DEC$ IF DEFINED (REGCOIL)
      ! REGCOIL files
      USE regcoil_variables
      USE regcoil_input_mod, ONLY: regcoil_write_input
      USE regcoil_validate_input
      USE regcoil_compute_lambda
      USE regcoil_init_plasma
      USE regcoil_init_coil_surface
      USE regcoil_read_bnorm
      USE regcoil_build_matrices
      USE regcoil_auto_regularization_solve
      USE regcoil_write_output
!DEC$ ENDIF

!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(inout) :: iflag
      LOGICAL, INTENT(inout)        :: lscreen

!-----------------------------------------------------------------------
!     Local Variables
!        iunit       File unit number
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Local Variables
!        istat         Error status
!        iunit         File unit number
      ! FOR REGCOIL
      INTEGER :: istat, iunit, m, n, ii, nummodes1, nummodes2

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
!      IF (iflag < 0) RETURN
      !lscreen = .true.
      IF (lscreen) then
         WRITE(6,'(a)') ' -------------  REGCOIL CALCULATION  ---------'
      ENDIF
      ! WRITE(6,'(a,a)') '<---- proc_string=', proc_string
      ! WRITE(6,'(a,i4.2)') ' -------------  REGCOIL: iflag=', iflag
!DEC$ IF DEFINED (REGCOIL)

      ! IF (lscreen) WRITE(6,'(a,a)') '<---- proc_string=', proc_string
      wout_filename = 'wout_'//TRIM(proc_string)//'.nc'
      ! STELLOPT (via lmdif->stellopt_fcn or similar) will modifiy the value of
      ! regcoil_winding_surface_separation, current_density, and/or the
      ! boundary coefficients. Here, the value of the REGCOIL variables
      ! are loaded with the new values, the correct mode of operation is
      ! determiend, control variables are set, and the regcoil-functions
      ! are called
      separation = regcoil_winding_surface_separation
      current_density_target = regcoil_current_density
     
      ! Loop over all of the spectral components of the winding surface
      ! and update the rc_*_stellopt  
      ! write(6,'(a)') '<----Looping over m and n'

      IF ((ANY(lregcoil_rcws_rbound_s_opt)) .or. (ANY(lregcoil_rcws_rbound_c_opt)) .or. &
          (ANY(lregcoil_rcws_zbound_s_opt)) .or. (ANY(lregcoil_rcws_zbound_c_opt)) ) THEN 
         nummodes1 = 0
         DO m = -my_mpol, my_mpol
             DO n = -my_ntor, my_ntor
                IF ( (regcoil_rcws_rbound_c(m,n) .ne. 0) .or. &
                     (regcoil_rcws_rbound_s(m,n) .ne. 0) .or. &
                     (regcoil_rcws_zbound_c(m,n) .ne. 0) .or. &
                     (regcoil_rcws_zbound_s(m,n) .ne. 0) ) THEN
                   nummodes1 = nummodes1 + 1
                END IF
             END do
         END do
         CALL safe_open(iunit, istat, TRIM('regcoil_nescout.'// &
                   TRIM(proc_string)), 'replace', 'formatted')
         !write(6,'(a)'), '<----JCSwrite_output'
         write(iunit,*), "Number of fourier modes in table"
         write(iunit,*), nummodes1
         write(iunit,*), "Table of fourier coefficients"
         write(iunit,*), "m,n,crc2,czs2,crs2,czc2"
         DO m = -my_mpol, my_mpol
             DO n = -my_ntor, my_ntor
                if ( (regcoil_rcws_rbound_c(m,n) .ne. 0) .or. &
                     (regcoil_rcws_rbound_s(m,n) .ne. 0) .or. &
                     (regcoil_rcws_zbound_c(m,n) .ne. 0) .or. &
                     (regcoil_rcws_zbound_s(m,n) .ne. 0) ) THEN
                   rc_rmnc_stellopt(m,n) = regcoil_rcws_rbound_c(m,n)
                   rc_rmns_stellopt(m,n) = regcoil_rcws_rbound_s(m,n)
                   rc_zmnc_stellopt(m,n) = regcoil_rcws_zbound_c(m,n)
                   rc_zmns_stellopt(m,n) = regcoil_rcws_zbound_s(m,n)
                   ! These are written in the same order as in a NESCIN
                   ! file: M N RC ZS RS ZC
                   write(iunit,*), m, n, &
                                   rc_rmnc_stellopt(m,n), rc_zmns_stellopt(m,n), &
                                   rc_rmns_stellopt(m,n), rc_zmnc_stellopt(m,n)
                END IF
              END DO
          END DO
          CLOSE(iunit)
      END IF

      ! regcoil will overwrite nlambda each time - need to restore it to
      ! the original value here
      nlambda = regcoil_nlambda

      ! This should be *almost* a duplicate of the main code from
      ! regcoil.f90
      ! write(6,'(a)') '<----Validate'

      ! check to make sure this doesn't muck up the winding surface
      call validate_input()
      ! check to make sure this doesn't muck up the winding surface
      ! write(6,'(a)') '<----Compute lambda'
      call compute_lambda(lscreen)

      ! Define the position vector and normal vector at each grid point for
      ! the surfaces:
      ! write(6,'(a)') '<----init_plasma'
      ! check to make sure this doesn't muck up the winding surface
      call init_plasma(lscreen)
      ! write(6,'(a)') '<----init coil surfs'
      IF ( (lregcoil_winding_surface_separation_opt) .and. &
           ((ANY(lregcoil_rcws_rbound_s_opt)) .or. (ANY(lregcoil_rcws_rbound_c_opt)) .or. &
            (ANY(lregcoil_rcws_zbound_s_opt)) .or. (ANY(lregcoil_rcws_zbound_c_opt))) ) THEN 
        write(6,'(a)') 'K====-----<<<<REGCOIL ERROR: Do not optimize both separation AND Fourier series simultaneously'
      END IF

      IF (lregcoil_winding_surface_separation_opt) then 
         call init_coil_surface(lscreen)
      END IF

      IF ((ANY(lregcoil_rcws_rbound_s_opt)) .or. (ANY(lregcoil_rcws_rbound_c_opt)) .or. &
          (ANY(lregcoil_rcws_zbound_s_opt)) .or. (ANY(lregcoil_rcws_zbound_c_opt)) ) THEN 
         !write(6,'(a)') '<----regcoil initupdate_nescin_coil_surface'
         call regcoil_initupdate_nescin_coil_surface(lscreen)
      END IF

      ! Initialize some of the vectors and matrices needed:
      !write(6,'(a)') '<----read bnorm'
      call read_bnorm(lscreen)
      ! write(6,'(a)') '<----build matrices'
      call build_matrices(lscreen)

      ! JCS: I disabled all options except for #5 (for now)
      ! As REGCOIL development continues, future cases can 
      ! be handled with case statements here.
      ! write(6,'(a)') '<----select a case'
      select case (general_option)
      !case (1)
      !   call solve()
      !case (2)
      !   call compute_diagnostics_for_nescout_potential()
      !case (3)
      !   call svd_scan()
      !case (4)
      !   call auto_regularization_solve()
      case (5)
         ! write(6,'(a)') '<----auto_reg solve'
         call auto_regularization_solve(lscreen)
         ! Now, the value we want should be in the variable
         ! 'chi2_B_target'. Normal termination of regcoil returns the
         ! achieved chi2_B (miniumum). If there is an 'error' (too high
         ! or too low of current), the chi2_B will contain the chi2_B
         ! that was achieved with infinite regularization 
         ! (spaced apart, straight-ish) coils
      case default
         print *,"Invalid general_option:",general_option
         stop
      END select


      !=================DEBUG SECTION==========================
      !  This section is for debugging purposes. If
      !  uncommented, the regcoil input files should be written to the
      !  working job directory and an output statement(s) will be displaed
      !  to the screen, regardless of the value of 'lscreen'
      !  - JCS
      !IF ((ANY(lregcoil_rcws_rbound_s_opt)) .or. (ANY(lregcoil_rcws_rbound_c_opt)) .or. &
      !    (ANY(lregcoil_rcws_zbound_s_opt)) .or. (ANY(lregcoil_rcws_zbound_c_opt)) ) THEN 
      ! write(6,'(a)') '<----REGCOIL DEBUG safe_open'
      !CALL safe_open(iunit, istat, TRIM('regcoil_nescout.'// &
      !          TRIM(proc_string)), 'replace', 'formatted')
      !write(6,'(a)'), '<----JCSwrite_output'
      !write(iunit,*), "Number of fourier modes in table"
      !write(iunit,*), nummodes1
      !write(iunit,*), "Table of fourier coefficients"
      !write(iunit,*), "m,n,crc2,czs2,crs2,czc2"

      !ii = 0
      !nummodes2 = 0
      !DO m = -mpol_rcws, mpol_rcws
      !    DO n = -ntor_rcws, ntor_rcws
      !       if ( (regcoil_rcws_rbound_c(m,n) .ne. 0) .or. &
      !            (regcoil_rcws_rbound_s(m,n) .ne. 0) .or. &
      !            (regcoil_rcws_zbound_c(m,n) .ne. 0) .or. &
      !            (regcoil_rcws_zbound_s(m,n) .ne. 0) ) THEN
      !         ii = ii+1
      !         rc_xm_stellopt(ii) = m
      !         rc_xn_stellopt(ii) = n
      !         rc_rmnc_stellopt(ii) = regcoil_rcws_rbound_c(m,n)
      !         rc_rmns_stellopt(ii) = regcoil_rcws_rbound_s(m,n)
      !         rc_zmnc_stellopt(ii) = regcoil_rcws_zbound_c(m,n)
      !         rc_zmns_stellopt(ii) = regcoil_rcws_zbound_s(m,n)
      !         if ( (rc_rmnc_stellopt(ii) .ne. 0) .or.  (rc_rmns_stellopt(ii) .ne. 0) .or. &
      !              (rc_zmnc_stellopt(ii) .ne. 0) .or.  (rc_zmns_stellopt(ii) .ne. 0) ) THEN
      !            nummodes2 = nummodes2 + 1
      !            write(iunit,*), m, n, &
      !                   rc_rmnc_stellopt(m,n), rc_rmns_stellopt(m,n), &
      !                   rc_zmnc_stellopt(m,n), rc_zmns_stellopt(m,n)
      !         END if
      !     END DO
      !END DO
      !if (nummodes1 .ne. nummodes2) then
      !  write(6,'(a)') '<----Hmmm....different number of modes???'
      !END if
      !END if

      !print *, chi2_B_target
      !print *,"REGCOIL complete. Total time=",totalTime,"sec."
       
      !=================END OF DEBUG SECTION==========================


      IF (lscreen) WRITE(6,'(a)') ' ---------------------------  REGCOIL CALCULATION DONE  ---------------------'
!DEC$ ENDIF
      RETURN

!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_regcoil_chi2_b