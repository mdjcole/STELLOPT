!-----------------------------------------------------------------------
!     Module:        beams3d_diagnostics
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          01/09/2014
!     Description:   This subroutine outputs a diagnostic text file
!                    of the run.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_diagnostics
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_lines
      USE beams3d_grid, ONLY: nr, nphi, nz, B_R, B_PHI, B_Z, raxis, &
                                 zaxis, phiaxis,vp_spl_s, plasma_mass
      USE beams3d_runtime, ONLY: id_string, npoinc, t_end, lbeam, lvac,&
                                 lvmec, charge_beams, &
                                 nbeams, beam, e_beams, charge_beams, &
                                 mass_beams, lverb, p_beams, MPI_BARRIER_ERR,&
                                 MPI_BCAST_ERR,nprocs_beams,handle_err, ldepo,&
                                 MPI_REDU_ERR, pi2, weight
      USE beams3d_physics_mod, ONLY: beams3d_SFLX, beams3d_COLLISIONS
      USE safe_open_mod, ONLY: safe_open
      USE EZspline
      USE mpi_params ! MPI
      USE mpi_sharmem
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          ier          Error Flag
!          iunit        File ID
!          ndist        Number of Vll divisions for dist function
!          ns           Number of flux divisions for current calculation
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: iunit, istat, i, j, k, l, m, n, sbeam, ebeam, ninj2, ier
      REAL(rprec) :: v1,v2, ddist, ninj
      LOGICAL, ALLOCATABLE     :: partmask(:), partmask2(:,:), partmask2t(:,:)
      INTEGER, ALLOCATABLE  :: int_mask(:), int_mask2(:,:)
      INTEGER, ALLOCATABLE  :: dist_func(:,:,:)
      REAL, ALLOCATABLE     :: real_mask(:), nlost(:)
      ! For 5D dist function
      LOGICAL :: lhelp
      INTEGER :: win_rho3d, win_efact, win_ifact, win_vdist2d, win_dVol, win_rho4d
      REAL(rprec) :: s1, s2, modb, ti, vp_temp, dr, dphi, dz, dvll, dvperp, drho
      REAL, ALLOCATABLE, DIMENSION(:)     :: rdistaxis,vllaxis,vperpaxis
      REAL, DIMENSION(:), POINTER         :: dVol
      REAL, DIMENSION(:,:), POINTER       :: vdist2d
      REAL, DIMENSION(:,:,:), POINTER     :: rho3d, efact3d, ifact3d
      REAL, DIMENSION(:,:,:,:), POINTER   :: rho4d, efact4d, ifact4d
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: help3d
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: help4d
      DOUBLE PRECISION, DIMENSION(3)      :: q

      INTEGER, PARAMETER :: ndist = 100
      DOUBLE PRECISION, PARAMETER :: electron_mass = 9.10938356D-31 !m_e
      DOUBLE PRECISION, PARAMETER :: sqrt_pi       = 1.7724538509   !pi^(1/2)
      DOUBLE PRECISION, PARAMETER :: e_charge      = 1.60217662E-19 !e_c
#if defined(MPI_OPT)
      INTEGER :: numprocs_local, mylocalid, mylocalmaster, mystart
      INTEGER :: MPI_COMM_LOCAL
#endif
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      IF (lverb) WRITE(6,'(A)')  '----- BEAM DIAGNOSTICS -----'


      ! DEALLOCATE stuff we do not need
      IF (ALLOCATED(Z_lines)) DEALLOCATE(Z_lines)
      IF (ALLOCATED(moment_lines)) DEALLOCATE(moment_lines)
      IF (ALLOCATED(U_lines)) DEALLOCATE(U_lines)
      IF (ALLOCATED(B_lines)) DEALLOCATE(B_lines)

      CALL FLUSH(6)

      mystart = mystart_save
      myend = myend_save

      ! Main Allocations
      IF (ALLOCATED(shine_through)) DEALLOCATE(shine_through)
      IF (ALLOCATED(shine_port)) DEALLOCATE(shine_port)
      IF (ALLOCATED(nlost)) DEALLOCATE(nlost)
      ALLOCATE(shine_through(nbeams))
      ALLOCATE(shine_port(nbeams))
      ALLOCATE(nlost(nbeams))
      ALLOCATE(partmask(mystart:myend))
      ALLOCATE(partmask2(0:npoinc,mystart:myend))
      ALLOCATE(partmask2t(0:npoinc,mystart:myend))
      ALLOCATE(int_mask(mystart:myend))
      ALLOCATE(int_mask2(0:npoinc,mystart:myend))
      ALLOCATE(real_mask(mystart:myend))
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= 0) CALL handle_err(MPI_BARRIER_ERR, 'beams3d_follow', ierr_mpi)
#endif

      ! Setup masking arrays
      FORALL(i=0:npoinc) int_mask2(i,mystart:myend) = beam(mystart:myend)              ! Index of beams
      WHERE(      ( (R_lines(:,mystart:myend)==0) .and. (PHI_lines(:,mystart:myend)==-1) )&
             .or. (neut_lines(:,mystart:myend)) ) int_mask2(:,mystart:myend) = 0  ! Mask the neutral and lost particles
      int_mask(mystart:myend) = 0
      IF (lbeam) int_mask(mystart:myend) = 3
      int_mask(mystart:myend)  = COUNT(neut_lines(:,mystart:myend),DIM=1)-1                  ! Starting index of every charge particle
      WHERE(int_mask < 0) int_mask = 0
      FORALL(j=mystart:myend) real_mask(j) = S_lines(int_mask(j),j) ! Starting points in s

      ! Do not need R_lines or PHI_lines after this point
      IF (ALLOCATED(R_lines)) DEALLOCATE(R_lines)
      IF (ALLOCATED(PHI_lines)) DEALLOCATE(PHI_lines)
      IF (ALLOCATED(neut_lines)) DEALLOCATE(neut_lines)

      ! Calculate distribution function
      ALLOCATE(dist_func(1:nbeams,1:ndist,0:npoinc))
      dist_func = 0
      ddist = 2*partvmax/ndist
      sbeam = MINVAL(beam(mystart:myend), DIM=1)
      ebeam = MAXVAL(beam(mystart:myend), DIM=1)
      DO k = 1, ndist
         v1 = -partvmax+(k-1)*ddist
         v2 = -partvmax+(k)*ddist
         partmask2(:,mystart:myend) = ((vll_lines(:,mystart:myend).ge.v1).and.(vll_lines(:,mystart:myend).lt.v2))
         DO i = sbeam, ebeam
            dist_func(i,k,0:npoinc) = COUNT(partmask2(:,mystart:myend).and.(int_mask2(:,mystart:myend)==i),DIM=2)
         END DO
      END DO

      ! Calculate shinethrough and loss
      shine_through = 0
      DO i = 1, nbeams
         nlost(i)          =      SUM(weight(mystart:myend), MASK = (end_state(mystart:myend) == 2 .and. (beam(mystart:myend)==i)))
         shine_through(i)  = 100.*SUM(weight(mystart:myend), MASK = (end_state(mystart:myend) == 3 .and. (beam(mystart:myend)==i)))/SUM(weight,MASK=(beam==i))
         shine_port(i)     = 100.*SUM(weight(mystart:myend), MASK = (end_state(mystart:myend) == 4 .and. (beam(mystart:myend)==i)))/SUM(weight,MASK=(beam==i))
         !shine_through(i) = 100.*COUNT(end_state(mystart:myend) == 3 .and. (beam(mystart:myend)==i),DIM=1)/COUNT(beam==i)
         !shine_port(i) = 100.*COUNT(end_state(mystart:myend) == 4 .and. (beam(mystart:myend)==i),DIM=1)/COUNT(beam==i)
         !nlost(i)         = COUNT(end_state(mystart:myend) == 2 .and. beam(mystart:myend) == i)
      END DO

#if defined(MPI_OPT)
      IF (myworkid == master) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE, dist_func,     nbeams*ndist*(npoinc+1), MPI_INTEGER,          MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE, shine_through, nbeams,                  MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE, shine_port,    nbeams,                  MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE, nlost,         nbeams,                  MPI_REAL,          MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
      ELSE
         CALL MPI_REDUCE(dist_func,     dist_func,     nbeams*ndist*(npoinc+1), MPI_INTEGER,          MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(shine_through, shine_through, nbeams,                  MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(shine_port,    shine_port,    nbeams,                  MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(nlost,         nlost,         nbeams,                  MPI_REAL,          MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
      END IF
#endif

      IF (myworkid == master) THEN
         ! Open the file
         iunit = 10
         CALL safe_open(iunit,istat,'beams3d_diag_'//TRIM(id_string)//'.txt','replace','formatted')
         ! Output number of beams
         WRITE(iunit,'(A,I5)') 'BEAMLINES: ',nbeams
         ! Screen Output
         DO i = 1, nbeams
            ! Output beam information
            WRITE(iunit,'(A)') ' BEAMLINE        ENERGY                 CHARGE                 MASS'
            WRITE(iunit,'((I5,3(2X,E22.12)))') i,E_BEAMS(i),CHARGE_BEAMS(i),MASS_BEAMS(i)

            ! Output beam losses
            ninj  = SUM(weight,MASK=(beam==i))
            ninj2  = COUNT(beam==i)
            WRITE(iunit,'(A)') ' Particles Launched  Particles Lost  Lost(%)  TIME_END'
            WRITE(iunit,'(6X,I10,11X,I5,7x,F5.1,6x,E22.12)') ninj2, NINT(ninj2*nlost(i)/ninj), 100.*nlost(i)/ninj, MAXVAL(t_end)
            WRITE(iunit,'(A)') ' '
            CALL FLUSH(iunit)

            ! Screen Output
            IF (lverb) THEN
               IF (i==1) WRITE(6,'(A)')  ' BEAMLINE     ENERGY [keV]   CHARGE [e]   MASS [Mp]   Particles [#]   Lost [%]  Shinethrough [%]  Port [%]'
               WRITE(6,'(I5,3(9X,I5),8X,I8,3(8X,F5.1))') i,NINT(E_BEAMS(i)*6.24150636309E15),NINT(CHARGE_BEAMS(i)*6.24150636309E18),&
                                         NINT(MASS_BEAMS(i)*5.97863320194E26), ninj2, 100.*nlost(i)/ninj, shine_through(i), shine_port(i)
               CALL FLUSH(6)
            END IF
            ! Write Distribution Function
            WRITE(iunit,'(A)') ' Parallel Velocity Distribution'
            WRITE(iunit,'(A,100(1X,E22.12))') 'VLL',(-partvmax+(j+0.5)*ddist,j=0,ndist-1)
            WRITE(iunit,'(A)') '==========================================='
            DO j = 0, npoinc
               WRITE(iunit,'(100(1X,I12))') dist_func(i,1:ndist,j)
            END DO
            WRITE(iunit,'(A)') '==========================================='
            WRITE(iunit,'(A)') ' '
            CALL FLUSH(iunit)
         END DO
         CLOSE(iunit)
      END IF

      DEALLOCATE(dist_func)
      DEALLOCATE(int_mask2,int_mask)
      DEALLOCATE(partmask2,partmask2t)
      DEALLOCATE(partmask,real_mask)

      ! Divide up Work
      mylocalid = myworkid
      numprocs_local = 1
#if defined(MPI_OPT)
      CALL MPI_COMM_DUP( MPI_COMM_SHARMEM, MPI_COMM_LOCAL, ierr_mpi)
      CALL MPI_COMM_RANK( MPI_COMM_LOCAL, mylocalid, ierr_mpi )              ! MPI
      CALL MPI_COMM_SIZE( MPI_COMM_LOCAL, numprocs_local, ierr_mpi )          ! MPI
#endif
      mylocalmaster = master

      lhelp = .false.
      IF (myworkid == master) lhelp = .true.
#if defined(MPI_OPT)
      CALL MPI_BCAST(lhelp, 1, MPI_LOGICAL, mylocalmaster, MPI_COMM_LOCAL, ierr_mpi)
#endif

      ! These diagnostics need Vp to be defined
      IF (lvmec .and. .not.lvac .and. .not.ldepo .and. lhelp) THEN

         ! Some helpers to avoid divisions
         dr     = 1.0/h1dist
         dphi   = 1.0/h2dist
         dz     = 1.0/h3dist
         dvll   = 1.0/h4dist
         dvperp = 1.0/h5dist
         drho   = 1.0/REAL(ndistns)

         ! ALLOCATE the profile arrays
         CALL mpialloc(dVol, ndistns, myid_sharmem, 0, MPI_COMM_LOCAL, win_dVol)
         CALL mpialloc(dense_prof,  nbeams, ndistns, myid_sharmem, 0, MPI_COMM_LOCAL, win_dense)
         CALL mpialloc(j_prof,      nbeams, ndistns, myid_sharmem, 0, MPI_COMM_LOCAL, win_jprof)
         CALL mpialloc(epower_prof, nbeams, ndistns, myid_sharmem, 0, MPI_COMM_LOCAL, win_ipower)
         CALL mpialloc(ipower_prof, nbeams, ndistns, myid_sharmem, 0, MPI_COMM_LOCAL, win_epower)
         CALL mpialloc(momll_prof,  nbeams, ndistns, myid_sharmem, 0, MPI_COMM_LOCAL, win_momll)
         CALL mpialloc(pperp_prof,  nbeams, ndistns, myid_sharmem, 0, MPI_COMM_LOCAL, win_pperp)
         CALL mpialloc(vdist2d,  ndist4, ndist5, myid_sharmem, 0, MPI_COMM_LOCAL, win_vdist2d)
         CALL mpialloc(efact4d,  nbeams, ndist1, ndist2, ndist3, myid_sharmem, 0, MPI_COMM_LOCAL, win_efact)
         CALL mpialloc(ifact4d,  nbeams, ndist1, ndist2, ndist3, myid_sharmem, 0, MPI_COMM_LOCAL, win_ifact)
         CALL mpialloc(rho4d,    nbeams, ndist1, ndist2, ndist3, myid_sharmem, 0, MPI_COMM_LOCAL, win_rho4d)

         ! Create the helper axis arrays
         ALLOCATE(vllaxis(ndist4),vperpaxis(ndist5),rdistaxis(ndist1))
         ALLOCATE(help4d(nbeams, ndist1,ndist2, ndist3))
         FORALL(k = 1:ndist1) rdistaxis(k) = rmin_dist+REAL(k-0.5)*dr
         FORALL(k = 1:ndist4) vllaxis(k) = -partvmax+REAL(k-0.5)*dvll
         FORALL(k = 1:ndist5) vperpaxis(k) = REAL(k-0.5)*dvperp

         ! Master needs to calculate dVol from spline
         IF (myworkid==master) PRINT *,'CALC: dVol'
         IF (myworkid == master) THEN
            DO k = 1, ndistns
               s1 = REAL(k-0.5)*drho ! Rho
               s2 = s1*s1
               CALL EZspline_interp(Vp_spl_s, s2, vp_temp, ier)
               dVol(k) = vp_temp*2*s1/REAL(ndistns)
            END DO
!            dVol = dVol
         END IF

         ! Calculte the helpers
         IF (myworkid==master) PRINT *,'CALC: S and DVel'
         fact_crit = SQRT(2*e_charge/mass_beams(1))*(0.75*sqrt_pi*sqrt(plasma_mass/electron_mass))**(1.0/3.0) ! Wesson pg 226 5.4.9
         mymass = mass_beams(1)
         myZ    = NINT(charge_beams(1)/e_charge)
         CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL, 1, ndist1*ndist2*ndist3, mystart, myend)
         DO l = mystart, myend
            i = MOD(l-1,ndist1)+1
            j = MOD(l-1,ndist1*ndist2)
            j = FLOOR(REAL(j) / REAL(ndist1))+1
            k = CEILING(REAL(l) / REAL(ndist1*ndist2))
            q(1) = rdistaxis(i)
            q(2) = REAL(j-0.5)*dphi
            q(3) = zmin_dist+REAL(k-0.5)*dz
            CALL beams3d_SFLX(q,s1)
            rho4d(:,i,j,k) = SQRT(s1)
            CALL beams3d_COLLISIONS(q,modb,ti,s1,s2)
            efact4d(:,i,j,k) = s1 ! tau_spit_inv
            ifact4d(:,i,j,k) = s2 ! vcrit
         END DO

         ! We need to correct for particle mass and charge
         IF (myworkid == master) THEN
            DO m = 1, nbeams
               s1 = mass_beams(1)/mass_beams(m)
               s2 = NINT(charge_beams(m)/charge_beams(1))
               efact4d(m,:,:,:) = efact4d(m,:,:,:)*s2*s2*s1
               ifact4d(m,:,:,:) = ifact4d(m,:,:,:)*sqrt(s1)
            END DO
            ifact4d = ifact4d*ifact4d*ifact4d*efact4d ! vcrit_cube*tau_spit_inv
            WRITE(327,*) rho4d(1,:,1,:)
            WRITE(328,*) efact4d(1,:,1,:)
            WRITE(329,*) ifact4d(1,:,1,:)
         END IF

         ! Create the velocity helpers
         IF (myworkid==master) PRINT *,'CALC: vdist2d'
         CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL, 1, ndist4*ndist5, mystart, myend)
         DO l = mystart, myend
            i = MOD(l-1,ndist4)+1
            j = MOD(l-1,ndist4*ndist5)
            j = FLOOR(REAL(j) / REAL(ndist4))+1
            vdist2d(i,j) = sqrt(vllaxis(i)*vllaxis(i)+vperpaxis(j)*vperpaxis(j))
         END DO

         ! Now calculate the rho density profile [part*m^-3]
         IF (myworkid==master) PRINT *,'CALC: dense_prof'
         CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL, 1, ndistns, mystart, myend)
         DO l = mystart, myend ! Edges
            s1 = REAL(l-1)*drho
            s2 = REAL(l)*drho
            help4d = SUM(SUM(dist5d_prof(:,:,:,:,:,:),DIM=6),DIM=5)
            WHERE ((rho4d<=s1) .or. (rho4d>s2)) help4d = 0
            dense_prof(:,l) = SUM(SUM(SUM(help4d,DIM=4),DIM=3),DIM=2)/dVol(l)
         END DO

         ! Now calculate the Current Profile [A*m^-2]
         ! Now calculate the Parallel momentum Profile [kg*m^-2*s^-1]
         j_prof = 0
         momll_prof = 0
         pperp_prof = 0
         IF (myworkid==master) PRINT *,'CALC: j, momll'
         CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL, 1, ndistns*ndist4, mystart, myend)
         DO l = mystart, myend ! Edges
            i = MOD(l-1,ndistns)+1
            j = MOD(l-1,ndistns*ndist4)
            j = FLOOR(REAL(j) / REAL(ndistns))+1
            s1 = REAL(i-1)*drho
            s2 = REAL(i)*drho
            help4d = SUM(dist5d_prof(:,:,:,:,j,:),DIM=5)
            WHERE ((rho4d<=s1) .or. (rho4d>s2)) help4d = 0
            j_prof(:,i) = j_prof(:,i) + SUM(SUM(SUM(help4d,DIM=4),DIM=3),DIM=2)*vllaxis(j)
         END DO
         !IF (myworkid == master) momll_prof = j_prof

         ! Now calculate the Perpendicular Pressure Profile [Pa]
         IF (myworkid==master) PRINT *,'CALC: pperp'
         CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL, 1, ndistns*ndist5, mystart, myend)
         DO l = mystart, myend ! Edges
            i = MOD(l-1,ndistns)+1
            j = MOD(l-1,ndistns*ndist5)
            j = FLOOR(REAL(j) / REAL(ndistns))+1
            s1 = REAL(i-1)*drho
            s2 = REAL(i)*drho
            help4d = SUM(dist5d_prof(:,:,:,:,:,j),DIM=5)
            WHERE ((rho4d<=s1) .or. (rho4d>s2)) help4d = 0
            pperp_prof(:,i) = pperp_prof(:,i) + SUM(SUM(SUM(help4d,DIM=4),DIM=3),DIM=2)*vperpaxis(j)*vperpaxis(j)
         END DO

         ! Now calculate the Heating Profile [W*m^-3]
         epower_prof = 0
         ipower_prof = 0
         IF (myworkid==master) PRINT *,'CALC: e_power, i_power'
         CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL, 1, ndistns*ndist4*ndist5, mystart, myend)
         DO l = mystart, myend
            i = MOD(l-1,ndistns)+1
            j = MOD(l-1,ndistns*ndist4)
            j = FLOOR(REAL(j) / REAL(ndistns))+1
            k = CEILING(REAL(l) / REAL(ndistns*ndist4))
            s1 = REAL(i-1)*drho
            s2 = REAL(i)*drho
            help4d = dist5d_prof(:,:,:,:,j,k)
            WHERE ((rho4d<=s1) .or. (rho4d>s2)) help4d = 0
            epower_prof(:,i) = epower_prof(:,i) + SUM(SUM(SUM(help4d*efact4d,DIM=4),DIM=3),DIM=2)*vdist2d(j,k)*vdist2d(j,k)
            ipower_prof(:,i) = ipower_prof(:,i) + SUM(SUM(SUM(help4d*ifact4d,DIM=4),DIM=3),DIM=2)/vdist2d(j,k)
         END DO

         ! Now finish the values
         CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL, 1, ndistns, mystart, myend)
         DO m = 1, nbeams
            pperp_prof(m,mystart:myend)  = pperp_prof(m,mystart:myend)  * mass_beams(m)   / dvol(mystart:myend) !p_perp
            momll_prof(m,mystart:myend)  = j_prof(m,mystart:myend)      * mass_beams(m)   / dvol(mystart:myend) ! mvll
            j_prof(m,mystart:myend)      = j_prof(m,mystart:myend)      * charge_beams(m) / dvol(mystart:myend)
            epower_prof(m,mystart:myend) = epower_prof(m,mystart:myend) * mass_beams(m)   / dvol(mystart:myend)
            ipower_prof(m,mystart:myend) = ipower_prof(m,mystart:myend) * mass_beams(m)   / dvol(mystart:myend)
         END DO

         IF (myworkid==master) PRINT *,'CALC: dist5d_norm'
         ! First normalize the 5D phase space density by dVolume
         !CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL, 1, ndist1, mystart, myend)
         !DO i = mystart, myend
         !   vp_temp = 1.0/(rdistaxis(i)*dr*dphi*dz) ! 1./R*dr*dphi*dz
         !   dist5d_prof(:,i,:,:,:,:) = dist5d_prof(:,i,:,:,:,:)*vp_temp
         !END DO

         ! Normalize to velocity space volume element
         !CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL, 1, ndist5, mystart, myend)
         !DO k = mystart, myend ! VPERP
         !   vp_temp = 1.0/(vperpaxis(k)*pi2*dvll*dvperp)
         !   dist5d_prof(:,:,:,:,:,k) = dist5d_prof(:,:,:,:,:,k)*vp_temp
         !END DO

         IF (myworkid==master) PRINT *,'CALC: DONE'
         ! DEALLOCATIONS
         CALL mpidealloc(rho4d,win_rho4d)
         CALL mpidealloc(efact4d,win_efact)
         CALL mpidealloc(ifact4d,win_ifact)
         CALL mpidealloc(vdist2d,win_vdist2d)
         CALL mpidealloc(dVol,win_dVol)
         DEALLOCATE(rdistaxis, vperpaxis, vllaxis)
         DEALLOCATE(help4d)

      END IF

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'beams3d_diagnostics',ierr_mpi)
#endif

      CALL beams3d_write('DIAG')


!-----------------------------------------------------------------------
!     End Subroutine
!----------------------------------------------------------------------- 
      END SUBROUTINE beams3d_diagnostics
