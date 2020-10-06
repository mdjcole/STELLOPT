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
                                 zaxis, phiaxis,vp_spl_s
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
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          ier          Error Flag
!          iunit        File ID
!          ndist        Number of Vll divisions for dist function
!          ns           Number of flux divisions for current calculation
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: iunit, istat, i, j, k, l, m, n, sbeam, ebeam, ninj2
      REAL(rprec) :: v1,v2, ddist, dvll, dvperp, ninj
      LOGICAL, ALLOCATABLE     :: partmask(:), partmask2(:,:), partmask2t(:,:)
      INTEGER, ALLOCATABLE  :: int_mask(:), int_mask2(:,:)
      INTEGER, ALLOCATABLE  :: dist_func(:,:,:)
      REAL, ALLOCATABLE     :: real_mask(:), nlost(:)
      ! For 5D dist function
      REAL(rprec) :: s1, s2, modb, ti, vp_temp
      REAL, ALLOCATABLE, DIMENSION(:)     :: rdistaxis,vllaxis,vperpaxis
      REAL, ALLOCATABLE, DIMENSION(:,:)   :: vdist2d
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: rho3d,help3d,efact3d,ifact3d
      DOUBLE PRECISION, DIMENSION(3)      :: q

      INTEGER, PARAMETER :: ndist = 100
#if defined(MPI_OPT)
      INTEGER :: mystart, mypace
      REAL(rprec), ALLOCATABLE :: buffer_mast(:,:), buffer_slav(:,:)
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

      ! These diagnostics need Vp to be defined
      IF (lvmec .and. .not.lvac .and. .not.ldepo .and. myworkid == master) THEN
         ! ALLOCATE the profile arrays
         ALLOCATE(dense_prof(nbeams,ndistns), j_prof(nbeams,ndistns), &
            epower_prof(nbeams,ndistns), ipower_prof(nbeams,ndistns), &
            momll_prof(nbeams,ndistns), pperp_prof(nbeams,ndistns))

         ! Create the helper axis arrays
         ALLOCATE(vllaxis(ndist4),vperpaxis(ndist5),rdistaxis(ndist1))
         FORALL(k = 1:ndist1) rdistaxis(k) = rmin_dist+REAL(k-0.5)/h1dist
         FORALL(k = 1:ndist4) vllaxis(k) = -partvmax+REAL(k-0.5)/h4dist
         FORALL(k = 1:ndist5) vperpaxis(k) = REAL(k-0.5)/h5dist

         ! Create the helper rho, ifact, and efact array
         ALLOCATE(rho3d(ndist1,ndist2,ndist3), &
            efact3d(ndist1,ndist2,ndist3), ifact3d(ndist1,ndist2,ndist3))
         DO l = 1, ndist1*ndist2*ndist3
            i = MOD(l-1,ndist1)+1
            j = MOD(l-1,ndist1*ndist2)
            j = FLOOR(REAL(j) / REAL(ndist1))+1
            k = CEILING(REAL(l) / REAL(ndist1*ndist2))
            q(1) = rdistaxis(i)
            q(2) = REAL(j-0.5)/h2dist
            q(3) = zmin_dist+REAL(k-0.5)/h3dist
            CALL beams3d_SFLX(q,s1)
            rho3d(i,j,k) = SQRT(s1)
            CALL beams3d_COLLISIONS(q,modb,ti,s1,s2)
            efact3d(i,j,k) = s1
            ifact3d(i,j,k) = s2*s2*s2*s1 ! vcrit_cube*tau_spit_inv
         END DO

         ! Create Velocity helper array
         ALLOCATE(vdist2d(ndist4,ndist5))
         DO l = 1, ndist4*ndist5
            i = MOD(l-1,ndist4)+1
            j = MOD(l-1,ndist4*ndist5)
            j = FLOOR(REAL(j) / REAL(ndist4))+1
            vdist2d = sqrt(vllaxis(i)*vllaxis(i)+vperpaxis(j)*vperpaxis(j))
         END DO

         ! First normalize the 5D phase space density by dVolume
         DO i = 1, ndist1
            vp_temp = (h1dist*h2dist*h3dist)/rdistaxis(i) ! 1./R*dr*dphi*dz
            dist5d_prof(:,i,:,:,:,:) = dist5d_prof(:,i,:,:,:,:)*vp_temp
!            ndot_prof(:,i)   =   ndot_prof(:,i)/vp_temp
         END DO

         ! Now calculate the rho density profile [part*m^-3]
         ALLOCATE(help3d(ndist1,ndist2, ndist3))
         DO l = 1, ndistns ! Edges
            s1 = REAL(l-1)/REAL(ndistns)
            s2 = REAL(l)/REAL(ndistns)
            DO m = 1, nbeams
               help3d = SUM(SUM(dist5d_prof(m,:,:,:,:,:),DIM=5),DIM=4)
               WHERE ((rho3d<=s1) .and. (rho3d>s2)) help3d = 0
               dense_prof(m,l) = SUM(SUM(SUM(help3d,DIM=3),DIM=2),DIM=1)
            END DO
         END DO

         ! Now calculate the Current Profile [A*m^-2]
         ! Now calculate the Parallel momentum Profile [kg*m^-2*s^-1]
         ! Now calculate the Perpendicular Pressure Profile [Pa]
         j_prof = 0
         momll_prof = 0
         DO m = 1, nbeams
            DO l = 1, ndistns ! Edges
               s1 = REAL(l-1)/REAL(ndistns)
               s2 = REAL(l)/REAL(ndistns)
               DO j = 1, ndist4
                  help3d = SUM(dist5d_prof(m,:,:,:,j,:),DIM=4)
                  WHERE ((rho3d<=s1) .and. (rho3d>s2)) help3d = 0
                  j_prof(m,l) = j_prof(m,l) + SUM(SUM(SUM(help3d,DIM=3),DIM=2),DIM=1)*vllaxis(j)
               END DO
               DO j = 1, ndist5
                  help3d = SUM(dist5d_prof(m,:,:,:,:,j),DIM=4)
                  WHERE ((rho3d<=s1) .and. (rho3d>s2)) help3d = 0
                  pperp_prof(m,l) = pperp_prof(m,l) + SUM(SUM(SUM(help3d,DIM=3),DIM=2),DIM=1)*vperpaxis(j)*vperpaxis(j)
               END DO
            END DO
            pperp_prof(m,:) = pperp_prof(m,:) * mass_beams(m) !p_perp
            momll_prof(m,:) = j_prof(m,:) * mass_beams(m) ! mvll
            j_prof(m,:) = j_prof(m,:) * charge_beams(m)
         END DO


         ! Now calculate the Heating Profile [W*m^-3]
         epower_prof = 0
         ipower_prof = 0
         DO m = 1, nbeams
            DO l = 1, ndistns ! Edges
               s1 = REAL(l-1)/REAL(ndistns)
               s2 = REAL(l)/REAL(ndistns)
               DO j = 1, ndist4
                  DO i = 1, ndist5
                     help3d = dist5d_prof(m,:,:,:,j,i)
                     WHERE ((rho3d<=s1) .and. (rho3d>s2)) help3d = 0
                     epower_prof(m,l) = epower_prof(m,l) + SUM(SUM(SUM(help3d*efact3d,DIM=3),DIM=2),DIM=1)*vdist2d(j,i)*mass_beams(m)*vdist2d(j,i)
                     ipower_prof(m,l) = ipower_prof(m,l) + SUM(SUM(SUM(help3d*ifact3d,DIM=3),DIM=2),DIM=1)*mass_beams(m)/vdist2d(j,i)
                  END DO
               END DO
            END DO
         END DO

         ! Normalize to velocity space volume element
         dvll = partvmax*2/ndist4 ! dVll
         dvperp = pi2*partvmax/ndist5 ! dVperp
         DO k = 1, ndist5 ! VPERP
            vp_temp = vperpaxis(k)*dvll*dvperp
            dist5d_prof(:,:,:,:,:,k) = dist5d_prof(:,:,:,:,:,k)/vp_temp
         END DO

         ! DEALLOCATIONS
         DEALLOCATE(rdistaxis, vperpaxis, vllaxis)
         DEALLOCATE(rho3d, efact3d, ifact3d, vdist2d, help3d)

      END IF

      CALL beams3d_write('DIAG')


!-----------------------------------------------------------------------
!     End Subroutine
!----------------------------------------------------------------------- 
      END SUBROUTINE beams3d_diagnostics
