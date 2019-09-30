subroutine lbfgsb_driver(fcn, m, n, x, l, u, nbd, dx_init, maxfev, ftol, &
                         m_lbfgsb, factr, pgtol, iprint, nprint, info, nfev )
!-----------------------------------------------------------------------
!     Subroutine:    lbfgsb_driver
!                  LBFGSB_DRIVER  in Fortran 90
!    --------------------------------------------------------------
!             CUSTOMIZED STELLOPT DRIVER FOR L-BFGS-B
!    --------------------------------------------------------------
!     Authors:       J.C. Schmitt
!     Date:          2018, 2019
!     Description:   This subroutine performs a BFGS quasi-newton
!                    optimization of user-supplied fcn beginning at x
!                    Uses L-BFGS-Bv3.0 (See license information below)
!
!    Input parameters
!    fcn             User-specified subroutine to compute objective function
!    m               Number of target functions.
!    n               Number of variables.
!    x               Vector of variables of length m
!    l               Upper bounds
!    u               Lower bounds
!    nbd             Defines the type of bounds (see lbfgsb.f documentation)
!    dx_init         The step-size 'dx' used for the Jacobian calculation
!                       To Do: Use the same stepsize varialbe as in Lev-Marq
!    maxfev          Maximum number of function evals.
!    ftol            Tolerance on norm of function value.
!    m_lbfgsb        The maximum number of variable metric corrections used to
!                       define the limited memory matrix
!    factr           On entry factr >= 0 is specified by the user.  The iteration
!                       will stop when
!         (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
!         where epsmch is the machine precision, which is automatically
!         generated by the code. Typical values for factr: 1.d+12 for
!         low accuracy; 1.d+7 for moderate accuracy; 1.d+1 for extremely
!         high accuracy.
!    pgtol          On entry pgtol >= 0 is specified by the user.  The iteration
!                       will stop when
!                 max{|proj g_i | i = 1, ..., n} <= pgtol
!                    where pg_i is the ith component of the projected gradient.   
!    iprint         Controls the frequency and type of output generated:
!        iprint<0    no output is generated;
!        iprint=0    print only one line at the last iteration;
!        0<iprint<99 print also f and |proj g| every iprint iterations;
!        iprint=99   print details of every iteration except n-vectors;
!        iprint=100  print also the changes of active set and final x;
!        iprint>100  print details of every iteration including x and g;
!       When iprint > 0, the file iterate.dat will be created to
!                        summarize the iteration.
!    nprint, info, nfev: Also passed in. To be updated
!

!-----------------------------------------------------------------------

!=======================================================================
!  L-BFGS-B Licencse information follows below. Also see the
!  LBFGSBv30_License.txt file
!=======================================================================
!  L-BFGS-B is released under the “New BSD License” (aka “Modified BSD License”        
!  or “3-clause license”)                                                              
!  Please read attached file License.txt                                               
!                                        
!       L-BFGS-B is a code for solving large nonlinear optimization
!            problems with simple bounds on the variables.
!
!       The code can also be used for unconstrained problems and is
!       as efficient for these problems as the earlier limited memory
!                         code L-BFGS.
!
!       This driver illustrates how to control the termination of the
!       run and how to design customized output.
!
!    References:
!
!       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
!       memory algorithm for bound constrained optimization'',
!       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
!
!       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
!       Subroutines for Large Scale Bound Constrained Optimization''
!       Tech. Report, NAM-11, EECS Department, Northwestern University,
!       1994.
!
!
!         (Postscript files of these papers are available via anonymous
!          ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
!
!                             *  *  *
!
!         February 2011   (latest revision)
!         Optimization Center at Northwestern University
!         Instituto Tecnologico Autonomo de Mexico
!
!         Jorge Nocedal and Jose Luis Morales
!
!    **************

! Other modules that are used by LBFGSB
  USE mpi_params
  USE safe_open_mod
  USE stel_kinds
  USE fdjac_mod, ONLY: FDJAC2_MP_QUEUE, FLAG_CLEANUP, &
                       FLAG_CLEANUP_LBFGSB, FLAG_SINGLETASK, &
                       flip, jac_order, h_order, jac_err, jac_index
  implicit none
 
!DEC$ IF DEFINED (MPI_OPT)
    INCLUDE 'mpif.h'
!DEC$ ENDIF

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  INTEGER :: m, n, maxfev, m_lbfgsb, iprint, nprint, info
  INTEGER, INTENT(inout) :: nfev
  INTEGER, DIMENSION(n) :: nbd
  REAL(rprec), INTENT(in) ::  dx_init, ftol, factr, pgtol
  REAL(rprec), DIMENSION(n) :: x, l, u

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  REAL(rprec) :: fvec(m)
  INTEGER :: iflag, istat, iunit, ikey, nvar, iter
  INTEGER :: ii, jj
  REAL(rprec) :: f            ! The function,
  REAL(rprec) :: g(n),  fjac(m,n)  ! gradient and the jacobian
  REAL(rprec) :: f_min 
  REAL(rprec), ALLOCATABLE, DIMENSION(:) :: x_min, fvec_min
  REAL(rprec), ALLOCATABLE, DIMENSION(:) :: f_array
  !real(dp), parameter    :: factr  = 0.0d0, pgtol  = 0.0d0
  character(len=60)      :: task, csave
  logical                :: lsave(4)
  integer                :: isave(44)
  real(dp)               :: dsave(29)
  integer,  allocatable  :: iwa(:)
  real(dp), allocatable  :: wa(:)
  real(dp)               :: t1, t2
  logical                :: f_success
  real(rprec) :: bigno = 1.0e10

!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
  EXTERNAL fcn
  REAL(rprec), EXTERNAL :: dpmpar, enorm

!DEC$ IF DEFINED (MPI_OPT)
    ! Get mpi parameters - check for mpi errors after each call 
    ! to help identify bugs, errors and other problems

    ! Assign the process rank to variable 'myid'
    CALL MPI_COMM_RANK (MPI_COMM_STEL, myid, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
    IF (ierr_mpi .ne. 0) STOP 'MPI_COMM_RANK error in LBFGSB_DRIVER'

    ! Assign the size of the group to 'numprocs'
    CALL MPI_COMM_SIZE (MPI_COMM_STEL, numprocs, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
    IF (ierr_mpi .ne. 0) STOP 'MPI_COMM_SIZE error in LBFGSB_DRIVER'
!DEC$ ELSE
    ! If not using mpi, assign a meaningful value to 'numprocs',
    !numprocs = 1
    ! 'myid' is set to the default value of 'master' in mpi_params
    ! if not using mpi, just stop
    stop
!DEC$ ENDIF

!DEC$ IF DEFINED (MPI_OPT)
    IF ((numprocs > n) .and. (myid .eq. master)) THEN
      WRITE (6, '(2A,2X,I5)'), &
             'K====Warning: more processors have been requested ', &
             'than the maximum (nvar) required = ', n
    END IF
!DEC$ ENDIF

  ! Check the input parameters for errors.
  IF ( (n<0) .or. (m<0) .or. (ftol<0) .or. (pgtol<0) &
    .or. (maxfev<0) ) THEN
     STOP "K====Error! LBFGSB_DRIVER called with improper input arguments."
  END IF

  ! Allocate memory for 
  !
  !      nbd(n)-in stellopt_optimize for now
  !      x(n), l(n), u(n) - in (or prior to) stellopt_optimize)
  !      g(n) -- need it
  !      allocate ( iwa(3*n) )
  !      allocate ( wa(2*m*n + 5*n + 11*m*m + 8*m) )

  ! Alllocate memory for fdjac2_mp_queue and LBFGSB
  !ALLOCATE (x_min(n), fvec_min(m), flip(n), jac_order(n), &
  !          f_array(n), h_order(n), jac_err(n), jac_index(n), &
  !          stat=istat)
  ! IF (istat .ne. 0) STOP 'K====ALLOCATION ERROR 1 IN LBFGSB_DRIVER'
  ALLOCATE (x_min(n), stat=istat)
  IF (istat .ne. 0) STOP 'K====ALLOCATION ERROR 2 IN LBFGSB_DRIVER'
  ALLOCATE (fvec_min(m), stat=istat)
  IF (istat .ne. 0) STOP 'K====ALLOCATION ERROR 3 IN LBFGSB_DRIVER'
  ALLOCATE (flip(n), stat=istat)
  IF (istat .ne. 0) STOP 'K====ALLOCATION ERROR 4 IN LBFGSB_DRIVER'
  ALLOCATE (jac_order(n), stat=istat)
  IF (istat .ne. 0) STOP 'K====ALLOCATION ERROR 5 IN LBFGSB_DRIVER'
  ALLOCATE (f_array(n), stat=istat)
  IF (istat .ne. 0) STOP 'K====ALLOCATION ERROR 6 IN LBFGSB_DRIVER'
  ALLOCATE (h_order(n), stat=istat)
  IF (istat .ne. 0) STOP 'K====ALLOCATION ERROR 7 IN LBFGSB_DRIVER'
  ALLOCATE (jac_err(n), stat=istat)
  IF (istat .ne. 0) STOP 'K====ALLOCATION ERROR 8 IN LBFGSB_DRIVER'
  ALLOCATE (jac_index(n), stat=istat)
  IF (istat .ne. 0) STOP 'K====ALLOCATION ERROR 9 IN LBFGSB_DRIVER'
  ALLOCATE (iwa(3*n), stat=istat)
  IF (istat .ne. 0) STOP 'K====ALLOCATION ERROR 10 IN LBFGSB_DRIVER'
  ALLOCATE (wa(2*m*n + 5*n + 11*m*m + 8*m), stat=istat)
  IF (istat .ne. 0) STOP 'K====ALLOCATION ERROR 11 IN LBFGSB_DRIVER'

  !  Initialize the newly allocated variables (at least the ones that need it)

  flip = .false. 
  jac_order = 0
  g = 0
  iwa = 0
  wa = 0

  ! Define to zero.
  f = 0

  !     Set up workers communicator (only master performs this task)
  !     for the initial run only
!DEC$ IF DEFINED (MPI_OPT)
    IF (myid .ne. master) THEN
       ikey = MPI_UNDEFINED
    ELSE
       ikey = WORKER_SPLIT_KEY+1  ! Mathces lmdif.f
    END IF
    CALL MPI_COMM_SPLIT(MPI_COMM_STEL, ikey, worker_id, &
                        MPI_COMM_WORKERS, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF

 
  !     We now define the starting point.
  !   x(1..n) is done. Either passed in or modified.
  !     We now write the heading of the LBFGSB output.

  IF (myid .eq. master) then
     WRITE (6, 500) numprocs
500    FORMAT (/,' Beginning LBFGS-B Finite Difference Iterations', /, &
             ' Number of Processors: ', i6, //, 40('='), / ,2x, &
             'Iteration', 3x, 'Processor', 7x, 'Chi-Sq', 7x, &
             /, 40('='))
  END IF


  !     We start the iteration by initializing task.
  ! 
  task = 'START'

  !        ------- the beginning of the loop ----------

  do while( task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. &
                task.eq.'START') 
    
    ! This is the call to the L-BFGS-B code.
    IF (myid .eq. master) THEN
      call setulb(n,m_lbfgsb,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint, &
                  csave,lsave,isave,dsave)
    END IF
!DEC$ IF DEFINED (MPI_OPT)
        CALL MPI_BCAST(task, 60, MPI_REAL8, master, &
                       MPI_COMM_STEL, ierr_mpi)
        IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF

    if (task(1:2) .eq. 'FG') then
      ! The minimization routine has returned to request the
      ! function f and gradient g values at the current x.

      ! Evaluate function at point and calculate its gradient.
      IF (myid .eq. master) THEN
        ! Set iflag to flag_singletask to tell fcn what to do
        iflag = FLAG_SINGLETASK
        f_success = .false.

        !  Compute function value f 
        CALL fcn (m, n, x, fvec, iflag, nfev)

        ! Increment nfev by 1 (Is this correct and necessary? JCS)
        ! nfev = nfev + 1

        IF (iflag .ne. 0) THEN
          WRITE(6,*) "<----Evaluation Failed!"
          !STOP "K====ERROR in LBFGSB_DRIVER!"
          do ii = 1,m
            fvec(ii) = bigno
          end do
        ELSE
          f_success = .true.
        END IF

        ! Calculate the Euclidean norm here- this is 'f'
        f = enorm(m,fvec)

        WRITE(6, '(2x,i6,8x,i3,7x,1es16.8,a,1es16.8,a)'), 0, myid, &
              f, '(', f**2, ')'

        ! Write useful information to 'xvec.dat'
        iunit = 12; istat = 0
        CALL safe_open(iunit,istat,'xvec.dat','unknown','formatted', &
                       ACCESS_IN='APPEND')
        ! Number of variables, 'n', followed by iteration count
        WRITE(iunit,'(2(2X,I12.12))') n, nfev
        ! The variables, 'x'
        WRITE(iunit,'(10ES22.12E3)') x(1:n)
        ! The function value, 'f'
        WRITE(iunit,'(ES22.12E3)') f
        CLOSE(iunit)

!DEC$ IF DEFINED (MPI_OPT)
          ! Mark the communicator for deallocation -- only required
          ! on first iteration (nfev .eq. 0)
          if (nfev .eq. 0) then
            CALL MPI_COMM_FREE(MPI_COMM_WORKERS, ierr_mpi)
            IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)      
          end if
!DEC$ ENDIF
      END IF ! End of 'Evaluate function at point and calculate its gradient'
      ! Workers jump to here.

!DEC$ IF DEFINED (MPI_OPT)
        CALL MPI_BCAST(f_success,1, MPI_REAL8, master, &
                       MPI_COMM_STEL, ierr_mpi)
        IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
        CALL MPI_BCAST(x,n, MPI_REAL8, master, &
                       MPI_COMM_STEL, ierr_mpi)
        IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF

      ! Cleanup after fcn call
      iflag = FLAG_CLEANUP
      ! What does master do, what does worker do?
      IF (myid .eq. master) iflag = FLAG_CLEANUP_LBFGSB
      CALL fcn (m, n, x, fvec, iflag, nfev)
      ! The master will do the bfgs cleanup.
      ! The workers will do 'regular' cleanup
 
      ! Increment nfev by 1 (Is this correct and necessary? JCS)
      nfev = nfev + 1


      !        Compute gradient g for the sample problem.
      !    Put a vector into g(1..n)
      !

!DEC$ IF DEFINED (MPI_OPT)
        ! The master will broadcast its value of iflag to all other workers
        CALL MPI_BCAST(iflag, 1, MPI_INTEGER, master, &
                       MPI_COMM_STEL, ierr_mpi)
        IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)

        ! If the value of iflag >= 0, then master will broadcast its 'fvec'
        ! array to all other workers (length of 'm')
        IF (iflag .ge. 0) CALL &
            MPI_BCAST(fvec, m, MPI_REAL8, master,  &
                      MPI_COMM_STEL, ierr_mpi)
        IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF


      ! Now that all of the processes have 'fvec', make sure they
      ! also have 'f'
      f = enorm(m, fvec)

!DEC$ IF DEFINED (MPI_OPT)
        ! Sync all of the processors at this point - otherwise
        CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)
        IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF

      !if (f_success) then
        CALL fdjac2_mp_queue(fcn, f, m, n, x, fvec, fjac, &
                             m, iflag, nfev, dx_init, f_min, x_min, &
                             fvec_min, f_array, .true.)
        nfev = nfev + n
        g = matmul(fvec, fjac) / f
      !else
      !  g(1:n) = bigno
      !end if

    else  !  if (task(1:2) .eq. 'FG') then
      if (task(1:5) .eq. 'NEW_X') then   
        ! The minimization routine has returned with a new iterate.
        ! At this point have the opportunity of stopping the iteration 
        ! or observing the values of certain parameters.
        ! Two examples of stopping tests are shown.
        ! Note: task(1:4) must be assigned the value 'STOP' to terminate  
        !  the iteration and ensure that the final results are
        !  printed in the default format. The rest of the character
        !  string TASK may be used to store other information.
        !  1) Terminate if the total number of f and g evaluations
        !     exceeds nfev.
        IF (myid .eq. master) THEN
          if (iprint .ge. 1) THEN
            write (6,*) 'Current X='
            write (6,'((1x,1p, 6(1x,d22.16)))') (x(ii),ii = 1,n)
            write (6,*) 'Current theta='
            write (6,'((1x,1p, 1x,d22.16))') (dsave(1))
            write (6,*) 'Previous f(x)='
            write (6,'((1x,1p, 1x,d22.16))') (dsave(2))
            write (6,*) 'factr*epsmch='
            write (6,'((1x,1p, 1x,d22.16))') (dsave(3))
            write (6,*) '2-norm of line search vector='
            write (6,'((1x,1p, 1x,d22.16))') (dsave(4))
            write (6,*) 'slope of line='
            write (6,'((1x,1p, 1x,d22.16))') (dsave(11))
            write (6,*) 'isave(34)='
            write (6,'((1x,1p, 1x,i11))') (isave(34))
            write (6,*) 'fev='
            write (6,'((1x,1p, 1x,i11))') (nfev)
            write (6,*) 'maxfev='
            write (6,'((1x,1p, 1x,i11))') (maxfev)
            write (6,*) 'lsave=',lsave
          end if

          if (isave(34) .ge.  maxfev)  &
              task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'

          !  2) Terminate if  |proj g|/(1+|f|) < 1.0d-10, where 
          !     "proj g" denoted the projected gradient

          if (dsave(13) .le. pgtol * (1.0d0 + abs(f))) &
              task='STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL'

          !  3) Terminate if  f()x  < factr*epsmch

          if (abs(f) .le. ftol ) &
              task='STOP: THE F IS SUFFICIENTLY SMALL'

          ! We now wish to print the following information at each
          ! iteration:
          !  1) the current iteration number, isave(30),
          !  2) the total number of f and g evaluations, isave(34),
          !  3) the value of the objective function f,
          !  4) the norm of the projected gradient,  dsve(13)

          !  See the comments at the end of driver1 for a description
          !  of the variables isave and dsave.
         
          write (6,'(2(a,i5,4x),a,1p,d12.5,4x,a,1p,d12.5)') 'Iterate', &
                 isave(30),'nfg =',isave(34),'f =',f,'|proj g| =',dsave(13)

          ! If the run is to be terminated, we print also the information
          ! contained in task as well as the final value of x.

          if (task(1:4) .eq. 'STOP') then
              write (6,*) task  
              if (n .le. 20) then
                write (6,*) 'Final X='
                write (6,'((1x,1p, 6(1x,d26.18)))') (x(ii),ii = 1,n)
              end if
          end if
          ! master broadcasts 'task' to the workers
!DEC$ IF DEFINED (MPI_OPT)
            CALL MPI_BCAST(task, 60, MPI_REAL8, master, &
                           MPI_COMM_STEL, ierr_mpi)
            IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF

        end if ! end of master activities
      end if !if (task(1:5) .eq. 'NEW_X') then   
    end if !  else  !  if (task(1:2) .eq. 'FG') then

  end do
!           ---------- the end of the loop -------------
 
!     If task is neither FG nor NEW_X we terminate execution.

! DEALLOCATE STUFF
  IF (ALLOCATED(x_min)) DEALLOCATE(x_min)
  IF (ALLOCATED(fvec_min)) DEALLOCATE(fvec_min)
  IF (ALLOCATED(flip)) DEALLOCATE(flip)
  IF (ALLOCATED(jac_order)) DEALLOCATE(jac_order)
  IF (ALLOCATED(f_array)) DEALLOCATE(f_array)
  IF (ALLOCATED(h_order)) DEALLOCATE(h_order)
  IF (ALLOCATED(jac_err)) DEALLOCATE(jac_err)
  IF (ALLOCATED(jac_index)) DEALLOCATE(jac_index)
  IF (ALLOCATED(iwa)) DEALLOCATE(iwa)
  IF (ALLOCATED(wa)) DEALLOCATE(wa)



! Exit message
  if (myid .eq. master) then
    write(*,"(A,I2,A)") "<----LBFGSB_DRIVER terminated" ! after ", iter, " iterations."
    write(*,"(A,I6)") "<----Function evaluations: ", nfev
  end if

  
  return
end subroutine lbfgsb_driver

!======================= The end of lbfgsb_driver ============================

