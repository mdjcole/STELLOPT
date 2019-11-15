!-----------------------------------------------------------------------
!     Module:        fieldlines_calc_surface_fit
!     Authors:       S. Lazerson (samuel.lazeson@ipp.mpg.de)
!     Date:          11/06/2019
!     Description:   This subroutine 
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_calc_surface_fit(fit_line)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------\
      USE fieldlines_runtime
      USE fieldlines_lines
!-----------------------------------------------------------------------
!     Input Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in) :: fit_line
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      LOGICAL, DIMENSION(:), ALLOCATABLE :: lmask
      INTEGER ::                  i, j, n1, n2, d1
      INTEGER, DIMENSION(:), ALLOCATABLE :: idex
      REAL(rprec) :: x_old, y_old
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: r, z, p, x, y, r0, z0, &
                                                th, dth, nth, tth, &
                                                xth, yth, dr, arr1, arr2
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: x2d, y2d, dx2d, dy2d, &
                                                  x_new, y_new

      ! Extract the surface and axis data
      i = nsteps
      ALLOCATE(r(i),z(i),p(i),r0(i),z0(i),x(i),y(i),th(i),dth(i))
      FORALL(j=0:i-1) r(j+1) = R_lines(fit_line,j)
      FORALL(j=0:i-1) z(j+1) = Z_lines(fit_line,j)
      FORALL(j=0:i-1) p(j+1) = PHI_lines(fit_line,j)
      FORALL(j=0:i-1) r0(j+1) = R_lines(1,j)
      FORALL(j=0:i-1) z0(j+1) = Z_lines(1,j)

!      DO j =1, nsteps
!         WRITE(325,*) r(j),z(j)
!      END DO

      ! Calc helpers
      x = r-r0
      y = z-z0
      th = atan2(y,x)
      dth(1:i-1)=th(2:i)-th(1:i-1)
      th(1) = 0
      DO i = 2, nsteps
         th(i) = th(i-1)+dth(i-1)
      END DO

!      DO j =1, nsteps
!         WRITE(326,*) x(j),y(j)
!      END DO

      ! Sort the data in 2D grid
      n1 = nsteps/npoinc
      n2 = npoinc
      ALLOCATE(lmask(n1))
      ALLOCATE(idex(n1))
      ALLOCATE(tth(n1),xth(n1),yth(n1),dr(n1),arr1(n1))
      ALLOCATE(x2d(n1,n2),y2d(n1,n2))
      DO i = 1, n2
         d1 = 1
         DO j = 1, nsteps, npoinc
            xth(d1) = x(j)
            yth(d1) = y(j)
            d1 = d1 + 1
         END DO
         tth=atan2(yth,xth)
         WHERE(tth<0) tth = tth + pi2
         CALL SORT(n1,tth,idex)
         DO j = 1, n1
            arr1(j) = xth(idex(j))
         END DO
         xth = arr1
         DO j = 1, n1
            arr1(j) = yth(idex(j))
         END DO
         yth = arr1
         x2d(1,i) = xth(1)
         y2d(1,i) = yth(1)
         x_old = xth(1)
         y_old = yth(1)
         lmask = .TRUE.
         lmask(1) = .FALSE.
         DO j = 2, n1
            dr = (xth-x_old)**2+(yth-y_old)**2
            d1 = MINLOC(dr,1,MASK=lmask)
            x2d(j,i) = xth(d1)
            y2d(j,i) = yth(d1)
            x_old = xth(d1)
            y_old = yth(d1)
            lmask(d1) = .false.
         END DO
      END DO

      DO j =1, n1
         WRITE(327,*) x2d(j,1),y2d(j,1)
      END DO

      ! Now we calculate a poloidal like angle
      ALLOCATE(dx2d(n1,n2),dy2d(n1,n2))
      DO i = 2, n1
         dx2d(i-1,:)=x2d(i,:)-x2d(i-1,:)
         dy2d(i-1,:)=y2d(i,:)-y2d(i-1,:)
      END DO
      dx2d = SQRT(dx2d*dx2d+dy2d*dy2d)
      dy2d(1,:) = 0
      DO i = 2, n1
         dy2d(i,:) = dy2d(i-1,:)+dx2d(i,:)
      END DO
      DO i = 1, n2
         dy2d(:,i) = dy2d(:,i)/MAXVAL(dy2d(:,i),1)
      END DO

      ! Now calculated new 2D grid
      ALLOCATE(nth(256),arr2(256))
      ALLOCATE(x_new(256,n2),y_new(256,n2))
      FORALL(i=1:256) nth(i) = DBLE(i-1)/DBLE(255)
      DO i = 1, n2-1
         xth = dy2d(:,i)
         yth = x2d(:,i)
         DO j = 1, n1
            WRITE(400+i,*) xth(i),yth(i)
         END DO
         CALL spline_it(n1,xth,yth,256,nth,arr1,0)
         x_new(:,i) = arr1
         yth = y2d(:,i)
         DO j = 1, n1
            WRITE(500+i,*) xth(i),yth(i)
         END DO
         CALL spline_it(n1,xth,yth,256,nth,arr1,0)
         y_new(:,i) = arr1
      END DO

      DO j =1, 256
         WRITE(327,*) x_new(j,1),y_new(j,1)
      END DO



      ! DEALLOCATION
      DEALLOCATE(x_new,y_new)
      DEALLOCATE(nth)
      DEALLOCATE(dx2d,dy2d)
      DEALLOCATE(x2d,y2d)
      DEALLOCATE(tth,xth,yth,dr,arr1)
      DEALLOCATE(idex)
      DEALLOCATE(lmask)
      DEALLOCATE(r,z,p,r0,z0,x,y,th,dth)      
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_calc_surface_fit
