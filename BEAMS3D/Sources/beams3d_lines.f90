!-----------------------------------------------------------------------
!     Module:        beams3d_lines
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/21/2012
!     Description:   This module contains the BEAMS3D field line
!                    variables.
!-----------------------------------------------------------------------
      MODULE beams3d_lines
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      
!-----------------------------------------------------------------------
!     Module Variables
!          myline    Dummy index
!          nparticles    Number of Particles
!          nsteps    Number of integration steps along fieldline
!          R_lines   Radial locations along fieldline [m] (npoinc per field period)
!          Z_lines   Vertical locations along field line [m]
!          PHI_lines Toroidal locations along field line [radians]
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL  ::  ltherm
      INTEGER  ::  ndist1, ndist2, ndist3, ndist4, ndist5, ndistns
      INTEGER  :: nparticles, nsteps, myline, mybeam, mytdex, myend, mystart_save,myend_save
      INTEGER  :: win_epower, win_ipower, win_ndot, win_dense, &
                  win_jprof, win_dist5d, win_momll, win_pperp
      REAL(rprec) :: xlast,ylast,zlast ! for storing position
      REAL(rprec) :: moment, mycharge, myZ, mymass, myv_neut(3), &
                     rand_prob, cum_prob, tau, next_t, &
                     partvmax, fact_crit, fact_pa, fact_vsound, &
                     partpmax, h1dist, h2dist, h3dist, h4dist, h5dist, &
                     rmin_dist, rmax_dist, zmin_dist, zmax_dist
      LOGICAL, ALLOCATABLE     :: neut_lines(:,:)
      INTEGER, ALLOCATABLE     :: end_state(:)
      REAL(rprec), ALLOCATABLE :: shine_through(:), shine_port(:), GFactor(:)
      REAL(rprec), DIMENSION(:,:), POINTER :: ndot_prof(:,:),epower_prof(:,:), &
                                  ipower_prof(:,:),j_prof(:,:), dense_prof(:,:), &
                                  momll_prof(:,:), pperp_prof(:,:)
      REAL(rprec), DIMENSION(:,:,:,:,:,:), POINTER :: dist5d_prof
      REAL(rprec), ALLOCATABLE :: R_lines(:,:),Z_lines(:,:),PHI_lines(:,:),vll_lines(:,:),moment_lines(:,:),&
                                  S_lines(:,:),U_lines(:,:),B_lines(:,:)

      END MODULE beams3d_lines
