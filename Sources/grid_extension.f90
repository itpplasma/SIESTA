!*******************************************************************************
!>  @file grid_extension.f90
!>  @brief Contains module @ref grid_extension.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  This file contains utilities for extending the vmec grid.
!*******************************************************************************
      MODULE grid_extension
      USE vmec_info, ONLY: jcurrumnc, jcurrvmnc
      USE stel_kinds
      USE stel_constants
      USE shared_data, ONLY: lasym, r1_i, z1_i, ru_i, zu_i, rv_i, zv_i, nsp,   &
                             lverbose
      USE island_params, ns=>ns_i, mpol=>mpol_i, ntor=>ntor_i, &
                         ntheta=>nu_i, nzeta=>nv_i, mnmax=>mnmax_i, nsv=>ns_v, &
                         ohs=>ohs_i, nfp=>nfp_i
      USE island_params, ONLY: phipf_i, chipf_i, presf_i, nsh
      USE v3_utilities, ONLY: assert
      USE siesta_namelist, ONLY: l_vessel
      USE descriptor_mod, ONLY: iam

      IMPLICIT NONE

!*******************************************************************************
!  grid extension module parameters.
!*******************************************************************************
!>  Linear extension type.
      INTEGER, PARAMETER   :: LIN = 1
!>  Quadtratic extension type.
      INTEGER, PARAMETER   :: QUAD = 2
!>  Cubic extension type.
      INTEGER, PARAMETER   :: CUB = 3
!>  Divide by two constant.
      REAL (dp), PARAMETER :: p5 = 0.5_dp

!*******************************************************************************
!  grid extension module variables.
!*******************************************************************************
!>
      INTEGER :: interp_type
      INTEGER :: nse

      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: ri, zi

      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: r1_vmec, z1_vmec
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: ru_vmec, zu_vmec
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: rv_vmec, zv_vmec
  
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: r1s_vmec, r1ss_vmec
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: z1s_vmec, z1ss_vmec
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: r1su_vmec, z1su_vmec

      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: rmnc_v, zmns_v
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: rmns_v, zmnc_v

      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: r_ext, z_ext
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: r1_ves, z1_ves

      REAL(dp)            :: s_edge, ming, maxg

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Extend the siesta grid out to the vacuum vessel.
!>
!>  The grid can be extended using a linear extension, quadratic, or cubic
!>  extension.
!>
!>  @params[in]  wout_file VMEC wout file.
!>  @params[in]  itype     Type of extension.
!>  @params[out] istat     Status of the extension.
!-------------------------------------------------------------------------------
      SUBROUTINE grid_extender(wout_file, itype, istat)
      USE fourier, ONLY: f_cos, f_sin, f_none, f_sum, f_du, f_dv
      USE v3_utilities, ONLY: assert
      USE vmec_info, ONLY: rmnc_i, rmns_i, zmns_i, zmnc_i,                     &
                           jcurrumnc, jcurrvmnc, jcurrumns, jcurrvmns
      USE descriptor_mod, ONLY: SIESTA_COMM
      USE siesta_namelist, ONLY: nsin_ext, ext_edge_fac
      USE shared_data, ONLY: lverbose
      USE island_params, ONLY: fourier_context
      USE metrics, ONLY: LoadGrid

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER(*), INTENT(in)                :: wout_file
      CHARACTER(*), INTENT(in)                :: itype
      INTEGER, INTENT(out)                    :: istat

!  Local Arguments.
      CHARACTER(LEN=10)                       :: chfit

      INTEGER                                 :: j
      INTEGER                                 :: m
      INTEGER                                 :: n
      INTEGER                                 :: u
      INTEGER                                 :: v
      INTEGER                                 :: i
      INTEGER                                 :: k
      INTEGER                                 :: js
      INTEGER                                 :: js1
      INTEGER                                 :: nsp1
      REAL(dp)                                :: rho
      REAL(dp)                                :: vol0
      REAL(dp)                                :: vol1
      REAL(dp)                                :: vols
      REAL(dp)                                :: rho1
      REAL(dp)                                :: rho2
      REAL(dp)                                :: rho3
      REAL(dp)                                :: vp1
      REAL(dp)                                :: s_ext
      REAL(dp)                                :: delv
      REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: r12
      REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: rs12
      REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: zs12
      REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: ru12
      REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: zu12
      REAL(dp), DIMENSION(:), ALLOCATABLE     :: vp
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: tmp
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: tmp2
      REAL(dp), DIMENSION(:), ALLOCATABLE     :: tmp1

      INTEGER, DIMENSION(:), ALLOCATABLE      :: vtor_modes
      REAL(dp)                                :: fac
      REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: ru_ves
      REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: zu_ves

      CLASS (fourier_class), POINTER          :: four => null()

!  Start of executable code.
      CALL init_extender

!  chfit = 'quad'
!  chfit = TRIM(itype)
      IF (INDEX(itype,'lin').EQ.1 .OR. INDEX(itype,'LIN').EQ.1) THEN
         interp_type = LIN; chfit = "linear"
      ELSE IF (INDEX(itype,'quad').EQ.1 .OR. INDEX(itype,'QUAD').EQ.1) THEN
         interp_type = QUAD; chfit = "quadratic"
      ELSE IF (INDEX(itype,'cub').EQ.1 .OR. INDEX(itype,'CUB').EQ.1) THEN
         interp_type = CUB; chfit = "cubic"
      ELSE
         CALL ASSERT(.FALSE., 'Unknown interpolation <TYPE> specified')
      END IF

      ALLOCATE(ri(ntheta,nzeta), zi(ntheta,nzeta))
      ALLOCATE(r1_ves(ntheta,nzeta), z1_ves(ntheta,nzeta))

!  FIXME: The vessel file could have a different periodicity than the equilbirum.
      ALLOCATE(vtor_modes(-ntor_v:ntor_v))
      
      DO n = 0, ntor_v
         vtor_modes(n) = n
      END DO
      four => fourier_class(mpol_v, ntor_v, ntheta, nzeta, nfp, lasym,         &
     &                      vtor_modes)
      DEALLOCATE(vtor_modes)

!  Set internal normalization
      rmnc_v = rmnc_v/four%orthonorm
      zmns_v = zmns_v/four%orthonorm

!  FIXME: What about non stellarator symmetric vessels?
      CALL four%toijsp(rmnc_v, r1_ves, f_none, f_cos)
      CALL four%toijsp(zmns_v, z1_ves, f_none, f_sin)
      ALLOCATE(ru_ves(ntheta,nzeta), zu_ves(ntheta,nzeta))
      CALL four%toijsp(rmnc_v, ru_ves, f_du, f_cos)
      CALL four%toijsp(zmns_v, zu_ves, f_du, f_sin)

      DEALLOCATE(four)

!  The quadratic blend reaches the wall with slope 2*(X_ves - X_vmec) - fac*X_s,
!  so a large fac compared with the plasma-wall gap folds the grid.
      fac = ext_edge_fac_limit(ext_edge_fac, r1_vmec, z1_vmec, ru_vmec,        &
                               zu_vmec, r1s_vmec, z1s_vmec, r1su_vmec,         &
                               z1su_vmec, r1_ves, z1_ves, ru_ves, zu_ves)
      CALL assert(fac .gt. zero, 'Vessel extension folds for every '  //       &
                  'EXT_EDGE_FAC; check the vessel parametrisation.')
      IF (fac .lt. ext_edge_fac .and. lverbose .and. iam .eq. 0) THEN
         WRITE (*,1001) ext_edge_fac, fac
      END IF
      r1s_vmec = fac*r1s_vmec
      z1s_vmec = fac*z1s_vmec
      DEALLOCATE(ru_ves, zu_ves, r1su_vmec, z1su_vmec)

!  Compute extended "s" coordinate for s>1, assuming Vp(s) = Vp(s=1)*s
!  where Vp(s=1) is obtained from extrapolating VMEC values. Results is
!    V(s) = Vp*(s^2-1)/2 + V0 => s_edge = SQRT(2(V1-V0)/Vp + 1)
      vol0 = get_volume(zero)    !rho = 0: s=1 vmec
      vol1 = get_volume(one)     !rho = 1: vessel

      vp1 = GetVP1()
      s_edge = SQRT(1 + 2*(vol1 - vol0)/vp1)

!  NOTE: at js=nse, s(js) == (js-1)*hs = s_edge
      nse = NINT(s_edge*ohs) + 1

      nsv = nse - ns
      s_edge = (nsv - 1)*hs_i + 1

!  Rounding discrepancy in volume due to integer nse
      delv = vol1 - (vol0 + p5*vp1*(s_edge*s_edge - 1))
      delv = delv/(s_edge-1)**2

      ALLOCATE(r_ext(ntheta,nzeta,2:nsv),                                      &
               z_ext(ntheta,nzeta,2:nsv),                                      &
               stat=istat)

      DO j = 2, nsv
         s_ext = 1 + (j-1)*hs_i    !at j=nsv, s_ext = s_edge

      !V(s) = Vp*(s^2-1)/2 + V0
         vols = vol0 + p5*vp1*(s_ext**2-1)   !matches vp exactly at s=1
         vols = vols + delv*(s_ext-1)**2     !matches vol1 exactly at s=s_edge
         vols = MIN(vols, vol1)
         CALL remap_radial(rho, vols)

         rho1 = 1 - rho
         IF (interp_type.EQ.LIN) THEN
            r_ext(:,:,j) = rho1*r1_vmec + rho*r1_ves
            z_ext(:,:,j) = rho1*z1_vmec + rho*z1_ves
         ELSE
            rho2 = rho*rho
            IF (interp_type .eq. QUAD) THEN
!  QUADRATIC (R' match at rho = 0)
               r_ext(:,:,j) = p5*r1s_vmec*(1 - rho2 - rho1*rho1)               &
                            + r1_vmec*(1-rho2) + r1_ves*rho2
               z_ext(:,:,j) = p5*z1s_vmec*(1 - rho2 - rho1*rho1)               &
                            + z1_vmec*(1-rho2) + z1_ves*rho2
            ELSE IF (interp_type .eq. CUB) THEN
               rho3 = rho*rho2
               CALL ASSERT(.FALSE., 'Unknown interpolation <TYPE> specified')
            END IF
         END IF
      END DO

      nse = nse - 1 !Not to double count s=1 surface
      nse = MIN(ns + nsin_ext, nse)

      IF (lverbose .and. iam .eq. 0) THEN
         WRITE (*,1000) nse, ns, nsin_ext
      END IF
1000  FORMAT('Using ',i4,' Total Surfaces.',/,                                 &
             'Number of VMEC surfaces ',i4,/,                                  &
             'Number of extended surfaces ',i4)
1001  FORMAT('EXT_EDGE_FAC reduced from ',f6.3,' to ',f6.3,                    &
             ' to keep the vessel extension from folding.')

!  Assign the new value of ns.
! REASSIGN THE NEW VALUE OF NS
      nsp = ns
      ns = nse
      nsh = ns - 1

      IF (nsp .NE. ns) THEN
         DEALLOCATE(r1_i)
         DEALLOCATE(z1_i)
         DEALLOCATE(ru_i)
         DEALLOCATE(zu_i)
         DEALLOCATE(rv_i)
         DEALLOCATE(zv_i)
      END IF

      IF (nsp .NE. ns) THEN
         nsp1 = nsp + 1
         ALLOCATE(tmp(0:mpol,-ntor:ntor,ns), tmp2(0:mpol,-ntor:ntor, ns))
         tmp(:,:,1:nsp) = rmnc_i
         tmp2(:,:,1:nsp) = zmns_i
         DEALLOCATE(rmnc_i, zmns_i)
         ALLOCATE(rmnc_i(0:mpol,-ntor:ntor,ns),                                &
                  zmns_i(0:mpol,-ntor:ntor,ns), stat=istat)
         rmnc_i(:,:,1:nsp) = tmp(:,:,1:nsp)
         zmns_i(:,:,1:nsp) = tmp2(:,:,1:nsp)
         CALL fourier_context%tomnsp(r_ext(:,:,2:nsin_ext + 1),                &
                                     rmnc_i(:,:,nsp1:ns), f_cos)
         CALL fourier_context%tomnsp(z_ext(:,:,2:nsin_ext + 1),                &
                                     zmns_i(:,:,nsp1:ns), f_sin)

! PUT PROFILES ONTO EXTENDED MESH
         tmp = 0
         tmp(:,:,1:nsp) = jcurrumnc
         DEALLOCATE(jcurrumnc)
         ALLOCATE(jcurrumnc(0:mpol,-ntor:ntor,ns))
         jcurrumnc = tmp

         tmp = 0
         tmp(:,:,1:nsp) = jcurrvmnc
         DEALLOCATE(jcurrvmnc)
         ALLOCATE(jcurrvmnc(0:mpol,-ntor:ntor,ns))
         jcurrvmnc = tmp

! DO THE SAME FOR THE phipf_i, chipf_i AND presf_i
! phipf, chipf will be computed in PCHELMS in the vacuum region
         ALLOCATE(tmp1(ns))
         tmp1 = 0
         tmp1(1:nsp) = phipf_i
         DEALLOCATE(phipf_i)
         ALLOCATE(phipf_i(ns))
         phipf_i = tmp1

         tmp1 = 0
         tmp1(1:nsp) = chipf_i
         DEALLOCATE(chipf_i)
         ALLOCATE(chipf_i(ns))
         chipf_i = tmp1

         tmp1(1:nsp) = presf_i
         tmp1(nsp + 1:) = presf_i(nsp)
         DEALLOCATE(presf_i)
         ALLOCATE(presf_i(ns))
         presf_i = tmp1
         DEALLOCATE(tmp1)

         IF (lasym) THEN
            tmp(:,:,1:nsp) = rmns_i
            tmp2(:,:,1:nsp) = zmnc_i
            DEALLOCATE(rmns_i, zmnc_i)
            ALLOCATE(rmns_i(0:mpol,-ntor:ntor,ns),                             &
                     zmnc_i(0:mpol,-ntor:ntor,ns), stat=istat)
            rmns_i(:,:,1:nsp) = tmp(:,:,1:nsp)
            zmnc_i(:,:,1:nsp) = tmp2(:,:,1:nsp)
            CALL fourier_context%tomnsp(r_ext(:,:,1:nsin_ext),                 &
     &                                  rmns_i(:,:,nsp1:ns), f_sin)
            CALL fourier_context%tomnsp(z1_i(:,:,nsp1:ns),                     &
     &                                  z_ext(:,:,1:nsin_ext), f_cos)

            tmp = 0
            tmp(:,:,1:nsp)=jcurrumns
            DEALLOCATE(jcurrumns)
            ALLOCATE(jcurrumns(0:mpol,-ntor:ntor,ns))
            jcurrumns=tmp

            tmp=0
            tmp(:,:,1:nsp)=jcurrvmns
            DEALLOCATE(jcurrvmns)
            ALLOCATE(jcurrvmns(0:mpol,-ntor:ntor,ns))
            jcurrvmns=tmp
         END IF

         DEALLOCATE(tmp, tmp2)

      END IF

#if defined(NOSKS)
      DO js=1, ns
         DO v=1, nzeta
            DO u=1, ntheta
               WRITE(1600,*) r1_i(u,v,js), z1_i(u,v,js)
               CALL FLUSH(1600)
            END DO
         END DO
      END DO

      LoadGrid(istat)

      ALLOCATE (vp(nse))
      CALL SurfAvgLocal(vp,sqrtg,1,nse)
      vp(1) = 0
      DO js=2, nse
         ming = MINVAL(sqrtg(:,:,js) )
         maxg = MAXVAL(sqrtg(:,:,js))
         WRITE(5000,'(i4,1p,3E14.4)') js, -vp(js), ming, maxg
      END DO
      CALL FLUSH(5000)
      DEALLOCATE (vp,r_ext, z_ext, ru_ext, zu_ext)
      DEALLOCATE (ri,zi)
#endif

      CALL LoadGrid(istat)

      DEALLOCATE(r1s_vmec)
      DEALLOCATE(z1s_vmec)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Initialize the grid extender.
!-------------------------------------------------------------------------------
      SUBROUTINE init_extender

      IMPLICIT NONE

!  Start of executable code.
      ALLOCATE(r1_vmec(ntheta,nzeta),z1_vmec(ntheta,nzeta))
      ALLOCATE(ru_vmec(ntheta,nzeta),zu_vmec(ntheta,nzeta))
      ALLOCATE(rv_vmec(ntheta,nzeta),zv_vmec(ntheta,nzeta))
      r1_vmec = r1_i(:,:,ns); z1_vmec = z1_i(:,:,ns)
      ru_vmec = ru_i(:,:,ns); zu_vmec(:,:) = zu_i(:,:,ns)
      rv_vmec = rv_i(:,:,ns); zv_vmec(:,:) = zv_i(:,:,ns)

!  First derivatives on s=1 surface, scaled by the edge factor in grid_extender
      ALLOCATE(r1s_vmec(ntheta,nzeta),z1s_vmec(ntheta,nzeta))
      ALLOCATE(r1su_vmec(ntheta,nzeta),z1su_vmec(ntheta,nzeta))
      r1s_vmec=(r1_vmec - r1_i(:,:,ns-1))*ohs
      z1s_vmec=(z1_vmec - z1_i(:,:,ns-1))*ohs
      r1su_vmec=(ru_vmec - ru_i(:,:,ns-1))*ohs
      z1su_vmec=(zu_vmec - zu_i(:,:,ns-1))*ohs

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Largest safe edge factor of the quadratic vessel extension.
!>
!>  The extension X(rho) = f*X_s*rho*(1 - rho) + X_0*(1 - rho^2) + X_w*rho^2
!>  matches the slope f*X_s at the plasma edge (rho = 0) and reaches the wall
!>  (rho = 1) with slope 2*(X_w - X_0) - f*X_s. Where f*X_s exceeds twice the
!>  plasma-wall gap the cross-section Jacobian X_rho x X_u changes sign and the
!>  grid folds. This returns the largest f = fac_max*k/nfac, k = nfac..1, for
!>  which the Jacobian keeps the sign of the plasma edge on 0 < rho <= 1, or 0
!>  if none does.
!>
!>  @param[in] fac_max Requested edge factor.
!>  @param[in] r0      R on the plasma edge.
!>  @param[in] z0      Z on the plasma edge.
!>  @param[in] r0u     dR/du on the plasma edge.
!>  @param[in] z0u     dZ/du on the plasma edge.
!>  @param[in] rs      dR/ds on the plasma edge.
!>  @param[in] zs      dZ/ds on the plasma edge.
!>  @param[in] rsu     d^2R/dsdu on the plasma edge.
!>  @param[in] zsu     d^2Z/dsdu on the plasma edge.
!>  @param[in] rw      R on the wall.
!>  @param[in] zw      Z on the wall.
!>  @param[in] rwu     dR/du on the wall.
!>  @param[in] zwu     dZ/du on the wall.
!>  @returns The edge factor.
!-------------------------------------------------------------------------------
      PURE FUNCTION ext_edge_fac_limit(fac_max, r0, z0, r0u, z0u, rs, zs,      &
                                       rsu, zsu, rw, zw, rwu, zwu)
      IMPLICIT NONE

!  Declare Arguments
      REAL(dp)                             :: ext_edge_fac_limit
      REAL(dp), INTENT(in)                 :: fac_max
      REAL(dp), DIMENSION(:,:), INTENT(in) :: r0, z0, r0u, z0u, rs, zs
      REAL(dp), DIMENSION(:,:), INTENT(in) :: rsu, zsu, rw, zw, rwu, zwu

!  Local Variables
      INTEGER                              :: k
      INTEGER                              :: i
      LOGICAL                              :: folds
      REAL(dp)                             :: f
      REAL(dp)                             :: rho
      REAL(dp)                             :: sgn
      REAL(dp), DIMENSION(SIZE(r0,1),SIZE(r0,2)) :: rr, zr, ru, zu

!  Local Parameters
      INTEGER, PARAMETER                   :: nfac = 50
      INTEGER, PARAMETER                   :: nrho = 100

!  Start of executable code.
      sgn = SIGN(1.0_dp, SUM(rs*z0u - r0u*zs))
      ext_edge_fac_limit = 0
      DO k = nfac, 1, -1
         f = fac_max*k/nfac
         folds = .false.
         DO i = 1, nrho
            rho = REAL(i, dp)/nrho
            rr = f*rs*(1 - 2*rho) + 2*rho*(rw - r0)
            zr = f*zs*(1 - 2*rho) + 2*rho*(zw - z0)
            ru = f*rsu*rho*(1 - rho) + r0u*(1 - rho*rho) + rwu*rho*rho
            zu = f*zsu*rho*(1 - rho) + z0u*(1 - rho*rho) + zwu*rho*rho
            IF (ANY(sgn*(rr*zu - ru*zr) .le. 0)) THEN
               folds = .true.
               EXIT
            END IF
         END DO
         IF (.not.folds) THEN
            ext_edge_fac_limit = f
            RETURN
         END IF
      END DO

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief UNKNOWN
!>
!>  @returns UNKNOWN
!-------------------------------------------------------------------------------
      FUNCTION GetVP1()

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp)                                ::  GetVP1

!  Local Variables
      INTEGER                                 :: js
      INTEGER                                 :: js1
      REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: r12
      REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: rs12
      REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: zs12
      REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: ru12
      REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: zu12
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: sqrtg_full
      REAL(dp), DIMENSION(:), ALLOCATABLE     :: vp

!  Start of executable code.
      ALLOCATE(vp(ns - 1:ns), sqrtg_full(ntheta,nzeta,ns - 1:ns))
      ALLOCATE(r12(ntheta,nzeta), ru12(ntheta,nzeta), rs12(ntheta,nzeta),      &
               zu12(ntheta,nzeta), zs12(ntheta,nzeta))

!  Guess assuming V(s) = V0*s**2:    vp1 = 2*vol0
!  NOTE: jacobian not computed yet

      DO js = ns - 1, ns
         js1 = js - 1
         r12 = p5*(r1_i(:,:,js) + r1_i(:,:,js1))
         rs12 = (r1_i(:,:,js) - r1_i(:,:,js1))*ohs
         zs12 = (z1_i(:,:,js) - z1_i(:,:,js1))*ohs
         ru12= p5*(ru_i(:,:,js) + ru_i(:,:,js1))
         zu12= p5*(zu_i(:,:,js) + zu_i(:,:,js1))
         sqrtg_full(:,:,js) = r12*(ru12*zs12 - rs12*zu12)
      END DO

      CALL SurfAvgLocal(vp, sqrtg_full, ns - 1, ns)
      vp = twopi*twopi*ABS(vp)
      GetVP1 = p5*(3*vp(ns) - vp(ns - 1))

      DEALLOCATE(vp, sqrtg_full, r12, ru12, rs12, zu12, zs12)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Find the new rho value correspond to the inital volume rho.
!>
!>  @params[inout] rho  Radial surface.
!>  @params[in]    vrho Volume inside rho.
!-------------------------------------------------------------------------------
      SUBROUTINE remap_radial(rho, vrho)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), INTENT(INOUT) :: rho
      REAL(dp), INTENT(IN)    :: vrho

!  Local variables
      INTEGER                 :: it
      REAL(dp)                :: rootl, rootr, vol, delv

!  Local parameters
      REAL(dp), PARAMETER     :: tol=1.E-6_dp

!  Start of executable code.
      rho = p5
      rootr = 1
      rootl = 0

      DO it = 1, 100
         vol = get_volume(rho)
         delv = vol-vrho
         IF (ABS(delv) .LT. tol*vrho) THEN
            EXIT
         END IF
         IF (delv .GT. 0) THEN
            rootr = rho
            rho = (rootl+rho)/2
         ELSE
            rootl = rho
            rho = rootl+(rootr-rho)/2
         END IF
      END DO
      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Get the volume at the surface rho.
!>
!>  Uses forumla volume(rho) = Integral[0.5*(R^2*Zu)] where integration is over
!>  u, v for a fixed rho. For simpplicity, we do the integral plane by plane.
!>  NOTE: This is not quite right.
!>
!>  @params[in] rho Radial surface.
!>  @returns The volume enclosed by surface rho.
!-------------------------------------------------------------------------------
      FUNCTION get_volume(rho)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp)             :: get_volume
      REAL(dp), INTENT(IN) :: rho

!  Local variables
      INTEGER              :: i
      INTEGER              :: j
      REAL(dp)             :: zu
      REAL(dp)             :: r1
      REAL(dp)             :: dnorm
      REAL(dp)             :: odelu

!  Start of executable code.
      get_volume = 0
      odelu = ntheta
      odelu = odelu/twopi
      dnorm = p5*twopi/(nzeta*odelu)
      IF (.NOT. lasym) THEN
         odelu = 2*odelu
      END IF
      CALL LoadSurfacePoints(rho)
      DO j = 1, nzeta
         DO i = 1, ntheta-1
            zu = (zi(i+1,j)-zi(i,j))*odelu
            r1 = p5*(ri(i+1,j)+ri(i,j))
            get_volume = get_volume + r1**2*zu
         END DO
      END DO

      IF (lasym) THEN
         DO j = 1, nzeta
            DO i = 1, ntheta-1
               zu = (zi(i+1,j)-zi(i,j))*odelu
               r1 = p5*(ri(i+1,j)+ri(i,j))
               get_volume = get_volume + r1**2*zu
            END DO
         END DO
      END IF

      get_volume = dnorm*get_volume

      END FUNCTION get_volume

!-------------------------------------------------------------------------------
!>  @brief Loads the theta and zeta values.
!>
!>  @param[in] rho Radial surface.
!-------------------------------------------------------------------------------
      SUBROUTINE LoadSurfacePoints(rho)
  
      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), INTENT(IN) :: rho

!  Local variables
      REAL(dp)             :: rho1
      REAL(dp)             :: rho2
      REAL(dp)             :: rho3

!  Start of executable code.
      rho1 = 1 - rho
      IF (interp_type .EQ. LIN) THEN
         ri = rho1*r1_vmec + rho*r1_ves
         zi = rho1*z1_vmec + rho*z1_ves
      ELSE
         rho2 = rho*rho
         IF (interp_type .EQ. QUAD) THEN
!  Quadratic (R' match at rho = 0)
            ri = p5*r1s_vmec*(1 - rho2 - rho1*rho1)                            &
     &         + r1_vmec*(1-rho2) + r1_ves*rho2
            zi = p5*z1s_vmec*(1 - rho2 - rho1*rho1)                            &
     &         + z1_vmec*(1-rho2) + z1_ves*rho2
         ELSE IF (interp_type .EQ. CUB) THEN
!  Cubic (R', R'' match at rho=0)
            rho3 = rho*rho2
            CALL ASSERT(.FALSE., 'Unknown interpolation <TYPE> specified')
         END IF
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Read the vessel file.
!>
!>  @returns Error status of the read.
!-------------------------------------------------------------------------------
      FUNCTION read_vessel_file()
      USE siesta_namelist, ONLY: vessel_file

      IMPLICIT NONE

!  Declare Arguments
      INTEGER   :: read_vessel_file

!  Local variables
      INTEGER   :: m
      INTEGER   :: n
      INTEGER   :: mn
      INTEGER   :: mnmax_v
      INTEGER   :: iou
      REAL (dp) :: rv1
      REAL (dp) :: zv1

!  Start of executable code.
      mpol_v = 0
      ntor_v = 0
      mnmax_v = 0

      iou = 36
      OPEN(UNIT=iou, FILE=vessel_file, STATUS="old", FORM="formatted",         &
           iostat=read_vessel_file)
      IF (read_vessel_file .NE. 0) THEN
         WRITE (*,1000) vessel_file
         RETURN
      END IF

      DO WHILE(.TRUE.)
         READ(iou,*,END=100) n, m, rv1, zv1
         mpol_v = MAX(m, mpol_v)
         ntor_v = MAX(ABS(n), ntor_v)
         mnmax_v = mnmax_v + 1
      END DO
  100 CONTINUE

      REWIND(UNIT=iou, iostat=read_vessel_file)
      IF (read_vessel_file .NE. 0 .and. lverbose) THEN
         WRITE (*,*) vessel_file
         RETURN
      END IF
      ALLOCATE(rmnc_v(0:mpol_v,-ntor_v:ntor_v),                                &
               zmns_v(0:mpol_v,-ntor_v:ntor_v))
      rmnc_v = 0
      zmns_v = 0
      IF (lasym) THEN
         ALLOCATE(rmns_v(0:mpol_v,-ntor_v:ntor_v),                             &
                  zmnc_v(0:mpol_v,-ntor_v:ntor_v))
         rmns_v = 0
         zmnc_v = 0
      END IF

      DO mn = 1, mnmax_v
         READ(iou,*,END=101) n, m, rv1, zv1
         CALL assert(m .ge. 0 .and. m .le. mpol_v,                             &
                     'm out of bounds: read_vessel_file')
         CALL ASSERT(n .ge. -ntor_v .and. n .le. ntor_v,                       &
                     'n out of bounds: read_vessel_file')
         rmnc_v(m,-n) = rv1
         zmns_v(m,-n) = zv1
      END DO
  101 CONTINUE
      CLOSE(UNIT=iou)

1000  FORMAT('Warning: Failed to open vessel file ',a)
1001  FORMAT('Warning: Failed to read vessel file ',a)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get surface average quantity.
!-------------------------------------------------------------------------------
      SUBROUTINE SurfAvgLocal(average, q3d, nsmin, nsmax)
      USE island_params, ONLY: fourier_context

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), INTENT(out)             :: average(nsmin:nsmax)
      REAL(dp), INTENT(in)              :: q3d(ntheta,nzeta,nsmin:nsmax)
      INTEGER, INTENT(in)               :: nsmin
      INTEGER, INTENT(in)               :: nsmax

!  Local Variables
      INTEGER                                 :: i
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: wint

!  Start of executable code.
      ALLOCATE(wint(ntheta,nzeta,nsmin:nsmax))

      DO i = 1, ntheta
         wint(i,:,:) = fourier_context%cosmui(i,0)
      END DO
      DO i = nsmin, nsmax
         average(i) = SUM(q3d(:,:,i)*wint(:,:,i))
      END DO

      DEALLOCATE(wint)

      END SUBROUTINE

!*******************************************************************************
!  UNIT TESTS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Test the edge factor limit of the vessel extension.
!>
!>  Circular plasma edge of minor radius a inside a concentric circular wall of
!>  radius a + d. The blend slope along the minor radius is f*a*(1 - 2*rho) +
!>  2*d*rho, so the grid folds at the wall exactly when f >= 2*d/a.
!>
!>  @returns True if the test passes.
!-------------------------------------------------------------------------------
      FUNCTION test_grid_extension()

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                         :: test_grid_extension

!  Local Variables
      INTEGER                         :: i
      REAL(dp), DIMENSION(32,1)       :: cu, su, r0, z0, r0u, z0u
      REAL(dp), DIMENSION(32,1)       :: rs, zs, rsu, zsu, rw, zw, rwu, zwu
      REAL(dp)                        :: fac

!  Local Parameters
      REAL(dp), PARAMETER             :: rmaj = 3.0_dp
      REAL(dp), PARAMETER             :: a = 0.5_dp
      REAL(dp), PARAMETER             :: d = 0.10125_dp

!  Start of executable code.
      DO i = 1, 32
         cu(i,1) = COS(twopi*(i - 1)/32)
         su(i,1) = SIN(twopi*(i - 1)/32)
      END DO
      r0 = rmaj + a*cu;        z0 = a*su
      r0u = -a*su;             z0u = a*cu
      rs = a*cu;               zs = a*su
      rsu = -a*su;             zsu = a*cu
      rw = rmaj + (a + d)*cu;  zw = (a + d)*su
      rwu = -(a + d)*su;       zwu = (a + d)*cu

!  The default 0.5 folds this grid; 2d/a = 0.405 lies between two steps.
      fac = ext_edge_fac_limit(0.5_dp, r0, z0, r0u, z0u, rs, zs, rsu, zsu,     &
                               rw, zw, rwu, zwu)
      test_grid_extension = ABS(fac - 0.40_dp) .lt. 1.0E-12_dp
      IF (.not.test_grid_extension) THEN
         WRITE (*,*) 'ext_edge_fac_limit(0.5): expected 0.40, got ', fac
      END IF

!  A factor below 2d/a is kept.
      fac = ext_edge_fac_limit(0.25_dp, r0, z0, r0u, z0u, rs, zs, rsu, zsu,    &
                               rw, zw, rwu, zwu)
      IF (fac .ne. 0.25_dp) THEN
         WRITE (*,*) 'ext_edge_fac_limit(0.25): expected 0.25, got ', fac
         test_grid_extension = .false.
      END IF

      END FUNCTION

      END MODULE
