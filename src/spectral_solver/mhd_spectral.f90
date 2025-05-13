module mhd_spectral
  use parameters
  use mhd_types
  use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_loc
  implicit none

 ! Mathematical constants
 real(kind=8), parameter :: pi    = acos(-1d0)
 real(kind=8), parameter :: two_pi = 2d0 * pi

  ! FFTW constants
  integer(c_int), parameter :: FFTW_FORWARD  = -1
  integer(c_int), parameter :: FFTW_BACKWARD = +1
  integer(c_int), parameter :: FFTW_ESTIMATE = 64

  ! C‑binding interfaces to FFTW3
  interface
    function fftw_plan_dft_3d(n0,n1,n2,inbuf,outbuf,sign,flags) bind(C,name="fftw_plan_dft_3d")
      import :: c_int, c_ptr
      integer(c_int), value :: n0,n1,n2
      type(c_ptr),   value :: inbuf, outbuf
      integer(c_int), value :: sign, flags
      type(c_ptr)          :: fftw_plan_dft_3d
    end function fftw_plan_dft_3d

    subroutine fftw_execute_dft(plan,inbuf,outbuf) bind(C,name="fftw_execute_dft")
      import :: c_ptr
      type(c_ptr), value :: plan, inbuf, outbuf
    end subroutine fftw_execute_dft

    subroutine fftw_destroy_plan(plan) bind(C,name="fftw_destroy_plan")
      import :: c_ptr
      type(c_ptr), value :: plan
    end subroutine fftw_destroy_plan
  end interface

  ! FFTW plans & buffer
  type(c_ptr)                        :: plan_fwd, plan_inv
  complex(kind=8), allocatable,target :: work(:,:,:)

  ! Real‑space arrays
  complex(kind=8), allocatable :: vx_r(:,:,:), vy_r(:,:,:), vz_r(:,:,:)
  complex(kind=8), allocatable :: bx_r(:,:,:), by_r(:,:,:), bz_r(:,:,:)

  ! Nonlinear RHS in real space
  complex(kind=8), allocatable :: nl_vx(:,:,:), nl_vy(:,:,:), nl_vz(:,:,:)
  complex(kind=8), allocatable :: nl_bx(:,:,:), nl_by(:,:,:), nl_bz(:,:,:)

contains

  subroutine setup_fftw_spectral()
    integer :: n(3), flags
    n     = [Nkx,Nky,Nkz]
    flags = FFTW_ESTIMATE
    allocate(work(n(1),n(2),n(3)))
    plan_fwd = fftw_plan_dft_3d(n(1),n(2),n(3), c_loc(work), c_loc(work), &
                                FFTW_FORWARD,  flags)
    plan_inv = fftw_plan_dft_3d(n(1),n(2),n(3), c_loc(work), c_loc(work), &
                                FFTW_BACKWARD, flags)
    ! Allocate real & nonlinear arrays
    allocate(vx_r(n(1),n(2),n(3)), vy_r(n(1),n(2),n(3)), vz_r(n(1),n(2),n(3)))
    allocate(bx_r(n(1),n(2),n(3)), by_r(n(1),n(2),n(3)), bz_r(n(1),n(2),n(3)))
    allocate(nl_vx(n(1),n(2),n(3)), nl_vy(n(1),n(2),n(3)), nl_vz(n(1),n(2),n(3)))
    allocate(nl_bx(n(1),n(2),n(3)), nl_by(n(1),n(2),n(3)), nl_bz(n(1),n(2),n(3)))
  end subroutine setup_fftw_spectral

  subroutine compute_nonlinear(state,rhs)
    type(mhd_state), intent(in)  :: state
    type(mhd_state), intent(out) :: rhs
    integer :: i,j,k,ip,im,jp,jm,kp,km
    real(kind=8) :: dx,dy,dz

    dx = Lx/Nkx; dy = Ly/Nky; dz = Lz/Nkz

    !— Inverse FFT k→real ——
    work = state%vx%data; call fftw_execute_dft(plan_inv,c_loc(work),c_loc(work)); vx_r = work
    work = state%vy%data; call fftw_execute_dft(plan_inv,c_loc(work),c_loc(work)); vy_r = work
    work = state%vz%data; call fftw_execute_dft(plan_inv,c_loc(work),c_loc(work)); vz_r = work
    work = state%bx%data; call fftw_execute_dft(plan_inv,c_loc(work),c_loc(work)); bx_r = work
    work = state%by%data; call fftw_execute_dft(plan_inv,c_loc(work),c_loc(work)); by_r = work
    work = state%bz%data; call fftw_execute_dft(plan_inv,c_loc(work),c_loc(work)); bz_r = work

    !— Compute nonlinear terms with central differences ——
    do i=1,Nkx
      ip = mod(i  ,Nkx)+1;  im = mod(i-2,Nkx)+1
      do j=1,Nky
        jp = mod(j  ,Nky)+1;  jm = mod(j-2,Nky)+1
        do k=1,Nkz
          kp = mod(k  ,Nkz)+1;  km = mod(k-2,Nkz)+1

          ! derivatives
          nl_vx(i,j,k) = - ( vx_r(i,j,k)*(vx_r(ip,j,k)-vx_r(im,j,k))/(2*dx) &
                           + vy_r(i,j,k)*(vx_r(i,jp,k)-vx_r(i,jm,k))/(2*dy) &
                           + vz_r(i,j,k)*(vx_r(i,j,kp)-vx_r(i,j,km))/(2*dz) )
          nl_vy(i,j,k) = - ( vx_r(i,j,k)*(vy_r(ip,j,k)-vy_r(im,j,k))/(2*dx) &
                           + vy_r(i,j,k)*(vy_r(i,jp,k)-vy_r(i,jm,k))/(2*dy) &
                           + vz_r(i,j,k)*(vy_r(i,j,kp)-vy_r(i,j,km))/(2*dz) )
          nl_vz(i,j,k) = - ( vx_r(i,j,k)*(vz_r(ip,j,k)-vz_r(im,j,k))/(2*dx) &
                           + vy_r(i,j,k)*(vz_r(i,jp,k)-vz_r(i,jm,k))/(2*dy) &
                           + vz_r(i,j,k)*(vz_r(i,j,kp)-vz_r(i,j,km))/(2*dz) )

          ! magnetic tension term: B·∇B
          nl_bx(i,j,k) =   bx_r(i,j,k)*(bx_r(ip,j,k)-bx_r(im,j,k))/(2*dx) &
                        + by_r(i,j,k)*(bx_r(i,jp,k)-bx_r(i,jm,k))/(2*dy) &
                        + bz_r(i,j,k)*(bx_r(i,j,kp)-bx_r(i,j,km))/(2*dz)

          nl_by(i,j,k) =   bx_r(i,j,k)*(by_r(ip,j,k)-by_r(im,j,k))/(2*dx) &
                        + by_r(i,j,k)*(by_r(i,jp,k)-by_r(i,jm,k))/(2*dy) &
                        + bz_r(i,j,k)*(by_r(i,j,kp)-by_r(i,j,km))/(2*dz)

          nl_bz(i,j,k) =   bx_r(i,j,k)*(bz_r(ip,j,k)-bz_r(im,j,k))/(2*dx) &
                        + by_r(i,j,k)*(bz_r(i,jp,k)-bz_r(i,jm,k))/(2*dy) &
                        + bz_r(i,j,k)*(bz_r(i,j,kp)-bz_r(i,j,km))/(2*dz)
        end do
      end do
    end do

    !— Forward FFT real→k & normalize ——
    work = nl_vx; call fftw_execute_dft(plan_fwd,c_loc(work),c_loc(work)); rhs%vx%data = work/(Nkx*Nky*Nkz)
    work = nl_vy; call fftw_execute_dft(plan_fwd,c_loc(work),c_loc(work)); rhs%vy%data = work/(Nkx*Nky*Nkz)
    work = nl_vz; call fftw_execute_dft(plan_fwd,c_loc(work),c_loc(work)); rhs%vz%data = work/(Nkx*Nky*Nkz)
    work = nl_bx; call fftw_execute_dft(plan_fwd,c_loc(work),c_loc(work)); rhs%bx%data = work/(Nkx*Nky*Nkz)
    work = nl_by; call fftw_execute_dft(plan_fwd,c_loc(work),c_loc(work)); rhs%by%data = work/(Nkx*Nky*Nkz)
    work = nl_bz; call fftw_execute_dft(plan_fwd,c_loc(work),c_loc(work)); rhs%bz%data = work/(Nkx*Nky*Nkz)

    ! copy density & energy through
    rhs%rho%data = state%rho%data
    rhs%E  %data = state%E  %data

    !— Dealias 2/3 rule ——
    call dealias(rhs)

    !— Divergence‐free projection ——
    call project(rhs)
  end subroutine compute_nonlinear

  !────────────────────────────────────────────────────────────────────────────
  subroutine destroy_fftw_spectral()
    call fftw_destroy_plan(plan_fwd)
    call fftw_destroy_plan(plan_inv)
    deallocate(work, vx_r,vy_r,vz_r, bx_r,by_r,bz_r)
    deallocate(nl_vx,nl_vy,nl_vz,nl_bx,nl_by,nl_bz)
  end subroutine destroy_fftw_spectral

  !— Dealias helper: zero modes beyond 2/3 kmax ——
  subroutine dealias(rhs)
    type(mhd_state), intent(inout) :: rhs
    integer :: i,j,k
    real(kind=8) :: kx,ky,kz, kxmax, kymax, kzmax
    kxmax = 2d0*pi*(Nkx/2)/Lx * 2d0/3d0
    kymax = 2d0*pi*(Nky/2)/Ly * 2d0/3d0
    kzmax = 2d0*pi*(Nkz/2)/Lz * 2d0/3d0
    do i=1,Nkx
      kx = 2d0*pi*merge(real(i-1,8), real(i-1,8)-real(Nkx,8), i>Nkx/2)/Lx
      do j=1,Nky
        ky = 2d0*pi*merge(real(j-1,8), real(j-1,8)-real(Nky,8), j>Nky/2)/Ly
        do k=1,Nkz
          kz = 2d0*pi*merge(real(k-1,8), real(k-1,8)-real(Nkz,8), k>Nkz/2)/Lz
          if (abs(kx)>kxmax .or. abs(ky)>kymax .or. abs(kz)>kzmax) then
            rhs%vx%data(i,j,k) = (0d0,0d0)
            rhs%vy%data(i,j,k) = (0d0,0d0)
            rhs%vz%data(i,j,k) = (0d0,0d0)
            rhs%bx%data(i,j,k) = (0d0,0d0)
            rhs%by%data(i,j,k) = (0d0,0d0)
            rhs%bz%data(i,j,k) = (0d0,0d0)
          endif
        end do
      end do
    end do
  end subroutine dealias

  !— Project onto divergence‐free subspace ——
  subroutine project(rhs)
    type(mhd_state), intent(inout) :: rhs
    integer :: i,j,k
    real(kind=8) :: kx,ky,kz,k2
    complex(kind=8) :: dotvk
    do i=1,Nkx
      kx = 2d0*pi*merge(real(i-1,8), real(i-1,8)-real(Nkx,8), i>Nkx/2)/Lx
      do j=1,Nky
        ky = 2d0*pi*merge(real(j-1,8), real(j-1,8)-real(Nky,8), j>Nky/2)/Ly
        do k=1,Nkz
          kz = 2d0*pi*merge(real(k-1,8), real(k-1,8)-real(Nkz,8), k>Nkz/2)/Lz
          k2 = kx*kx + ky*ky + kz*kz
          if (k2>0d0) then
            dotvk = kx*rhs%vx%data(i,j,k) + ky*rhs%vy%data(i,j,k) + kz*rhs%vz%data(i,j,k)
            rhs%vx%data(i,j,k) = rhs%vx%data(i,j,k) - dotvk*kx/k2
            rhs%vy%data(i,j,k) = rhs%vy%data(i,j,k) - dotvk*ky/k2
            rhs%vz%data(i,j,k) = rhs%vz%data(i,j,k) - dotvk*kz/k2

            dotvk = kx*rhs%bx%data(i,j,k) + ky*rhs%by%data(i,j,k) + kz*rhs%bz%data(i,j,k)
            rhs%bx%data(i,j,k) = rhs%bx%data(i,j,k) - dotvk*kx/k2
            rhs%by%data(i,j,k) = rhs%by%data(i,j,k) - dotvk*ky/k2
            rhs%bz%data(i,j,k) = rhs%bz%data(i,j,k) - dotvk*kz/k2
          endif
        end do
      end do
    end do
  end subroutine project

end module mhd_spectral
