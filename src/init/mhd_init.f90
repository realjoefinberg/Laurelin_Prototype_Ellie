module mhd_init
  use parameters
  use mhd_types
  implicit none

  ! Real‑space arrays for initialization
  complex(kind=8), allocatable :: vx_real(:,:,:), vy_real(:,:,:)
  complex(kind=8), allocatable :: bx_real(:,:,:), bz_real(:,:,:)

  real(kind=8), parameter :: two_pi = 2d0 * acos(-1d0)

contains

  !──────────────────────────────────────────────────────────────────────────
  subroutine setup_fftw_init()
    ! No-op stub for FFTW plan setup in initialization
  end subroutine setup_fftw_init
  !──────────────────────────────────────────────────────────────────────────

  subroutine initialize(state)
    type(mhd_state), intent(out) :: state
    integer :: ix, iy, iz
    real(kind=8) :: x, y, z

    ! 1) Allocate k‑space fields
    call allocate_fields(state)

    ! 2) Allocate real‑space temporaries
    allocate(vx_real(Nkx,Nky,Nkz))
    allocate(vy_real(Nkx,Nky,Nkz))
    allocate(bx_real(Nkx,Nky,Nkz))
    allocate(bz_real(Nkx,Nky,Nkz))

    ! 3) Fill Orszag–Tang vortex in real space
    do ix = 1, Nkx
      x = (ix-1) * Lx / Nkx
      do iy = 1, Nky
        y = (iy-1) * Ly / Nky
        do iz = 1, Nkz
          z = (iz-1) * Lz / Nkz

          vx_real(ix,iy,iz) = cmplx(-sin(two_pi * y / Ly), 0d0)
          vy_real(ix,iy,iz) = cmplx( sin(two_pi * x / Lx), 0d0)
          bx_real(ix,iy,iz) = cmplx( sin(2*two_pi * y / Ly), 0d0)
          bz_real(ix,iy,iz) = cmplx( sin(    two_pi * x / Lx), 0d0)
        end do
      end do
    end do

    ! 4) Uniform density & zero energy in k‑space
    state%rho%data = rho0
    state%E  %data = (0d0, 0d0)

    ! 5) **BYPASS FFTW**: copy real→k directly
    state%vx%data = vx_real
    state%vy%data = vy_real
    state%bx%data = bx_real
    state%bz%data = bz_real

    ! 6) Free real‑space arrays
    deallocate(vx_real, vy_real, bx_real, bz_real)
  end subroutine initialize

end module mhd_init
