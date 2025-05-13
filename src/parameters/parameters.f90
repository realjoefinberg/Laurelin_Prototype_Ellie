module parameters
  implicit none
  ! grid sizes
  integer, parameter :: Nkx = 64, Nky = 64, Nkz = 64
  real(kind=8), parameter :: Lx = 1d0, Ly = 1d0, Lz = 1d0
  real(kind=8), parameter :: dealiasing = 2d0/3d0
  ! time stepping
  real(kind=8), parameter :: CFL = 0.5d0
  real(kind=8) :: dt
  ! physics constants
  real(kind=8), parameter :: mu0 = 4d0*acos(-1d0)*1e-7
  real(kind=8), parameter :: rho0 = 1d0
contains
  subroutine compute_dt(max_speed)
    real(kind=8), intent(in) :: max_speed
    dt = CFL * min(Lx/Nkx, Ly/Nky, Lz/Nkz) / max_speed
  end subroutine compute_dt
end module parameters
