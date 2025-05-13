module mhd_time_integration
  use parameters
  use mhd_types
  use mhd_spectral
  implicit none

contains

  subroutine advance_rk3(state)
    type(mhd_state), intent(inout) :: state
    type(mhd_state) :: k1, k2, k3
    integer :: nx, ny, nz

    ! Determine grid dims once
    nx = Nkx; ny = Nky; nz = Nkz

    ! Allocate the three RHS states
    call allocate_fields(k1)
    call allocate_fields(k2)
    call allocate_fields(k3)

    !── Stage 1 ───────────────────────────────────────────────────────────────
    call compute_nonlinear(state, k1)
    ! state = state + dt*k1
    state%rho%data = state%rho%data + dt * k1%rho%data
    state%vx %data = state%vx %data + dt * k1%vx %data
    state%vy %data = state%vy %data + dt * k1%vy %data
    state%vz %data = state%vz %data + dt * k1%vz %data
    state%bx %data = state%bx %data + dt * k1%bx %data
    state%by %data = state%by %data + dt * k1%by %data
    state%bz %data = state%bz %data + dt * k1%bz %data
    state%E  %data = state%E  %data + dt * k1%E  %data

    !── Stage 2 ───────────────────────────────────────────────────────────────
    call compute_nonlinear(state, k2)
    ! state = 3/4*state + 1/4*(state + dt*k2)
    state%rho%data = 0.75d0*state%rho%data + 0.25d0*(state%rho%data + dt*k2%rho%data)
    state%vx %data = 0.75d0*state%vx %data + 0.25d0*(state%vx %data + dt*k2%vx %data)
    state%vy %data = 0.75d0*state%vy %data + 0.25d0*(state%vy %data + dt*k2%vy %data)
    state%vz %data = 0.75d0*state%vz %data + 0.25d0*(state%vz %data + dt*k2%vz %data)
    state%bx %data = 0.75d0*state%bx %data + 0.25d0*(state%bx %data + dt*k2%bx %data)
    state%by %data = 0.75d0*state%by %data + 0.25d0*(state%by %data + dt*k2%by %data)
    state%bz %data = 0.75d0*state%bz %data + 0.25d0*(state%bz %data + dt*k2%bz %data)
    state%E  %data = 0.75d0*state%E  %data + 0.25d0*(state%E  %data + dt*k2%E  %data)

    !── Stage 3 ───────────────────────────────────────────────────────────────
    call compute_nonlinear(state, k3)
    ! state = 1/3*state + 2/3*(state + dt*k3)
    state%rho%data = (1d0/3d0)*state%rho%data + (2d0/3d0)*(state%rho%data + dt*k3%rho%data)
    state%vx %data = (1d0/3d0)*state%vx %data + (2d0/3d0)*(state%vx %data + dt*k3%vx %data)
    state%vy %data = (1d0/3d0)*state%vy %data + (2d0/3d0)*(state%vy %data + dt*k3%vy %data)
    state%vz %data = (1d0/3d0)*state%vz %data + (2d0/3d0)*(state%vz %data + dt*k3%vz %data)
    state%bx %data = (1d0/3d0)*state%bx %data + (2d0/3d0)*(state%bx %data + dt*k3%bx %data)
    state%by %data = (1d0/3d0)*state%by %data + (2d0/3d0)*(state%by %data + dt*k3%by %data)
    state%bz %data = (1d0/3d0)*state%bz %data + (2d0/3d0)*(state%bz %data + dt*k3%bz %data)
    state%E  %data = (1d0/3d0)*state%E  %data + (2d0/3d0)*(state%E  %data + dt*k3%E  %data)
  end subroutine advance_rk3

end module mhd_time_integration