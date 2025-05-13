program mhd3d
  !── Import modules & symbols ─────────────────────────────────────────────
  use parameters,           only: compute_dt, Lx, Ly, Lz, Nkx, Nky, Nkz, mu0
  use mhd_types,            only: mhd_state
  use mhd_init,             only: setup_fftw_init, initialize
  use mhd_spectral,         only: setup_fftw_spectral, compute_nonlinear
  use mhd_time_integration, only: advance_rk3
  use mhd_io,               only: write_state

  implicit none

  !── Declarations ────────────────────────────────────────────────────────
  type(mhd_state) :: state
  integer         :: nsteps, i, csvUnit
  real(kind=8)    :: cellVol, Ekin, Emag

  !── Simulation settings ─────────────────────────────────────────────────
  nsteps = 1000   ! or read from command‑line / config

  !── Compute cell volume once ────────────────────────────────────────────
  cellVol = (Lx / real(Nkx,8)) * (Ly / real(Nky,8)) * (Lz / real(Nkz,8))

  !── Open CSV for energy diagnostics ─────────────────────────────────────
  open(newunit=csvUnit, file="energy.csv", status="replace", action="write")
  write(csvUnit,'("step, E_kinetic, E_magnetic")')

  !── Build FFTW plans in init & spectral modules ─────────────────────────
  call setup_fftw_init()
  call setup_fftw_spectral()

  !── Initialize the field & compute dt ───────────────────────────────────
  call initialize(state)
  call compute_dt(1d0)
  !── Main RK3 time‑stepping loop with diagnostics ────────────────────────
  integer, parameter :: write_stride = 10

  do i = 1, nsteps
    call advance_rk3(state)

    !── Compute energies each step ──────────────────────────────────────
    Ekin = 0.5d0 * sum( real(state%rho%data) * &
             ( abs(state%vx%data)**2 + abs(state%vy%data)**2 + abs(state%vz%data)**2 ) ) * cellVol

    Emag = 0.5d0/mu0 * sum( abs(state%bx%data)**2 + abs(state%by%data)**2 + abs(state%bz%data)**2 ) * cellVol

    !── Write energies every step to CSV ───────────────────────────────
    write(csvUnit,'(I8,",",F20.8,",",F20.8)') i, Ekin, Emag

    !── Write full state only every write_stride steps ────────────────
    if (mod(i, write_stride) == 0) then
      call write_state(state, i)
    end if

  end do


  close(csvUnit)
end program mhd3d
