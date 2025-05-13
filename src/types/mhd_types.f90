module mhd_types
  use parameters
  implicit none

  !––– Field & state definitions –––
  type :: field3d
    complex(kind=8), allocatable :: data(:,:,:)
  end type field3d

  type :: mhd_state
    type(field3d) :: rho, vx, vy, vz, bx, by, bz, E
  end type mhd_state

  !––– Overload “+” for mhd_state –––
  interface operator(+)
    module procedure add_state
  end interface operator(+)

  !––– Overload “*” (real * state) –––
  interface operator(*)
    module procedure scale_state
  end interface operator(*)

contains

  !––– Allocate each component array –––
  subroutine allocate_fields(state)
    type(mhd_state), intent(inout) :: state
    integer :: nx, ny, nz
    nx = Nkx; ny = Nky; nz = Nkz
    allocate(state%rho%data (nx,ny,nz))
    allocate(state%vx%data  (nx,ny,nz))
    allocate(state%vy%data  (nx,ny,nz))
    allocate(state%vz%data  (nx,ny,nz))
    allocate(state%bx%data  (nx,ny,nz))
    allocate(state%by%data  (nx,ny,nz))
    allocate(state%bz%data  (nx,ny,nz))
    allocate(state%E %data  (nx,ny,nz))
  end subroutine allocate_fields

  !––– Addition of two states –––
  function add_state(a, b) result(c)
    type(mhd_state), intent(in) :: a, b
    type(mhd_state)            :: c

    call allocate_fields(c)

    c%rho%data = a%rho%data + b%rho%data
    c%vx %data = a%vx %data + b%vx %data
    c%vy %data = a%vy %data + b%vy %data
    c%vz %data = a%vz %data + b%vz %data
    c%bx %data = a%bx %data + b%bx %data
    c%by %data = a%by %data + b%by %data
    c%bz %data = a%bz %data + b%bz %data
    c%E  %data = a%E  %data + b%E  %data
  end function add_state

  !––– Scale a state by a real –––
  function scale_state(alpha, a) result(c)
    real(kind=8), intent(in)      :: alpha
    type(mhd_state), intent(in)   :: a
    type(mhd_state)               :: c

    call allocate_fields(c)

    c%rho%data = alpha * a%rho%data
    c%vx %data = alpha * a%vx %data
    c%vy %data = alpha * a%vy %data
    c%vz %data = alpha * a%vz %data
    c%bx %data = alpha * a%bx %data
    c%by %data = alpha * a%by %data
    c%bz %data = alpha * a%bz %data
    c%E  %data = alpha * a%E  %data
  end function scale_state

end module mhd_types
