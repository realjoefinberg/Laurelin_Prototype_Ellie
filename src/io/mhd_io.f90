module mhd_io
  use iso_fortran_env, only: real64
  use parameters
  use mhd_types
  use hdf5                ! Fortran HDF5 module
  implicit none

contains

  subroutine write_state(state, step)
    type(mhd_state), intent(in) :: state
    integer,        intent(in) :: step

    ! HDF5 handles
    integer(hid_t)         :: file_id, dset_id, space_id
    integer                :: error
    character(len=32)      :: fname
    integer(hsize_t), dimension(3) :: dims

    ! Build filename "state_0001.h5"
    write(fname, '("state_", I4.4, ".h5")') step

    ! Dimensions of the 3D grid
    dims = (/ Nkx, Nky, Nkz /)

    ! Initialize HDF5 library
    call h5open_f(error)

    ! Create (or overwrite) the file
    call h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, error)

    ! Create a simple dataspace
    call h5screate_simple_f(3, dims, space_id, error)

    ! --- Write density ---
    call h5dcreate_f(file_id, "rho", H5T_NATIVE_DOUBLE, space_id, &
                     dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,                         &
                    real(state%rho%data), dims, error)
    call h5dclose_f(dset_id, error)

    ! --- Write velocity components ---
    call h5dcreate_f(file_id, "vx", H5T_NATIVE_DOUBLE, space_id, dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, real(state%vx %data), dims, error)
    call h5dclose_f(dset_id, error)

    call h5dcreate_f(file_id, "vy", H5T_NATIVE_DOUBLE, space_id, dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, real(state%vy %data), dims, error)
    call h5dclose_f(dset_id, error)

    call h5dcreate_f(file_id, "vz", H5T_NATIVE_DOUBLE, space_id, dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, real(state%vz %data), dims, error)
    call h5dclose_f(dset_id, error)

    ! --- Write magnetic field components ---
    call h5dcreate_f(file_id, "bx", H5T_NATIVE_DOUBLE, space_id, dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, real(state%bx %data), dims, error)
    call h5dclose_f(dset_id, error)

    call h5dcreate_f(file_id, "by", H5T_NATIVE_DOUBLE, space_id, dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, real(state%by %data), dims, error)
    call h5dclose_f(dset_id, error)

    call h5dcreate_f(file_id, "bz", H5T_NATIVE_DOUBLE, space_id, dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, real(state%bz %data), dims, error)
    call h5dclose_f(dset_id, error)

    ! --- Write energy ---
    call h5dcreate_f(file_id, "E", H5T_NATIVE_DOUBLE, space_id, dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, real(state%E  %data), dims, error)
    call h5dclose_f(dset_id, error)

    ! Close dataspace and file
    call h5sclose_f(space_id, error)
    call h5fclose_f(file_id, error)
    call h5close_f(error)
  end subroutine write_state

end module mhd_io