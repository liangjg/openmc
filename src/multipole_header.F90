module multipole_header

  use hdf5

  use constants
  use dict_header,      only: DictIntInt
  use error,            only: fatal_error
  use hdf5_interface

  implicit none

  !========================================================================
  ! Multipole related constants

  integer, parameter :: MAX_CF_ORDER = 10

  ! Constants that determine which value to access
  integer, parameter :: MP_EA = 1       ! Pole

  ! Reich-Moore indices
  integer, parameter :: MP_RT = 2, &    ! Residue total
                        MP_RA = 3, &    ! Residue absorption
                        MP_RF = 4       ! Residue fission

  ! Polynomial fit indices
  integer, parameter :: CF_RT = 1, &    ! Total
                        CF_RA = 2, &    ! Absorption
                        CF_RF = 3       ! Fission

!===============================================================================
! MULTIPOLE contains all the components needed for the windowed multipole
! temperature dependent cross section libraries for the resolved resonance
! region.
!===============================================================================

  type MultipoleArray

    logical                 :: fissionable     ! Is this isotope fissionable?
    real(8)                 :: sqrtAWR         ! Square root of the atomic
    real(8)                 :: start_E         ! Start energy for the windows
    real(8)                 :: end_E           ! End energy for the windows
    integer                 :: n_windows
    real(8), allocatable    :: w_grid(:)       ! Contains the fitting function.
    integer, allocatable    :: grid_index(:)   ! Contains the fitting function.
    integer, allocatable    :: mp_offsets(:)   ! Contains the index of the pole at
    integer, allocatable    :: cf_offsets(:)   ! Contains the index of the curvefit at
    complex(8), allocatable :: mp_data(:,:)    ! Poles and residues
    real(8), allocatable    :: curvefit(:,:)   ! Contains the fitting function.

  contains

    procedure :: from_hdf5 => multipole_from_hdf5

  end type MultipoleArray

contains

!===============================================================================
! FROM_HDF5 loads multipole data from an HDF5 file.
!===============================================================================

  subroutine multipole_from_hdf5(this, filename)
    class(MultipoleArray), intent(inout) :: this
    character(len=*),      intent(in)    :: filename

    character(len=10) :: version
    integer :: i, n_poles, n_residue_types
    integer(HSIZE_T) :: dims_1d(1), dims_2d(2)
    integer(HID_T) :: file_id
    integer(HID_T) :: group_id
    integer(HID_T) :: dset

    ! Open file for reading and move into the /isotope group
    file_id = file_open(filename, 'r', parallel=.true.)
    group_id = open_group(file_id, "/nuclide")

    ! Check the file version number.
    call read_dataset(version, file_id, "version")
    if (version /= VERSION_MULTIPOLE) call fatal_error("The current multipole&
         & format version is " // trim(VERSION_MULTIPOLE) // " but the file "&
         // trim(filename) // " uses version " // trim(version) // ".")

    ! Read scalar values.
    call read_dataset(this % sqrtAWR, group_id, "sqrtAWR")

    ! read w_grid and determine n_windows
    dset = open_dataset(group_id, "w_grid")
    call get_shape(dset, dims_1d)
    this % n_windows = int(dims_1d(1), 4) - 1
    allocate(this % w_grid(this % n_windows+1))
    call read_dataset(this % w_grid, dset)
    call close_dataset(dset)

    this % start_E = this % w_grid(1)
    this % end_E = this % w_grid(size(this % w_grid))

    allocate(this % mp_offsets(this % n_windows+1))
    call read_dataset(this % mp_offsets, group_id, "mp_offsets")
    allocate(this % cf_offsets(this % n_windows+1))
    call read_dataset(this % cf_offsets, group_id, "cf_offsets")

    ! Read the multipole and curvefit array
    dset = open_dataset(group_id, "mp_data")
    call get_shape(dset, dims_2d)
    n_residue_types = int(dims_2d(1), 4) - 1
    n_poles = int(dims_2d(2), 4)
    allocate(this % mp_data(n_residue_types+1, n_poles))
    call read_dataset(this % mp_data, dset)
    call close_dataset(dset)

    ! Check to see if this data includes fission residues.
    this % fissionable = (n_residue_types == 3)

    ! Read the "curvefit" array.
    dset = open_dataset(group_id, "curvefit")
    call get_shape(dset, dims_2d)
    if (n_residue_types /= int(dims_2d(1), 4)) call fatal_error("reaction &
         &types in curvefit is inconsistent with mp_data")
    allocate(this % curvefit(dims_2d(1), dims_2d(2)))
    call read_dataset(this % curvefit, dset)
    call close_dataset(dset)

    ! Close the group and file.
    call close_group(group_id)
    call file_close(file_id)
  end subroutine multipole_from_hdf5
end module multipole_header
