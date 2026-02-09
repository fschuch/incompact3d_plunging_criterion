! Stub module to provide compatibility for legacy USE var statements
! This module re-exports variables from the variables module and adds
! missing array declarations

module var
  
  use variables
  use decomp_2d, only : mytype, xsize, ysize, zsize
  
  ! Re-export all from variables module
  
  ! Working arrays for velocity, pressure, and other fields
  ! These are allocated dynamically during initialization
  real(mytype), allocatable, dimension(:,:,:) :: ux1, uy1, uz1  ! velocity in x-pencil
  real(mytype), allocatable, dimension(:,:,:) :: ux2, uy2, uz2  ! velocity in y-pencil
  real(mytype), allocatable, dimension(:,:,:) :: ux3, uy3, uz3  ! velocity in z-pencil
  real(mytype), allocatable, dimension(:,:,:) :: pre1, ep1      ! pressure and volume fraction in x-pencil
  real(mytype), allocatable, dimension(:,:,:) :: ta1, tb1, tc1, td1, te1, tf1, tg1, th1, ti1  ! temp arrays x-pencil
  real(mytype), allocatable, dimension(:,:,:) :: ta2, tb2, tc2, td2, te2, tf2, tg2, th2, ti2, tj2  ! temp arrays y-pencil
  real(mytype), allocatable, dimension(:,:,:) :: ta3, tb3, tc3, td3, te3, tf3, tg3, th3, ti3  ! temp arrays z-pencil
  real(mytype), allocatable, dimension(:,:,:) :: di1, di2, di3  ! work arrays
  real(mytype), allocatable, dimension(:,:,:) :: uvisu  ! visualization arrays
  
  ! Scalar fields (passive scalars, concentration, etc.)
  real(mytype), allocatable, dimension(:,:,:,:) :: phi1, phi2, phi3  ! scalar fields in different pencils
  real(mytype), allocatable, dimension(:,:,:,:) :: phis1, phiss1  ! scalar fields derivatives
  
  ! Gradient and other field arrays
  real(mytype), allocatable, dimension(:,:,:) :: gx1, gy1, gz1  ! gradients in x-pencil
  real(mytype), allocatable, dimension(:,:,:) :: hx1, hy1, hz1  ! helper arrays in x-pencil
  real(mytype), allocatable, dimension(:,:,:) :: px1, py1, pz1  ! pressure gradients
  real(mytype), allocatable, dimension(:,:,:) :: pp3  ! pressure in z-pencil
  
  ! Dissipation and other derived quantities
  real(mytype), allocatable, dimension(:,:,:) :: diss1, dv3  ! dissipation arrays
  real(mytype), allocatable, dimension(:,:,:) :: nut1  ! turbulent viscosity
  
  ! Mesh size arrays
  integer, allocatable, dimension(:) :: nxmsize, nymsize, nzmsize
  
end module var
