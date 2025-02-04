module solid_volume

  use, intrinsic :: iso_fortran_env

  implicit none

  type :: point
    real(REAL64) :: x
    real(REAL64) :: y
    real(REAL64) :: z
  end type point

  type :: plane
    real(REAL64) :: a
    real(REAL64) :: b
    real(REAL64) :: c
    real(REAL64) :: d
  end type plane

contains

!==============================================================================

subroutine plane_by_3_points(p1,p2,p3,pln)

  type(point), intent(in) :: p1
  type(point), intent(in) :: p2
  type(point), intent(in) :: p3
  type(plane), intent(out) :: pln

  type(point) :: v1 ! p2 - p1
  type(point) :: v2 ! p3 - p1
  real(REAL64) :: alpha_
  real(REAL64) :: beta_
  real(REAL64) :: gamma_

  v1%x = p2%x - p1%x
  v1%y = p2%y - p1%y
  v1%z = p2%z - p1%z

  v2%x = p3%x - p1%x
  v2%y = p3%y - p1%y
  v2%z = p3%z - p1%z

  alpha_ = v1%y*v2%z - v1%z*v2%y
  beta_  = v1%x*v2%z - v1%z*v2%x
  gamma_ = v1%x*v2%y - v1%y*v2%x

  pln%a = alpha_
  pln%b = -beta_
  pln%c = gamma_
  pln%d = -alpha_*p1%x + beta_*p1%y - gamma_*p1%z

end subroutine plane_by_3_points

!==============================================================================

real(REAL64) function plane_point_distance(pln,p)

  type(plane), intent(in) :: pln
  type(point), intent(in) :: p

  real(REAL64) :: n
  real(REAL64) :: d

  n = abs(pln%a*p%x + pln%b*p%y + pln%c*p%z + pln%d)
  d = sqrt(pln%a**2 + pln%b**2 + pln%c**2)

  plane_point_distance = n/d

end function plane_point_distance

!==============================================================================

real(REAL64) function triangle_area(p1,p2,p3)

  type(point), intent(in) :: p1
  type(point), intent(in) :: p2
  type(point), intent(in) :: p3

  type(point) :: v1 ! p2 - p1
  type(point) :: v2 ! p3 - p1
  real(REAL64) :: alpha_
  real(REAL64) :: beta_
  real(REAL64) :: gamma_

  v1%x = p2%x - p1%x
  v1%y = p2%y - p1%y
  v1%z = p2%z - p1%z

  v2%x = p3%x - p1%x
  v2%y = p3%y - p1%y
  v2%z = p3%z - p1%z

  alpha_ = v1%y*v2%z - v1%z*v2%y
  beta_  = v1%x*v2%z - v1%z*v2%x
  gamma_ = v1%x*v2%y - v1%y*v2%x

  triangle_area = sqrt(alpha_**2 + beta_**2 + gamma_**2) / 2.0_REAL64

end function triangle_area

!==============================================================================

real(REAL64) function triangular_pyramid_volume(p1,p2,p3,p4)

  type(point), intent(in) :: p1
  type(point), intent(in) :: p2
  type(point), intent(in) :: p3
  type(point), intent(in) :: p4

  type(plane) :: pln
  real(REAL64) :: area
  real(REAL64) :: height

  area = triangle_area(p1,p2,p3)
  call plane_by_3_points(p1,p2,p3,pln)
  height = plane_point_distance(pln,p4)
  triangular_pyramid_volume = area * height / 3.0_REAL64

end function triangular_pyramid_volume

!==============================================================================

real(REAL64) function octahedron_volume(p1,p2,p3,p4,p5,p6)

!             .p1.
!           .'/  \'.
!         .' /    \ '.
!       p2--/------\--p5
!       | \/        \/ |
!       | /\        /\ |
!       p3------------p4
!         '. \    / .'
!           '.\  /.'
!             'p6'

  type(point), intent(in) :: p1
  type(point), intent(in) :: p2
  type(point), intent(in) :: p3
  type(point), intent(in) :: p4
  type(point), intent(in) :: p5
  type(point), intent(in) :: p6

  real(REAL64) :: v1
  real(REAL64) :: v2
  real(REAL64) :: v3
  real(REAL64) :: v4

  v1 = triangular_pyramid_volume(p2,p3,p4,p1)
  v2 = triangular_pyramid_volume(p2,p3,p4,p6)
  v3 = triangular_pyramid_volume(p2,p4,p5,p1)
  v4 = triangular_pyramid_volume(p2,p4,p5,p6)

  octahedron_volume = v1 + v2 + v3 + v4

end function octahedron_volume

!==============================================================================

end module solid_volume
