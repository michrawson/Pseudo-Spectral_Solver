function eone ( x )

!*****************************************************************************80
!
!! EONE evaluates the exponential integral E1(X).
!
!  Discussion:
!
!    This routine computes approximate values for the
!    exponential integral E1(x), where x is real.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) EONE, the value of the function.
!
  implicit none

  real ( kind = 8 ) eone
  integer ( kind = 4 ) jint
  real ( kind = 8 ) result
  real ( kind = 8 ) x

  jint = 2
  call calcei ( x, result, jint )
  eone = result

  return
end
