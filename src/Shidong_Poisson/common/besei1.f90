function besei1 ( x )

!*****************************************************************************80
!
!! BESEI1 evaluates the exponentially scaled Bessel I1(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for the
!    modified Bessel function of the first kind of order one
!    multiplied by EXP(-ABS(X)).
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
!    Output, real ( kind = 8 ) BESEI1, the value of the function.
!
  implicit none

  real ( kind = 8 ) besei1
  integer ( kind = 4 ) jint
  real ( kind = 8 ) result
  real ( kind = 8 ) x

  jint = 2
  call calci1 ( x, result, jint )
  besei1 = result

  return
end
