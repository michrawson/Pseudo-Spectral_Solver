function besei0 ( x )

!*****************************************************************************80
!
!! BESEI0 evaluates the exponentially scaled Bessel I0(X) function.
!
!  Discussion:
!
!    This routine computes approximate values for the modified Bessel
!    function of the first kind of order zero multiplied by EXP(-ABS(X)).
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
!    Output, real ( kind = 8 ) BESEI0, the value of the function.
!
  implicit none

  real ( kind = 8 ) besei0
  integer ( kind = 4 ) jint
  real ( kind = 8 ) result
  real ( kind = 8 ) x

  jint = 2
  call calci0 ( x, result, jint )
  besei0 = result

  return
end
