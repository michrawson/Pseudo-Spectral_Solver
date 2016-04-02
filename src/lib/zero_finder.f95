module zero_finder
implicit none

contains

    real ( kind = 8 ) function triginterp_caller(p_xi)
        implicit none

        integer, parameter           :: n = 128
        real ( kind = 8 ), intent(in)                :: p_xi

        real ( kind = 8 )                 :: p_x(n)=0
        real ( kind = 8 )                 :: p_y(n)=0
        real ( kind = 8 )                   :: pi
        integer                             :: j

        pi = 4.*atan(1.)

        p_x = [ (j, j = 1,n) ]
        p_x = p_x/n*2*pi
        p_y = cos(p_x)
        triginterp_caller = triginterp(p_xi,p_x,p_y)

    end function triginterp_caller

    real ( kind = 8 ) function triginterp(p_xi,p_x,p_y)
        integer, parameter          :: n = 128

        real ( kind = 8 ), intent(in)                :: p_xi
        real ( kind = 8 ), intent(in)                   :: p_x(n)
        real ( kind = 8 ), intent(in)                   :: p_y(n)

        real ( kind = 8 )                   :: h,scale,P,pi,x(n),xi
        integer                             :: k

        pi = 4.*atan(1.)
        h = 2./N
        scale = (p_x(2)-p_x(1)) / h
        x = p_x/scale
        xi = p_xi/scale
        P = 0
        do k = 1, N
            if (xi-x(k) == 0) then
                P = P + p_y(k)
            else
                P = P + p_y(k) * sin(N*pi*(xi-x(k))/2) / (N*tan(pi*(xi-x(k))/2))
            end if
        end do
        triginterp = P
    end function triginterp

end module zero_finder
