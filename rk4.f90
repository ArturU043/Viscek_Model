subroutine rk4(x,y,dx,n,deriv)
! 4th order Runge-Kutta. See Numerical Recipes p. 701ff
! f90 & double precision
implicit none
integer           , intent(in)    :: n
real              , intent(in)    :: x, dx
real, dimension(n), intent(inout) :: y
real                              :: ddx
real, dimension(n)                :: yp, k1, k2, k3, k4

ddx = 0.5*dx
call deriv(x,y,k1,n)      ; yp = y + ddx*k1
call deriv(x+ddx,yp,k2,n) ; yp = y + ddx*k2
call deriv(x+ddx,yp,k3,n) ; yp = y +  dx*k3
call deriv(x+dx,yp,k4,n)  ; y  = y +  dx*( k1 + 2.0*k2 + 2.0*k3 + k4 )/6.0
end subroutine rk4
