program main
    use, intrinsic :: iso_fortran_env, only: wp => real128
    use findiffmod
    implicit none
    real(wp), parameter :: a=1._wp, b=1000._wp
    real(wp) :: x(2), &
              & truth(2,2), &
              & truth2(2,2,2), &
              & fd(2,2), &
              & fd2(2,2,2)

    integer i, j

    call random_number(x)
    x = x * 10._wp
    truth = deriv(x)
    truth2 = dderiv(x)
    fd = findiff_multiscale(func, x, [1.e-2_wp, 1.e-6_wp], 9)
    fd2 = findiffhes_multiscale(deriv, x, [1.e-2_wp, 1.e-6_wp], 9)
    do i=1,2
    print *, truth(i,:) - fd(i,:)
    end do
    print *, "Hess"
    print *, "~~"
    do i=1,2
    do j = 1,2
    print *, truth2(i,j,:)
    end do
    print *, ""
    end do
    print *, "FD Hess"
    print *, "~~"
    do i=1,2
    do j = 1,2
    print *, fd2(i,j,:)
    end do
    print *, ""
    end do
    print *, "Hess Err"
    print *, "~~"
    do i=1,2
    do j = 1,2
    print *, truth2(i,j,:) - fd2(i,j,:)
    end do
    print *, ""
    end do

    contains
        function func(x) result(res)
            real(wp), intent(in) :: x(:)
            real(wp)             :: res(size(x))
            res(1) = cos(x(1)*a)*cos(x(2)*b)
            res(2) = sin(x(1)*a)*sin(x(2)*b)
        end function
        function deriv(x) result(res)
            real(wp), intent(in) :: x(:)
            real(wp)             :: res(size(x),size(x))
            res(1,1) = -sin(x(1)*a)*cos(x(2)*b)*a
            res(1,2) = -cos(x(1)*a)*sin(x(2)*b)*b
            res(2,1) =  cos(x(1)*a)*sin(x(2)*b)*a
            res(2,2) =  sin(x(1)*a)*cos(x(2)*b)*b
        end function
        function dderiv(x) result(res)
            real(wp), intent(in) :: x(:)
            real(wp)             :: res(size(x),size(x),size(x))
            res(1,1,1) = -cos(x(1)*a)*cos(x(2)*b)*a**2
            res(1,1,2) =  sin(x(1)*a)*sin(x(2)*b)*(a*b)
            res(1,2,1) =  sin(x(1)*a)*sin(x(2)*b)*(a*b)
            res(1,2,2) = -cos(x(1)*a)*cos(x(2)*b)*b**2
            res(2,1,1) = -sin(x(1)*a)*sin(x(2)*b)*a**2
            res(2,1,2) =  cos(x(1)*a)*cos(x(2)*b)*(a*b)
            res(2,2,1) =  cos(x(1)*a)*cos(x(2)*b)*(a*b)
            res(2,2,2) = -sin(x(1)*a)*sin(x(2)*b)*b**2
        end function
        ! function rosen(x) result(res)
        !     real(wp), intent(in) :: x(:)
        !     real(wp)             :: res(size(x))
        !     res(1) = (a-x(1))**2 + b*(x(2)-x(1)**2)**2
        !     res(2) = (a-x(2))**2 + b*(x(1)-x(2)**2)**2
        ! end function
        ! function deriv(x) result(res)
        !     real(wp), intent(in) :: x(:)
        !     real(wp)             :: res(size(x),size(x))
        !     res(1,1) = -2._wp*a - 4._wp*b*x(1)*(-x(1)**2 + x(2)) + 2._wp*x(1)
        !     res(1,2) =  b*(-2._wp*x(1)**2 + 2*x(2))
        !     res(2,1) = b*(-2._wp*x(2)**2 + 2*x(1))
        !     res(2,2) = -2._wp*a - 4._wp*b*x(2)*(-x(2)**2 + x(1)) + 2._wp*x(2)
        ! end function
        ! function dderiv(x) result(res)
        !     real(wp), intent(in) :: x(:)
        !     real(wp)             :: res(size(x),size(x),size(x))
        !     res(1,1,1) = 2*(4*b*x(1)**2 + 2*b*(x(1)**2 - x(2)) + 1)
        !     res(1,1,2) = -4*b*x(1)

        !     res(1,2,1) = -4*b*x(1)
        !     res(1,2,2) = 2*b

        !     res(2,1,1) = 2*b
        !     res(2,1,2) = -4*b*x(2)

        !     res(2,2,1) = -4*b*x(2)
        !     res(2,2,2) = 2*(4*b*x(2)**2 + 2*b*(x(2)**2 - x(1)) + 1)

        ! end function
end program main
