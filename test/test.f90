program main
    use, intrinsic :: iso_fortran_env, only: wp => real64
    use findiffmod
    implicit none
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
    fd = findiff(func, x, 1.e-2_wp, 9)
    do i=1,2
    print *, truth(i,:) - fd(i,:)
    end do
    print *, "Hess Err"
    print *, "~~"
    fd2 = findiffhes(deriv, x, 1.e-2_wp, 9)
    do i=1,2
    do j = 1,2
    print *, truth2(i,j,:) - fd2(i,j,:)
    end do
    print *, ""
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

    contains
        function func(x) result(res)
            real(wp), intent(in) :: x(:)
            real(wp)             :: res(size(x))

            res(1) = x(1)**2 * sin(x(2)) + x(2)
            res(2) = exp(x(1))
        end function
        function deriv(x) result(res)
            real(wp), intent(in) :: x(:)
            real(wp)             :: res(size(x),size(x))
            res(1,1) = 2*x(1) * sin(x(2))
            res(1,2) = x(1)**2 * cos(x(2)) + 1._wp

            res(2,1) = exp(x(1))
            res(2,2) = 0._wp
        end function
        function dderiv(x) result(res)
            real(wp), intent(in) :: x(:)
            real(wp)             :: res(size(x),size(x),size(x))
            res(1,1,1) = 2 * sin(x(2))
            res(1,1,2) = 2*x(1) * cos(x(2))

            res(1,2,1) = 2*x(1) * cos(x(2))
            res(1,2,2) = -x(1)**2 * sin(x(2))

            res(2,1,1) = exp(x(1))
            res(2,1,2) =  0._wp

            res(2,2,1) =  0._wp
            res(2,2,2) =  0._wp
        end function

end program main
