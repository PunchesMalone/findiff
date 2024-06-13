program main
    use, intrinsic :: iso_fortran_env, only: wp => real64
    use findiffmod
    implicit none
    character(len=19) :: fname
    character(len=2)  :: row(6), col(6), ord(4), page(6)
    real(wp) :: x(6), &
              & truth(4,20,6,6), &
              & truth2(4,20,6,6,6), &
              & fd(4,20,6,6), &
              & fd2(4,20,6,6,6), err(4,20,6,6), err2(4,20,6,6,6), precs(20)

    integer i,j,k,l, order(4)

    row = ["r1","r2","r3","r4","r5","r6"]
    col = ["c1","c2","c3","c4","c5","c6"]
    page = ["p1","p2","p3","p4","p5","p6"]
    ord = ["o3","o5","o7","o9"]
    order = [3, 5, 7, 9]
    precs = [(10**real(i,wp), i=4,-15,-1)]
    call random_number(x)
    x = x*10._wp
    do i=1,4
        do j=1,20
           truth(i,j,:,:) = jac(x)
           truth2(i,j,:,:,:) = hes(x)
           fd(i,j,:,:) = findiff(acc,x,precs(j),order(i))
           fd2(i,j,:,:,:) = findiffhes(jac,x,precs(j),order(i))
        end do
    end do
    err = abs(fd - truth)
    err2 = abs(fd2 - truth2)

    do i=1,6
    do j=1,6
    do k=1,4
        fname ="outs/err"//row(i)//col(j)//ord(k)
        call print_to_file(trim(fname),err(k,:,i,j))
        
    end do
    end do
    end do
    do i=1,6
    do j=1,6
    do k=1,6
    do l=1,4
        fname ="outs/err2"//row(i)//col(j)//page(k)//ord(l)
        call print_to_file(trim(fname),err2(l,:,i,j,k))
    end do
    end do
    end do
    end do


contains

    ! BEGIN AUTOCODE OUTPUT FOR ACC
    function acc(x) result(res)
        implicit none
        real(wp), intent(in) :: x (:)
        real(wp)             :: &
                              & x0

        real(wp) :: res(size(x))

        x0  =  (x(1)**2 + x(2)**2 + x(3)**2)**(-3.0_wp/2.0_wp)

        res(1) =  x(4)
        res(2) =  x(5)
        res(3) =  x(6)
        res(4) =  -x(1)*x0
        res(5) =  -x(2)*x0
        res(6) =  -x(3)*x0
    end function acc
    ! END AUTOCODE OUTPUT FOR ACC

    ! BEGIN AUTOCODE OUTPUT FOR JAC
    function jac(x) result(res)
        implicit none
        real(wp), intent(in) :: x (:)
        real(wp)             :: &
                              & x0, x1, x2, x3, x4, x5, x6, x7, & 
                              & x8, x9
        real(wp), dimension(36) :: inter
        real(wp)  :: res(size(x),size(x))

        x0  =  x(1)**2
        x1  =  x(2)**2
        x2  =  x(3)**2
        x3  =  x0 + x1 + x2
        x4  =  -1/x3**(3.0_wp/2.0_wp)
        x5  =  3/x3**(5.0_wp/2.0_wp)
        x6  =  x(1)*x5
        x7  =  x(2)*x6
        x8  =  x(3)*x6
        x9  =  x(2)*x(3)*x5

        inter(1) =  0
        inter(2) =  0
        inter(3) =  0
        inter(4) =  x0*x5 + x4
        inter(5) =  x7
        inter(6) =  x8
        inter(7) =  0
        inter(8) =  0
        inter(9) =  0
        inter(10) =  x7
        inter(11) =  x1*x5 + x4
        inter(12) =  x9
        inter(13) =  0
        inter(14) =  0
        inter(15) =  0
        inter(16) =  x8
        inter(17) =  x9
        inter(18) =  x2*x5 + x4
        inter(19) =  1
        inter(20) =  0
        inter(21) =  0
        inter(22) =  0
        inter(23) =  0
        inter(24) =  0
        inter(25) =  0
        inter(26) =  1
        inter(27) =  0
        inter(28) =  0
        inter(29) =  0
        inter(30) =  0
        inter(31) =  0
        inter(32) =  0
        inter(33) =  1
        inter(34) =  0
        inter(35) =  0
        inter(36) =  0
        res = reshape(inter,[6,6])
    end function jac
    ! END AUTOCODE OUTPUT FOR JAC

    ! BEGIN AUTOCODE OUTPUT FOR HES
    function hes(x) result(res)
        implicit none
        real(wp), intent(in) :: x (:)
        real(wp)             :: &
                              & x0, x1, x2, x3, x4, x5, x6, x7, & 
                              & x8, x9, x10, x11, x12, x13, x14, x15, & 
                              & x16, x17, x18, x19, x20, x21, x22, x23, & 
                              & x24, x25
        real(wp), dimension(216) :: inter
        real(wp), dimension(6,6,6) :: res

        x0  =  x(1)**2
        x1  =  x(2)**2
        x2  =  x(3)**2
        x3  =  x0 + x1 + x2
        x4  =  5/x3
        x5  =  x0*x4
        x6  =  3/x3**(5.0_wp/2.0_wp)
        x7  =  x(1)*x6
        x8  =  x5 - 1
        x9  =  x6*x8
        x10  =  -x6*x8
        x11  =  x(2)*x10
        x12  =  x1*x4
        x13  =  x12 - 1
        x14  =  -x13
        x15  =  x14*x7
        x16  =  -15*x(1)*x(2)*x(3)/x3**(7.0_wp/2.0_wp)
        x17  =  x(3)*x10
        x18  =  x2*x4
        x19  =  x18 - 1
        x20  =  -x19
        x21  =  x20*x7
        x22  =  x(2)*x6
        x23  =  x(3)*x6
        x24  =  x14*x23
        x25  =  x20*x22

        inter(1) =  0
        inter(2) =  0
        inter(3) =  0
        inter(4) =  x7*(3 - x5)
        inter(5) =  -x(2)*x9
        inter(6) =  -x(3)*x9
        inter(7) =  0
        inter(8) =  0
        inter(9) =  0
        inter(10) =  x11
        inter(11) =  x15
        inter(12) =  x16
        inter(13) =  0
        inter(14) =  0
        inter(15) =  0
        inter(16) =  x17
        inter(17) =  x16
        inter(18) =  x21
        inter(19) =  0
        inter(20) =  0
        inter(21) =  0
        inter(22) =  0
        inter(23) =  0
        inter(24) =  0
        inter(25) =  0
        inter(26) =  0
        inter(27) =  0
        inter(28) =  0
        inter(29) =  0
        inter(30) =  0
        inter(31) =  0
        inter(32) =  0
        inter(33) =  0
        inter(34) =  0
        inter(35) =  0
        inter(36) =  0
        inter(37) =  0
        inter(38) =  0
        inter(39) =  0
        inter(40) =  x11
        inter(41) =  x15
        inter(42) =  x16
        inter(43) =  0
        inter(44) =  0
        inter(45) =  0
        inter(46) =  -x13*x7
        inter(47) =  x22*(3 - x12)
        inter(48) =  -x13*x23
        inter(49) =  0
        inter(50) =  0
        inter(51) =  0
        inter(52) =  x16
        inter(53) =  x24
        inter(54) =  x25
        inter(55) =  0
        inter(56) =  0
        inter(57) =  0
        inter(58) =  0
        inter(59) =  0
        inter(60) =  0
        inter(61) =  0
        inter(62) =  0
        inter(63) =  0
        inter(64) =  0
        inter(65) =  0
        inter(66) =  0
        inter(67) =  0
        inter(68) =  0
        inter(69) =  0
        inter(70) =  0
        inter(71) =  0
        inter(72) =  0
        inter(73) =  0
        inter(74) =  0
        inter(75) =  0
        inter(76) =  x17
        inter(77) =  x16
        inter(78) =  x21
        inter(79) =  0
        inter(80) =  0
        inter(81) =  0
        inter(82) =  x16
        inter(83) =  x24
        inter(84) =  x25
        inter(85) =  0
        inter(86) =  0
        inter(87) =  0
        inter(88) =  -x19*x7
        inter(89) =  -x19*x22
        inter(90) =  x23*(3 - x18)
        inter(91) =  0
        inter(92) =  0
        inter(93) =  0
        inter(94) =  0
        inter(95) =  0
        inter(96) =  0
        inter(97) =  0
        inter(98) =  0
        inter(99) =  0
        inter(100) =  0
        inter(101) =  0
        inter(102) =  0
        inter(103) =  0
        inter(104) =  0
        inter(105) =  0
        inter(106) =  0
        inter(107) =  0
        inter(108) =  0
        inter(109) =  0
        inter(110) =  0
        inter(111) =  0
        inter(112) =  0
        inter(113) =  0
        inter(114) =  0
        inter(115) =  0
        inter(116) =  0
        inter(117) =  0
        inter(118) =  0
        inter(119) =  0
        inter(120) =  0
        inter(121) =  0
        inter(122) =  0
        inter(123) =  0
        inter(124) =  0
        inter(125) =  0
        inter(126) =  0
        inter(127) =  0
        inter(128) =  0
        inter(129) =  0
        inter(130) =  0
        inter(131) =  0
        inter(132) =  0
        inter(133) =  0
        inter(134) =  0
        inter(135) =  0
        inter(136) =  0
        inter(137) =  0
        inter(138) =  0
        inter(139) =  0
        inter(140) =  0
        inter(141) =  0
        inter(142) =  0
        inter(143) =  0
        inter(144) =  0
        inter(145) =  0
        inter(146) =  0
        inter(147) =  0
        inter(148) =  0
        inter(149) =  0
        inter(150) =  0
        inter(151) =  0
        inter(152) =  0
        inter(153) =  0
        inter(154) =  0
        inter(155) =  0
        inter(156) =  0
        inter(157) =  0
        inter(158) =  0
        inter(159) =  0
        inter(160) =  0
        inter(161) =  0
        inter(162) =  0
        inter(163) =  0
        inter(164) =  0
        inter(165) =  0
        inter(166) =  0
        inter(167) =  0
        inter(168) =  0
        inter(169) =  0
        inter(170) =  0
        inter(171) =  0
        inter(172) =  0
        inter(173) =  0
        inter(174) =  0
        inter(175) =  0
        inter(176) =  0
        inter(177) =  0
        inter(178) =  0
        inter(179) =  0
        inter(180) =  0
        inter(181) =  0
        inter(182) =  0
        inter(183) =  0
        inter(184) =  0
        inter(185) =  0
        inter(186) =  0
        inter(187) =  0
        inter(188) =  0
        inter(189) =  0
        inter(190) =  0
        inter(191) =  0
        inter(192) =  0
        inter(193) =  0
        inter(194) =  0
        inter(195) =  0
        inter(196) =  0
        inter(197) =  0
        inter(198) =  0
        inter(199) =  0
        inter(200) =  0
        inter(201) =  0
        inter(202) =  0
        inter(203) =  0
        inter(204) =  0
        inter(205) =  0
        inter(206) =  0
        inter(207) =  0
        inter(208) =  0
        inter(209) =  0
        inter(210) =  0
        inter(211) =  0
        inter(212) =  0
        inter(213) =  0
        inter(214) =  0
        inter(215) =  0
        inter(216) =  0
        res = reshape(inter,[6,6,6])
    end function hes
    ! END AUTOCODE OUTPUT FOR HES
    subroutine print_to_file(fname, var)
        integer io,j
        real(wp), intent(in) :: var(:)
        character(len=*) :: fname
        open(newunit=io, file=trim(adjustl(fname))//".txt")
        do j = 1,size(var)
        write(io,*) var(j)
        end do
        close(io)
    end subroutine print_to_file
end program main
