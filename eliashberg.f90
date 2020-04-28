module constants
implicit none
    real, parameter :: RY_TO_K = 157887.6633481157
end module constants

module utils
implicit none
contains

    ! Integrate y(x) using the trapesium rule
    function integrate(x, y) result(integral)
        real, intent(in)  :: x(:), y(:)
        real              :: integral
        integer           :: i

        integral = 0
        do i=1, size(x)-1
            integral = integral + (x(i+1) - x(i)) * (y(i+1)+y(i)) * 0.5
        enddo

    end function integrate

    ! Calculate the Allen-Dynes approximateion to Tc
    function allen_dynes(w_a2f, a2f, mu) result (tc_ad)
        use                  constants
        real, intent(in)  :: w_a2f(:), a2f(:), mu
        real, allocatable :: to_int(:)
        real              :: lambda, w_log, w_rms 
        real              :: g1, g2, f1, f2       
        real              :: tc_ad

        write(*,*) "Calculating Allen-Dynes Tc"

        allocate(to_int(size(a2f)))
        to_int = 2*a2f/w_a2f
        lambda = integrate(w_a2f, to_int)
        to_int = log(w_a2f)*a2f/w_a2f
        w_log  = exp(2.0*integrate(w_a2f, to_int)/lambda)
        to_int = w_a2f*a2f
        w_rms  = (2.0*integrate(w_a2f, to_int)/lambda)**0.5
        g1     = 2.46*(1.0 + 3.8*mu)
        g2     = 1.82*(1.0 + 6.3*mu)*(w_rms/w_log)
        f1     = (1.0 + (lambda/g1)**(3.0/2.0))**(1.0/3.0)
        f2     = 1.0 + (w_rms/w_log - 1.0)*(lambda**2.0)/(lambda**2.0 + g2**2.0)
        tc_ad  = exp(-1.04*(1.0 + lambda)/(lambda - mu - 0.62*lambda*mu))
        tc_ad  = tc_ad * f1 * f2 * RY_TO_K * w_log / 1.2
        deallocate(to_int)

        write(*,"(A16, A16, A16, A16)") "Lambda", "W_log", "W_rms", "Tc_ad"
        write(*,"(F16.6, F16.6, F16.6, F16.6)") lambda, w_log, w_rms, tc_ad

    end function allen_dynes

    ! Read the eliashberg function from 
    ! the file "a2f" with the format:
    !     w_1 a2f(w_1)
    !     w_2 a2f(w_2)
    !     ...
    subroutine read_a2f(w_a2f, a2f)
        real, allocatable, intent(inout) :: w_a2f(:), a2f(:)
        real                             :: tmp1, tmp2
        integer                          :: n_a2f, i, ierr

        open(unit=10, file="a2f")

        ! Read the number of lines in the a2f file
        n_a2f = 0
        do

            read(10,*,iostat=ierr) tmp1, tmp2
            if (ierr > 0) then
                write(*,*) "Error reading a2f file! line", (n_a2f+1)
                stop
            endif
            if (ierr < 0) then
                exit ! End of file reached
            endif
            n_a2f = n_a2f + 1

        end do
        rewind(10)

        ! Read the eliashberg function from
        ! the a2f file. 
        allocate(w_a2f(n_a2f))
        allocate(a2f(n_a2f))
        do i=1, n_a2f

            read(10,*) w_a2f(i), a2f(i)

            ! Ignore negative (unstable) phonon modes 
            ! and negative eliashberg function values.
            if (w_a2f(i) < 0) a2f(i) = 0
            if (a2f(i)   < 0) a2f(i) = 0

        end do
        close(10)

    end subroutine read_a2f

end module utils

! Calculate the critical temperature of a superconductor with
! the Eliashberg function stored in the file a2f, and the value
! of the coulomb pseudopotential stored in the file mu. Units
! of Ry for energies are assumed throughout.
program eliashberg
    use                  utils
    implicit none
    real              :: mu                   ! Coulomb pseudopotential
    real              :: tc_ad                ! Allen-Dynes critical temperature
    real, allocatable :: w_a2f(:), a2f(:)     ! The eliashberg function a2f(w_a2f)

    call read_a2f(w_a2f, a2f)

    ! Read the value of the coulomb pseudopotential
    ! mu* from the mu file
    open(unit=10, file="mu")
    read(10,*) mu
    close(10)

    ! Calculate the Allen-Dynes approximation to Tc
    tc_ad = allen_dynes(w_a2f, a2f, mu)

end program eliashberg
