module mod_usr
  use mod_hd
  implicit none
contains

  !! Some AMRVAC bookkeeping
  subroutine usr_init()
    usr_init_one_grid => kh_init
    usr_special_bc    => kh_boundaries
    call set_coordinate_system('Cartesian')
    call hd_activate()
  end subroutine usr_init

  !! Setting the initial condition
  subroutine kh_init(ixG^L,ix^L,w,x)
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S,1:ndim) ! The coordinates
    double precision, intent(inout) :: w(ixG^S,1:nw)   ! The variables (rho, momentum, pressure)
    double precision :: rho, pint, uinf, delta0, cn    ! Problem parameters
    logical          :: first = .true.

    ! ==================
    ! Problem parameters
    ! ==================
    ! Value of the density field (which is uniform)
    rho = 1.0d0

    ! Pressure at the interface (took this from `amrvac/tests/hd/Kelvin_Helmholtz_2D/mod_usr.t`)
    pint = 2.5d0

    ! Some parameters from the paper
    uinf = 1.0
    cn = 1.0d-3
    delta0 = 1.0 / 28.0

    ! ==================
    ! Setting the fields
    ! ==================
    ! Initial density field
    w(ixG^S, rho_) = rho

    ! Initial pressure field
    w(ixG^S,p_) = pint

    ! Initial horizontal momentum
    w(ixG^S,mom(1)) = uinf * tanh((2*x(ixG^S,2)-1)/delta0) &
        + cn * uinf * (-(2*x(ixG^S,2)-1) / delta0**2) * exp(-((x(ixG^S,2)-0.5)/delta0)**2) &
        * (cos(8*dpi*x(ixG^S,1)) + cos(20*dpi*x(ixG^S,1)))

    ! Initial vertical momentum
    w(ixG^S,mom(2)) = cn * uinf * exp(-((x(ixG^S,2)-0.5)/delta0)**2) &
        * (8*dpi*sin(8*dpi*x(ixG^S,1)) + 20*dpi*sin(20*dpi*x(ixG^S,1)))

    call hd_to_conserved(ixG^L,ix^L,w,x)
  end subroutine kh_init

  !! Setting the boundary condition
  subroutine kh_boundaries(qt,ixI^L,ixO^L,iB,w,x)
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer :: ix^D

    select case(iB)
    case(3)
        ! Bottom boundary: Free-slip
        w(ixOmin1:ixOmax1, ixOmin2:ixOmax2, mom(1)) = &
            w(ixOmin1:ixOmax1, ixOmax2+nghostcells:ixOmax2+1:-1, mom(1)) &
            / w(ixOmin1:ixOmax1, ixOmax2+nghostcells:ixOmax2+1:-1, rho_)
        w(ixOmin1:ixOmax1, ixOmin2:ixOmax2, mom(2)) = &
            -w(ixOmin1:ixOmax1, ixOmax2+nghostcells:ixOmax2+1:-1, mom(2)) &
            / w(ixOmin1:ixOmax1, ixOmax2+nghostcells:ixOmax2+1:-1, rho_)
    case(4)
        ! Top boundary: Free-slip
        w(ixOmin1:ixOmax1, ixOmin2:ixOmax2, mom(1)) = &
            w(ixOmin1:ixOmax1, ixOmin2-1:ixOmin2-nghostcells:-1, mom(1)) &
            / w(ixOmin1:ixOmax1, ixOmin2-1:ixOmin2-nghostcells:-1, rho_)
        w(ixOmin1:ixOmax1, ixOmin2:ixOmax2, mom(2)) = &
            -w(ixOmin1:ixOmax1, ixOmin2-1:ixOmin2-nghostcells:-1, mom(2)) &
            / w(ixOmin1:ixOmax1, ixOmin2-1:ixOmin2-nghostcells:-1, rho_)
    case default
       call mpistop("Special boundary is not defined for this region")
    end select

    call hd_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine kh_boundaries

end module mod_usr
