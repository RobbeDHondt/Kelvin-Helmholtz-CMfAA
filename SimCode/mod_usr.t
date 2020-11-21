!> See http://amrvac.org/mod__usr__methods_8t_source.html for an overview
!> of user-definable methods.
module mod_usr
  use mod_hd
  use mod_viscosity

  implicit none

  integer :: omega_, grad_omega(2)
contains

  !! Some AMRVAC bookkeeping
  subroutine usr_init()
    usr_init_one_grid => kh_init
    usr_special_bc    => kh_boundaries
    usr_modify_output => set_output_vars
    usr_refine_grid   => specialrefine_grid

    call set_coordinate_system('Cartesian')
    call hd_activate()

    omega_        = var_set_extravar("omega", "omega")
    grad_omega(1) = var_set_extravar("grad_omega_x", "grad_omega_x")
    grad_omega(2) = var_set_extravar("grad_omega_y", "grad_omega_y")
  end subroutine usr_init

  !! Setting the initial condition
  subroutine kh_init(ixG^L,ix^L,w,x)
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S,1:ndim) ! The coordinates
    double precision, intent(inout) :: w(ixG^S,1:nw)   ! The variables (rho, momentum, pressure)
    double precision :: rho, uinf, delta0, cn, pint    ! Problem parameters
    logical          :: first = .true.

    ! ==================
    ! Problem parameters
    ! ==================
    ! Value of the density field (which is uniform)
    rho = 1.0d0

    ! Value of pressure at interface
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

    ! Initial pressure
    if (hd_energy) then
        w(ixG^S, p_) = pint
    end if

    ! Initial horizontal momentum
    w(ixG^S,mom(1)) = uinf * tanh((2.0d0*x(ixG^S,2)-1.0d0)/delta0) &
        + cn * uinf * (-(2.0d0*x(ixG^S,2)-1.0d0) / delta0**2.0d0) * exp(-((x(ixG^S,2)-0.5)/delta0)**2.0d0) &
        * (cos(8.0d0*dpi*x(ixG^S,1)) + cos(20.0d0*dpi*x(ixG^S,1)))

    ! Initial vertical momentum
    w(ixG^S,mom(2)) = cn * uinf * exp(-((x(ixG^S,2)-0.5)/delta0)**2.0d0) &
        * (8.0d0*dpi*sin(8.0d0*dpi*x(ixG^S,1)) + 20.0d0*dpi*sin(20.0d0*dpi*x(ixG^S,1)))

    call hd_to_conserved(ixG^L,ix^L,w,x)
  end subroutine kh_init

  !! Setting the boundary condition
  subroutine kh_boundaries(qt,ixI^L,ixO^L,iB,w,x)
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer :: ix^D
    double precision :: rho, pint

    ! Value for the Dirichlet boundary
    rho  = 1.0d0
    pint = 2.5d0

    select case(iB)
    case(3)
        ! Bottom boundary: Free-slip
        ! w(ixO^S, mom(1)) = w(ixO^LIM1, ..., mom(1)) / w(ixO^LIM1, ..., mom(2))
        w(ixOmin1:ixOmax1, ixOmin2:ixOmax2, mom(1)) = &
            w(ixOmin1:ixOmax1, ixOmax2+nghostcells:ixOmax2+1:-1, mom(1)) &
            / w(ixOmin1:ixOmax1, ixOmax2+nghostcells:ixOmax2+1:-1, rho_)
        w(ixOmin1:ixOmax1, ixOmin2:ixOmax2, mom(2)) = &
            -w(ixOmin1:ixOmax1, ixOmax2+nghostcells:ixOmax2+1:-1, mom(2)) &
            / w(ixOmin1:ixOmax1, ixOmax2+nghostcells:ixOmax2+1:-1, rho_)
        ! density: wall
        w(ixO^S,rho_)=rho
        !pressure: constant
        if (hd_energy) then
            w(ixO^S,p_)=pint
        end if
    case(4)
        ! Top boundary: Free-slip
        w(ixOmin1:ixOmax1, ixOmin2:ixOmax2, mom(1)) = &
            w(ixOmin1:ixOmax1, ixOmin2-1:ixOmin2-nghostcells:-1, mom(1)) &
            / w(ixOmin1:ixOmax1, ixOmin2-1:ixOmin2-nghostcells:-1, rho_)
        w(ixOmin1:ixOmax1, ixOmin2:ixOmax2, mom(2)) = &
            -w(ixOmin1:ixOmax1, ixOmin2-1:ixOmin2-nghostcells:-1, mom(2)) &
            / w(ixOmin1:ixOmax1, ixOmin2-1:ixOmin2-nghostcells:-1, rho_)
        ! density: wall
        w(ixO^S,rho_)=rho
        !pressure: constant
        if (hd_energy) then
            w(ixO^S,p_)=pint
        end if
    case default
       call mpistop("Special boundary is not defined for this region")
    end select

    call hd_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine kh_boundaries

  subroutine set_output_vars(ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    use mod_radiative_cooling    
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,nw)

    double precision :: drho(ixI^S), vrot(ixI^S), tmp(ixI^S)
    double precision :: wlocal(ixI^S,1:nw), domega_x(ixI^S), domega_y(ixI^S)
    integer          :: idims

    ! Make a copy for local computations
    wlocal(ixI^S,1:nw) = w(ixI^S,1:nw)

    ! ================
    ! Output vorticity
    ! ================
    vrot(ixO^S) = zero

    ! x-direction: compute d(v_y)/dx
    idims = 1
    tmp(ixI^S) = wlocal(ixI^S,mom(2)) / wlocal(ixI^S,rho_) ! = velocity v_y
    call gradient(tmp, ixI^L, ixO^L, idims, drho) ! drho = d(tmp)/dx
    vrot(ixO^S) = vrot(ixO^S) + drho(ixO^S)

    ! y-direction: compute d(v_x)/dy
    idims = 2
    tmp(ixI^S) = wlocal(ixI^S,mom(1)) / wlocal(ixI^S,rho_) ! = velocity v_x
    call gradient(tmp, ixI^L, ixO^L, idims, drho) ! drho = d(tmp)/dy
    vrot(ixO^S) = vrot(ixO^S) - drho(ixO^S)

    ! Write the result: omega = d(v_y)/dx - d(v_x)/dy
    w(ixO^S,omega_) = vrot(ixO^S)

    ! ======================
    ! Output grad(vorticity)
    ! ======================
    tmp(ixI^S) = w(ixI^S,omega_) !vrot(ixO^S)
    call gradient(tmp, ixI^L, ixO^L, 1, domega_x) ! x-direction
    call gradient(tmp, ixI^L, ixO^L, 2, domega_y) ! y-direction
    w(ixO^S, grad_omega(1)) = domega_x(ixO^S)
    w(ixO^S, grad_omega(2)) = domega_y(ixO^S)

  end subroutine set_output_vars 

  subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen
    double precision:: R(ixG^S)

    ! Prohibit refinement in the low vorticity parts. Coarsen if possible.
    if (all(abs(w(ix^S, omega_)) <= 10)) then
        refine  = -1
        coarsen =  1
    end if
    ! Enforce refinement to highest level in the high vorticity parts
    if (any(abs(w(ix^S, omega_)) > 10)) then
        refine = 1
    end if

  end subroutine specialrefine_grid

end module mod_usr
