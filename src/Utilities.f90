module Utilities

    use Model_Parameters

    implicit none

contains

    function labor_tax_single(y) result(res)
        real(8) :: y
        real(8) :: res
        
        res = y - after_tax_labor_inc_single(y)
    end function labor_tax_single
    
    function labor_tax_married(y) result(res)
        real(8) :: y
        real(8) :: res
        
        res = y - after_tax_labor_inc_married(y)
    end function labor_tax_married

    function after_tax_labor_inc_single(y) result(res)
        real(8) :: y
        real(8) :: res
        
        if (y < Deduct_cutoff) then
            res = y
        else
            res = Deduct_cutoff + thetas(1)*(y - Deduct_cutoff)**(1.0d0-thetas(2))
        end if
    end function after_tax_labor_inc_single
    
    function after_tax_labor_inc_married(y) result(res)
        real(8) :: y
        real(8) :: res
        
        if (y < Deduct_cutoff_Mar) then
            res = y
        else
            res = Deduct_cutoff_Mar + theta(1)*(y - Deduct_cutoff_Mar)**(1.0d0-theta(2))
        end if
    end function after_tax_labor_inc_married   

    function wage(gender,aval,x,uval)
        integer :: i,gender
        real(8) :: aval, x, uval
        real(8) :: xdum, w1, w2
        real(8) :: wage

        xdum=min(x,(T-1)*(1d0))
        !xdum=x
        
        if(gender==1) then
            wage = w*exp(gamma0+gamma(gender,1)*xdum + gamma(gender,2)*xdum**2d0 + gamma(gender,3)*xdum**3d0 + aval + uval)
        else
            wage = w*exp(gamma0f+gamma(gender,1)*xdum + gamma(gender,2)*xdum**2d0 + gamma(gender,3)*xdum**3d0 + aval + uval)
        end if
        
    end function wage

    function tSS_employer(y)
        real(8) :: y
        real(8) :: tSS_employer

        tSS_employer = t_employer

    end function tSS_employer

    function tSS_employee(y)
        real(8) :: y
        real(8) :: tSS_employee

        tSS_employee = t_employee

    end function tSS_employee

    function Uc(c)
        real(8) :: c, Uc

        Uc = log(c)
        !Uc = c**(1d0-sigma)/(1d0-sigma)
    end function Uc

    function Ul(lm,lf)
        real(8) :: lm,lf, Ul

        Ul = -chim*(lm**(1d0+etam))/(1d0+etam)-chif*(lf**(1d0+etaf))/(1d0+etaf)
    end function Ul

    function dUc(c)
        real(8) :: c, dUc

        dUc = c**(-1d0)
    end function dUc


    function LinInterp(x,xGrid,fVals,nx)
        real(8) :: LinInterp
        real(8) :: x
        integer :: nx
        real(8) :: xGrid(nx), fVals(nx)
        integer :: i
        real(8) :: b

        !        i = 2
        !        do while (x > xGrid(i) .and. i < nx)
        !            i = i + 1
        !        end do
        do i = 1, nx-1
            if (xGrid(i+1) >= x) exit
        end do
        i = min(i+1,nx)

        b = (x - xGrid(i-1))/(xGrid(i) - xGrid(i-1))

        LinInterp = (1d0-b)*fVals(i-1) + b*fVals(i)

        if(b>1.001d0) Then
            LinInterp=fVals(nx)+((x - xGrid(nx))/(xGrid(nx) - xGrid(nx-1)))*(fVals(nx)-fVals(nx-1))
        end if

    end function LinInterp

    subroutine MakeGrid(nx,grid,lbound,ubound,convex_param)
        integer :: nx
        real(8) :: grid(nx)
        integer :: i
        real(8) :: ubound, lbound, convex_param
        real(8) :: dx

        dx = 1d0/dble(nx-1)

        grid(1) = 0d0

        do i = 2, nx
            grid(i) = grid(i-1) + dx
        end do

        grid = grid**convex_param

        grid = grid*(ubound-lbound) + lbound

    end subroutine MakeGrid

    function locate(xx,x,n) result(loc)
        integer :: loc
        integer :: n
        !real(8), dimension (:), allocatable :: xx
        real(8) :: xx(n)
        real(8) :: x
        integer :: jl, ju, jm

        jl = 1
        ju = n

        do
            if (ju-jl <= 1) exit
            jm=(ju+jl)/2
            if (x >= xx(jm)) then
                jl=jm
            else
                ju=jm
            end if
        end do
        if (x == xx(1)) then
            loc=1
        else if (x == xx(n)) then
            loc=n-1
        else
            loc=jl
        end if

    end function locate

    function bilin_interp(xx,yy,V,nx,ny,arg) result(val)
        real(8) :: val
        integer :: nx, ny
        !real(8), dimension (:), allocatable :: xx, yy
        real(8) :: xx(nx)
        real(8) :: yy(ny)
        !real(8), dimension (:,:), allocatable :: V
        real(8) :: V(nx,ny)
        real(8) :: arg(2)
        real(8) :: x, y
        integer :: jx, jy
        real(8) :: dx, dy

        x = arg(1)
        y = arg(2)

        jx = locate(xx,x,nx)
        jy = locate(yy,y,ny)   
        dx = (x - xx(jx))/(xx(jx+1) - xx(jx))
        dy = (y - yy(jy))/(yy(jy+1) - yy(jy))        

        val = V(jx,jy)*(1d0-dx)*(1d0-dy) + V(min(jx+1,nx),jy)*dx*(1d0-dy) + V(jx,min(jy+1,ny))*(1d0-dx)*dy + V(min(jx+1,nx),min(jy+1,ny))*dx*dy

    end function bilin_interp

    function trilin_interp(xx,yy,zz,V,nx,ny,nz,arg) result(val)
        real(8) :: val
        integer :: nx, ny, nz
        !real(8), dimension (:), allocatable :: xx, yy, zz
        real(8) :: xx(nx), yy(ny), zz(nz)
        !real(8), dimension (:,:,:), allocatable :: V
        real(8) :: V(nx,ny,nz)
        real(8) :: arg(3)
        real(8) :: x, y, z
        integer :: jx, jy, jz
        real(8) :: dx, dy, dz
        real(8) :: c00, c10, c01, c11
        real(8) :: c0, c1

        x = arg(1)
        y = arg(2)
        z = arg(3)

        jx = locate(xx,x,nx)
        jy = locate(yy,y,ny)
        jz = locate(zz,z,nz)

        dx = (x - xx(jx))/(xx(jx+1) - xx(jx))
        dy = (y - yy(jy))/(yy(jy+1) - yy(jy))
        dz = (z - zz(jz))/(zz(jz+1) - zz(jz))

        c00 = V(jx,jy,jz)*(1-dx)                     + V(min(jx+1,nx),jy,jz)*dx
        c10 = V(jx,min(jy+1,ny),jz)*(1-dx)           + V(min(jx+1,nx),min(jy+1,ny),jz)*dx
        c01 = V(jx,jy,min(jz+1,nz))*(1-dx)           + V(min(jx+1,nx),jy,min(jz+1,nz))*dx
        c11 = V(jx,min(jy+1,ny),min(jz+1,nz))*(1-dx) + V(min(jx+1,nx),min(jy+1,ny),min(jz+1,nz))*dx

        c0 = c00*(1-dy) + c10*dy
        c1 = c01*(1-dy) + c11*dy

        val = c0*(1-dz) + c1*dz

    end function trilin_interp

    
     function cov_xy(x,y,n) result(val)
        integer,intent(in) :: n
        real(8),intent(in) :: x(n), y(n)
        real(8) :: val

        val = sum(x*y)/dble(n) - sum(x)/dble(n)*sum(y)/dble(n)

    end function cov_xy

    function corr_xy(x,y,n) result(val)
        integer,intent(in) :: n
        real(8),intent(in) :: x(n), y(n)
        real(8) :: val
        real(8) :: cov
        real(8) :: sx, sy

        cov = cov_xy(x,y,n)
        sx = sqrt(cov_xy(x,x,n))
        sy = sqrt(cov_xy(y,y,n))

        val = cov/sx/sy

    end function corr_xy
    
end module Utilities
