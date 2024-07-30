module PolicyFunctions_obj
    
    use Model_Parameters
    
    implicit none
    
    integer, parameter :: nx = nk
    integer, parameter :: ny = nexp
    integer, parameter :: nz = nexp
    integer, parameter :: kxx = 4 !3
    integer, parameter :: kyy = 4 !3
    integer, parameter :: kzz = 4 !3
    
    type policy_fn_2d
        real(8) :: coefs(nx, ny)
        integer :: inbvx
        integer :: inbvy
        integer :: iloy
        real(8) :: tx(nx+kxx)
        real(8) :: ty(ny+kyy)
    contains
        procedure :: eval => eval_2d
        procedure :: set => set_2d
    end type policy_fn_2d
    
    type policy_fn_3d
        real(8) :: coefs(nx, ny, nz)
        real(8) :: tx(nx+kxx)
        real(8) :: ty(ny+kyy)
        real(8) :: tz(nz+kzz)
        real(8) :: xgrid(nx)
        real(8) :: ygrid(ny)
        real(8) :: zgrid(nz)
    contains
        procedure :: set => set_3d
        procedure :: eval => eval_3d
    end type policy_fn_3d
    
    type (policy_fn_3d), dimension(na, nu, na, nu, nfc, nfcm, 2), target :: pol_v_aux
    type (policy_fn_3d), dimension(na, nu, na, nu, nfc, nfcm, 2), target :: pol_ev_aux
    type (policy_fn_3d), dimension(na, nu, na, nu, T, nfc, nfcm, 2, 2, 2), target :: pol_v_mar_lfp
    type (policy_fn_2d), dimension(2, na, nu, T, nfc, 2) :: pol_v_sing_lfp
    type (policy_fn_2d), dimension(2,2,na,nu,Tret,nfc), target :: pol_vs_ret
    type (policy_fn_3d), dimension(2,2,na,nu,na,nu,Tret,nfc,nfcm), target :: pol_vm_ret
    
    
contains

    subroutine update_ev_aux(i_age)
        use PolicyFunctions, only: exp_grid, ev_aux, v_aux, k_grid
        use GlobParams, only: MEN, WOMEN
        implicit none
        integer :: i_age
        integer :: i, ia, iu, ifc
        real(8), pointer :: exp_grid_ptr(:) 
        integer :: iam, iaf, ium, iuf, ifcm, ifcf
        integer :: iexp1, iexp2
        real(8) :: ev_sing_tmp(nk,nexp), ev_mar_tmp(nk,nexp,nexp)
        real(8), pointer :: v_mar_tmp(:,:,:)
        integer :: i_m, i_f, jj
        
        exp_grid_ptr => exp_grid(:,i_age)     
        
        do iam = 1, na
            do iaf = 1, na
                do ium = 1, nu
                    do iuf = 1, nu
                        do ifcm = 1, nfc
                            do ifcf = 1, nfc                                
                                do jj = MEN, WOMEN
                                    v_mar_tmp => ev_aux(:,:,:,iam,ium,iaf,iuf,ifcm,ifcf,jj)
                                    call pol_ev_aux(iam, ium, iaf, iuf, ifcm, ifcf, jj)%set(v_mar_tmp, k_grid, exp_grid_ptr, exp_grid_ptr)
                                    !v_mar_tmp => v_aux(:,:,:,iam,ium,iaf,iuf,ifcm,ifcf,jj)
                                    !call pol_v_aux(iam, ium, iaf, iuf, ifcm, ifcf, jj)%set(v_mar_tmp, k_grid, exp_grid_ptr, exp_grid_ptr)
                                end do             
                            end do
                        end do
                    end do
                end do
            end do
        end do
                
    end subroutine update_ev_aux   

    subroutine update_lfp_policies(i_age)
        use PolicyFunctions, only: exp_grid, evs, evm, ev, k_grid, v, vs, v_lfp, vs_lfp
        use GlobParams, only: LFP_M1, LFP_M0, LFP_F1, LFP_F0, LFP_0, LFP_1, MEN, WOMEN
        !use pyplot_module, only : pyplot, wp => pyplot_wp
        
        integer :: i_age
        integer :: i, ia, iu, ifc
        real(8), pointer :: exp_grid_ptr(:) 
        integer :: iam, iaf, ium, iuf, ifcm, ifcf
        integer :: iexp1, iexp2
        real(8) :: ev_sing_tmp(nk,nexp), ev_mar_tmp(nk,nexp,nexp)
        real(8) :: v_mar_tmp(nk,nexp,nexp), v_sing_tmp(nk,nexp)
        !type(pyplot) :: plt   !! pytplot handler
        character(len=*), parameter :: testdir = "Plots/"
        integer :: istat !! status code
        character(len=2) :: age_ch
        integer :: i_m, i_f, ii
        
        write(age_ch, '(i2)') i_age
        
        exp_grid_ptr => exp_grid(:,i_age)
        do i = 1, 2
            do ia = 1, na
                do iu = 1, nu
                    do ifc = 1, nfc
                        !ev_sing_tmp = beta*OmegaActive(i_age-1)*( (1d0-Probm(i_age-1))*evs(i,:,:,ia,iu,i_age,ifc) + Probm(i_age-1)*evm(i,:,:,ia,iu,i_age,ifc) )
                        !call pol_ev_single(i, ia, iu, ifc)%set(ev_sing_tmp, k_grid, exp_grid_ptr)
                        !v_sing_tmp = Vs(i,:,:,ia,iu,i_age,ifc)
                        !call pol_v_sing(i, ia, iu, i_age, ifc)%set(v_sing_tmp, k_grid, exp_grid_ptr)
                        do ii = LFP_1, LFP_0
                            v_sing_tmp = Vs_lfp(i,:,:,ia,iu,i_age,ifc,ii)
                            call pol_v_sing_lfp(i, ia, iu, i_age, ifc, ii)%set(v_sing_tmp, k_grid, exp_grid_ptr)                            
                        end do
                    end do
                end do
            end do
        end do          
        
        
        do iam = 1, na
            do iaf = 1, na
                do ium = 1, nu
                    do iuf = 1, nu
                        do ifcm = 1, nfc
                            do ifcf = 1, nfc
                                !ev_mar_tmp = beta*OmegaActive(i_age-1)*ev(:,:,:,iam,ium,iaf,iuf,i_age,ifcm,ifcf)
                                !ev_mar_tmp = ev(:,:,:,iam,ium,iaf,iuf,i_age,ifcm,ifcf)
                                !call pol_ev(iam, ium, iaf, iuf, ifcm, ifcf)%set(ev_mar_tmp, k_grid, exp_grid_ptr, exp_grid_ptr)
                                !v_mar_tmp = v(:,:,:,iam,ium,iaf,iuf,i_age,ifcm,ifcf)
                                !call pol_v_mar(iam, ium, iaf, iuf, i_age, ifcm, ifcf)%set(v_mar_tmp, k_grid, exp_grid_ptr, exp_grid_ptr)
                                
                                do i_m = LFP_M1, LFP_M0
                                    do i_f = LFP_F1, LFP_F0
                                        do ii = MEN, WOMEN
                                            v_mar_tmp = v_lfp(:,:,:,iam,ium,iaf,iuf,i_age,ifcm,ifcf,i_m,i_f,ii)
                                            call pol_v_mar_lfp(iam, ium, iaf, iuf, i_age, ifcm, ifcf, i_m, i_f,ii)%set(v_mar_tmp, k_grid, exp_grid_ptr, exp_grid_ptr)
                                        end do
                                    end do
                                end do
                            
                            end do
                        end do
                    end do
                end do
            end do
        end do
        
        !call plt%initialize(grid=.true.,xlabel='Savings',figsize=[20,10],&
        !                    title='Married',legend=.true.,axis_equal=.true.,&
        !                    tight_layout=.true.)  
        !iam = 2
        !ium = 2
        !iaf = 2
        !iuf = 2
        !ifcm = 1
        !ifcf = 1
        !call plt%add_plot(k_grid,v(:,iexp1,iexp2,iam,ium,iaf,iuf,i_age,ifcm,ifcf),label='iexp=1',linestyle='b-o',markersize=5,linewidth=2,istat=istat)
        !call plt%savefig(testdir//'Vmarried'//age_ch//'.png', pyfile=testdir//'plottest.py',istat=istat)   
        
    end subroutine update_lfp_policies    
    
     subroutine update_ret_policies(i_age, make_plot)
    use PolicyFunctions, only: exp_grid, k_grid, ev_ret, evs_ret, vs_ret, v_ret
    !use fplot_core
        integer :: i_age
        integer :: ik2, j, ifcm, iaf, iuf, ium, iam, ifc
        integer :: ia, iu
        real(8), pointer :: exp_grid_ptr(:) 
        real(8) :: ev_sing_tmp(nk,nexp), ev_mar_tmp(nk,nexp,nexp), v_sing_tmp(nk,nexp), v_mar_tmp(nk,nexp,nexp)
        real(8) :: Omega
        !type(plot_2d) :: plt
        !type(plot_data_2d) :: d1, d2
        !class(plot_axis), pointer :: xAxis, yAxis
        !type(legend), pointer :: leg  
        logical :: make_plot

        exp_grid_ptr => exp_grid(:,T+i_age)
        !exp_grid_ptr => exp_grid(:,i_age+1)
        !if (i_age == T) then
        !    Omega = OmegaActive(T)
        !else
        !    Omega = OmegaRet(i_age - T) 
        !end if
        
        do ik2=1,2
        do j=1,2
        do ifcm=1,nfc
        do ifc=1,nfc
        do iaf=1,na
        do iam=1,na
        do iuf=1,nu
        do ium=1,nu
            v_mar_tmp = v_ret(ik2,j,:,:,:,iam,ium,iaf,iuf,i_age,ifc,ifcm)
            call pol_vm_ret(ik2,j,iam,ium,iaf,iuf,i_age,ifc,ifcm)%set(v_mar_tmp, k_grid, exp_grid_ptr, exp_grid_ptr)                
        end do
        end do
        end do
        end do
        end do
        end do    
        end do
        end do        
        
        
        do ik2=1,2
            do j=1,2
                do ia=1,na
                    do iu=1,nu
                        do ifc=1,nfc
                            v_sing_tmp = vs_ret(ik2,j,:,:,ia,iu,i_age,ifc)
                            call pol_vs_ret(ik2,j,ia,iu,i_age,ifc)%set(v_sing_tmp, k_grid, exp_grid_ptr)                            
                        end do
                    end do
                end do
            end do
        end do        
        
    end subroutine update_ret_policies   

    subroutine set_2d(self, fvals, xgrid, ygrid)
        use bspline_sub_module, only: db2ink
        class(policy_fn_2d) :: self
        real(8) :: xgrid(:), ygrid(:)
        real(8), intent(in) :: fvals(nx, ny)
        integer :: iflag
        integer :: iknot
        
        iknot = 0
        call db2ink(xgrid, nx, ygrid, ny, fvals, &
            kxx, kyy, iknot, self%tx, self%ty, self%coefs, iflag) 
    end subroutine set_2d
    
    subroutine set_3d(self, fvals, xgrid, ygrid, zgrid)
        use bspline_sub_module, only: db3ink
        class(policy_fn_3d) :: self
        real(8) :: xgrid(:), ygrid(:), zgrid(:)
        real(8), intent(in) :: fvals(nx, ny, nz)
        integer :: iflag, iknot
        
        iknot = 0
        call db3ink(xgrid, nx, ygrid, ny, zgrid, nz, &
                    fvals, kxx, kyy, kzz, iknot, self%tx, self%ty, self%tz, &
                    self%coefs, iflag)      
        !self%coefs = fvals
        !self%xgrid = xgrid
        !self%ygrid = ygrid
        !self%zgrid = zgrid
    end subroutine set_3d

    function eval_2d(self, x)
        use bspline_sub_module, only: db2val
        class(policy_fn_2d) :: self
        real(8), intent(in) :: x(2)
        real(8) :: eval_2d
        integer :: iflag
        real(8) :: w1(kyy)
        real(8) :: w0(3*max(kxx,kyy) )
        integer :: inbvx, inbvy, iloy
        
        inbvx = 1
        inbvy = 1
        iloy = 1        
        call db2val(x(1), x(2), 0, 0, &
                    self%tx, self%ty, nx, ny, kxx, kyy, &
                    self%coefs, eval_2d, iflag, &
                    inbvx, inbvy, iloy, w1, w0, extrap=.true.)         
    end function eval_2d
    
    function eval_3d(self, x)
        use bspline_sub_module, only: db3val
        use Utilities, only: trilin_interp
        class(policy_fn_3d) :: self
        real(8), intent(in) :: x(3)
        real(8) :: eval_3d
        integer :: iflag
        integer :: inbvx, inbvy, inbvz, iloy, iloz
        real(8) :: ww2(kyy,kzz), ww1(kzz), ww0(3*max(kxx,kyy,kzz))

        inbvx = 1
        inbvy = 1
        inbvz = 1
        iloy = 1        
        iloz = 1        
        
        call db3val(x(1), x(2), x(3), 0, 0, 0, &
                    self%tx, self%ty, self%tz, &
                    nx, ny, nz, kxx, kyy, kzz, &
                    self%coefs, eval_3d, iflag,&
                    inbvx, inbvy, inbvz, iloy, iloz, &
                    ww2, ww1, ww0, extrap=.true.)         
        
        !eval_3d = trilin_interp(self%xgrid, self%ygrid, self%zgrid, self%coefs, nx, ny, nz, x)
        
    end function eval_3d
    
end module PolicyFunctions_obj