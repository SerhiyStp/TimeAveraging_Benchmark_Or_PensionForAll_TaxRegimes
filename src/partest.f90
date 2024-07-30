subroutine partest(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use Utilities
    use bspline_sub_module
    use GlobParams
    use PolicyFunctions_obj, only: pol_v_aux

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixm,iam,ium,iaf,iuf,iu2,j,ik2,ifc,ifcm,ik
    real(8) :: dum, dum3
    real(8) :: vnext
    integer :: iflag
    integer :: idx, idy, idz, iloy, iloz
    integer :: inbvx, inbvy, inbvz
    real(8) :: ww2(ky,kz),ww1(kz),ww0(3*max(kx,ky,kz))
    real(8) :: exp_m_prime, exp_f_prime

    idx=0
    idy=0
    idz=0
    inbvx=1
    inbvy=1
    inbvz=1
    iloy=1
    iloz=1    
    
    !Assigning the grid points
    dum3=((counter*1d0)/(nexp*nu*na*1d0))-0.00001d0
    ik=int(dum3)+1
    dum3=(((counter-(ik-1)*nexp*nu*na)*1d0)/(nu*na*1d0))-0.00001d0
    ix=int(dum3)+1
    dum3=(((counter-(ik-1)*nexp*nu*na-(ix-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
    ium=int(dum3)+1
    iam=counter-(ik-1)*nexp*na*nu-(ix-1)*na*nu-(ium-1)*na


    j=2
    exp_f_prime = exp_grid(ix,T-it)
    evm(j,ik,ix,iam,ium,T-it,:)=0d0
    do ifc=1,nfc
        do iu2=1,nu
            do ik2=1,nk
                do ixm = 1, nexp
                    exp_m_prime = exp_grid(ixm,T-it)
                    do iaf=1,na
                        do iuf=1,nu
                            do ifcm=1,nfcm
                                dum=k_grid(ik)+k_grid(ik2)
                                
                                vnext = pol_v_aux(iaf,iuf,iam,iu2,ifc,ifcm,WOMEN)%eval([dum,exp_f_prime,exp_m_prime])
                                evm(j,ik,ix,iam,ium,T-it,ifc) = evm(j,ik,ix,iam,ium,T-it,ifc) + &
                                                                trans_u(WOMEN,ium,iu2)*ability_prob(iam,iaf)*mpartner(ik2,ixm,iaf,iuf,T-it,ifcm)*vnext                                
                                
                                !if(dum<k_grid(nk)-0.001d0) then
                                !    call db3val(dum,exp_grid(ix,T-it),exp_grid(ixm,T-it),idx,idy,idz,&
                                !            tx,ty(:,T-it),tz(:,T-it),&
                                !            nk,nexp,nexp,kx,ky,kz,&
                                !            v_bspl(iaf,iuf,iam,iu2,ifc,ifcm)%coefs,vnext,iflag,&
                                !            inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.true.)
                                !    evm(j,ik,ix,iam,ium,T-it,ifc)=evm(j,ik,ix,iam,ium,T-it,ifc)+trans_u(2,ium,iu2)*ability_prob(iam,iaf)*mpartner(ik2,ixm,iaf,iuf,T-it,ifcm)*vnext
                                !else
                                !    evm(j,ik,ix,iam,ium,T-it,ifc)=evm(j,ik,ix,iam,ium,T-it,ifc) &
                                !        +trans_u(2,ium,iu2)*ability_prob(iam,iaf)*mpartner(ik2,ixm,iaf,iuf,T-it,ifcm)*&
                                !        LinInterp(dum,k_grid,v(:,ix,ixm,iaf,iuf,iam,iu2,T-it,ifc,ifcm),nk)
                                !end if
                                

                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do

end subroutine partest
    
subroutine partest2(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use Utilities
    use bspline_sub_module
    use GlobParams
    use PolicyFunctions_obj, only: pol_v_aux

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixm,iam,ium,iaf,iuf,iu2,iu3,j,ik2,ifc,ifcm,ik
    real(8) :: dum, dum3
    real(8) :: vnext, vnext_test
    integer :: iflag
    integer :: idx, idy, idz, iloy, iloz
    integer :: inbvx, inbvy, inbvz
    real(8) :: ww2(ky,kz),ww1(kz),ww0(3*max(kx,ky,kz))
    real(8) :: exp_m_prime, exp_f_prime

    idx=0
    idy=0
    idz=0
    inbvx=1
    inbvy=1
    inbvz=1
    iloy=1
    iloz=1    

    !Assigning the grid points
    dum3=((counter*1d0)/(nexp*nu*na*1d0))-0.00001d0
    ik=int(dum3)+1
    dum3=(((counter-(ik-1)*nexp*nu*na)*1d0)/(nu*na*1d0))-0.00001d0
    ixm=int(dum3)+1
    dum3=(((counter-(ik-1)*nexp*nu*na-(ixm-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
    ium=int(dum3)+1
    iam=counter-(ik-1)*nexp*na*nu-(ixm-1)*na*nu-(ium-1)*na


    j=1
    exp_m_prime = exp_grid(ixm,T-it)
    evm(j,ik,ixm,iam,ium,T-it,:)=0d0
    do ifcm=1,nfcm
        do iu2=1,nu
            do ik2=1,nk
                do ix=1,nexp
                    exp_f_prime = exp_grid(ix,T-it)
                    do iaf=1,na
                        do iuf=1,nu
                            do ifc=1,nfc
                                dum=k_grid(ik)+k_grid(ik2)
                                
                                vnext = pol_v_aux(iam,iu2,iaf,iuf,ifc,ifcm,MEN)%eval([dum,exp_f_prime,exp_m_prime])
                                evm(j,ik,ixm,iam,ium,T-it,ifcm) = evm(j,ik,ixm,iam,ium,T-it,ifcm)&
                                         + trans_u(MEN,ium,iu2)*ability_prob(iam,iaf)*fpartner(ik2,ix,iaf,iuf,T-it,ifc)*vnext   
                                
                                !if(dum<k_grid(nk)-0.001d0) then
                                !    call db3val(dum,exp_grid(ix,T-it),exp_grid(ixm,T-it),idx,idy,idz,&
                                !            tx,ty(:,T-it),tz(:,T-it),&
                                !            nk,nexp,nexp,kx,ky,kz,&
                                !            v_bspl(iam,iu2,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                !            inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.true.)
                                !    evm(j,ik,ixm,iam,ium,T-it,ifcm)=evm(j,ik,ixm,iam,ium,T-it,ifcm)&
                                !        +trans_u(1,ium,iu2)*ability_prob(iam,iaf)*fpartner(ik2,ix,iaf,iuf,T-it,ifc)*vnext
                                !else
                                !    evm(j,ik,ixm,iam,ium,T-it,ifcm)=evm(j,ik,ixm,iam,ium,T-it,ifcm)&
                                !        +trans_u(1,ium,iu2)*ability_prob(iam,iaf)*fpartner(ik2,ix,iaf,iuf,T-it,ifc)*&
                                !        LinInterp(dum,k_grid,v(:,ix,ixm,iam,iu2,iaf,iuf,T-it,ifc,ifcm),nk)
                                !end if
                            
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do

end subroutine partest2
    
subroutine partest3(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use Utilities
    use bspline_sub_module
    use PolicyFunctions_obj, only: pol_v_aux

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,iam,iaf,iuf,iu2,iu3,j,ik2,ifc,ium,ifcm
    real(8) :: dum3
    integer :: iflag
    integer :: iknot
    real(8), pointer :: exp_grid_ptr(:)
    real(8), pointer :: v_mar_tmp(:,:,:)

    !Assigning the grid points
    dum3=((counter*1d0)/(na*nfc*1d0))-0.00001d0
    ium=int(dum3)+1
    dum3=(((counter-(ium-1)*nfc*na)*1d0)/(na*1d0))-0.00001d0
    ifc=int(dum3)+1
    iam=counter-(ium-1)*nfc*na-(ifc-1)*na

    exp_grid_ptr => exp_grid(:,T-it)
    iknot = 0
    do ifcm=1,nfcm    
        do iaf=1,na
            do iuf=1,nu
                call db3ink(k_grid,nk,exp_grid_ptr,nexp,exp_grid_ptr,nexp,&
                            ev(:,:,:,iam,ium,iaf,iuf,T-it,ifc,ifcm),&
                            kx,ky,kz,iknot,tx,ty(:,T-it),tz(:,T-it),&
                            ev_bspl(iam,ium,iaf,iuf,ifc,ifcm)%coefs,iflag) 
                call db3ink(k_grid,nk,exp_grid_ptr,nexp,exp_grid_ptr,nexp,&
                            v(:,:,:,iam,ium,iaf,iuf,T-it,ifc,ifcm),&
                            kx,ky,kz,iknot,tx,ty(:,T-it),tz(:,T-it),&
                            v_bspl(iam,ium,iaf,iuf,ifc,ifcm)%coefs,iflag) 
                
                do j = 1, 2
                    v_mar_tmp => v_aux(:,:,:,iam,ium,iaf,iuf,ifc,ifcm,j)
                    call pol_v_aux(iam, ium, iaf, iuf, ifc, ifcm, j)%set(v_mar_tmp, k_grid, exp_grid_ptr, exp_grid_ptr)   
                end do       
                
            end do
        end do
    end do

end subroutine partest3

subroutine partest5(counter)
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    use bspline_sub_module

    implicit none

    integer, INTENT(IN) :: counter
    integer :: iam, j, ifc, ium
    real(8) :: dum3
    integer :: iknot
    integer :: iflag
    real(8), pointer:: exp_grid_ptr(:)

    !Assigning the grid points
    dum3=((counter*1d0)/(na*nfc*1d0))-0.00001d0
    ium=int(dum3)+1
    dum3=(((counter-(ium-1)*nfc*na)*1d0)/(na*1d0))-0.00001d0
    ifc=int(dum3)+1
    iam=counter-(ium-1)*nfc*na-(ifc-1)*na

    exp_grid_ptr => exp_grid(:,T-it)
    iknot = 0

    do j=1,2
        call db2ink(k_grid, nk, exp_grid_ptr, nexp, evs(j,:,:,iam,ium,T-it,ifc), &
            kx, ky, iknot, tx, ty(:,T-it), evs_bspl(j,iam,ium,ifc)%coefs, iflag)
        call db2ink(k_grid, nk, exp_grid_ptr, nexp, vs(j,:,:,iam,ium,T-it,ifc), &
            kx, ky, iknot, tx, ty(:,T-it), vs_bspl(j,iam,ium,ifc)%coefs, iflag)
    end do

end subroutine partest5   
    
subroutine partest6(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use Utilities
    use bspline_sub_module

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixd,iam,iaf,iuf,iu2,iu3,j,ik2,ifc,ium
    real(8) :: dum3
    integer :: iknot
    integer :: iflag
    real(8), pointer:: exp_grid_ptr(:)

    !Assigning the grid points
    dum3=((counter*1d0)/(na*nfc*1d0))-0.00001d0
    ium=int(dum3)+1
    dum3=(((counter-(ium-1)*nfc*na)*1d0)/(na*1d0))-0.00001d0
    ifc=int(dum3)+1
    iam=counter-(ium-1)*nfc*na-(ifc-1)*na

    iknot = 0
    exp_grid_ptr => exp_grid(:,T-it)

    do j=1,2
        call db2ink(k_grid, nk, exp_grid_ptr, nexp, evm(j,:,:,iam,ium,T-it,ifc), &
                    kx, ky, iknot, tx, ty(:,T-it), evm_bspl(j,iam,ium,ifc)%coefs, iflag)
    end do

end subroutine partest6
    
subroutine partest7(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use bspline_sub_module
    use GlobParams

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixm,ixd,iam,ium,iaf,iuf,iu2,iu3,j,ik2,j2,ifc,ifcm,ik
    real(8) :: vnext
    real(8) :: dum3
    integer :: iflag
    real(8) :: w1_d2(ky) 
    real(8) :: w0_d2(3*max(kx,ky)) 
    integer :: idx, idy, idz, iloy, iloz
    integer :: inbvx, inbvy, inbvz
    
    idx=0
    idy=0
    idz=0
    inbvx=1
    inbvy=1
    inbvz=1
    iloy=1
    iloz=1

    !Assigning the grid points
    dum3=((counter*1d0)/(nexp*nu*na*1d0))-0.00001d0
    ik=int(dum3)+1
    dum3=(((counter-(ik-1)*nexp*nu*na)*1d0)/(nu*na*1d0))-0.00001d0
    ix=int(dum3)+1
    dum3=(((counter-(ik-1)*nexp*nu*na-(ix-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
    ium=int(dum3)+1
    iam=counter-(ik-1)*nexp*na*nu-(ix-1)*na*nu-(ium-1)*na

    if(T-it>1) then
        ev(ik,ix,:,iam,ium,:,:,T-it,:,:)=0d0
        ev_aux(ik,ix,:,iam,ium,:,:,:,:,:) = 0d0
        dum3=K_grid(ik)/2d0
        do ifc=1,nfc
        do ifcm=1,nfcm
        do ixm = 1, nexp
        do iuf = 1, nu
        do iaf = 1, na
            do iu2=1,nu
                do iu3=1,nu
                    ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm) =ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)+(1d0-Probd(T-it-1))*trans_u(1,ium,iu2)*trans_u(2,iuf,iu3)*V(ik,ix,ixm,iam,iu2,iaf,iu3,T-it,ifc,ifcm)
                    ev_aux(ik,ix,ixm,iam,ium,iaf,iuf,ifc,ifcm,MEN) = ev_aux(ik,ix,ixm,iam,ium,iaf,iuf,ifc,ifcm,MEN)+(1d0-Probd(T-it-1))*trans_u(1,ium,iu2)*trans_u(2,iuf,iu3)*V_aux(ik,ix,ixm,iam,iu2,iaf,iu3,ifc,ifcm,MEN)
                    ev_aux(ik,ix,ixm,iam,ium,iaf,iuf,ifc,ifcm,WOMEN) = ev_aux(ik,ix,ixm,iam,ium,iaf,iuf,ifc,ifcm,WOMEN)+(1d0-Probd(T-it-1))*trans_u(1,ium,iu2)*trans_u(2,iuf,iu3)*V_aux(ik,ix,ixm,iam,iu2,iaf,iu3,ifc,ifcm,WOMEN)                    
                end do
            end do
            do iu3=1,nu
                call db2val(dum3,exp_grid(ix,T-it),idx,idy,&
                    tx,ty(:,T-it),nk,nexp,kx,ky,&
                    vs_bspl(2,iaf,iu3,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)
                ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm) =ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)+Probd(T-it-1)*trans_u(2,iuf,iu3)*0.5d0*vnext
                ev_aux(ik,ix,ixm,iam,ium,iaf,iuf,ifc,ifcm,WOMEN) = ev_aux(ik,ix,ixm,iam,ium,iaf,iuf,ifc,ifcm,WOMEN)+Probd(T-it-1)*trans_u(WOMEN,iuf,iu3)*vnext                
            end do
            do iu2=1,nu
                call db2val(dum3,exp_grid(ixm,T-it),idx,idy,&
                    tx,ty(:,T-it),nk,nexp,kx,ky,&
                    vs_bspl(1,iam,iu2,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)
                ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm) =ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)+Probd(T-it-1)*trans_u(1,ium,iu2)*0.5d0*vnext
                ev_aux(ik,ix,ixm,iam,ium,iaf,iuf,ifc,ifcm,MEN) = ev_aux(ik,ix,ixm,iam,ium,iaf,iuf,ifc,ifcm,MEN)+Probd(T-it-1)*trans_u(MEN,ium,iu2)*vnext    
            end do
        end do
        end do
        end do
        end do
        end do
    end if

end subroutine partest7
       
subroutine partest8(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ik, ium, iam, j, ifc, ix, iu2
    real(8) :: dum3

    !Assigning the grid points
    dum3=((counter*1d0)/(nu*na*1d0))-0.00001d0
    ik=int(dum3)+1
    dum3=(((counter-(ik-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
    ium=int(dum3)+1
    iam=counter-(ik-1)*na*nu-(ium-1)*na

    evs(:,ik,:,iam,ium,T-it,:)=0d0

    do j=1,2
        do ifc=1,nfc
            do ix = 1, nexp
                do iu2=1,nu
                    evs(j,ik,ix,iam,ium,T-it,ifc)=evs(j,ik,ix,iam,ium,T-it,ifc) &
                        +trans_u(j,ium,iu2)*Vs(j,ik,ix,iam,iu2,T-it,ifc)
                end do
            end do
        end do
    end do        

end subroutine partest8
    
subroutine partest9(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixm,ixd,iam,ium,iaf,iuf,iu2,iu3,j,ik2,ifc,ifcm,ik
    real(8) :: dum3
    real(8) :: vnext

    !Assigning the grid points
    dum3=((counter*1d0)/(nu*na*1d0))-0.00001d0
    ik=int(dum3)+1
    dum3=(((counter-(ik-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
    ium=int(dum3)+1
    iam=counter-(ik-1)*na*nu-(ium-1)*na

    !Singles
    evs_ret(2,:,ik,:,iam,ium,Tret-it+1,:)=0d0
    do j=1,2
        do ifc=1,nfc
            do ix = 1, nexp
                do iu2=1,nu
                    evs_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc) &
                        =evs_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc) &
                        +trans_u(j,ium,iu2)*vs_ret(2,j,ik,ix,iam,iu2,Tret-it+1,ifc)
                end do
            end do
        end do
    end do

    !Woman works
    ev_ret(2,1,ik,:,:,:,:,iam,ium,Tret-it+1,:,:)=0d0
    do ifc=1,nfc
        do ix = 1, nexp
            do ixm = 1, nexp
                do iaf = 1, na
                    do iu2=1,nu
                        ev_ret(2,1,ik,ix,ixm,iaf,:,iam,ium,Tret-it+1,ifc,:) &
                            =ev_ret(2,1,ik,ix,ixm,iaf,:,iam,ium,Tret-it+1,ifc,:) &
                            +trans_u(2,ium,iu2)*v_ret(2,1,ik,ix,ixm,iaf,:,iam,iu2,Tret-it+1,ifc,:)
                    end do
                end do
            end do
        end do
    end do

    !Man works
    ev_ret(1,2,ik,:,:,iam,ium,:,:,Tret-it+1,:,:)=0d0
    do ifc=1,nfcm
        do ix = 1, nexp
            do ixm = 1, nexp
                do iaf = 1, na
                    do iu2=1,nu
                        ev_ret(1,2,ik,ixm,ix,iam,ium,iaf,:,Tret-it+1,:,ifc) &
                            =ev_ret(1,2,ik,ixm,ix,iam,ium,iaf,:,Tret-it+1,:,ifc) &
                            +trans_u(1,ium,iu2)*v_ret(1,2,ik,ixm,ix,iam,iu2,iaf,:,Tret-it+1,:,ifc)
                    end do
                end do
            end do
        end do
    end do

    !Both spouses work
    ev_ret(2,2,ik,:,:,iam,ium,:,:,Tret-it+1,:,:)=0d0
    do ifc=1,nfc
        do ifcm=1,nfcm
            do ix = 1, nexp
                do ixm = 1, nexp
                    do iuf = 1, nu
                        do iaf = 1, na
                            do iu2=1,nu
                                do iu3=1,nu
                                    ev_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm) &
                                        =ev_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm) &
                                        +trans_u(1,ium,iu2)*trans_u(2,iuf,iu3) &
                                        *v_ret(2,2,ik,ix,ixm,iam,iu2,iaf,iu3,Tret-it+1,ifc,ifcm)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do

end subroutine partest9
    
subroutine partest10(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use Utilities
    use bspline_sub_module

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixd,iam,iaf,iuf,iu2,iu3,j,ik2,ifc,ifcm,ium
    real(8) :: dum3 
    integer :: iknot
    integer :: iflag
    real(8), pointer:: exp_grid_ptr(:)

    !Assigning the grid points
    dum3=((counter*1d0)/(na*nfc*1d0))-0.00001d0
    ium=int(dum3)+1
    dum3=(((counter-(ium-1)*nfc*na)*1d0)/(na*1d0))-0.00001d0
    ifc=int(dum3)+1
    iam=counter-(ium-1)*nfc*na-(ifc-1)*na
    
    exp_grid_ptr => exp_grid(:,T+Tret+1-it)

    iknot = 0

    do ik2=1,2
        do j=1,2
            call db2ink(k_grid, nk, exp_grid_ptr, nexp, evs_ret(ik2,j,:,:,iam,ium,Tret-it+1,ifc), &
                        kx, ky, iknot, tx, ty(:,T+Tret+1-it), evs_ret_bspl(ik2,j,iam,ium,ifc)%coefs, iflag)
            call db2ink(k_grid, nk, exp_grid_ptr, nexp, edcs_ret(ik2,j,:,:,iam,ium,Tret-it+1,ifc), &
                        kx, ky, iknot, tx, ty(:,T+Tret+1-it), edcs_ret_bspl(ik2,j,iam,ium,ifc)%coefs, iflag)             
        end do
    end do

    do ik2=1,2
    do j=1,2
    do ifcm=1,nfcm
    do iaf=1,na
        do iuf=1,nu
            call db3ink(k_grid,nk,exp_grid_ptr,nexp,exp_grid_ptr,nexp,&
                        ev_ret(ik2,j,:,:,:,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm),&
                        kx,ky,kz,iknot,tx,ty(:,T+Tret+1-it),tz(:,T+Tret+1-it),&
                        ev_ret_bspl(ik2,j,iam,ium,iaf,iuf,ifc,ifcm)%coefs,iflag) 
            call db3ink(k_grid,nk,exp_grid_ptr,nexp,exp_grid_ptr,nexp,&
                        edc_ret(ik2,j,:,:,:,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm),&
                        kx,ky,kz,iknot,tx,ty(:,T+Tret+1-it),tz(:,T+Tret+1-it),&
                        edc_ret_bspl(ik2,j,iam,ium,iaf,iuf,ifc,ifcm)%coefs,iflag)         
        end do
    end do
    end do    
    end do
    end do

end subroutine partest10
