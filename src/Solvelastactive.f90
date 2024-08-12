subroutine Solvelastactive(counter)

    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    use bspline_sub_module
    use GlobParams

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ik,ix,ixm,ixd,iam,ium,iaf,iuf,iu2,iu3,j,ifc,ifcm
    real(8) :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8) :: ces,cus,kes,kus,nes,nus,ves,vus
    real(8), dimension (:), allocatable :: Expdum
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y,r_ret_next
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2,pnt2(3),pnt1(2)
    real(8) :: vnext, exp_grid_dum(nexp), INTERP2D(nk,nexp), INTERP3D(nk,nexp,nexp)
    real(8) :: dd1, dd2, dd3, dd4, dd5

    integer :: iflag 
    real(8) :: vnext_test
    integer :: idx, idy, idz, iloy, iloz
    integer :: inbvx, inbvy, inbvz
    real(8) :: ww2(ky,kz),ww1(kz),ww0(3*max(kx,ky,kz))
    real(8) :: w1_d2(ky) 
    real(8) :: w0_d2(3*max(kx,ky)) 
    real(8) :: dum2_test, nonlab_inc, aftertax_lab_inc
    real(8) :: lfpm_e, lfpf_e, lfpm_u, lfpf_u

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

    exp_grid_dum=exp_grid(:,T+1)
    r_ret_next=((1d0+r)/OmegaRet2(1))-1d0

    ce=0d0
    cu=0d0
    ke=0d0
    ku=0d0
    nem=0d0
    nef=0d0
    num=0d0
    nuf=0d0
    ve=0d0
    vu=0d0

    ces=0d0
    cus=0d0
    kes=0d0
    kus=0d0
    nes=0d0
    nus=0d0
    ves=0d0
    vus=0d0


    ! Last period of active life
    !if employed
    do ifc=1,nfc
        do ifcm=1,nfcm
            do iaf = 1, na
                do iuf = 1, nu
                    do ixm = 1, nexp
                        wagem = (1d0/(1d0+exp(kappa*(T-agestart))))*wage(1,a(1,iam),exp_grid(ixm,T),u(1,ium))/(1d0+t_employer)
                        wagef = (1d0/(1d0+exp(kappa*(T-agestart))))*wage(2,a(2,iaf),exp_grid(ix,T),u(2,iuf))/(1d0+t_employer)
                        P1=0.001d0
                        nonlab_inc = (k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum
                        aftertax_lab_inc = after_tax_labor_inc_married(wagem+wagef)
                        P4 = min((k_grid(nk)-0.001d0)/(1d0+tc), (nonlab_inc + (1d0-t_const)*(aftertax_lab_inc - (wagem+wagef)*t_employee))/(1d0+tc))                   

                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                            pnt2 = (/P2, wagem, wagef/)
                            dum4 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborm, nc, nw, nw, pnt2),1d-10),1d0)
                            dum5 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborf, nc, nw, nw, pnt2),1d-10),1d0)
                            y=dum4*wagem+dum5*wagef

                            dum2 = (nonlab_inc + (1d0-t_const)*(after_tax_labor_inc_married(y)-y*t_employee)-P2*(1d0+tc))/(1d0+mu)

                            if(dum2<0.0001d0) then
                                V2=-999999999d0
                            elseif(dum2>k_grid(nk)-0.001d0) then
                                V2=Uc(P2)+Ul(dum4,dum5)-fc(1,ifc)-fcm(1,ifcm)
                                pnt2=(/dum2, exp_grid(ix,T)+1d0, exp_grid(ixm,T)+1d0/)
                                INTERP3D=ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm)
                                vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                V2=V2+beta*OmegaActive(T)*vnext
                            else
                                V2=Uc(P2)+Ul(dum4,dum5)-fc(1,ifc)-fcm(1,ifcm)
                                call db3val(dum2,exp_grid(ix,T)+1d0, exp_grid(ixm,T)+1d0,idx,idy,idz,&
                                        tx,ty(:,T+1),tz(:,T+1),&
                                        nk,nexp,nexp,kx,ky,kz,&
                                        ev_ret_bspl(2,2,iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                        inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.false.)
                                V2=V2+beta*OmegaActive(T)*vnext
                            end if

                            pnt2 = (/P3, wagem, wagef/)
                            dum4 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborm, nc, nw, nw, pnt2),1d-10),1d0)
                            dum5 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborf, nc, nw, nw, pnt2),1d-10),1d0)
                            y=dum4*wagem+dum5*wagef

                            dum2 = (nonlab_inc + (1d0-t_const)*(after_tax_labor_inc_married(y)-y*t_employee)-P3*(1d0+tc))/(1d0+mu)
                            
                            if(dum2<0.0001d0) then
                                V3=-999999999d0
                            elseif(dum2>k_grid(nk)-0.001d0) then
                                V3=Uc(P3)+Ul(dum4,dum5)-fc(1,ifc)-fcm(1,ifcm)
                                pnt2=(/dum2, exp_grid(ix,T)+1d0, exp_grid(ixm,T)+1d0/)
                                INTERP3D=ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm)
                                vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                V3=V3+beta*OmegaActive(T)*vnext
                            else
                                V3=Uc(P3)+Ul(dum4,dum5)-fc(1,ifc)-fcm(1,ifcm)
                                call db3val(dum2,exp_grid(ix,T)+1d0, exp_grid(ixm,T)+1d0,idx,idy,idz,&
                                        tx,ty(:,T+1),tz(:,T+1),&
                                        nk,nexp,nexp,kx,ky,kz,&
                                        ev_ret_bspl(2,2,iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                        inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.false.)
                                V3=V3+beta*OmegaActive(T)*vnext
                            end if

                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-6) exit
                        end do

                        if(dum4<minhours) then
                            V2=-999999999d0
                        end if
                        if(dum5<minhours) then
                            V2=-999999999d0
                        end if
                        
                        Ve=V2
                        ke=dum2
                        nem=dum4
                        nef=dum5
                        ce=P2
                        
                        lfpm_e = 1d0
                        lfpf_e = 1d0
                        c_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M1,LFP_F1)=ce
                        k_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M1,LFP_F1)=ke
                        nm_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M1,LFP_F1)=dum4
                        nf_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M1,LFP_F1)=dum5   
                        v_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M1,LFP_F1,MEN)=V2
                        v_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M1,LFP_F1,WOMEN)=V2                        
                        

                        !If female unemployed
                        P1=0.001d0
                        nonlab_inc = (k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+Unemp_benefit
                        aftertax_lab_inc = after_tax_labor_inc_married(wagem)
                        P4 = min((k_grid(nk)-0.001d0)/(1d0+tc), (nonlab_inc + (1d0-t_const)*(aftertax_lab_inc - wagem*t_employee))/(1d0+tc))                   

                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                            pnt1 = (/P2, wagem/)
                            dum4 = min(max(bilin_interp(c_grid, wage_grid, labormwork, nc, nw, pnt1),0d0),1d0)
                            y=dum4*wagem
                            
                            dum2 = (nonlab_inc + (1d0-t_const)*(after_tax_labor_inc_married(y)-y*t_employee)-P2*(1d0+tc))/(1d0+mu)

                            if(dum2<0.0001d0) then
                                V2=-999999999d0
                            elseif(dum2>k_grid(nk)-0.001d0) then
                                V2=Uc(P2)+Ul(dum4,0d0)-fcm(1,ifcm)
                                pnt2=(/dum2, exp_grid(ix,T)*(1d0-deltaexp), exp_grid(ixm,T)+1d0/)
                                INTERP3D=ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm)
                                vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                V2=V2+beta*OmegaActive(T)*vnext
                            else
                                V2=Uc(P2)+Ul(dum4,0d0)-fcm(1,ifcm)
                                !vnext = D_BS3VL(dum2, exp_grid(ix,T)*(1d0-deltaexp), exp_grid(ixm,T)+1d0, KORDER, EXPORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T+1),EXP_KNOT(:,T+1), nk, nexp, nexp, ev_ret_spln_coefs(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm))
                                call db3val(dum2,exp_grid(ix,T)*(1d0-deltaexp),exp_grid(ixm,T)+1d0,idx,idy,idz,&
                                        tx,ty(:,T+1),tz(:,T+1),&
                                        nk,nexp,nexp,kx,ky,kz,&
                                        ev_ret_bspl(2,2,iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                        inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.false.)
                                V2=V2+beta*OmegaActive(T)*vnext
                            end if

                            pnt1 = (/P3, wagem/)
                            dum4 = min(max(bilin_interp(c_grid, wage_grid, labormwork, nc, nw, pnt1),0d0),1d0)
                            y=dum4*wagem

                            dum2 = (nonlab_inc + (1d0-t_const)*(after_tax_labor_inc_married(y)-y*t_employee)-P3*(1d0+tc))/(1d0+mu)                           

                            if(dum2<0.0001d0) then
                                V3=-999999999d0
                            elseif(dum2>k_grid(nk)-0.001d0) then
                                V3=Uc(P3)+Ul(dum4,0d0)-fcm(1,ifcm)
                                pnt2=(/dum2, exp_grid(ix,T)*(1d0-deltaexp), exp_grid(ixm,T)+1d0/)
                                INTERP3D=ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm)
                                vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                V3=V3+beta*OmegaActive(T)*vnext
                            else
                                V3=Uc(P3)+Ul(dum4,0d0)-fcm(1,ifcm)
                                !vnext = D_BS3VL(dum2, exp_grid(ix,T)*(1d0-deltaexp), exp_grid(ixm,T)+1d0, KORDER, EXPORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T+1),EXP_KNOT(:,T+1), nk, nexp, nexp, ev_ret_spln_coefs(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm))
                                call db3val(dum2,exp_grid(ix,T)*(1d0-deltaexp),exp_grid(ixm,T)+1d0,idx,idy,idz,&
                                        tx,ty(:,T+1),tz(:,T+1),&
                                        nk,nexp,nexp,kx,ky,kz,&
                                        ev_ret_bspl(2,2,iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                        inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.false.)
                                V3=V3+beta*OmegaActive(T)*vnext
                            end if

                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-6) exit
                        end do

                        if(dum4<minhours) then
                            V2=-999999999d0
                        end if
                        
                        Vu=V2
                        ku=dum2
                        num=dum4
                        nuf=0d0
                        cu=P2
                        
                        lfpm_u=1d0
                        lfpf_u=0d0
                        c_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M1,LFP_F0)=cu
                        k_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M1,LFP_F0)=ku  
                        nm_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M1,LFP_F0)=dum4
                        nf_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M1,LFP_F0)=0d0  
                        v_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M1,LFP_F0,MEN)=V2
                        v_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M1,LFP_F0,WOMEN)=V2                        
                        
                        
                        !If male unemployed
                        P1=0.001d0
                        nonlab_inc = (k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+Unemp_benefit
                        aftertax_lab_inc = after_tax_labor_inc_married(wagef)
                        P4 = min((k_grid(nk)-0.001d0)/(1d0+tc), (nonlab_inc + (1d0-t_const)*(aftertax_lab_inc - (wagef)*t_employee))/(1d0+tc))                 

                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                            pnt1 = (/P2, wagef/)
                            dum4 = min(max(bilin_interp(c_grid, wage_grid, laborfwork, nc, nw, pnt1),0d0),1d0)
                            y=dum4*wagef

                            dum2 = (nonlab_inc + (1d0-t_const)*(after_tax_labor_inc_married(y)-y*t_employee)-P2*(1d0+tc))/(1d0+mu)                           

                            if(dum2<0.0001d0) then
                                V2=-999999999d0
                            elseif(dum2>k_grid(nk)-0.001d0) then
                                V2=Uc(P2)+Ul(0d0,dum4)-fc(1,ifc)
                                pnt2=(/dum2, exp_grid(ix,T)+1d0, exp_grid(ixm,T)*(1d0-deltaexp)/)
                                INTERP3D=ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm)
                                vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                V2=V2+beta*OmegaActive(T)*vnext
                            else
                                V2=Uc(P2)+Ul(0d0,dum4)-fc(1,ifc)
                                call db3val(dum2,exp_grid(ix,T)+1d0,exp_grid(ixm,T)*(1d0-deltaexp),idx,idy,idz,&
                                        tx,ty(:,T+1),tz(:,T+1),&
                                        nk,nexp,nexp,kx,ky,kz,&
                                        ev_ret_bspl(2,2,iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                        inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.false.)
                                V2=V2+beta*OmegaActive(T)*vnext
                            end if

                            pnt1 = (/P3, wagef/)
                            dum4 = min(max(bilin_interp(c_grid, wage_grid, laborfwork, nc, nw, pnt1),0d0),1d0)
                            y=dum4*wagef

                            dum2 = (nonlab_inc + (1d0-t_const)*(after_tax_labor_inc_married(y)-y*t_employee)-P3*(1d0+tc))/(1d0+mu)                           

                            if(dum2<0.0001d0) then
                                V3=-999999999d0
                            elseif(dum2>k_grid(nk)-0.001d0) then
                                V3=Uc(P3)+Ul(0d0,dum4)-fc(1,ifc)
                                pnt2=(/dum2, exp_grid(ix,T)+1d0, exp_grid(ixm,T)*(1d0-deltaexp)/)
                                INTERP3D=ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm)
                                vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                V3=V3+beta*OmegaActive(T)*vnext
                            else
                                V3=Uc(P3)+Ul(0d0,dum4)-fc(1,ifc)
                                !vnext = D_BS3VL(dum2, exp_grid(ix,T)+1d0, exp_grid(ixm,T)*(1d0-deltaexp), KORDER, EXPORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T+1),EXP_KNOT(:,T+1), nk, nexp, nexp, ev_ret_spln_coefs(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm))
                                call db3val(dum2,exp_grid(ix,T)+1d0,exp_grid(ixm,T)*(1d0-deltaexp),idx,idy,idz,&
                                        tx,ty(:,T+1),tz(:,T+1),&
                                        nk,nexp,nexp,kx,ky,kz,&
                                        ev_ret_bspl(2,2,iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                        inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.false.)
                                V3=V3+beta*OmegaActive(T)*vnext
                            end if

                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-6) exit
                        end do

                        if(dum4<minhours) then
                            V2=-999999999d0
                        end if
                        
                        if (V2 >= vu) then
                            Vu=V2
                            ku=dum2
                            num=0d0
                            nuf=dum4
                            cu=P2
                            lfpm_u=0d0
                            lfpf_u=1d0                              
                        end if
                        c_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M0,LFP_F1)=P2
                        k_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M0,LFP_F1)=dum2 
                        nm_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M0,LFP_F1)=0d0
                        nf_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M0,LFP_F1)=dum4    
                        v_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M0,LFP_F1,MEN)=V2
                        v_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M0,LFP_F1,WOMEN)=V2                        

                        
                        !If both spouses unemployed
                        nonlab_inc = (k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk)) + lumpsum + 2d0*Unemp_benefit + after_tax_labor_inc_married(0d0)
                        P1=0.001d0
                        P4 = min( k_grid(nk)-0.001d0, nonlab_inc )/ (1d0+tc)
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                            dum2 = (nonlab_inc - P2*(1d0+tc))/(1d0+mu)
                            if(dum2<0.0001d0) then
                                V2=-999999999d0
                            elseif(dum2>k_grid(nk)-0.001d0) then
                                V2=Uc(P2)+Ul(0d0,0d0)
                                pnt2=(/dum2, exp_grid(ix,T)*(1d0-deltaexp), exp_grid(ixm,T)*(1d0-deltaexp)/)
                                INTERP3D=ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm)
                                vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                V2=V2+beta*OmegaActive(T)*vnext
                            else
                                V2=Uc(P2)+Ul(0d0,0d0)
                                call db3val(dum2,exp_grid(ix,T)*(1d0-deltaexp),exp_grid(ixm,T)*(1d0-deltaexp),idx,idy,idz,&
                                        tx,ty(:,T+1),tz(:,T+1),&
                                        nk,nexp,nexp,kx,ky,kz,&
                                        ev_ret_bspl(2,2,iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                        inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.false.)
                                V2=V2+beta*OmegaActive(T)*vnext
                            end if

                            dum2 = (nonlab_inc - P3*(1d0+tc))/(1d0+mu)
                            if(dum2<0.0001d0) then
                                V3=-999999999d0
                            elseif(dum2>k_grid(nk)-0.001d0) then
                                V3=Uc(P3)+Ul(0d0,0d0)
                                pnt2=(/dum2, exp_grid(ix,T)*(1d0-deltaexp), exp_grid(ixm,T)*(1d0-deltaexp)/)
                                INTERP3D=ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm)
                                vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                V3=V3+beta*OmegaActive(T)*vnext
                            else
                                V3=Uc(P3)+Ul(0d0,0d0)
                                call db3val(dum2,exp_grid(ix,T)*(1d0-deltaexp),exp_grid(ixm,T)*(1d0-deltaexp),idx,idy,idz,&
                                        tx,ty(:,T+1),tz(:,T+1),&
                                        nk,nexp,nexp,kx,ky,kz,&
                                        ev_ret_bspl(2,2,iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                        inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.false.)
                                V3=V3+beta*OmegaActive(T)*vnext
                            end if

                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-6) exit
                        end do

                        if (V2 >= vu) then
                            Vu=V2
                            ku=dum2
                            num=0d0
                            nuf=0d0
                            cu=P2
                            lfpm_u=0d0
                            lfpf_u=0d0                             
                        end if
                        
                        c_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M0,LFP_F0)=P2
                        k_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M0,LFP_F0)=dum2    
                        nm_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M0,LFP_F0)=0d0
                        nf_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M0,LFP_F0)=0d0   
                        v_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M0,LFP_F0,MEN)=V2
                        v_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm,LFP_M0,LFP_F0,WOMEN)=V2                        

                        if (ve >= vu) then
                            v(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=ve
                            c(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=ce
                            k(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=ke
                            nm(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=nem
                            nf(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=nef
                            
                            lfpm(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=lfpm_e
                            lfpf(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=lfpf_e                              
                            
                        else
                            v(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=vu
                            c(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=cu
                            k(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=ku
                            nm(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=num
                            nf(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=nuf
                            
                            lfpm(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=lfpm_u
                            lfpf(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=lfpf_u                               
                        end if

                        v_aux(ik,ix,ixm,iam,ium,iaf,iuf,ifc,ifcm,:) = v(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)


                    end do
                end do
            end do
        end do
    end do

    !Singles    
    j=2
    do ifc=1,nfc
        wagef = (1d0/(1d0+exp(kappa*(T-agestart))))*wage(2,a(2,iam),exp_grid(ix,T),u(2,ium))/(1d0+t_employer)
        P1=0.001d0
        nonlab_inc = (k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0
        aftertax_lab_inc = after_tax_labor_inc_single(wagef)
        P4 = min((k_grid(nk)-0.001d0)/(1d0+tc), (nonlab_inc + (1d0-t_const)*(aftertax_lab_inc - (wagef)*t_employee))/(1d0+tc))
        do
            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

            pnt1 = (/P2, wagef/)
            dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglef, nc, nw, pnt1),0d0),1d0)
            y=dum4*wagef

            dum2 = (nonlab_inc + (1d0-t_const)*(after_tax_labor_inc_single(y)-y*t_employee)-P2*(1d0+tc))/(1d0+mu)        

            if(dum2<0.0001d0) then
                V2=-999999999d0
            elseif(dum2>k_grid(nk)) then
                pnt1=(/dum2, exp_grid(ix,T+1)+1d0/)
                INTERP2D=evs_ret(2,j,:,:,iam,ium,1,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V2=Uc(P2)-chifs*(dum4**(1d0+etaf))/(1d0+etaf)-fc(2,ifc)
                V2=V2+beta*OmegaActive(T)*vnext
            else
                V2=Uc(P2)-chifs*(dum4**(1d0+etaf))/(1d0+etaf)-fc(2,ifc)
                call db2val(dum2,exp_grid(ix,T)+1d0,idx,idy,&
                    tx,ty(:,T+1),nk,nexp,kx,ky,&
                    evs_ret_bspl(2,j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.false.)
                V2=V2+beta*OmegaActive(T)*vnext
            end if

            pnt1 = (/P3, wagef/)
            dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglef, nc, nw, pnt1),0d0),1d0)
            y=dum4*wagef

            dum2 = (nonlab_inc + (1d0-t_const)*(after_tax_labor_inc_single(y)-y*t_employee)-P3*(1d0+tc))/(1d0+mu)            

            if(dum2<0.0001d0) then
                V3=-999999999d0
            elseif(dum2>k_grid(nk)) then
                pnt1=(/dum2, exp_grid(ix,T+1)+1d0/)
                INTERP2D=evs_ret(2,j,:,:,iam,ium,1,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V3=Uc(P3)-chifs*(dum4**(1d0+etaf))/(1d0+etaf)-fc(2,ifc)
                V3=V3+beta*OmegaActive(T)*vnext
            else
                call db2val(dum2,exp_grid(ix,T)+1d0,idx,idy,&
                    tx,ty(:,T+1),nk,nexp,kx,ky,&
                    evs_ret_bspl(2,j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.false.)
                V3=Uc(P3)-chifs*(dum4**(1d0+etaf))/(1d0+etaf)-fc(2,ifc)
                V3=V3+beta*OmegaActive(T)*vnext
            end if

            if (V2 < V3) then
                P1=P2
            else
                P4=P3
            end if
            if((P4-P1)<1d-6) exit
        end do

        if(dum4<minhours) then
            V2=-999999999d0
        end if
                    
        Ves=V2
        kes=dum2
        nes=dum4
        ces=P2
        
        cs_lfp(j,ik,ix,iam,ium,T,ifc,LFP_1)=P2
        ks_lfp(j,ik,ix,iam,ium,T,ifc,LFP_1)=dum2    
        vs_lfp(j,ik,ix,iam,ium,T,ifc,LFP_1)=V2
        ns_lfp(j,ik,ix,iam,ium,T,ifc,LFP_1)=dum4

        !If single female unemployed
        nonlab_inc = (k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk)) + lumpsum*0.5d0 + Unemp_benefit + after_tax_labor_inc_single(0d0)
        P1=0.001d0
        P4 = min( k_grid(nk)-0.001d0, nonlab_inc )/(1d0+tc)
        do
            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

            dum2 = (nonlab_inc - P2*(1d0+tc))/(1d0+mu)

            if(dum2<0.0001d0) then
                V2=-999999999d0
            elseif(dum2>k_grid(nk)) then
                pnt1=(/dum2, exp_grid(ix,T+1)*(1d0-deltaexp)/)
                INTERP2D=evs_ret(2,j,:,:,iam,ium,1,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V2=Uc(P2)
                V2=V2+beta*OmegaActive(T)*vnext
            else
                call db2val(dum2,exp_grid(ix,T)*(1d0-deltaexp),idx,idy,&
                    tx,ty(:,T+1),nk,nexp,kx,ky,&
                    evs_ret_bspl(2,j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.false.)
                V2=Uc(P2)
                V2=V2+beta*OmegaActive(T)*vnext
            end if

            dum2 = (nonlab_inc - P3*(1d0+tc))/(1d0+mu)
            if(dum2<0.0001d0) then
                V3=-999999999d0
            elseif(dum2>k_grid(nk)) then
                pnt1=(/dum2, exp_grid(ix,T+1)*(1d0-deltaexp)/)
                INTERP2D=evs_ret(2,j,:,:,iam,ium,1,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V3=Uc(P3)
                V3=V3+beta*OmegaActive(T)*vnext
            else
                call db2val(dum2,exp_grid(ix,T)*(1d0-deltaexp),idx,idy,&
                    tx,ty(:,T+1),nk,nexp,kx,ky,&
                    evs_ret_bspl(2,j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.false.)
                V3=Uc(P3)
                V3=V3+beta*OmegaActive(T)*vnext
            end if

            if (V2 < V3) then
                P1=P2
            else
                P4=P3
            end if
            if((P4-P1)<1d-6) exit
        end do

        Vus=V2
        kus=dum2
        nus=0d0
        cus=P2
        
        cs_lfp(j,ik,ix,iam,ium,T,ifc,LFP_0)=P2
        ks_lfp(j,ik,ix,iam,ium,T,ifc,LFP_0)=dum2  
        vs_lfp(j,ik,ix,iam,ium,T,ifc,LFP_0)=V2  
        ns_lfp(j,ik,ix,iam,ium,T,ifc,LFP_0)=0d0

        if (ves >= vus) then
            vs(j,ik,ix,iam,ium,T,ifc)=ves
            cs(j,ik,ix,iam,ium,T,ifc)=ces
            Uprimes(j,ik,ix,iam,ium,T,ifc)=dUc(ces)
            ks(j,ik,ix,iam,ium,T,ifc)=kes
            ns(j,ik,ix,iam,ium,T,ifc)=nes
            
            ns(j,ik,ix,iam,ium,T,ifc)=nes
            lfps(j,ik,ix,iam,ium,T,ifc)=1d0            
        else
            vs(j,ik,ix,iam,ium,T,ifc)=vus
            cs(j,ik,ix,iam,ium,T,ifc)=cus
            Uprimes(j,ik,ix,iam,ium,T,ifc)=dUc(cus)
            ks(j,ik,ix,iam,ium,T,ifc)=kus
            ns(j,ik,ix,iam,ium,T,ifc)=nus
            
            ns(j,ik,ix,iam,ium,T,ifc)=nus
            lfps(j,ik,ix,iam,ium,T,ifc)=0d0            
        end if
    end do


    !Men   
    j=1
    do ifc=1,nfcm
        wagem = (1d0/(1d0+exp(kappa*(T-agestart))))*wage(1,a(1,iam),exp_grid(ix,T),u(1,ium))/(1d0+t_employer)
        P1=0.001d0
        nonlab_inc = (k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0
        aftertax_lab_inc = after_tax_labor_inc_single(wagem)
        P4 = min((k_grid(nk)-0.001d0)/(1d0+tc), (nonlab_inc + (1d0-t_const)*(aftertax_lab_inc - (wagem)*t_employee))/(1d0+tc))
        do
            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

            pnt1 = (/P2, wagem/)
            dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglem, nc, nw, pnt1),0d0),1d0)
            y=dum4*wagem

            dum2 = (nonlab_inc + (1d0-t_const)*(after_tax_labor_inc_single(y)-y*t_employee)-P2*(1d0+tc))/(1d0+mu)            

            if(dum2<0.0001d0) then
                V2=-999999999d0
            elseif(dum2>k_grid(nk)) then
                pnt1=(/dum2, exp_grid(ix,T+1)+1d0/)
                INTERP2D=evs_ret(2,j,:,:,iam,ium,1,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V2=Uc(P2)-chims*(dum4**(1d0+etam))/(1d0+etam)-fcm(2,ifc)
                V2=V2+beta*OmegaActive(T)*vnext
            else
                V2=Uc(P2)-chims*(dum4**(1d0+etam))/(1d0+etam)-fcm(2,ifc)
                call db2val(dum2,exp_grid(ix,T)+1d0,idx,idy,&
                    tx,ty(:,T+1),nk,nexp,kx,ky,&
                    evs_ret_bspl(2,j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.false.)
                V2=V2+beta*OmegaActive(T)*vnext
            end if

            pnt1 = (/P3, wagem/)
            dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglem, nc, nw, pnt1),0d0),1d0)
            y=dum4*wagem

            dum2 = (nonlab_inc + (1d0-t_const)*(after_tax_labor_inc_single(y)-y*t_employee)-P3*(1d0+tc))/(1d0+mu)             

            if(dum2<0.0001d0) then
                V3=-999999999d0
            elseif(dum2>k_grid(nk)) then
                pnt1=(/dum2, exp_grid(ix,T+1)+1d0/)
                INTERP2D=evs_ret(2,j,:,:,iam,ium,1,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V3=Uc(P3)-chims*(dum4**(1d0+etam))/(1d0+etam)-fcm(2,ifc)
                V3=V3+beta*OmegaActive(T)*vnext
            else
                !vnext = D_BS2VL(dum2, exp_grid(ix,T)+1d0, KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T+1), nk, nexp,evs_ret_spln_coefs(2,j,:,:,iam,ium,1,ifc))
                call db2val(dum2,exp_grid(ix,T)+1d0,idx,idy,&
                    tx,ty(:,T+1),nk,nexp,kx,ky,&
                    evs_ret_bspl(2,j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.false.)
                V3=Uc(P3)-chims*(dum4**(1d0+etam))/(1d0+etam)-fcm(2,ifc)
                V3=V3+beta*OmegaActive(T)*vnext
            end if

            if (V2 < V3) then
                P1=P2
            else
                P4=P3
            end if
            if((P4-P1)<1d-6) exit
        end do

        if(dum4<minhours) then
            V2=-999999999d0
        end if
        
        Ves=V2
        kes=dum2
        nes=dum4
        ces=P2
        
        cs_lfp(j,ik,ix,iam,ium,T,ifc,LFP_1)=P2
        ks_lfp(j,ik,ix,iam,ium,T,ifc,LFP_1)=dum2  
        vs_lfp(j,ik,ix,iam,ium,T,ifc,LFP_1)=V2     
        ns_lfp(j,ik,ix,iam,ium,T,ifc,LFP_1)=dum4


        !If single male unemployed
        nonlab_inc = (k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk)) + lumpsum*0.5d0 + Unemp_benefit + after_tax_labor_inc_single(0d0)
        P1=0.001d0
        P4 = min( k_grid(nk)-0.001d0, nonlab_inc )/(1d0+tc)
        do
            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

            dum2 = (nonlab_inc - P2*(1d0+tc))/(1d0+mu)

            if(dum2<0.0001d0) then
                V2=-999999999d0
            elseif(dum2>k_grid(nk)) then
                pnt1=(/dum2, exp_grid(ix,T+1)*(1d0-deltaexp)/)
                INTERP2D=evs_ret(2,j,:,:,iam,ium,1,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V2=Uc(P2)
                V2=V2+beta*OmegaActive(T)*vnext
            else
                call db2val(dum2,exp_grid(ix,T)*(1d0-deltaexp),idx,idy,&
                    tx,ty(:,T+1),nk,nexp,kx,ky,&
                    evs_ret_bspl(2,j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.false.)
                V2=Uc(P2)
                V2=V2+beta*OmegaActive(T)*vnext
            end if

            dum2 = (nonlab_inc - P3*(1d0+tc))/(1d0+mu)

            if(dum2<0.0001d0) then
                V3=-999999999d0
            elseif(dum2>k_grid(nk)) then
                pnt1=(/dum2, exp_grid(ix,T+1)*(1d0-deltaexp)/)
                INTERP2D=evs_ret(2,j,:,:,iam,ium,1,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V3=Uc(P3)
                V3=V3+beta*OmegaActive(T)*vnext
            else
                call db2val(dum2,exp_grid(ix,T)*(1d0-deltaexp),idx,idy,&
                    tx,ty(:,T+1),nk,nexp,kx,ky,&
                    evs_ret_bspl(2,j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.false.)
                V3=Uc(P3)
                V3=V3+beta*OmegaActive(T)*vnext
            end if

            if (V2 < V3) then
                P1=P2
            else
                P4=P3
            end if
            if((P4-P1)<1d-6) exit
        end do

        Vus=V2
        kus=dum2
        nus=0d0
        cus=P2
        
        cs_lfp(j,ik,ix,iam,ium,T,ifc,LFP_0)=P2
        ks_lfp(j,ik,ix,iam,ium,T,ifc,LFP_0)=dum2  
        vs_lfp(j,ik,ix,iam,ium,T,ifc,LFP_0)=V2  
        ns_lfp(j,ik,ix,iam,ium,T,ifc,LFP_0)=0d0

        if (ves >= vus) then
            vs(j,ik,ix,iam,ium,T,ifc)=ves
            cs(j,ik,ix,iam,ium,T,ifc)=ces
            Uprimes(j,ik,ix,iam,ium,T,ifc)=dUc(ces)
            ks(j,ik,ix,iam,ium,T,ifc)=kes
            ns(j,ik,ix,iam,ium,T,ifc)=nes
            lfps(j,ik,ix,iam,ium,T,ifc)=1d0
        else
            vs(j,ik,ix,iam,ium,T,ifc)=vus
            cs(j,ik,ix,iam,ium,T,ifc)=cus
            Uprimes(j,ik,ix,iam,ium,T,ifc)=dUc(cus)
            ks(j,ik,ix,iam,ium,T,ifc)=kus
            ns(j,ik,ix,iam,ium,T,ifc)=nus
            lfps(j,ik,ix,iam,ium,T,ifc)=0d0
        end if
    end do


end subroutine Solvelastactive

