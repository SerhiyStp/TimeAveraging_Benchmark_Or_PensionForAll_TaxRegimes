subroutine SolveActiveLife(counter)
    !This subroutine computes optimal policies durin active worklife.
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    use bspline_sub_module
    use PolicyFunctions_obj, only: policy_fn_3d, pol_ev_aux
    use GlobParams

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixm,ixd,iam,ium,iaf,iuf,iu2,iu3,tprint,ikd,j,ifc,ifcm,ik
    real(8) :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8) :: ces,cus,kes,kus,nes,nus,ves,vus
    real(8) :: dum3,dum4,dum5,dum6,y
    real(8) :: P1,P2,P3,P4,V2,V3,dum2,pnt2(3),pnt1(2)
    real(8) :: vnext, exp_grid_dum(nexp), INTERP2D(nk,nexp), INTERP3D(nk,nexp,nexp)

    real(8) :: vnext_test
    integer :: idx, idy, idz, iloy, iloz
    integer :: inbvx, inbvy, inbvz
    integer :: iflag
    real(8) :: ww2(ky,kz),ww1(kz),ww0(3*max(kx,ky,kz))
    real(8) :: w1_d2(ky) 
    real(8) :: w0_d2(3*max(kx,ky)) 
    real(8) :: nonlab_inc, aftertax_lab_inc
    real(8) :: lfpm_e, lfpf_e, lfpm_u, lfpf_u
    real(8) :: lfp_e, lfp_u    
    real(8) :: Ucur, vnext_m, vnext_f, ve_m_aux, ve_f_aux, vu_m_aux, vu_f_aux
    real(8) :: vtmp_m_aux, vtmp_f_aux
    type(policy_fn_3d), pointer :: ev_aux_ptr
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

    exp_grid_dum=exp_grid(:,T+1-it)


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

    !if both spouses employed

    do ifc=1,nfc
        do ifcm=1,nfcm
            do iaf = 1, na
                do iuf = 1, nu
                    do ixm = 1, nexp

                        !Finding optimal capital by golden search
                        wagem = (1d0/(1d0+exp(kappa*(T-(it+agestart)))))*wage(1,a(1,iam),exp_grid(ixm,T-it),u(1,ium))/(1d0+t_employer)
                        wagef = (1d0/(1d0+exp(kappa*(T-(it+agestart)))))*wage(2,a(2,iaf),exp_grid(ix,T-it),u(2,iuf))/(1d0+t_employer)
                        P1=0.001d0
                        nonlab_inc = (k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum
                        aftertax_lab_inc = after_tax_labor_inc_married(wagem+wagef)
                        P4 = min((k_grid(nk)-0.001d0)/(1d0+tc), (nonlab_inc + (1d0-t_const)*(aftertax_lab_inc - (wagem+wagef)*t_employee))/(1d0+tc))                   
                        
                        exp_m_prime = exp_grid(ixm,T-it)+1d0
                        exp_f_prime = exp_grid(ix,T-it)+1d0
                        
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                            pnt2 = (/P2, wagem, wagef/)
                            dum4 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborm, nc, nw, nw, pnt2),1d-10),1d0)
                            dum5 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborf, nc, nw, nw, pnt2),1d-10),1d0)
                            y=dum4*wagem+dum5*wagef

                            dum2 = (nonlab_inc + (1d0-t_const)*(after_tax_labor_inc_married(y)-y*t_employee) - P2*(1d0+tc))/(1d0+mu)

                            if(dum2<0.0001d0) then
                                V2=-999999999d0
                            elseif(dum2>k_grid(nk)-0.001d0) then
                                V2=Uc(P2)+Ul(dum4,dum5)-fc(1,ifc)-fcm(1,ifcm)
                                pnt2=(/dum2, exp_grid(ix,T-it)+1d0, exp_grid(ixm,T-it)+1d0/)
                                INTERP3D=ev(:,:,:,iam,ium,iaf,iuf,T+1-it,ifc,ifcm)
                                vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                V2=V2+beta*OmegaActive(T-it)*vnext
                            else
                                V2=Uc(P2)+Ul(dum4,dum5)-fc(1,ifc)-fcm(1,ifcm)
                                call db3val(dum2,exp_grid(ix,T-it)+1d0,exp_grid(ixm,T-it)+1d0,idx,idy,idz,&
                                        tx,ty(:,T+1-it),tz(:,T+1-it),&
                                        nk,nexp,nexp,kx,ky,kz,&
                                        ev_bspl(iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                        inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.true.)
                                V2=V2+beta*OmegaActive(T-it)*vnext

                            end if

                            pnt2 = (/P3, wagem, wagef/)
                            dum4 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborm, nc, nw, nw, pnt2),1d-10),1d0)
                            dum5 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborf, nc, nw, nw, pnt2),1d-10),1d0)
                            y=dum4*wagem+dum5*wagef

                            dum2 = (nonlab_inc + (1d0-t_const)*(after_tax_labor_inc_married(y)-y*t_employee) - P3*(1d0+tc))/(1d0+mu)

                            if(dum2<0.0001d0) then
                                V3=-999999999d0
                            elseif(dum2>k_grid(nk)-0.001d0) then
                                V3=Uc(P3)+Ul(dum4,dum5)-fc(1,ifc)-fcm(1,ifcm)
                                pnt2=(/dum2, exp_grid(ix,T-it)+1d0, exp_grid(ixm,T-it)+1d0/)
                                INTERP3D=ev(:,:,:,iam,ium,iaf,iuf,T+1-it,ifc,ifcm)
                                vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                V3=V3+beta*OmegaActive(T-it)*vnext
                            else
                                V3=Uc(P3)+Ul(dum4,dum5)-fc(1,ifc)-fcm(1,ifcm)
                                call db3val(dum2,exp_grid(ix,T-it)+1d0,exp_grid(ixm,T-it)+1d0,idx,idy,idz,&
                                        tx,ty(:,T+1-it),tz(:,T+1-it),&
                                        nk,nexp,nexp,kx,ky,kz,&
                                        ev_bspl(iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                        inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.true.)
                                V3=V3+beta*OmegaActive(T-it)*vnext
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
                        
                        pnt2 = [dum2, exp_f_prime, exp_m_prime]
                        associate (pl_ev_aux => pol_ev_aux(iam, ium, iaf, iuf, ifc, ifcm, MEN))
                            vnext_m = pl_ev_aux%eval(pnt2)
                        end associate
                        associate (pl_ev_aux => pol_ev_aux(iam, ium, iaf, iuf, ifc, ifcm, WOMEN))
                            vnext_f = pl_ev_aux%eval(pnt2)
                        end associate
                        Ucur = Uc(ce)+Ul(nem,nef)-fc(1,ifc)-fcm(1,ifcm)
                        ve_m_aux = Ucur + beta*OmegaActive(T-it)*vnext_m
                        ve_f_aux = Ucur + beta*OmegaActive(T-it)*vnext_f
                        
                        lfpm_e = 1d0
                        lfpf_e = 1d0
                        c_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M1,LFP_F1)=ce
                        k_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M1,LFP_F1)=ke
                        nm_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M1,LFP_F1)=dum4
                        nf_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M1,LFP_F1)=dum5
                        v_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M1,LFP_F1,MEN)=ve_m_aux !V2
                        v_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M1,LFP_F1,WOMEN)=ve_f_aux                        
                        

                        !If female unemployed
                        P1=0.001d0
                        nonlab_inc = (k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+Unemp_benefit
                        aftertax_lab_inc = after_tax_labor_inc_married(wagem)
                        P4 = min((k_grid(nk)-0.001d0)/(1d0+tc), (nonlab_inc + (1d0-t_const)*(aftertax_lab_inc - (wagem)*t_employee))/(1d0+tc))                  

                        exp_f_prime = exp_grid(ix,T-it)*(1d0-deltaexp)
                        exp_m_prime = exp_grid(ixm,T-it)+1d0                        
                        
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
                                pnt2=(/dum2, exp_grid(ix,T-it)*(1d0-deltaexp), exp_grid(ixm,T-it)+1d0/)
                                INTERP3D=ev(:,:,:,iam,ium,iaf,iuf,T+1-it,ifc,ifcm)
                                vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                V2=V2+beta*OmegaActive(T-it)*vnext
                            else
                                V2=Uc(P2)+Ul(dum4,0d0)-fcm(1,ifcm)
                                call db3val(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),exp_grid(ixm,T-it)+1d0,idx,idy,idz,&
                                        tx,ty(:,T+1-it),tz(:,T+1-it),&
                                        nk,nexp,nexp,kx,ky,kz,&
                                        ev_bspl(iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                        inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.true.)
                                V2=V2+beta*OmegaActive(T-it)*vnext
                            end if

                            pnt1 = (/P3, wagem/)
                            dum4 = min(max(bilin_interp(c_grid, wage_grid, labormwork, nc, nw, pnt1),0d0),1d0)
                            y=dum4*wagem

                            dum2 = (nonlab_inc + (1d0-t_const)*(after_tax_labor_inc_married(y)-y*t_employee)-P3*(1d0+tc))/(1d0+mu)

                            if(dum2<0.0001d0) then
                                V3=-999999999d0
                            elseif(dum2>k_grid(nk)-0.001d0) then
                                V3=Uc(P3)+Ul(dum4,0d0)-fcm(1,ifcm)
                                dum2=max(dum2,0d0)
                                pnt2=(/dum2, exp_grid(ix,T-it)*(1d0-deltaexp), exp_grid(ixm,T-it)+1d0/)
                                INTERP3D=ev(:,:,:,iam,ium,iaf,iuf,T+1-it,ifc,ifcm)
                                vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                V3=V3+beta*OmegaActive(T-it)*vnext
                            else
                                V3=Uc(P3)+Ul(dum4,0d0)-fcm(1,ifcm)
                                call db3val(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),exp_grid(ixm,T-it)+1d0,idx,idy,idz,&
                                        tx,ty(:,T+1-it),tz(:,T+1-it),&
                                        nk,nexp,nexp,kx,ky,kz,&
                                        ev_bspl(iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                        inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.true.)
                                V3=V3+beta*OmegaActive(T-it)*vnext
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
                        
                        pnt2 = [dum2, exp_f_prime, exp_m_prime]
                        associate (pl_ev_aux => pol_ev_aux(iam, ium, iaf, iuf, ifc, ifcm, MEN))
                            vnext_m = pl_ev_aux%eval(pnt2)
                        end associate         
                        associate (pl_ev_aux => pol_ev_aux(iam, ium, iaf, iuf, ifc, ifcm, WOMEN))
                            vnext_f = pl_ev_aux%eval(pnt2)
                        end associate
                        Ucur = Uc(P2)+Ul(dum4,0d0)-fcm(1,ifcm)
                        vu_m_aux = Ucur + beta*OmegaActive(T-it)*vnext_m
                        vu_f_aux = Ucur + beta*OmegaActive(T-it)*vnext_f                        
                        
                        lfpm_u=1d0
                        lfpf_u=0d0
                        c_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M1,LFP_F0)=cu
                        k_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M1,LFP_F0)=ku   
                        nm_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M1,LFP_F0)=dum4
                        nf_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M1,LFP_F0)=0d0   
                        v_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M1,LFP_F0,MEN)=vu_m_aux !V2
                        v_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M1,LFP_F0,WOMEN)=vu_f_aux                        
                        
                        

                        !If male unemployed
                        P1=0.001d0
                        nonlab_inc = (k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+Unemp_benefit
                        aftertax_lab_inc = after_tax_labor_inc_married(wagef)
                        P4 = min((k_grid(nk)-0.001d0)/(1d0+tc), (nonlab_inc + (1d0-t_const)*(aftertax_lab_inc - (wagef)*t_employee))/(1d0+tc))                  

                        exp_f_prime = exp_grid(ix,T-it)+1d0
                        exp_m_prime = exp_grid(ixm,T-it)*(1d0-deltaexp)                                                  
                                 
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
                                pnt2=(/dum2, exp_grid(ix,T-it)+1d0, exp_grid(ixm,T-it)*(1d0-deltaexp)/)
                                INTERP3D=ev(:,:,:,iam,ium,iaf,iuf,T+1-it,ifc,ifcm)
                                vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                V2=V2+beta*OmegaActive(T-it)*vnext
                            else
                                V2=Uc(P2)+Ul(0d0,dum4)-fc(1,ifc)
                                call db3val(dum2,exp_grid(ix,T-it)+1d0,exp_grid(ixm,T-it)*(1d0-deltaexp),idx,idy,idz,&
                                        tx,ty(:,T+1-it),tz(:,T+1-it),&
                                        nk,nexp,nexp,kx,ky,kz,&
                                        ev_bspl(iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                        inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.true.)
                                V2=V2+beta*OmegaActive(T-it)*vnext
                            end if

                            pnt1 = (/P3, wagef/)
                            dum4 = min(max(bilin_interp(c_grid, wage_grid, laborfwork, nc, nw, pnt1),0d0),1d0)
                            y=dum4*wagef

                            dum2 = (nonlab_inc + (1d0-t_const)*(after_tax_labor_inc_married(y)-y*t_employee)-P3*(1d0+tc))/(1d0+mu)

                            if(dum2<0.0001d0) then
                                V3=-999999999d0
                            elseif(dum2>k_grid(nk)-0.001d0) then
                                V3=Uc(P3)+Ul(0d0,dum4)-fc(1,ifc)
                                dum2=max(dum2,0d0)
                                pnt2=(/dum2, exp_grid(ix,T-it)+1d0, exp_grid(ixm,T-it)*(1d0-deltaexp)/)
                                INTERP3D=ev(:,:,:,iam,ium,iaf,iuf,T+1-it,ifc,ifcm)
                                vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                V3=V3+beta*OmegaActive(T-it)*vnext
                            else
                                V3=Uc(P3)+Ul(0d0,dum4)-fc(1,ifc)
                                call db3val(dum2,exp_grid(ix,T-it)+1d0,exp_grid(ixm,T-it)*(1d0-deltaexp),idx,idy,idz,&
                                        tx,ty(:,T+1-it),tz(:,T+1-it),&
                                        nk,nexp,nexp,kx,ky,kz,&
                                        ev_bspl(iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                        inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.true.)
                                V3=V3+beta*OmegaActive(T-it)*vnext
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
                        
                        pnt2 = [dum2, exp_f_prime, exp_m_prime]
                        associate (pl_ev_aux => pol_ev_aux(iam, ium, iaf, iuf, ifc, ifcm, MEN))
                            vnext_m = pl_ev_aux%eval(pnt2)
                        end associate         
                        associate (pl_ev_aux => pol_ev_aux(iam, ium, iaf, iuf, ifc, ifcm, WOMEN))
                            vnext_f = pl_ev_aux%eval(pnt2)
                        end associate
                        Ucur = Uc(P2)+Ul(0d0,dum4)-fc(1,ifc)                            
                            
                        vtmp_m_aux = Ucur + beta*OmegaActive(T-it)*vnext_m
                        vtmp_f_aux = Ucur + beta*OmegaActive(T-it)*vnext_f                          
                        
                        

                        if (V2 >= vu) then
                            Vu=V2
                            ku=dum2
                            num=0d0
                            nuf=dum4
                            cu=P2
                            
                            lfpm_u=0d0
                            lfpf_u=1d0
                            
                            vu_m_aux = vtmp_m_aux
                            vu_f_aux = vtmp_f_aux                              
                            
                        end if
                        
                        c_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M0,LFP_F1)=P2
                        k_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M0,LFP_F1)=dum2        
                        nm_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M0,LFP_F1)=0d0
                        nf_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M0,LFP_F1)=dum4    
                        v_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M0,LFP_F1,MEN)=vtmp_m_aux !V2
                        v_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M0,LFP_F1,WOMEN)=vtmp_f_aux                        
                        
                        

                        !If both spouses unemployed
                        nonlab_inc = (k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk)) + lumpsum + 2d0*Unemp_benefit + after_tax_labor_inc_married(0d0)
                        P1=0.001d0
                        P4 = min( k_grid(nk)-0.001d0, nonlab_inc )/ (1d0+tc)
                        
                        exp_f_prime = exp_grid(ix,T-it)*(1d0-deltaexp)
                        exp_m_prime = exp_grid(ixm,T-it)*(1d0-deltaexp)                          
                        
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                            dum2 = (nonlab_inc - P2*(1d0+tc))/(1d0+mu)
                            if(dum2<0.0001d0) then
                                V2=-999999999d0
                            elseif(dum2>k_grid(nk)-0.001d0) then
                                V2=Uc(P2)+Ul(0d0,0d0)
                                pnt2=(/dum2, exp_grid(ix,T-it)*(1d0-deltaexp), exp_grid(ixm,T-it)*(1d0-deltaexp)/)
                                INTERP3D=ev(:,:,:,iam,ium,iaf,iuf,T+1-it,ifc,ifcm)
                                vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                V2=V2+beta*OmegaActive(T-it)*vnext
                            else
                                V2=Uc(P2)+Ul(0d0,0d0)
                                call db3val(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),exp_grid(ixm,T-it)*(1d0-deltaexp),idx,idy,idz,&
                                        tx,ty(:,T+1-it),tz(:,T+1-it),&
                                        nk,nexp,nexp,kx,ky,kz,&
                                        ev_bspl(iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                        inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.true.)
                                V2=V2+beta*OmegaActive(T-it)*vnext
                            end if

                            dum2 = (nonlab_inc - P3*(1d0+tc))/(1d0+mu)
                            if(dum2<0.0001d0) then
                                V3=-999999999d0
                            elseif(dum2>k_grid(nk)-0.001d0) then
                                V3=Uc(P3)+Ul(0d0,0d0)
                                dum2=max(dum2,0d0)
                                pnt2=(/dum2, exp_grid(ix,T-it)*(1d0-deltaexp), exp_grid(ixm,T-it)*(1d0-deltaexp)/)
                                INTERP3D=ev(:,:,:,iam,ium,iaf,iuf,T+1-it,ifc,ifcm)
                                vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                V3=V3+beta*OmegaActive(T-it)*vnext
                            else
                                V3=Uc(P3)+Ul(0d0,0d0)
                                call db3val(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),exp_grid(ixm,T-it)*(1d0-deltaexp),idx,idy,idz,&
                                        tx,ty(:,T+1-it),tz(:,T+1-it),&
                                        nk,nexp,nexp,kx,ky,kz,&
                                        ev_bspl(iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                        inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.true.)
                                V3=V3+beta*OmegaActive(T-it)*vnext
                            end if

                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-6) exit
                        end do

                        
                        pnt2 = [dum2, exp_f_prime, exp_m_prime]
                        associate (pl_ev_aux => pol_ev_aux(iam, ium, iaf, iuf, ifc, ifcm, MEN))
                            vnext_m = pl_ev_aux%eval(pnt2)
                        end associate         
                        associate (pl_ev_aux => pol_ev_aux(iam, ium, iaf, iuf, ifc, ifcm, WOMEN))
                            vnext_f = pl_ev_aux%eval(pnt2)
                        end associate
                        Ucur = Uc(P2)+Ul(0d0,0d0)
                        vtmp_m_aux = Ucur + beta*OmegaActive(T-it)*vnext_m
                        vtmp_f_aux = Ucur + beta*OmegaActive(T-it)*vnext_f                          
                        
                        
                        if (V2 >= vu) then
                            Vu=V2
                            ku=dum2
                            num=0d0
                            nuf=0d0
                            cu=P2
                            
                            lfpm_u=0d0
                            lfpf_u=0d0
                            
                            vu_m_aux = vtmp_m_aux
                            vu_f_aux = vtmp_f_aux                              
                        end if
                        
                        c_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M0,LFP_F0)=P2
                        k_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M0,LFP_F0)=dum2    
                        nm_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M0,LFP_F0)=0d0
                        nf_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M0,LFP_F0)=0d0   
                        v_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M0,LFP_F0,MEN)=vtmp_m_aux !V2
                        v_lfp(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm,LFP_M0,LFP_F0,WOMEN)=vtmp_f_aux                        
                        


                        if (ve >= vu) then
                            v(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)=ve
                            c(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)=ce
                            k(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)=ke
                            nm(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)=nem
                            nf(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)=nef
                            
                            v_aux(ik,ix,ixm,iam,ium,iaf,iuf,ifc,ifcm,MEN)=ve_m_aux
                            v_aux(ik,ix,ixm,iam,ium,iaf,iuf,ifc,ifcm,WOMEN)=ve_f_aux
                            lfpm(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)=lfpm_e
                            lfpf(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)=lfpf_e                                                        
                        else
                            v(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)=vu
                            c(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)=cu
                            k(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)=ku
                            nm(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)=num
                            nf(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)=nuf
                            
                            v_aux(ik,ix,ixm,iam,ium,iaf,iuf,ifc,ifcm,MEN)=vu_m_aux
                            v_aux(ik,ix,ixm,iam,ium,iaf,iuf,ifc,ifcm,WOMEN)=vu_f_aux                             
                            lfpm(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)=lfpm_u
                            lfpf(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)=lfpf_u                                                         
                        end if

                    end do
                end do
            end do
        end do
    end do


    !Singles

    j=2

    do ifc=1,nfc
        !Print *,'ix is',ix
        !Finding optimal capital by golden search
        wagef = (1d0/(1d0+exp(kappa*(T-(it+agestart)))))*wage(2,a(2,iam),exp_grid(ix,T-it),u(2,ium))/(1d0+t_employer)
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
            elseif(dum2>k_grid(nk)-0.001d0) then
                V2=Uc(P2)-chifs*(dum4**(1d0+etaf))/(1d0+etaf)-fc(2,ifc)
                pnt1=(/dum2, exp_grid(ix,T-it)+1d0/)
                INTERP2D=evs(j,:,:,iam,ium,T+1-it,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext
                INTERP2D=evm(j,:,:,iam,ium,T+1-it,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V2=V2+beta*OmegaActive(T-it)*Probm(T-it)*vnext
            else
                V2=Uc(P2)-chifs*(dum4**(1d0+etaf))/(1d0+etaf)-fc(2,ifc)
                call db2val(dum2,exp_grid(ix,T-it)+1d0,idx,idy,&
                    tx,ty(:,T-it+1),nk,nexp,kx,ky,&
                    evs_bspl(j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)
                V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext
                call db2val(dum2,exp_grid(ix,T-it)+1d0,idx,idy,&
                    tx,ty(:,T-it+1),nk,nexp,kx,ky,&
                    evm_bspl(j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)
                V2=V2+beta*OmegaActive(T-it)*Probm(T-it)*vnext
            end if

            pnt1 = (/P3, wagef/)
            dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglef, nc, nw, pnt1),0d0),1d0)
            y=dum4*wagef

            dum2 = (nonlab_inc + (1d0-t_const)*(after_tax_labor_inc_single(y)-y*t_employee)-P3*(1d0+tc))/(1d0+mu)

            if(dum2<0.0001d0) then
                V3=-999999999d0
            elseif(dum2>k_grid(nk)-0.001d0) then
                V3=Uc(P3)-chifs*(dum4**(1d0+etaf))/(1d0+etaf)-fc(2,ifc)
                pnt1=(/dum2, exp_grid(ix,T-it)+1d0/)
                INTERP2D=evs(j,:,:,iam,ium,T+1-it,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V3=V3+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext
                INTERP2D=evm(j,:,:,iam,ium,T+1-it,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V3=V3+beta*OmegaActive(T-it)*Probm(T-it)*vnext
            else
                V3=Uc(P3)-chifs*(dum4**(1d0+etaf))/(1d0+etaf)-fc(2,ifc)
                call db2val(dum2,exp_grid(ix,T-it)+1d0,idx,idy,&
                    tx,ty(:,T-it+1),nk,nexp,kx,ky,&
                    evs_bspl(j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)
                V3=V3+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext
                call db2val(dum2,exp_grid(ix,T-it)+1d0,idx,idy,&
                    tx,ty(:,T-it+1),nk,nexp,kx,ky,&
                    evm_bspl(j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)
                V3=V3+beta*OmegaActive(T-it)*Probm(T-it)*vnext
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
                        
        cs_lfp(j,ik,ix,iam,ium,T-it,ifc,LFP_1)=P2
        ks_lfp(j,ik,ix,iam,ium,T-it,ifc,LFP_1)=dum2
        vs_lfp(j,ik,ix,iam,ium,T-it,ifc,LFP_1)=V2
        ns_lfp(j,ik,ix,iam,ium,T-it,ifc,LFP_1)=dum4


        !If female unemployed
        nonlab_inc = (k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk)) + lumpsum*0.5d0 + Unemp_benefit + after_tax_labor_inc_single(0d0)
        P1=0.001d0
        P4 = min( k_grid(nk)-0.001d0, nonlab_inc )/(1d0+tc)
        do
            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

            dum2 = (nonlab_inc - P2*(1d0+tc))/(1d0+mu)

            if(dum2<0.0001d0) then
                V2=-999999999d0
            elseif(dum2>k_grid(nk)-0.001d0) then
                V2=Uc(P2)
                pnt1=(/dum2, exp_grid(ix,T-it)*(1d0-deltaexp)/)
                INTERP2D=evs(j,:,:,iam,ium,T+1-it,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext
                INTERP2D=evm(j,:,:,iam,ium,T+1-it,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V2=V2+beta*OmegaActive(T-it)*Probm(T-it)*vnext
            else
                V2=Uc(P2)
                call db2val(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),idx,idy,&
                    tx,ty(:,T-it+1),nk,nexp,kx,ky,&
                    evs_bspl(j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)
                V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext
                call db2val(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),idx,idy,&
                    tx,ty(:,T-it+1),nk,nexp,kx,ky,&
                    evm_bspl(j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)
                V2=V2+beta*OmegaActive(T-it)*Probm(T-it)*vnext
            end if

            dum2 = (nonlab_inc - P3*(1d0+tc))/(1d0+mu)
            if(dum2<0.0001d0) then
                V3=-999999999d0
            elseif(dum2>k_grid(nk)-0.001d0) then
                V3=Uc(P3)
                pnt1=(/dum2, exp_grid(ix,T-it)*(1d0-deltaexp)/)
                INTERP2D=evs(j,:,:,iam,ium,T+1-it,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V3=V3+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext
                INTERP2D=evm(j,:,:,iam,ium,T+1-it,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V3=V3+beta*OmegaActive(T-it)*Probm(T-it)*vnext
            else
                V3=Uc(P3)
                call db2val(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),idx,idy,&
                    tx,ty(:,T-it+1),nk,nexp,kx,ky,&
                    evs_bspl(j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)
                V3=V3+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext
                call db2val(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),idx,idy,&
                    tx,ty(:,T-it+1),nk,nexp,kx,ky,&
                    evm_bspl(j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)
                V3=V3+beta*OmegaActive(T-it)*Probm(T-it)*vnext                      
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

        cs_lfp(j,ik,ix,iam,ium,T-it,ifc,LFP_0)=P2
        ks_lfp(j,ik,ix,iam,ium,T-it,ifc,LFP_0)=dum2
        vs_lfp(j,ik,ix,iam,ium,T-it,ifc,LFP_0)=V2      
        ns_lfp(j,ik,ix,iam,ium,T-it,ifc,LFP_0)=0d0
        

        if (ves >= vus) then
            vs(j,ik,ix,iam,ium,T-it,ifc)=ves
            cs(j,ik,ix,iam,ium,T-it,ifc)=ces
            ks(j,ik,ix,iam,ium,T-it,ifc)=kes
            ns(j,ik,ix,iam,ium,T-it,ifc)=nes
            lfps(j,ik,ix,iam,ium,T-it,ifc)=1d0            
        else
            vs(j,ik,ix,iam,ium,T-it,ifc)=vus
            cs(j,ik,ix,iam,ium,T-it,ifc)=cus
            ks(j,ik,ix,iam,ium,T-it,ifc)=kus
            ns(j,ik,ix,iam,ium,T-it,ifc)=nus
            lfps(j,ik,ix,iam,ium,T-it,ifc)=0d0
        end if
    end do

    !Men   
    j=1
    do ifc=1,nfcm
        wagem = (1d0/(1d0+exp(kappa*(T-(it+agestart)))))*wage(1,a(1,iam),exp_grid(ix,T-it),u(1,ium))/(1d0+t_employer)
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
            elseif(dum2>k_grid(nk)-0.001d0) then
                V2=Uc(P2)-chims*(dum4**(1d0+etam))/(1d0+etam)-fcm(2,ifc)
                pnt1=(/dum2, exp_grid(ix,T-it)+1d0/)
                INTERP2D=evs(j,:,:,iam,ium,T+1-it,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext
                INTERP2D=evm(j,:,:,iam,ium,T+1-it,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V2=V2+beta*OmegaActive(T-it)*Probm(T-it)*vnext
            else
                V2=Uc(P2)-chims*(dum4**(1d0+etam))/(1d0+etam)-fcm(2,ifc)
                call db2val(dum2,exp_grid(ix,T-it)+1d0,idx,idy,&
                    tx,ty(:,T-it+1),nk,nexp,kx,ky,&
                    evs_bspl(j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)
                V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext
                call db2val(dum2,exp_grid(ix,T-it)+1d0,idx,idy,&
                    tx,ty(:,T-it+1),nk,nexp,kx,ky,&
                    evm_bspl(j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)
                V2=V2+beta*OmegaActive(T-it)*Probm(T-it)*vnext
            end if

            pnt1 = (/P3, wagem/)
            dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglem, nc, nw, pnt1),0d0),1d0)
            y=dum4*wagem

            dum2 = (nonlab_inc + (1d0-t_const)*(after_tax_labor_inc_single(y)-y*t_employee)-P3*(1d0+tc))/(1d0+mu)

            if(dum2<0.0001d0) then
                V3=-999999999d0
            elseif(dum2>k_grid(nk)-0.001d0) then
                V3=Uc(P3)-chims*(dum4**(1d0+etam))/(1d0+etam)-fcm(2,ifc)
                pnt1=(/dum2, exp_grid(ix,T-it)+1d0/)
                INTERP2D=evs(j,:,:,iam,ium,T+1-it,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V3=V3+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext
                INTERP2D=evm(j,:,:,iam,ium,T+1-it,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V3=V3+beta*OmegaActive(T-it)*Probm(T-it)*vnext
            else
                V3=Uc(P3)-chims*(dum4**(1d0+etam))/(1d0+etam)-fcm(2,ifc)
                call db2val(dum2,exp_grid(ix,T-it)+1d0,idx,idy,&
                    tx,ty(:,T-it+1),nk,nexp,kx,ky,&
                    evs_bspl(j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)
                V3=V3+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext
                call db2val(dum2,exp_grid(ix,T-it)+1d0,idx,idy,&
                    tx,ty(:,T-it+1),nk,nexp,kx,ky,&
                    evm_bspl(j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)
                V3=V3+beta*OmegaActive(T-it)*Probm(T-it)*vnext
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
        
        cs_lfp(j,ik,ix,iam,ium,T-it,ifc,LFP_1)=P2
        ks_lfp(j,ik,ix,iam,ium,T-it,ifc,LFP_1)=dum2
        vs_lfp(j,ik,ix,iam,ium,T-it,ifc,LFP_1)=V2
        ns_lfp(j,ik,ix,iam,ium,T-it,ifc,LFP_1)=dum4

        !If male unemployed
        nonlab_inc = (k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk)) + lumpsum*0.5d0 + Unemp_benefit + after_tax_labor_inc_single(0d0)
        P1=0.001d0
        P4 = min( k_grid(nk)-0.001d0, nonlab_inc )/(1d0+tc)
        do
            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
            dum2 = (nonlab_inc - P2*(1d0+tc))/(1d0+mu)
            if(dum2<0.0001d0) then
                V2=-999999999d0
            elseif(dum2>k_grid(nk)-0.001d0) then
                V2=Uc(P2)
                pnt1=(/dum2, exp_grid(ix,T-it)*(1d0-deltaexp)/)
                INTERP2D=evs(j,:,:,iam,ium,T+1-it,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext
                INTERP2D=evm(j,:,:,iam,ium,T+1-it,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V2=V2+beta*OmegaActive(T-it)*Probm(T-it)*vnext
            else
                V2=Uc(P2)
                call db2val(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),idx,idy,&
                    tx,ty(:,T-it+1),nk,nexp,kx,ky,&
                    evs_bspl(j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)
                V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext
                call db2val(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),idx,idy,&
                    tx,ty(:,T-it+1),nk,nexp,kx,ky,&
                    evm_bspl(j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)
                V2=V2+beta*OmegaActive(T-it)*Probm(T-it)*vnext
            end if

            dum2 = (nonlab_inc - P3*(1d0+tc))/(1d0+mu)
            if(dum2<0.0001d0) then
                V3=-999999999d0
            elseif(dum2>k_grid(nk)-0.001d0) then
                V3=Uc(P3)
                pnt1=(/dum2, exp_grid(ix,T-it)*(1d0-deltaexp)/)
                INTERP2D=evs(j,:,:,iam,ium,T+1-it,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V3=V3+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext
                INTERP2D=evm(j,:,:,iam,ium,T+1-it,ifc)
                vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                V3=V3+beta*OmegaActive(T-it)*Probm(T-it)*vnext
            else
                V3=Uc(P3)
                call db2val(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),idx,idy,&
                    tx,ty(:,T-it+1),nk,nexp,kx,ky,&
                    evs_bspl(j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)
                V3=V3+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext
                call db2val(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),idx,idy,&
                    tx,ty(:,T-it+1),nk,nexp,kx,ky,&
                    evm_bspl(j,iam,ium,ifc)%coefs,vnext,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)
                V3=V3+beta*OmegaActive(T-it)*Probm(T-it)*vnext                      
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
        
        cs_lfp(j,ik,ix,iam,ium,T-it,ifc,LFP_0)=P2
        ks_lfp(j,ik,ix,iam,ium,T-it,ifc,LFP_0)=dum2   
        vs_lfp(j,ik,ix,iam,ium,T-it,ifc,LFP_0)=V2    
        ns_lfp(j,ik,ix,iam,ium,T-it,ifc,LFP_0)=0d0


        if (ves >= vus) then
            vs(j,ik,ix,iam,ium,T-it,ifc)=ves
            cs(j,ik,ix,iam,ium,T-it,ifc)=ces
            ks(j,ik,ix,iam,ium,T-it,ifc)=kes
            ns(j,ik,ix,iam,ium,T-it,ifc)=nes
            lfps(j,ik,ix,iam,ium,T-it,ifc)=1d0
        else
            vs(j,ik,ix,iam,ium,T-it,ifc)=vus
            cs(j,ik,ix,iam,ium,T-it,ifc)=cus
            ks(j,ik,ix,iam,ium,T-it,ifc)=kus
            ns(j,ik,ix,iam,ium,T-it,ifc)=nus
            lfps(j,ik,ix,iam,ium,T-it,ifc)=0d0
        end if
    end do

end subroutine SolveActiveLife
