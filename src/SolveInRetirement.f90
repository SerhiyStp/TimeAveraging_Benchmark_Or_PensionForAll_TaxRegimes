subroutine SolveInRetirement(counter)
    ! Both men and women are retired (absorbing state)
    use Model_Parameters
    use PolicyFunctions
    use Utilities
    use bspline_sub_module
    use GlobParams

    implicit none

    integer, INTENT(IN) :: counter
    real(8) :: c2, v2,dum2,dum3, Psi_pension, Psi_pensionf
    real(8) :: d1, d2,P1,P2,P3,P4,v3
    integer :: ix,ixm,ixd,iam,ium,iaf,iuf,iu2,iu3,tprint,ikd,j,ifc,ifcm,ik
    integer :: it2
    real(8) :: test1, test2,r_ret, r_ret_next
    real(8) :: vnext_p2, vnext_p3
    real(8) :: vnext, exp_grid_dum(nexp), INTERP2D(nk,nexp),pnt1(2),INTERP3D(nk,nexp,nexp),pnt2(3)
    integer :: iflag 
    real(8) :: vnext_test
    integer :: idx, idy, idz, iloy, iloz
    integer :: inbvx, inbvy, inbvz
    real(8) :: ww2(ky,kz),ww1(kz),ww0(3*max(kx,ky,kz))
    real(8) :: w1_d2(ky) 
    real(8) :: w0_d2(3*max(kx,ky)) 

    idx=0
    idy=0
    idz=0
    inbvx=1
    inbvy=1
    inbvz=1
    iloy=1
    iloz=1

    !Assigning the grid points
    dum3=((counter*1d0)/(nexp*na*1d0))-0.00001d0
    ik=int(dum3)+1
    dum3=(((counter-(ik-1)*nexp*na)*1d0)/(na*1d0))-0.00001d0
    ix=int(dum3)+1
    iam=counter-(ik-1)*nexp*na-(ix-1)*na

    r_ret=((1d0+r)/OmegaRet2(Tret+1-it))-1d0
    if(it>1) then
        r_ret_next=((1d0+r)/OmegaRet2(Tret+2-it))-1d0
    end if

    if(it>1) then
        exp_grid_dum=exp_grid(:,T+Tret+2-it)
    end if

    !Married    
    do iaf = 1, na
        do ixm = 1, nexp
            v_ret(1,1,ik,ix,ixm,iam,:,iaf,:,Tret-it+1,:,:)=-999999999d0
            v_ret2(1,1,ik,ix,ixm,iam,:,iaf,:,Tret-it+1,:,:)=-999999999d0                 
        end do
    end do

    !Single
    do j=1,2
        vs_ret(1,j,ik,ix,iam,:,Tret-it+1,:)=-999999999d0
        vs_ret2(1,j,ik,ix,iam,:,Tret-it+1,:)=-999999999d0
    end do

end subroutine SolveInRetirement

    
subroutine SolveInRetirement2(counter)
    ! Woman or Man can choose to work
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    use bspline_sub_module
    use GlobParams

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixm,ixd,iam,ium,iaf,iuf,iu2,iu3,tprint,ikd,j,ifc,ik
    real(8) :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu,Psi_pension,Psi_pensionf
    real(8) :: ces,cus,kes,kus,nes,nus,Inces,ves,vus,r_ret,r_ret_next
    real(8) :: dum3,dum4,dum5,dum6,y
    real(8) :: P1,P2,P3,P4,V2,V3,dum2,pnt2(3),pnt1(2)
    real(8) :: vnext, exp_grid_dum(nexp), INTERP2D(nk,nexp), INTERP3D(nk,nexp,nexp)

    integer :: iflag 
    real(8) :: vnext_test
    integer :: idx, idy, idz, iloy, iloz
    integer :: inbvx, inbvy, inbvz
    real(8) :: ww2(ky,kz),ww1(kz),ww0(3*max(kx,ky,kz))
    real(8) :: w1_d2(ky) 
    real(8) :: w0_d2(3*max(kx,ky)) 
    real(8) :: nonlab_inc, aftertax_lab_inc, dum2_test, P4_test
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

    if(it>1) then
        exp_grid_dum=exp_grid(:,T+Tret+2-it)
    end if


    r_ret=((1d0+r)/OmegaRet2(Tret+1-it))-1d0
    if(it>1) then
        r_ret_next=((1d0+r)/OmegaRet2(Tret+2-it))-1d0
    end if

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
    Inces=0d0
    ves=0d0
    vus=0d0

    !Married woman is working

    do ifc=1,nfc
        do iaf = 1, na
            do ixm = 1, nexp

                v_ret2(2,1,ik,ix,ixm,iaf,:,iam,ium,Tret-it+1,ifc,:)=-999999999d0 
                v_ret(2,1,ik,ix,ixm,iaf,:,iam,ium,Tret-it+1,ifc,:)=-999999999d0

            end do
        end do
    end do

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

    !Married man is working

    do ifc=1,nfcm
        do iaf = 1, na
            do ixm = 1, nexp

                v_ret2(1,2,ik,ixm,ix,iam,ium,iaf,:,Tret-it+1,:,ifc)=-999999999d0
                v_ret(1,2,ik,ixm,ix,iam,ium,iaf,:,Tret-it+1,:,ifc)=-999999999d0

            end do
        end do
    end do

    !Single women
    j=2
    do ifc=1,nfc
        wagef = (1d0/(1d0+exp(kappa*(T+Tret+1-it-agestart))))*wage(2,a(2,iam),exp_grid(ix,T+Tret+1-it),u(2,ium))/(1d0+t_employer)
        Psi_pensionf=(psi0+psi1*av_earnings(2,2,iam)*min(1d0,exp_grid(ix,T+Tret+1-it)/35d0))*(1d0-t_const)

        if(it==1) then
            ! Solve the very last period problem:           
            kes = 0d0
            ces = ((k_grid(ik)+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk)) + Psi_pensionf + lumpsum*0.5d0)/(1d0+tc)
            nes = 0d0
            ves = Uc(ces)
            
            lfps_ret(1,j,ik,ix,iam,:,Tret,:) = 0d0
            cs_ret_lfp(2,j,ik,ix,iam,:,Tret,:,LFP_0) = ces
            ks_ret_lfp(2,j,ik,ix,iam,:,Tret,:,LFP_0) = kes
            ns_ret_lfp(2,j,ik,ix,iam,:,Tret,:,LFP_0) = 0d0               
        else
            ! Solve the Tret-1 to 1st period of retirement:
            !Finding optimal capital by golden search
            P1=0.001d0
            nonlab_inc = (k_grid(ik)+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk))+lumpsum*0.5d0 + Psi_pensionf
            aftertax_lab_inc = after_tax_labor_inc_single(wagef)
            P4 = min((k_grid(nk)-0.001d0)/(1d0+tc), (nonlab_inc + (1d0-t_const)*(aftertax_lab_inc - wagef*t_employee))/(1d0+tc))
            do
                P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                pnt1 = (/P2, wagef/)
                dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglef, nc, nw, pnt1),0d0),1d0)
                y=dum4*wagef
                aftertax_lab_inc = after_tax_labor_inc_single(y)
                dum2 = (nonlab_inc + (1d0-t_const)*(aftertax_lab_inc-y*t_employee)-P2*(1d0+tc))/(1d0+mu) 
                if(dum2<0.0001d0) then
                    V2=-999999999d0
                elseif(dum2>k_grid(nk)-0.001d0) then
                    V2=Uc(P2)-chifs*(dum4**(1d0+etaf))/(1d0+etaf)-fc(2,ifc)

                    pnt1=(/dum2, exp_grid(ix,T+Tret+1-it)+1d0/)
                    INTERP2D=evs_ret(2,j,:,:,iam,ium,Tret+2-it,ifc)
                    vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    V2=V2+beta*OmegaRet(Tret-it+1)*vnext

                else
                    V2=Uc(P2)-chifs*(dum4**(1d0+etaf))/(1d0+etaf)-fc(2,ifc)
                    call db2val(dum2,exp_grid(ix,T+Tret+1-it)+1d0,idx,idy,&
                        tx,ty(:,T+Tret+2-it),nk,nexp,kx,ky,&
                        evs_ret_bspl(2,j,iam,ium,ifc)%coefs,vnext,iflag,&
                        inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.false.)
                    !vnext = D_BS2VL(dum2, exp_grid(ix,T+Tret+1-it), KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T+Tret+2-it), nk, nexp,evs_ret_spln_coefs(2,j,:,:,iam,ium,Tret+2-it,ifc))
                    V2=V2+beta*OmegaRet(Tret-it+1)*vnext

                end if

                pnt1 = (/P3, wagef/)
                dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglef, nc, nw, pnt1),0d0),1d0)
                y=dum4*wagef

                aftertax_lab_inc = after_tax_labor_inc_single(y)
                dum2 = (nonlab_inc + (1d0-t_const)*(aftertax_lab_inc-y*t_employee)-P3*(1d0+tc))/(1d0+mu) 

                if(dum2<0.0001d0) then
                    V3=-999999999d0
                elseif(dum2>k_grid(nk)-0.001d0) then
                    V3=Uc(P3)-chifs*(dum4**(1d0+etaf))/(1d0+etaf)-fc(2,ifc)

                    pnt1=(/dum2, exp_grid(ix,T+Tret+1-it)+1d0/)
                    INTERP2D=evs_ret(2,j,:,:,iam,ium,Tret+2-it,ifc)
                    vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    V3=V3+beta*OmegaRet(Tret-it+1)*vnext

                else
                    V3=Uc(P3)-chifs*(dum4**(1d0+etaf))/(1d0+etaf)-fc(2,ifc)
                    call db2val(dum2,exp_grid(ix,T+Tret+1-it)+1d0,idx,idy,&
                        tx,ty(:,T+Tret+2-it),nk,nexp,kx,ky,&
                        evs_ret_bspl(2,j,iam,ium,ifc)%coefs,vnext,iflag,&
                        inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.false.)
                    V3=V3+beta*OmegaRet(Tret-it+1)*vnext

                end if
                if (V2 < V3) then
                    P1=P2
                else
                    P4=P3
                end if
                if((P4-P1)<1d-8) exit
            end do
            
            if(dum4<minhours) then
                V2=-999999999d0
            end if
            
            Ves=V2
            kes=dum2
            nes=dum4
            ces=P2
            !Inces=dum4*wagef
            cs_ret_lfp(2,j,ik,ix,iam,ium,Tret-it+1,ifc,LFP_1) = P2
            ks_ret_lfp(2,j,ik,ix,iam,ium,Tret-it+1,ifc,LFP_1) = dum2
            ns_ret_lfp(2,j,ik,ix,iam,ium,Tret-it+1,ifc,LFP_1) = dum4               

            !Working 0 hours but not retiring

            P1=0.001d0
            P4=min((k_grid(nk)-0.001d0)/(1d0+tc),((k_grid(ik)+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk))+lumpsum*0.5d0+Psi_pensionf+Unemp_benefit)/(1d0+tc))
            do
                P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                dum2=((k_grid(ik)+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk))+lumpsum*0.5d0+Psi_pensionf+Unemp_benefit-P2*(1d0+tc))/(1d0+mu)
                if(dum2<0.0001d0) then
                    V2=-999999999d0
                elseif(dum2>k_grid(nk)-0.001d0) then
                    V2=Uc(P2)

                    pnt1=(/dum2, exp_grid(ix,T+Tret+1-it)/)
                    INTERP2D=evs_ret(2,j,:,:,iam,ium,Tret+2-it,ifc)
                    vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    V2=V2+beta*OmegaRet(Tret-it+1)*vnext

                else
                    V2=Uc(P2)
                    call db2val(dum2,exp_grid(ix,T+Tret+1-it),idx,idy,&
                        tx,ty(:,T+Tret+2-it),nk,nexp,kx,ky,&
                        evs_ret_bspl(2,j,iam,ium,ifc)%coefs,vnext,iflag,&
                        inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.false.)
                    V2=V2+beta*OmegaRet(Tret-it+1)*vnext

                end if

                dum2=((k_grid(ik)+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk))+lumpsum*0.5d0+Psi_pensionf+Unemp_benefit-P3*(1d0+tc))/(1d0+mu)
                if(dum2<0.0001d0) then
                    V3=-999999999d0
                elseif(dum2>k_grid(nk)-0.001d0) then
                    V3=Uc(P3)

                    pnt1=(/dum2, exp_grid(ix,T+Tret+1-it)/)
                    INTERP2D=evs_ret(2,j,:,:,iam,ium,Tret+2-it,ifc)
                    vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    V3=V3+beta*OmegaRet(Tret-it+1)*vnext

                else
                    V3=Uc(P3)
                    call db2val(dum2,exp_grid(ix,T+Tret+1-it),idx,idy,&
                        tx,ty(:,T+Tret+2-it),nk,nexp,kx,ky,&
                        evs_ret_bspl(2,j,iam,ium,ifc)%coefs,vnext,iflag,&
                        inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.false.)
                    V3=V3+beta*OmegaRet(Tret-it+1)*vnext

                end if
                if (V2 < V3) then
                    P1=P2
                else
                    P4=P3
                end if
                if((P4-P1)<1d-8) exit
            end do

            if (V2 >= ves) then
                Ves=V2
                kes=dum2
                nes=0d0
                Inces=0d0
                ces=P2
            end if
            
            cs_ret_lfp(2,j,ik,ix,iam,ium,Tret-it+1,ifc,LFP_0) = P2
            ks_ret_lfp(2,j,ik,ix,iam,ium,Tret-it+1,ifc,LFP_0) = dum2
            ns_ret_lfp(2,j,ik,ix,iam,ium,Tret-it+1,ifc,LFP_0) = 0d0              

        end if

        vs_ret2(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=ves
        cs_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=ces
        Uprimes_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=dUc(ces)
        ks_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=kes
        ns_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=nes
        Incomes_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=Inces

        vus=vs_ret(1,j,ik,ix,iam,ium,Tret-it+1,ifc)
        if (ves >= vus) then
            vs_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=ves
            rets(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=0d0
            lfps_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=1d0
        else
            vs_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=vus
            rets(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=1d0
            lfps_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=0d0
        end if

        if(it==1) then
            rets(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=1d0
            lfps_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=0d0
        end if

    end do

    !Single men
    j=1
    do ifc=1,nfcm

        !Pension depends on expected wage conditional on ability and experience
        Psi_pension=(psi0+psi1*av_earnings(1,2,iam)*min(1d0,exp_grid(ix,T+Tret+1-it)/35d0))*(1d0-t_const)
        wagem = (1d0/(1d0+exp(kappa*(T+Tret+1-it-agestart))))*wage(1,a(1,iam),exp_grid(ix,T+Tret+1-it),u(1,ium))/(1d0+t_employer)

        if(it==1) then
            ! Solve the very last period problem:           
            kes = 0d0
            ces = ((k_grid(ik)+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk)) + Psi_pension+lumpsum*0.5d0)/(1d0+tc)
            nes = 0d0
            ves = Uc(ces)
            
            cs_ret_lfp(2,j,ik,ix,iam,:,Tret,:,LFP_0) = ces
            ks_ret_lfp(2,j,ik,ix,iam,:,Tret,:,LFP_0) = kes
            ns_ret_lfp(2,j,ik,ix,iam,:,Tret,:,LFP_0) = 0d0                
            
        else
            ! Solve the Tret-1 to 1st period of retirement:
            !Finding optimal capital by golden search
            P1=0.001d0
            nonlab_inc = (k_grid(ik)+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk))+lumpsum*0.5d0+Psi_pension
            aftertax_lab_inc = after_tax_labor_inc_single(wagem)
            P4 = min((k_grid(nk)-0.001d0)/(1d0+tc), (nonlab_inc + (1d0-t_const)*(aftertax_lab_inc - wagem*t_employee))/(1d0+tc))
            do
                P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                pnt1 = (/P2, wagem/)
                dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglem, nc, nw, pnt1),0d0),1d0)
                y=dum4*wagem
                aftertax_lab_inc = after_tax_labor_inc_single(y)
                dum2 = (nonlab_inc + (1d0-t_const)*(aftertax_lab_inc-y*t_employee)-P2*(1d0+tc))/(1d0+mu) 
                if(dum2<0.0001d0) then
                    V2=-999999999d0
                elseif(dum2>k_grid(nk)-0.001d0) then
                    V2=Uc(P2)-chims*(dum4**(1d0+etam))/(1d0+etam)-fcm(2,ifc)

                    pnt1=(/dum2, exp_grid(ix,T+Tret+1-it)+1d0/)
                    INTERP2D=evs_ret(2,j,:,:,iam,ium,Tret+2-it,ifc)
                    vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    V2=V2+beta*OmegaRet(Tret-it+1)*vnext

                else
                    V2=Uc(P2)-chims*(dum4**(1d0+etam))/(1d0+etam)-fcm(2,ifc)
                    call db2val(dum2,exp_grid(ix,T+Tret+1-it)+1d0,idx,idy,&
                        tx,ty(:,T+Tret+2-it),nk,nexp,kx,ky,&
                        evs_ret_bspl(2,j,iam,ium,ifc)%coefs,vnext,iflag,&
                        inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.false.)
                    !vnext = D_BS2VL(dum2, exp_grid(ix,T+Tret+1-it), KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T+Tret+2-it), nk, nexp,evs_ret_spln_coefs(2,j,:,:,iam,ium,Tret+2-it,ifc))
                    V2=V2+beta*OmegaRet(Tret-it+1)*vnext

                end if

                pnt1 = (/P3, wagem/)
                dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglem, nc, nw, pnt1),0d0),1d0)
                y=dum4*wagem

                aftertax_lab_inc = after_tax_labor_inc_single(y)
                dum2 = (nonlab_inc + (1d0-t_const)*(aftertax_lab_inc-y*t_employee)-P3*(1d0+tc))/(1d0+mu) 

                if(dum2<0.0001d0) then
                    V3=-999999999d0
                elseif(dum2>k_grid(nk)-0.001d0) then
                    V3=Uc(P3)-chims*(dum4**(1d0+etam))/(1d0+etam)-fcm(2,ifc)

                    pnt1=(/dum2, exp_grid(ix,T+Tret+1-it)+1d0/)
                    INTERP2D=evs_ret(2,j,:,:,iam,ium,Tret+2-it,ifc)
                    vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    V3=V3+beta*OmegaRet(Tret-it+1)*vnext

                else
                    V3=Uc(P3)-chims*(dum4**(1d0+etam))/(1d0+etam)-fcm(2,ifc)
                    call db2val(dum2,exp_grid(ix,T+Tret+1-it)+1d0,idx,idy,&
                        tx,ty(:,T+Tret+2-it),nk,nexp,kx,ky,&
                        evs_ret_bspl(2,j,iam,ium,ifc)%coefs,vnext,iflag,&
                        inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.false.)
                    V3=V3+beta*OmegaRet(Tret-it+1)*vnext

                end if
                if (V2 < V3) then
                    P1=P2
                else
                    P4=P3
                end if
                if((P4-P1)<1d-8) exit
            end do
            
            if(dum4<minhours) then
                V2=-999999999d0
            end if
            
            Ves=V2
            kes=dum2
            nes=dum4
            ces=P2
            
            cs_ret_lfp(2,j,ik,ix,iam,ium,Tret-it+1,ifc,LFP_1) = P2
            ks_ret_lfp(2,j,ik,ix,iam,ium,Tret-it+1,ifc,LFP_1) = dum2
            ns_ret_lfp(2,j,ik,ix,iam,ium,Tret-it+1,ifc,LFP_1) = dum4                 

            !Working 0 hours but not retiring

            P1=0.001d0
            P4=min((k_grid(nk)-0.001d0)/(1d0+tc),((k_grid(ik)+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk))+lumpsum*0.5d0+Psi_pension+Unemp_benefit)/(1d0+tc))
            do
                P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                dum2=((k_grid(ik)+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk))+lumpsum*0.5d0+Psi_pension+Unemp_benefit-P2*(1d0+tc))/(1d0+mu)
                if(dum2<0.0001d0) then
                    V2=-999999999d0
                elseif(dum2>k_grid(nk)-0.001d0) then
                    V2=Uc(P2)

                    pnt1=(/dum2, exp_grid(ix,T+Tret+1-it)/)
                    INTERP2D=evs_ret(2,j,:,:,iam,ium,Tret+2-it,ifc)
                    vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    V2=V2+beta*OmegaRet(Tret-it+1)*vnext

                else
                    V2=Uc(P2)
                    call db2val(dum2,exp_grid(ix,T+Tret+1-it),idx,idy,&
                        tx,ty(:,T+Tret+2-it),nk,nexp,kx,ky,&
                        evs_ret_bspl(2,j,iam,ium,ifc)%coefs,vnext,iflag,&
                        inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.false.)
                    V2=V2+beta*OmegaRet(Tret-it+1)*vnext

                end if

                dum2=((k_grid(ik)+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk))+lumpsum*0.5d0+Psi_pension+Unemp_benefit-P3*(1d0+tc))/(1d0+mu)
                if(dum2<0.0001d0) then
                    V3=-999999999d0
                elseif(dum2>k_grid(nk)-0.001d0) then
                    V3=Uc(P3)

                    pnt1=(/dum2, exp_grid(ix,T+Tret+1-it)/)
                    INTERP2D=evs_ret(2,j,:,:,iam,ium,Tret+2-it,ifc)
                    vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    V3=V3+beta*OmegaRet(Tret-it+1)*vnext

                else
                    V3=Uc(P3)
                    call db2val(dum2,exp_grid(ix,T+Tret+1-it),idx,idy,&
                        tx,ty(:,T+Tret+2-it),nk,nexp,kx,ky,&
                        evs_ret_bspl(2,j,iam,ium,ifc)%coefs,vnext,iflag,&
                        inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.false.)
                    V3=V3+beta*OmegaRet(Tret-it+1)*vnext

                end if
                if (V2 < V3) then
                    P1=P2
                else
                    P4=P3
                end if
                if((P4-P1)<1d-8) exit
            end do

            if (V2 >= ves) then
                Ves=V2
                kes=dum2
                nes=0d0
                ces=P2
            end if
            
            cs_ret_lfp(2,j,ik,ix,iam,ium,Tret-it+1,ifc,LFP_0) = P2
            ks_ret_lfp(2,j,ik,ix,iam,ium,Tret-it+1,ifc,LFP_0) = dum2
            ns_ret_lfp(2,j,ik,ix,iam,ium,Tret-it+1,ifc,LFP_0) = 0d0               

        end if

        vs_ret2(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=ves
        cs_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=ces
        Uprimes_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=dUc(ces)
        ks_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=kes
        ns_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=nes

        vus=vs_ret(1,j,ik,ix,iam,ium,Tret-it+1,ifc)

        if (ves >= vus) then
            vs_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=ves
            rets(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=0d0
            lfps_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=1d0
        else
            vs_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=vus
            rets(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=1d0
            lfps_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=0d0
        end if

        if(it==1) then
            rets(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=1d0
            lfps_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=0d0
        end if

    end do

end subroutine SolveInRetirement2


subroutine SolveInRetirement3(counter)
    ! Woman and Man can is working
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    use bspline_sub_module
    use GlobParams

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixm,ixd,iam,ium,iaf,iuf,iu2,iu3,tprint,ikd,j,ifc,ifcm,ik
    real(8) :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8) :: ces,cus,kes,kus,nes,nus,ves,vus,r_ret,r_ret_next
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y, Psi_pension, Psi_pensionf
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2,pnt2(3),pnt1(2)
    real(8) :: vnext, exp_grid_dum(nexp), INTERP2D(nk,nexp), INTERP3D(nk,nexp,nexp)

    integer :: iflag
    real(8) :: vnext_test
    integer :: idx, idy, idz, iloy, iloz
    integer :: inbvx, inbvy, inbvz
    real(8) :: ww2(ky,kz),ww1(kz),ww0(3*max(kx,ky,kz))
    real(8) :: w1_d2(ky) 
    real(8) :: w0_d2(3*max(kx,ky)) 
    real(8) :: nonlab_inc, aftertax_lab_inc, dum2_test, P4_test
    real(8) :: lfpm_loc, lfpf_loc
    
    !Assigning the grid points
    dum3=((counter*1d0)/(nexp*nu*na*1d0))-0.00001d0
    ik=int(dum3)+1
    dum3=(((counter-(ik-1)*nexp*nu*na)*1d0)/(nu*na*1d0))-0.00001d0
    ix=int(dum3)+1
    dum3=(((counter-(ik-1)*nexp*nu*na-(ix-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
    ium=int(dum3)+1
    iam=counter-(ik-1)*nexp*na*nu-(ix-1)*na*nu-(ium-1)*na

    idx=0
    idy=0
    idz=0
    inbvx=1
    inbvy=1
    inbvz=1
    iloy=1
    iloz=1

    if(it>1) then
        exp_grid_dum=exp_grid(:,T+Tret+2-it)
    end if


    r_ret=((1d0+r)/OmegaRet2(Tret+1-it))-1d0
    if(it>1) then
        r_ret_next=((1d0+r)/OmegaRet2(Tret+2-it))-1d0
    end if

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

    !Married man and woman is working

    do ifc=1,nfc
        do ifcm=1,nfcm
            do iaf = 1, na
                do iuf = 1, nu
                    do ixm = 1, nexp

                        !Pension depends on expected wage conditional on ability and experience
                        Psi_pension=(psi0+psi1*av_earnings(1,1,iam)*min(1d0,exp_grid(ixm,T+Tret+1-it)/35d0))*(1d0-t_const)
                        Psi_pensionf=(psi0+psi1*av_earnings(2,1,iaf)*min(1d0,exp_grid(ix,T+Tret+1-it)/35d0))*(1d0-t_const)                

                        wagem = (1d0/(1d0+exp(kappa*(T+Tret+1-it-agestart))))*wage(1,a(1,iam),exp_grid(ixm,T+Tret+1-it),u(1,ium))/(1d0+t_employer)
                        wagef = (1d0/(1d0+exp(kappa*(T+Tret+1-it-agestart))))*wage(2,a(2,iaf),exp_grid(ix,T+Tret+1-it),u(2,iuf))/(1d0+t_employer)

                        if(it==1) then
                            ! Solve the very last period problem:           
                            ke = 0d0
                            ce = ((k_grid(ik) + Gamma_redistr)*(1d0 + r_ret*(1d0-tk)) + Psi_pension+Psi_pensionf+lumpsum)/(1d0+tc)
                            nef = 0d0
                            nem = 0d0
                            ve = Uc(ce)
                            
                            c_ret_lfp(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm,LFP_M0,LFP_F0)=ce
                            k_ret_lfp(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm,LFP_M0,LFP_F0)=ke
                            nm_ret_lfp(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm,LFP_M0,LFP_F0)=0d0
                            nf_ret_lfp(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm,LFP_M0,LFP_F0)=0d0                                
                            
                        else
                            ! Solve the Tret-1 to 1st period of retirement:
                            !Finding optimal capital by golden search
                            P1=0.001d0
                            nonlab_inc = (k_grid(ik)+Gamma_redistr)*(1d0+r_ret*(1d0-tk))+lumpsum+Psi_pension+Psi_pensionf
                            aftertax_lab_inc = after_tax_labor_inc_married(wagem+wagef)
                            P4 = min((k_grid(nk)-0.001d0)/(1d0+tc), (nonlab_inc + (1d0-t_const)*(aftertax_lab_inc -(wagem+wagef)*t_employee))/(1d0+tc))
                            do
                                P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                                P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                                pnt2 = (/P2, wagem, wagef/)
                                dum4 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborm, nc, nw, nw, pnt2),1d-10),1d0)
                                dum5 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborf, nc, nw, nw, pnt2),1d-10),1d0)
                                y=dum4*wagem+dum5*wagef
                                aftertax_lab_inc = after_tax_labor_inc_married(y)
                                dum2 = (nonlab_inc + (1d0-t_const)*(aftertax_lab_inc-y*t_employee)-P2*(1d0+tc))/(1d0+mu) 
                                if(dum2<0.0001d0) then
                                    V2=-999999999d0
                                elseif(dum2>k_grid(nk)-0.001d0) then
                                    V2=Uc(P2)+Ul(dum4,dum5)-fc(1,ifc)-fcm(1,ifcm)

                                    pnt2=(/dum2, exp_grid(ix,T+Tret-it)+1d0, exp_grid(ixm,T+Tret-it)+1d0/)
                                    INTERP3D=ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,Tret+2-it,ifc,ifcm)
                                    vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                    V2=V2+beta*OmegaRet(Tret-it+1)*vnext
                                else
                                    V2=Uc(P2)+Ul(dum4,dum5)-fc(1,ifc)-fcm(1,ifcm)
                                    call db3val(dum2,exp_grid(ix,T+Tret+1-it)+1d0,exp_grid(ixm,T+Tret+1-it)+1d0,idx,idy,idz,&
                                            tx,ty(:,T+Tret-it+2),tz(:,T+Tret-it+2),&
                                            nk,nexp,nexp,kx,ky,kz,&
                                            ev_ret_bspl(2,2,iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                            inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.false.)
                                    V2=V2+beta*OmegaRet(Tret-it+1)*vnext

                                end if

                                pnt2 = (/P3, wagem, wagef/)
                                dum4 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborm, nc, nw, nw, pnt2),1d-10),1d0)
                                dum5 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborf, nc, nw, nw, pnt2),1d-10),1d0)
                                y=dum4*wagem+dum5*wagef
                                aftertax_lab_inc = after_tax_labor_inc_married(y)
                                dum2 = (nonlab_inc + (1d0-t_const)*(aftertax_lab_inc-y*t_employee)-P3*(1d0+tc))/(1d0+mu) 
                                if(dum2<0.0001d0) then
                                    V3=-999999999d0
                                elseif(dum2>k_grid(nk)-0.001d0) then
                                    V3=Uc(P3)+Ul(dum4,dum5)-fc(1,ifc)-fcm(1,ifcm)
                                    pnt2=(/dum2, exp_grid(ix,T+Tret-it)+1d0, exp_grid(ixm,T+Tret-it)+1d0/)
                                    INTERP3D=ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,Tret+2-it,ifc,ifcm)
                                    vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                    V3=V3+beta*OmegaRet(Tret-it+1)*vnext
                                else
                                    V3=Uc(P3)+Ul(dum4,dum5)-fc(1,ifc)-fcm(1,ifcm)
                                    call db3val(dum2,exp_grid(ix,T+Tret+1-it)+1d0,exp_grid(ixm,T+Tret+1-it)+1d0,idx,idy,idz,&
                                            tx,ty(:,T+Tret-it+2),tz(:,T+Tret-it+2),&
                                            nk,nexp,nexp,kx,ky,kz,&
                                            ev_ret_bspl(2,2,iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                            inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.false.)
                                    V3=V3+beta*OmegaRet(Tret-it+1)*vnext

                                end if


                                if (V2 < V3) then
                                    P1=P2
                                else
                                    P4=P3
                                end if
                                if((P4-P1)<1d-8) exit
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
                            
                            lfpm_loc=1d0
                            lfpm_loc=1d0
                            
                            c_ret_lfp(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm,LFP_M1,LFP_F1)=P2
                            k_ret_lfp(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm,LFP_M1,LFP_F1)=dum2
                            nm_ret_lfp(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm,LFP_M1,LFP_F1)=dum4
                            nf_ret_lfp(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm,LFP_M1,LFP_F1)=dum5
                            

                            !Female not working but not retiring

                            ! Solve the Tret-1 to 1st period of retirement:
                            !Finding optimal capital by golden search
                            P1=0.001d0
                            P4=min((k_grid(nk)-0.001d0)/(1d0+tc),((k_grid(ik)+Gamma_redistr)*(1d0+r_ret*(1d0-tk))+lumpsum+Psi_pension+Psi_pensionf+Unemp_benefit+(wagem)*(1d0-tax_labor(wagem)-tSS_employee(wagem)))/(1d0+tc))
                            do
                                P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                                P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                                pnt1 = (/P2, wagem/)
                                dum4 = min(max(bilin_interp(c_grid, wage_grid, labormwork, nc, nw, pnt1),0d0),1d0)
                                y=dum4*wagem
                                dum2=((k_grid(ik)+Gamma_redistr)*(1d0+r_ret*(1d0-tk))+lumpsum+Psi_pension+Psi_pensionf+Unemp_benefit+y*(1d0-t_const)*(1d0-tax_labor(y)-tSS_employee(y))-P2*(1d0+tc))/(1d0+mu)
                                if(dum2<0.0001d0) then
                                    V2=-999999999d0
                                elseif(dum2>k_grid(nk)-0.001d0) then
                                    V2=Uc(P2)+Ul(dum4,0d0)-fcm(1,ifcm)

                                    pnt2=(/dum2, exp_grid(ix,T+Tret-it), exp_grid(ixm,T+Tret-it)+1d0/)
                                    INTERP3D=ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,Tret+2-it,ifc,ifcm)
                                    vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                    V2=V2+beta*OmegaRet(Tret-it+1)*vnext
                                else
                                    V2=Uc(P2)+Ul(dum4,0d0)-fcm(1,ifcm)
                                    call db3val(dum2,exp_grid(ix,T+Tret+1-it),exp_grid(ixm,T+Tret+1-it)+1d0,idx,idy,idz,&
                                            tx,ty(:,T+Tret-it+2),tz(:,T+Tret-it+2),&
                                            nk,nexp,nexp,kx,ky,kz,&
                                            ev_ret_bspl(2,2,iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                            inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.false.)

                                    !vnext = D_BS3VL(dum2, exp_grid(ix,T+Tret+1-it), exp_grid(ixm,T+Tret+1-it), KORDER, EXPORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T+Tret-it+2),EXP_KNOT(:,T+Tret-it+2), nk, nexp, nexp, ev_ret_spln_coefs(2,2,:,:,:,iam,ium,iaf,iuf,Tret+2-it,ifc,ifcm))
                                    V2=V2+beta*OmegaRet(Tret-it+1)*vnext

                                end if

                                pnt1 = (/P3, wagem/)
                                dum4 = min(max(bilin_interp(c_grid, wage_grid, labormwork, nc, nw, pnt1),0d0),1d0)
                                y=dum4*wagem
                                dum2=((k_grid(ik)+Gamma_redistr)*(1d0+r_ret*(1d0-tk))+lumpsum+Psi_pension+Psi_pensionf+Unemp_benefit+y*(1d0-t_const)*(1d0-tax_labor(y)-tSS_employee(y))-P3*(1d0+tc))/(1d0+mu)
                                if(dum2<0.0001d0) then
                                    V3=-999999999d0
                                elseif(dum2>k_grid(nk)-0.001d0) then
                                    V3=Uc(P3)+Ul(dum4,0d0)-fcm(1,ifcm)
                                    !V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*D_CSVAL(dum2,BREAK,evs_spln_coefs(j,:,:,ix+1,iam,ium,T+1-it,ifc))

                                    pnt2=(/dum2, exp_grid(ix,T+Tret-it), exp_grid(ixm,T+Tret-it)+1d0/)
                                    INTERP3D=ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,Tret+2-it,ifc,ifcm)
                                    vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                    V3=V3+beta*OmegaRet(Tret-it+1)*vnext
                                else
                                    V3=Uc(P3)+Ul(dum4,0d0)-fcm(1,ifcm)
                                    call db3val(dum2,exp_grid(ix,T+Tret+1-it),exp_grid(ixm,T+Tret+1-it)+1d0,idx,idy,idz,&
                                            tx,ty(:,T+Tret-it+2),tz(:,T+Tret-it+2),&
                                            nk,nexp,nexp,kx,ky,kz,&
                                            ev_ret_bspl(2,2,iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                            inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.false.)
                                    !vnext = D_BS3VL(dum2, exp_grid(ix,T+Tret+1-it), exp_grid(ixm,T+Tret+1-it), KORDER, EXPORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T+Tret-it+2),EXP_KNOT(:,T+Tret-it+2), nk, nexp, nexp, ev_ret_spln_coefs(2,2,:,:,:,iam,ium,iaf,iuf,Tret+2-it,ifc,ifcm))
                                    V3=V3+beta*OmegaRet(Tret-it+1)*vnext

                                end if


                                if (V2 < V3) then
                                    P1=P2
                                else
                                    P4=P3
                                end if
                                if((P4-P1)<1d-8) exit
                            end do

                            if(dum4<minhours) then
                                V2=-999999999d0
                            end if
                            
                            if (V2 >= ve) then
                                Ve=V2
                                ke=dum2
                                nem=dum4
                                nef=0d0
                                ce=P2
                                lfpm_loc=1d0
                                lfpf_loc=0d0                                   
                            end if
                            
                            c_ret_lfp(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm,LFP_M1,LFP_F0)=P2
                            k_ret_lfp(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm,LFP_M1,LFP_F0)=dum2
                            nm_ret_lfp(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm,LFP_M1,LFP_F0)=dum4
                            nf_ret_lfp(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm,LFP_M1,LFP_F0)=0d0                           
                            

                            !Male not working but not retiring

                            P1=0.001d0
                            nonlab_inc = (k_grid(ik)+Gamma_redistr)*(1d0+r_ret*(1d0-tk))+lumpsum+Psi_pension+Psi_pensionf+Unemp_benefit
                            aftertax_lab_inc = after_tax_labor_inc_married(wagef)
                            P4 = min((k_grid(nk)-0.001d0)/(1d0+tc), (nonlab_inc + (1d0-t_const)*(aftertax_lab_inc -(wagef)*t_employee))/(1d0+tc))
                            do
                                P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                                P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                                pnt1 = (/P2, wagef/)
                                dum5 = min(max(bilin_interp(c_grid, wage_grid, laborfwork, nc, nw, pnt1),0d0),1d0)
                                y=dum5*wagef
                                aftertax_lab_inc = after_tax_labor_inc_married(y)
                                dum2 = (nonlab_inc + (1d0-t_const)*(aftertax_lab_inc-y*t_employee)-P2*(1d0+tc))/(1d0+mu) 
                                if(dum2<0.0001d0) then
                                    V2=-999999999d0
                                elseif(dum2>k_grid(nk)-0.001d0) then
                                    V2=Uc(P2)+Ul(0d0,dum5)-fc(1,ifc)

                                    pnt2=(/dum2, exp_grid(ix,T+Tret-it)+1d0, exp_grid(ixm,T+Tret-it)/)
                                    INTERP3D=ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,Tret+2-it,ifc,ifcm)
                                    vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                    V2=V2+beta*OmegaRet(Tret-it+1)*vnext
                                else
                                    V2=Uc(P2)+Ul(0d0,dum5)-fc(1,ifc)
                                    call db3val(dum2,exp_grid(ix,T+Tret+1-it)+1d0,exp_grid(ixm,T+Tret+1-it),idx,idy,idz,&
                                            tx,ty(:,T+Tret-it+2),tz(:,T+Tret-it+2),&
                                            nk,nexp,nexp,kx,ky,kz,&
                                            ev_ret_bspl(2,2,iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                            inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.false.)
                                    !vnext = D_BS3VL(dum2, exp_grid(ix,T+Tret+1-it), exp_grid(ixm,T+Tret+1-it), KORDER, EXPORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T+Tret-it+2),EXP_KNOT(:,T+Tret-it+2), nk, nexp, nexp, ev_ret_spln_coefs(2,2,:,:,:,iam,ium,iaf,iuf,Tret+2-it,ifc,ifcm))
                                    V2=V2+beta*OmegaRet(Tret-it+1)*vnext

                                end if

                                pnt1 = (/P3, wagef/)
                                dum5 = min(max(bilin_interp(c_grid, wage_grid, laborfwork, nc, nw, pnt1),0d0),1d0)
                                y=dum5*wagef

                                aftertax_lab_inc = after_tax_labor_inc_married(y)
                                dum2 = (nonlab_inc + (1d0-t_const)*(aftertax_lab_inc-y*t_employee)-P3*(1d0+tc))/(1d0+mu) 

                                if(dum2<0.0001d0) then
                                    V3=-999999999d0
                                elseif(dum2>k_grid(nk)-0.001d0) then
                                    V3=Uc(P3)+Ul(0d0,dum5)-fc(1,ifc)
                                    pnt2=(/dum2, exp_grid(ix,T+Tret-it)+1d0, exp_grid(ixm,T+Tret-it)/)
                                    INTERP3D=ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,Tret+2-it,ifc,ifcm)
                                    vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                    V3=V3+beta*OmegaRet(Tret-it+1)*vnext
                                else
                                    V3=Uc(P3)+Ul(0d0,dum5)-fc(1,ifc)
                                    call db3val(dum2,exp_grid(ix,T+Tret+1-it)+1d0,exp_grid(ixm,T+Tret+1-it),idx,idy,idz,&
                                            tx,ty(:,T+Tret-it+2),tz(:,T+Tret-it+2),&
                                            nk,nexp,nexp,kx,ky,kz,&
                                            ev_ret_bspl(2,2,iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                            inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.false.)
                                    V3=V3+beta*OmegaRet(Tret-it+1)*vnext

                                end if


                                if (V2 < V3) then
                                    P1=P2
                                else
                                    P4=P3
                                end if
                                if((P4-P1)<1d-8) exit
                            end do

                            if(dum5<minhours) then
                                V2=-999999999d0
                            end if
                            
                            if (V2 >= ve) then
                                Ve=V2
                                ke=dum2
                                nem=0d0
                                nef=dum5
                                ce=P2
                                lfpm_loc=0d0
                                lfpf_loc=1d0                                  
                            end if
                            
                            c_ret_lfp(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm,LFP_M0,LFP_F1)=P2
                            k_ret_lfp(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm,LFP_M0,LFP_F1)=dum2
                            nm_ret_lfp(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm,LFP_M0,LFP_F1)=0d0
                            nf_ret_lfp(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm,LFP_M0,LFP_F1)=dum5                             
                            

                            !Both spouses do not work but do not retire

                            P1=0.001d0
                            P4=min((k_grid(nk)-0.001d0)/(1d0+tc),((k_grid(ik)+Gamma_redistr)*(1d0+r_ret*(1d0-tk))+lumpsum+Psi_pension+Psi_pensionf+2d0*Unemp_benefit)/(1d0+tc))
                            do
                                P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                                P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                                dum2=((k_grid(ik)+Gamma_redistr)*(1d0+r_ret*(1d0-tk))+lumpsum+Psi_pension+Psi_pensionf+2d0*Unemp_benefit-P2*(1d0+tc))/(1d0+mu)
                                if(dum2<0.0001d0) then
                                    V2=-999999999d0
                                elseif(dum2>k_grid(nk)-0.001d0) then
                                    V2=Uc(P2)

                                    pnt2=(/dum2, exp_grid(ix,T+Tret-it), exp_grid(ixm,T+Tret-it)/)
                                    INTERP3D=ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,Tret+2-it,ifc,ifcm)
                                    vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                    V2=V2+beta*OmegaRet(Tret-it+1)*vnext
                                else
                                    V2=Uc(P2)
                                    call db3val(dum2,exp_grid(ix,T+Tret+1-it),exp_grid(ixm,T+Tret+1-it),idx,idy,idz,&
                                            tx,ty(:,T+Tret-it+2),tz(:,T+Tret-it+2),&
                                            nk,nexp,nexp,kx,ky,kz,&
                                            ev_ret_bspl(2,2,iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                            inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.false.)

                                    V2=V2+beta*OmegaRet(Tret-it+1)*vnext

                                end if

                                dum2=((k_grid(ik)+Gamma_redistr)*(1d0+r_ret*(1d0-tk))+lumpsum+Psi_pension+Psi_pensionf+2d0*Unemp_benefit-P3*(1d0+tc))/(1d0+mu)
                                if(dum2<0.0001d0) then
                                    V3=-999999999d0
                                elseif(dum2>k_grid(nk)-0.001d0) then
                                    V3=Uc(P3)
                                    pnt2=(/dum2, exp_grid(ix,T+Tret-it), exp_grid(ixm,T+Tret-it)/)
                                    INTERP3D=ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,Tret+2-it,ifc,ifcm)
                                    vnext = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                                    V3=V3+beta*OmegaRet(Tret-it+1)*vnext
                                else
                                    V3=Uc(P3)
                                    call db3val(dum2,exp_grid(ix,T+Tret+1-it),exp_grid(ixm,T+Tret+1-it),idx,idy,idz,&
                                            tx,ty(:,T+Tret-it+2),tz(:,T+Tret-it+2),&
                                            nk,nexp,nexp,kx,ky,kz,&
                                            ev_ret_bspl(2,2,iam,ium,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                            inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.false.)
                                    V3=V3+beta*OmegaRet(Tret-it+1)*vnext
                                end if


                                if (V2 < V3) then
                                    P1=P2
                                else
                                    P4=P3
                                end if
                                if((P4-P1)<1d-8) exit
                            end do

                            if (V2 >= ve) then
                                Ve=V2
                                ke=dum2
                                nem=0d0
                                nef=0d0
                                ce=P2
                                lfpm_loc=0d0
                                lfpf_loc=0d0                                   
                            end if
                            c_ret_lfp(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm,LFP_M0,LFP_F0)=P2
                            k_ret_lfp(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm,LFP_M0,LFP_F0)=dum2
                            nm_ret_lfp(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm,LFP_M0,LFP_F0)=0d0
                            nf_ret_lfp(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm,LFP_M0,LFP_F0)=0d0                            


                        end if

                        v_ret2(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=ve
                        c_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=ce
                        k_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=ke
                        nf_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=nef
                        nm_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=nem

                        vu=v_ret(1,1,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)
                        if (ve >= vu) then
                            v_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=ve
                            retf(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=0d0
                            retm(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=0d0
                        else
                            v_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=vu
                            retf(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=1d0
                            retm(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=1d0
                        end if



                        if (v_ret(2,1,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm) >= v_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)) then
                            v_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=v_ret(2,1,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)
                            retf(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=0d0
                            retm(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=1d0
                        end if

                        if (v_ret(1,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm) >= v_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)) then
                            v_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=v_ret(1,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)
                            retf(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=1d0
                            retm(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=0d0
                        end if

                        if(it==1) then
                            retf(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=1d0
                            retm(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=1d0
                            lfpm_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=0d0
                            lfpf_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=0d0                                                                                    
                        end if

                    end do
                end do
            end do
        end do
    end do


end subroutine SolveInRetirement3
