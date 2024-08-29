program Laffer

    use Utilities
    use Model_Parameters
    use PolicyFunctions
    use SolveRetirement
    use Tauchen
    use PolicyFunctions_obj, only: update_ev_aux, update_ret_policies, update_lfp_policies
    use Simulations_mod
    use Summarize_Simulation
    use LaborSupply
    use, intrinsic :: iso_fortran_env, only: output_unit

    implicit none
    integer :: ik,tprint,it2,it3,it4,it6,it7,it8,it9,ium,iam,iuf,iaf,ix,j,iu2,ik2,ifc,counter,iter_ratio
    real(8) :: dum,dum2,dum3,dum4,dum5,dum6,dum7,epsilon_ratio=1d0,epsilon_ratio_old=1d0,step_ratio=0.05
    real(8) :: P1, P2, P4
    integer :: i_k, i
    integer :: istat
    logical :: solve_lifecycle
    integer :: short_testing
    integer :: iu_simres
    integer :: get_paths
    character(len=15) :: model_version_str
    character(len=15) :: tax_reg_str
    character(len=4) :: deduct_str

    !call OMP_SET_NUM_THREADS(106)
    
    solve_lifecycle = .true.
    call Init_Model_TaxRegime_version()
    
    
    if (pension_for_all == 0) then
        model_version_str = 'Benchmark'
    else
        model_version_str = 'PensionForAll'
    end if
    if (tax_regime == 1) then
        tax_reg_str = '_bench'
    else if (tax_regime == 2) then
        tax_reg_str = '_nit_w_deduct'
    else if (tax_regime == 3) then
        tax_reg_str = '_ubi'
    else if (tax_regime == 4) then
        tax_reg_str = '_eitc'
    else if (tax_regime == 5) then
        write(deduct_str, '(f4.2)') Deduct_Cutoff
        tax_reg_str = '_flattax_'//TRIM(deduct_str)
    else if (tax_regime == 6) then
        tax_reg_str = '_nit1'
    else if (tax_regime == 7) then
        tax_reg_str = '_nit2'        
    end if

    open(newunit=iu_simres, file="Laffer_Results_"//TRIM(model_version_str)//TRIM(tax_reg_str)//".txt")
    call Initialize(iu_simres)

    iter_ratio=1

    !Compute optimal policies in retirement
    t_const = 0.0d0
    write(iu_simres,*) "========================="
    write(iu_simres, "(a, f10.6)") "t_const = ", t_const
    write(iu_simres,*) "========================="


    tax_level_scale=1.0d0
    tau_L = 0.25d0    
    
    epsilon3=1d0
    write(iu_simres,*) '********epsilon3 is',epsilon3
    write(output_unit,*) '********epsilon3 is',epsilon3
    

    if (tax_regime == 1 .or. tax_regime == 3 .or. tax_regime == 5) then
        P2=tax_level_scale
        P4=tax_level_scale+0.15d0
        P1=tax_level_scale-0.15d0
        write(iu_simres,*) 'tax_level_scale is',tax_level_scale
        write(output_unit,*) 'tax_level_scale is',tax_level_scale
    else
        P2=tau_L
        P4=tau_L+0.15d0
        P1=tau_L-0.15d0        
    end if
    write(iu_simres,*) '             '
    write(output_unit,*) '             '

    short_testing = 0
    do while(abs(epsilon3)>0.0001d0)


        thetas = (/0.81773322*tax_level_scale, 0.11060017*tax_prog_scale /)
        theta = (/0.93124354*tax_level_scale, 0.15002363*tax_prog_scale /)
        if (tax_regime == 5 .or. tax_regime == 3) then
            ! Flat tax no deduction
            theta = thetas
        end if
        epsilon_ratio=1d0

        do while(abs(epsilon_ratio)>0.003d0)

            epsilon=1d0
            epsilon2=1d0
            epsilon5=1d0

            iter=0

            !do while(abs(epsilon3)>0.001d0)

            iter=iter+1

            !dum2= 1.5d0*wage(1,a(1,5),dble(T),u(1,5))/(1d0+t_employer)
            dum2= 1.5d0*wage(1,a(1,na),dble(T),u(1,nu))/(1d0+t_employer)
            call MakeGrid(nw,wage_grid,0.01d0,dum2,2d0)
            
            !!$OMP PARALLEL PRIVATE(ik)
            !!$OMP DO SCHEDULE(DYNAMIC)
            !do ik=1,nc
            !    call lsupply(ik)
            !end do
            !!$OMP END DO    
            !!$OMP END PARALLEL
            
            call Solve_Hours(tax_regime)
            
            if (solve_lifecycle == .true.) then
            
            do it=1,Tret

                print *, 't is',T+Tret+1-it

                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter=1,nk*nexp*na
                    call SolveInRetirement(counter, pension_for_all)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL

                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter=1,nk*nexp*na*nu
                    call SolveInRetirement2(counter, pension_for_all)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL

                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter=1,nk*nexp*na*nu
                    call SolveInRetirement3(counter, pension_for_all)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL

                ev_ret(1,1,:,:,:,:,:,:,:,Tret-it+1,:,:)=v_ret(1,1,:,:,:,:,:,:,:,Tret-it+1,:,:)
                edc_ret(1,1,:,:,:,:,:,:,:,Tret-it+1,:,:)=Uprime_ret(1,1,:,:,:,:,:,:,:,Tret-it+1,:,:)
                evs_ret(1,:,:,:,:,:,Tret-it+1,:)=vs_ret(1,:,:,:,:,:,Tret-it+1,:)
                edcs_ret(1,:,:,:,:,:,Tret-it+1,:)=Uprimes_ret(1,:,:,:,:,:,Tret-it+1,:)

                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter=1,nk*na*nu
                    call partest9(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL

                ! Compute spline coefficients:
                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter = 1, nu*na*nfc
                    call partest10(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL

            end do

            !Compute optimal policies at age 64
            !$OMP PARALLEL PRIVATE(counter)
            !$OMP DO SCHEDULE(DYNAMIC)
            do counter=1,nk*nexp*na*nu
                call Solvelastactive(counter)
            end do
            !$OMP END DO    
            !$OMP END PARALLEL

            it=0
            print *, 't is',T-it

            !$OMP PARALLEL PRIVATE(counter)
            !$OMP DO SCHEDULE(DYNAMIC)
            do counter=1,nk*na*nu
                call partest8(counter)
            end do
            !$OMP END DO    
            !$OMP END PARALLEL

            !$OMP PARALLEL PRIVATE(counter)
            !$OMP DO SCHEDULE(DYNAMIC)
            do counter = 1, nu*na*nfc
                call partest5(counter)
            end do
            !$OMP END DO    
            !$OMP END PARALLEL

            !$OMP PARALLEL PRIVATE(counter)
            !$OMP DO SCHEDULE(DYNAMIC)
            do counter=1,nk*nexp*na*nu
                call partest7(counter)
            end do
            !$OMP END DO    
            !$OMP END PARALLEL


            ! Compute spline coefficients:

            !$OMP PARALLEL PRIVATE(counter)
            !$OMP DO SCHEDULE(DYNAMIC)
            do counter = 1, nu*na*nfc
                call partest3(counter)
            end do
            !$OMP END DO    
            !$OMP END PARALLEL

            ! Compute spline coefficients:


            !$OMP PARALLEL PRIVATE(counter)
            !$OMP DO SCHEDULE(DYNAMIC)
            do counter=1,nk*nexp*na*nu
                call partest(counter)
            end do
            !$OMP END DO    
            !$OMP END PARALLEL


            !$OMP PARALLEL PRIVATE(counter)
            !$OMP DO SCHEDULE(DYNAMIC)
            do counter=1,nk*nexp*na*nu
                call partest2(counter)
            end do
            !$OMP END DO    
            !$OMP END PARALLEL

            !$OMP PARALLEL PRIVATE(counter)
            !$OMP DO SCHEDULE(DYNAMIC)
            do counter = 1, nu*na*nfc
                call partest6(counter)
            end do
            !$OMP END DO    
            !$OMP END PARALLEL


            !Compute optimal policies for age 2-63
            do it=1,T-2

                print *, 't is',T-it
                call update_lfp_policies(i_age=T-it+1)
                call update_ev_aux(i_age=T-it+1)                

                !CALL SYSTEM_CLOCK(t1)
                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter=1,nk*nexp*na*nu
                    call SolveActiveLife(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL

                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter=1,nk*na*nu
                    call partest8(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL


                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter = 1, nu*na*nfc
                    call partest5(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL

                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter = 1, nk*nexp*na*nu
                    call partest7(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL   

                ! Compute spline coefficients:

                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter = 1, nfc*na*nu
                    call partest3(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL



                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter=1,nk*nexp*na*nu
                    call partest(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL


                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter=1,nk*nexp*na*nu
                    call partest2(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL

                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter = 1, nfc*na*nu
                    call partest6(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL

            end do


            it=T-1               

            print *, 't is',T-it
            call update_lfp_policies(i_age=T-it+1)
            call update_ev_aux(i_age=T-it+1)             

            !$OMP PARALLEL PRIVATE(counter)
            !$OMP DO SCHEDULE(DYNAMIC)
            do counter=1,nk*na*nu
                call Solvefirstactive(counter)
            end do
            !$OMP END DO    
            !$OMP END PARALLEL
            
            call update_lfp_policies(i_age=1)
            
            !open(71,file='saved_policies.dat')
            !write(71,*) lfpm, lfpf, k_lfp, c_lfp, nm_lfp, nf_lfp, v_lfp, lfps, ks_lfp, cs_lfp, ns, Vs_lfp, V_ret2, k_ret, c_ret, nm_ret, nf_ret, v_ret, Vs_ret2, ks_ret, cs_ret, ns_ret, vs_ret
            !close(71)
            
            else
                open(71,file='saved_policies.dat')
                read(71,*) lfpm, lfpf, k_lfp, c_lfp, nm_lfp, nf_lfp, v_lfp, lfps, ks_lfp, cs_lfp, ns, Vs_lfp, V_ret2, k_ret, c_ret, nm_ret, nf_ret, v_ret, Vs_ret2, ks_ret, cs_ret, ns_ret, vs_ret
                close(71)    
                do it = 1, T
                    call update_lfp_policies(i_age=it)    
                end do
            end if

            print *, 'Simulations: '

            !$OMP PARALLEL PRIVATE(ik)
            !$OMP DO SCHEDULE(DYNAMIC)
            do ik=1,nsim2
                call Simulation(ik)
            end do
            !$OMP END DO    
            !$OMP END PARALLEL

            !write(output_unit,*) 'Simulations: '
            call Statistics(file_id=output_unit)
            !epsilon3=0d0
            !end do

            epsilon_ratio=ratiodum-ratio

            print *, 'epsilon_ratio is',epsilon_ratio

            ! --------------- Prepare for next iteration

            iter_ratio=iter_ratio+1
            ratio=ratio+epsilon_ratio*0.1d0
            print *, "-----------------"
            print *, "iter_ratio = ", iter_ratio
            print *, "ratio = ", ratio

            w=(1d0-alpha)*ratio**alpha
            r=alpha*ratio**(alpha-1d0)-delta
            !end if

            !epsilon_ratio=0d0
        end do

        !STOP

        !t_const = t_const + 0.01d0

        ! --------------- Prepare for next iteration
        write(output_unit,*) '             '
        write(output_unit,*) 'epsilon_ratio loop converged'        
        write(output_unit,*) '             '
        write(iu_simres,*) '             '
        write(iu_simres,*) 'epsilon_ratio loop converged'        
        write(iu_simres,*) '             '
        call Statistics(file_id=iu_simres)
        if (short_testing == 1) STOP

        write(iu_simres,*) '             '
        write(iu_simres,*) '********epsilon3 is',epsilon3        
        if (tax_regime == 1 .or. tax_regime == 3 .or. tax_regime == 5) then
        ! Baseline, UBI and FlatTax tax regimes
            if (epsilon3 < 0d0) then
                P4=P2
            else
                P1=P2
            end if
            P2 = (P4+P1)/2d0            
            tax_level_scale = P2
            write(iu_simres,*) 'tax_level_scale is',tax_level_scale
            tau_L = 1d0-0.81773322*tax_level_scale
            write(iu_simres,*) 'tau_L is',tau_L
        else
            if (epsilon3 > 0d0) then
                P4=P2
            else
                P1=P2
            end if
            P2 = (P4+P1)/2d0              
            tau_L = P2
            write(iu_simres,*) 'tau_L is',tau_L
        end if

    end do


    close(iu_simres)

    get_paths = 1
    if (get_paths == 1) then
        call generate_paths()
    end if
        
contains

    subroutine Init_Model_TaxRegime_version()
    
        implicit none
        
        !if (pension_for_all == 1) then
        !    results_folder = 'pensionforall/'    
        !    write(iunit, *) 'Computing pensionforall model '
        !else
        !    results_folder = 'benchmark/'   
        !    write(iunit, *) 'Computing benchmark model '
        !end if
        
        if (tax_regime == 1) then
            !tax_folder = 'benchmark/'
            Deduct_Cutoff = 0d0
            Deduct_Cutoff_Mar = Deduct_Cutoff
            tax_prog_scale = 1.0d0
            I_ubi = 0
            after_tax_labor_inc_single => after_tax_labor_inc_single_base
            after_tax_labor_inc_married => after_tax_labor_inc_married_base
            !write(iunit, *) 'with benchmark tax system'
        else if (tax_regime == 2) then
            !tax_folder = 'nit_w_deduct/'
            Deduct_Cutoff = 0.1d0
            Deduct_Cutoff_Mar = 2.0d0*Deduct_Cutoff
            yhat = 0.50d0
            yhat_mar = 2d0*yhat            
            I_ubi = 0
            after_tax_labor_inc_single => after_tax_labor_inc_single_nit_w_deduct
            after_tax_labor_inc_married => after_tax_labor_inc_married_nit_w_deduct
            !write(iunit, *) 'with NIT with deduction tax system'
        else if (tax_regime == 3) then
            !tax_folder = 'ubi/'
            Deduct_Cutoff = 0d0
            Deduct_Cutoff_Mar = Deduct_Cutoff
            tax_prog_scale = 1.0d0
            I_ubi = 1
            after_tax_labor_inc_single => after_tax_labor_inc_single_base
            after_tax_labor_inc_married => after_tax_labor_inc_married_base  
            !write(iunit, *) 'with UBI tax system'
        else if (tax_regime == 4) then
            !tax_folder = 'eitc/'
            Deduct_Cutoff = 0.1d0
            Deduct_Cutoff_Mar = 2.0d0*Deduct_Cutoff
            yhat = 0.50d0
            yhat_mar = 2d0*yhat
            b_nit = s_nit*yhat
            b_nit_mar = s_nit*yhat_mar   
            I_ubi = 0
            after_tax_labor_inc_single => after_tax_labor_inc_single_eitc
            after_tax_labor_inc_married => after_tax_labor_inc_married_eitc  
            !write(iunit, *) 'with UBI tax system'    
        else if (tax_regime == 5) then
            theta = thetas
            !tax_folder = 'flattax/'
            Deduct_Cutoff = 0.0d0 ! There will be 3 versions of this, 0.0, 0.2 and 0.4
            Deduct_Cutoff_Mar = 2.0d0*Deduct_Cutoff 
            tax_prog_scale = 0.001d0
            I_ubi = 0
            after_tax_labor_inc_single => after_tax_labor_inc_single_base
            after_tax_labor_inc_married => after_tax_labor_inc_married_base             
            !write(iunit, *) 'with FlatTax tax system' 
        else if (tax_regime == 6) then
            !tax_folder = 'nit'
            Deduct_Cutoff = 0.0d0 ! There will be 3 versions of this
            Deduct_Cutoff_Mar = 2.0d0*Deduct_Cutoff 
            yhat = 0.50d0
            yhat_mar = 2d0*yhat
            b_nit = s_nit*yhat
            b_nit_mar = s_nit*yhat_mar
            !tax_prog_scale = 0.001d0
            I_ubi = 0
            after_tax_labor_inc_single => after_tax_labor_inc_single_nit
            after_tax_labor_inc_married => after_tax_labor_inc_married_nit            
            !write(iunit, *) 'with NIT tax system' 
        else if (tax_regime == 7) then
            !tax_folder = 'nit'
            Deduct_Cutoff = 0.0d0 ! There will be 3 versions of this
            Deduct_Cutoff_Mar = 2.0d0*Deduct_Cutoff 
            yhat = 1.00d0
            yhat_mar = 2d0*yhat
            b_nit = s_nit*yhat
            b_nit_mar = s_nit*yhat_mar
            !tax_prog_scale = 0.001d0
            I_ubi = 0
            after_tax_labor_inc_single => after_tax_labor_inc_single_nit
            after_tax_labor_inc_married => after_tax_labor_inc_married_nit            
            !write(iunit, *) 'with NIT tax system' 
        end if
    
    end subroutine Init_Model_TaxRegime_version


    subroutine Initialize(iunit)

        !USE ANORDF_INT
        implicit none
        
        integer, intent(in) :: iunit
        integer :: iu_tmp
        character(len=4) :: deduct_str
        
        if (pension_for_all == 1) then
            results_folder = 'pensionforall/'    
            write(iunit, *) 'Computing pensionforall model '
        else
            results_folder = 'benchmark/'   
            write(iunit, *) 'Computing benchmark model '
        end if
        
        
        if (tax_regime == 1) then
            tax_folder = 'benchmark/'
            !Deduct_Cutoff = 0d0
            !Deduct_Cutoff_Mar = Deduct_Cutoff
            !tax_prog_scale = 1.0d0
            !I_ubi = 0
            !after_tax_labor_inc_single => after_tax_labor_inc_single_base
            !after_tax_labor_inc_married => after_tax_labor_inc_married_base
            write(iunit, *) 'with benchmark tax system'
        else if (tax_regime == 2) then
            tax_folder = 'nit_w_deduct/'
            !Deduct_Cutoff = 0.1d0
            !Deduct_Cutoff_Mar = 2.0d0*Deduct_Cutoff
            !I_ubi = 0
            !after_tax_labor_inc_single => after_tax_labor_inc_single_nit_w_deduct
            !after_tax_labor_inc_married => after_tax_labor_inc_married_nit_w_deduct
            !write(iunit, *) 'with NIT with deduction tax system'
        else if (tax_regime == 3) then
            tax_folder = 'ubi/'
            !Deduct_Cutoff = 0d0
            !Deduct_Cutoff_Mar = Deduct_Cutoff
            !tax_prog_scale = 1.0d0
            !I_ubi = 1
            !after_tax_labor_inc_single => after_tax_labor_inc_single_base
            !after_tax_labor_inc_married => after_tax_labor_inc_married_base  
            write(iunit, *) 'with UBI tax system'
        else if (tax_regime == 4) then
            tax_folder = 'eitc/'
            !Deduct_Cutoff = 0.1d0
            !Deduct_Cutoff_Mar = 2.0d0*Deduct_Cutoff
            !after_tax_labor_inc_single => after_tax_labor_inc_single_eitc
            !after_tax_labor_inc_married => after_tax_labor_inc_married_eitc  
            write(iunit, *) 'with UBI tax system'    
        else if (tax_regime == 5) then
            !theta = thetas
            write(deduct_str, '(f4.2)') Deduct_Cutoff
            tax_folder = 'flattax_'//TRIM(deduct_str)//'/'
            !Deduct_Cutoff = 0.0d0 ! There will be 3 versions of this
            !Deduct_Cutoff_Mar = 2.0d0*Deduct_Cutoff 
            !tax_prog_scale = 0.001d0
            !I_ubi = 0
            !after_tax_labor_inc_single => after_tax_labor_inc_single_base
            !after_tax_labor_inc_married => after_tax_labor_inc_married_base             
            write(iunit, *) 'with FlatTax tax system and deduction = '//deduct_str 
        else if (tax_regime == 6) then
            tax_folder = 'nit1/'
            !Deduct_Cutoff = 0.0d0 ! There will be 3 versions of this
            !Deduct_Cutoff_Mar = 2.0d0*Deduct_Cutoff 
            !!tax_prog_scale = 0.001d0
            !I_ubi = 0
            !after_tax_labor_inc_single => after_tax_labor_inc_single_nit
            !after_tax_labor_inc_married => after_tax_labor_inc_married_nit            
            write(iunit, *) 'with NIT1 tax system' 
        else if (tax_regime == 7) then
            tax_folder = 'nit2/'
            !Deduct_Cutoff = 0.0d0 ! There will be 3 versions of this
            !Deduct_Cutoff_Mar = 2.0d0*Deduct_Cutoff 
            !!tax_prog_scale = 0.001d0
            !I_ubi = 0
            !after_tax_labor_inc_single => after_tax_labor_inc_single_nit
            !after_tax_labor_inc_married => after_tax_labor_inc_married_nit            
            write(iunit, *) 'with NIT2 tax system'             
        end if
        
        
        open(newunit=iu_tmp, file=trim(results_folder)//trim(tax_folder)//'test.txt')
        write(iu_tmp, '(i0)') testing 
        write(iu_tmp, '(i0)') pension_for_all
        write(iu_tmp, '(i0)') tax_regime
        write(iu_tmp, '(i0)') I_ubi
        close(iu_tmp)        

        allocate(v(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm))
        allocate(ev(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm))
        allocate(c(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm))
        allocate(k(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm))
        allocate(nm(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm))
        allocate(nf(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm))

        allocate(v_aux(nk,nexp,nexp,na,nu,na,nu,nfc,nfcm,2))
        allocate(ev_aux(nk,nexp,nexp,na,nu,na,nu,nfc,nfcm,2))

        allocate(lfpm(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm))
        allocate(lfpf(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm))
        allocate(c_lfp(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm,2,2))
        allocate(v_lfp(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm,2,2,2))
        allocate(k_lfp(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm,2,2))  
        allocate(nm_lfp(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm,2,2))
        allocate(nf_lfp(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm,2,2))   

        allocate(lfps(2,nk,nexp,na,nu,T,nfc))
        allocate(cs_lfp(2,nk,nexp,na,nu,T,nfc,2))
        allocate(ns_lfp(2,nk,nexp,na,nu,T,nfc,2))
        allocate(ks_lfp(2,nk,nexp,na,nu,T,nfc,2))  
        allocate(vs_lfp(2,nk,nexp,na,nu,T,nfc,2))                 

        allocate(vs(2,nk,nexp,na,nu,T,nfc))
        allocate(evs(2,nk,nexp,na,nu,T,nfc))
        allocate(evm(2,nk,nexp,na,nu,T,nfc))
        allocate(cs(2,nk,nexp,na,nu,T,nfc))
        allocate(edcs(2,nk,nexp,na,nu,T,nfc))
        allocate(Uprimes(2,nk,nexp,na,nu,T,nfc))
        allocate(ks(2,nk,nexp,na,nu,T,nfc))
        allocate(ns(2,nk,nexp,na,nu,T,nfc))

        allocate(a_couple(na,na))
        allocate(a(2,na))
        allocate(Prob_a(2,na))

        allocate(fc(2,nfc))
        allocate(Prob_fc(2,nfc))
        allocate(fcm(2,nfcm))
        allocate(Prob_fcm(2,nfcm))
        allocate(u(2,nu))
        allocate(Prob_u(2,nu))
        allocate(trans_u(2,nu,nu))
        allocate(trans_a(2,na,na))
        allocate(trans_fc(2,nfc,nfc))
        allocate(trans_fcm(2,nfcm,nfcm))
        allocate(OmegaRet(Tret))
        allocate(OmegaRet2(Tret))
        allocate(OmegaActive(T))
        allocate(Probm(T))
        allocate(Probd(T))
        !allocate(Share_single(T))
        allocate(WeightRet(Tret))
        allocate(WeightActive(T))

        allocate(fpartner(nk,nexp,na,nu,T,nfc))
        allocate(mpartner(nk,nexp,na,nu,T,nfcm))
        allocate(fpartnerdum(nk,nexp,na,nu,T,nfc))
        allocate(mpartnerdum(nk,nexp,na,nu,T,nfcm))
        allocate(fpartnerdum2(nk,nexp,na,nu,T,nfc))
        allocate(mpartnerdum2(nk,nexp,na,nu,T,nfcm))
        allocate(ability_prob(na,na))
        allocate(av_earnings(2,2,na))
        allocate(laborm(nc,nw,nw))
        allocate(laborf(nc,nw,nw))
        allocate(labormwork(nc,nw))
        allocate(laborfwork(nc,nw))
        allocate(laborsinglem(nc,nw))
        allocate(laborsinglef(nc,nw))

        allocate(c_grid(nc))
        allocate(wage_grid(nw))
        allocate(k_grid(nk))
        allocate(exp_grid(nexp,T+Tret))

        allocate(c_ret(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfcm))
        allocate(edc_ret(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfcm))
        allocate(edc_ret_spln_coefs(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfcm))
        allocate(Uprime_ret(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfcm))
        allocate(v_ret(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfcm))
        allocate(v_ret2(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfcm))
        allocate(ev_ret(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfcm))
        allocate(k_ret(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfcm))
        allocate(nm_ret(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfcm))
        allocate(nf_ret(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfcm))
        allocate(retm(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfcm))
        allocate(retf(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfcm))
        allocate(vs_ret2(2,2,nk,nexp,na,nu,Tret,nfc))
        allocate(vs_ret(2,2,nk,nexp,na,nu,Tret,nfc))
        allocate(evs_ret(2,2,nk,nexp,na,nu,Tret,nfc))
        allocate(evs_ret_spln_coefs(2,2,nk,nexp,na,nu,Tret,nfc))
        allocate(cs_ret(2,2,nk,nexp,na,nu,Tret,nfc))
        allocate(Eulers_ret(2,2,nk,nexp,na,nu,Tret,nfc))
        allocate(edcs_ret(2,2,nk,nexp,na,nu,Tret,nfc))
        allocate(edcs_ret_spln_coefs(2,2,nk,nexp,na,nu,Tret,nfc))
        allocate(Uprimes_ret(2,2,nk,nexp,na,nu,Tret,nfc))
        allocate(ks_ret(2,2,nk,nexp,na,nu,Tret,nfc))
        allocate(ns_ret(2,2,nk,nexp,na,nu,Tret,nfc))
        allocate(Incomes_ret(2,2,nk,nexp,na,nu,Tret,nfc))
        allocate(rets(2,2,nk,nexp,na,nu,Tret,nfc))
        allocate(break(nk))
        
        allocate(lfpm_ret(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfcm))
        allocate(lfpf_ret(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfcm))
        allocate(k_ret_lfp(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfcm,2,2))
        allocate(c_ret_lfp(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfcm,2,2))
        allocate(nm_ret_lfp(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfcm,2,2))
        allocate(nf_ret_lfp(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfcm,2,2))
        allocate(lfps_ret(2,2,nk,nexp,na,nu,Tret,nfc))
        allocate(cs_ret_lfp(2,2,nk,nexp,na,nu,Tret,nfc,2))
        allocate(ks_ret_lfp(2,2,nk,nexp,na,nu,Tret,nfc,2))
        allocate(ns_ret_lfp(2,2,nk,nexp,na,nu,Tret,nfc,2))
        

        open(1, file='singledist.txt')

        read (1, *) fpartner,mpartner

        close(1)

        open(1, file='abilityprob.txt')

        read (1, *) ability_prob

        close(1)

        open(1, file='av_earnings.txt')

        read (1, *) av_earnings

        close(1)

        trans_u = 0d0
        prob_u=0d0
        trans_a = 0d0
        prob_a=0d0
        gamma(1,:) = (/ 0.0605927d0, -0.0010648d0, 0.0000093d0 /)
        gamma(2,:) = (/ 0.0784408d0, -0.0025596d0, 0.0000256d0 /)
        gamma0=0.0676d0
        gamma0f=-0.0685d0

        theta(:) = (/ 0.93124354*tax_level_scale, 0.15002363*tax_prog_scale /)
        thetas(:) = (/ 0.81773322*tax_level_scale, 0.11060017*tax_prog_scale /)
        AE = 1d0
        !Unemp_benefit=0.201795*AE
        !ratio=4.325d0
        Unemp_benefit=0d0
        r=alpha*ratio**(alpha-1d0)-delta
        w=(1d0-alpha)*ratio**alpha
        call MakeGrid(nk,k_grid,0d0,100d0,3d0)
        call MakeGrid(nc,c_grid,0.01d0,100d0,3d0)

        exp_grid=0d0
        do it2=2,T
            call MakeGrid(nexp,exp_grid(:,it2),0d0,1d0*(it2-1),1d0)
        end do

        do it2=T+1,T+Tret
            call MakeGrid(nexp,exp_grid(:,it2),0d0,1d0*(it2-1),1d0)
        end do


        exp_grid(:,1)=exp_grid(:,2)

        !Filling in US divorce and marriage probabilities

        probd(1)=0.11908938
        probd(2)=0.09947980
        probd(3)=0.08335590
        probd(4)=0.07024181
        probd(5)=0.05970394
        probd(6)=0.05134912

        probm(1)=0.08301004
        probm(2)=0.09522406
        probm(3)=0.10423554
        probm(4)=0.11041682
        probm(5)=0.11411954
        probm(6)=0.11567485       

        OmegaActive=1d0

        OmegaRet(1)=1d0-0.014319d0
        OmegaRet(2)=1d0-0.015540d0    

        if (testing == 0) then
            probd(7)=0.04482276
            probd(8)=0.03980704
            probd(9)=0.03601902
            probd(10)=0.03320882
            probd(11)=0.03115778
            probd(12)=0.02967661
            probd(13)=0.02860353
            probd(14)=0.02780249
            probd(15)=0.02716124
            probd(16)=0.02658955
            probd(17)=0.02601734
            probd(18)=0.02539286
            probd(19)=0.02468080
            probd(20)=0.02386051
            probd(21)=0.02292411
            probd(22)=0.02187467
            probd(23)=0.02072434
            probd(24)=0.01949255
            probd(25)=0.01820413
            probd(26)=0.01688749
            probd(27)=0.01557275
            probd(28)=0.01428994
            probd(29)=0.01306712
            probd(30)=0.01192853
            probd(31)=0.01089280
            probd(32)=0.00997105
            probd(33)=0.00916507
            probd(34)=0.00846549
            probd(35)=0.00784992
            probd(36)=0.00728109
            probd(37)=0.00670507
            probd(38)=0.00604934
            probd(39)=0.00522102
            probd(40)=0.00410500
            probd(41)=0.00256208
            probd(42)=0.00042715
            probd(43)=0.00000000
            probd(44)=0.00000000
            probd(45)=0.00000000

            probm(7)=0.11539366
            probm(8)=0.11356688
            probm(9)=0.11046570
            probm(10)=0.10634176
            probm(11)=0.10142749
            probm(12)=0.09593630
            probm(13)=0.09006280
            probm(14)=0.08398313
            probm(15)=0.07785512
            probm(16)=0.07181858
            probm(17)=0.06599554
            probm(18)=0.06049049
            probm(19)=0.05539063
            probm(20)=0.05076610
            probm(21)=0.04667024
            probm(22)=0.04313985
            probm(23)=0.04019541
            probm(24)=0.03784131
            probm(25)=0.03606616
            probm(26)=0.03484296
            probm(27)=0.03412940
            probm(28)=0.03386807
            probm(29)=0.03398674
            probm(30)=0.03439857
            probm(31)=0.03500238
            probm(32)=0.03568288
            probm(33)=0.03631093
            probm(34)=0.03674376
            probm(35)=0.03682526
            probm(36)=0.03638617
            probm(37)=0.03524438
            probm(38)=0.03320513
            probm(39)=0.03006129
            probm(40)=0.02559356
            probm(41)=0.01957079
            probm(42)=0.01175015
            probm(43)=0.00187740
            probm(44)=0.00000000
            probm(45)=0.00000000

            OmegaRet(3)=1d0-0.016920d0
            OmegaRet(4)=1d0-0.018448d0
            OmegaRet(5)=1d0-0.020170d0         
            OmegaRet(6)=1d0-0.022022d0
            OmegaRet(7)=1d0-0.023973d0
            OmegaRet(8)=1d0-0.026203d0
            OmegaRet(9)=1d0-0.028771d0
            OmegaRet(10)=1d0-0.031629d0
            OmegaRet(11)=1d0-0.034611d0
            OmegaRet(12)=1d0-0.037710d0
            OmegaRet(13)=1d0-0.041264d0
            OmegaRet(14)=1d0-0.045405d0
            OmegaRet(15)=1d0-0.050128d0
            OmegaRet(16)=1d0-0.055339d0
            OmegaRet(17)=1d0-0.061005d0
            OmegaRet(18)=1d0-0.067396d0
            OmegaRet(19)=1d0-0.074476d0
            OmegaRet(20)=1d0-0.082272d0
            OmegaRet(21)=1d0-0.091816d0
            OmegaRet(22)=1d0-0.101898d0
            OmegaRet(23)=1d0-0.112870d0
            OmegaRet(24)=1d0-0.124763d0
            OmegaRet(25)=1d0-0.137597d0
            OmegaRet(26)=1d0-0.151383d0
            OmegaRet(27)=1d0-0.166117d0
            OmegaRet(28)=1d0-0.181778d0
            OmegaRet(29)=1d0-0.198331d0
            OmegaRet(30)=1d0-0.215721d0
            OmegaRet(31)=1d0-0.233874d0
            OmegaRet(32)=1d0-0.252699d0
            OmegaRet(33)=1d0-0.272086d0
            OmegaRet(34)=1d0-0.291912d0
            OmegaRet(35)=1d0-0.312040d0
            OmegaRet(36)=1d0-1d0
        end if

        !OmegaRet=1d0

        OmegaRet2(1)=1d0
        do it2=1,Tret-1
            OmegaRet2(it2+1)=OmegaRet(it2)
        end do


        call tauchen_hans(sigma_am,rho_am,na,a(1,:),trans_a(1,:,:),prob_a(1,:))
        call tauchen_hans(sigma_um,rho_um,nu,u(1,:),trans_u(1,:,:),prob_u(1,:))

        call tauchen_hans(sigma_af,rho_af,na,a(2,:),trans_a(2,:,:),prob_a(2,:))
        call tauchen_hans(sigma_uf,rho_uf,nu,u(2,:),trans_u(2,:,:),prob_u(2,:))

        fc(1,:)=exp(mu_fcm)
        fc(2,:)=exp(mu_fcs)
        fcm(1,:)=exp(mu_fcmm)
        fcm(2,:)=exp(mu_fcsm)
        
        call InitSimulation()


        
    end subroutine Initialize  

    subroutine generate_paths()
        implicit none
        integer :: it2, it3, it4
        real(8) :: dum2, dum3
        integer :: gen_earnings_distr
        
        gen_earnings_distr = 0
        
        print *, 'av_earnings(1,1,:) is',av_earnings(1,1,:)
        print *, 'av_earnings(1,2,:) is',av_earnings(1,2,:)
        print *, 'av_earnings(2,1,:) is',av_earnings(2,1,:)
        print *, 'av_earnings(2,2,:) is',av_earnings(2,2,:)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WELFARE BY AGE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Welfare of everyone

        open(1,file=trim(results_folder)//trim(tax_folder)//'welfarepath_all.txt')


        do it2=1,T
            dum2=0.0d0
            dum3=0.0d0
            do it4=1,nsim2
                do it3=1,nsim
                    dum3=dum3+2d0
                    dum2=dum2+Sim1f(it4,it3,it2,16)+(beta**(it2-1))*Sim1f(it4,it3,it2,11)
                    dum2=dum2+Sim1m(it4,it3,it2,16)+(beta**(it2-1))*Sim1m(it4,it3,it2,11)
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F12.5,F12.5)') (19+it2)*1d0, dum2
        end do

        !Welfare of everyone by age after 65


        do it4=1,Tret

            dum2=0d0
            dum3=0d0

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+(SimR1m(it2,it,it4,19)+(beta**(45+(it4-1)))*SimR1m(it2,it,it4,11))*WeightRet(it4)
                    dum2=dum2+(SimR1f(it2,it,it4,19)+(beta**(45+(it4-1)))*SimR1f(it2,it,it4,11))*WeightRet(it4)
                    dum3=dum3+2d0*WeightRet(it4)
                end do
            end do



            dum2=dum2/dum3

            write (1,'(F12.5,F12.5)') (64+it4)*1d0, dum2

        end do

        close(1)



        !Single female welfare

        open(1,file=trim(results_folder)//trim(tax_folder)//'welfarepathsingle_female.txt')


        do it2=1,T
            dum2=0.0d0
            dum3=0.0d0
            do it4=1,nsim2
                do it3=1,nsim
                    if(Sim1f(it4,it3,it2,10)<0.5d0) then
                        dum3=dum3+1d0
                        dum2=dum2+Sim1f(it4,it3,it2,16)+(beta**(it2-1))*Sim1f(it4,it3,it2,11)
                    end if
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F12.5,F12.5)') (19+it2)*1d0, dum2
        end do

        !Single female welfare by age after 65


        do it4=1,Tret

            dum2=0d0
            dum3=0d0

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,T,10)<0.5) then
                        dum2=dum2+(SimR1f(it2,it,it4,19)+(beta**(45+(it4-1)))*SimR1f(it2,it,it4,11))*WeightRet(it4)
                        dum3=dum3+1d0*WeightRet(it4)
                    end if
                end do
            end do



            dum2=dum2/dum3

            write (1,'(F12.5,F12.5)') (64+it4)*1d0, dum2

        end do

        close(1)

        !Single female welfare ability 1

        open(1,file=trim(results_folder)//trim(tax_folder)//'welfarepathsingle_female_a1.txt')


        do it2=1,T
            dum2=0.0d0
            dum3=0.0d0
            do it4=1,nsim2
                do it3=1,nsim
                    if(Sim1f(it4,it3,it2,10)<0.5d0) then
                        if(exp1f(it4,it3,it2,2)<2) then
                            dum3=dum3+1d0
                            dum2=dum2+Sim1f(it4,it3,it2,16)+(beta**(it2-1))*Sim1f(it4,it3,it2,11)
                        end if
                    end if
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F12.5,F12.5)') (19+it2)*1d0, dum2
        end do

        !Single female welfare by age after 65


        do it4=1,Tret

            dum2=0d0
            dum3=0d0

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,T,10)<0.5) then
                        if(exp1f(it2,it,T,2)<2) then
                            dum2=dum2+(SimR1f(it2,it,it4,19)+(beta**(45+(it4-1)))*SimR1f(it2,it,it4,11))*WeightRet(it4)
                            dum3=dum3+1d0*WeightRet(it4)
                        end if
                    end if
                end do
            end do



            dum2=dum2/dum3

            write (1,'(F12.5,F12.5)') (64+it4)*1d0, dum2

        end do

        close(1)

        !Single female welfare ability 1 and 2

        open(1,file=trim(results_folder)//trim(tax_folder)//'welfarepathsingle_female_a12.txt')


        do it2=1,T
            dum2=0.0d0
            dum3=0.0d0
            do it4=1,nsim2
                do it3=1,nsim
                    if(Sim1f(it4,it3,it2,10)<0.5d0) then
                        if(exp1f(it4,it3,it2,2)<3) then
                            dum3=dum3+1d0
                            dum2=dum2+Sim1f(it4,it3,it2,16)+(beta**(it2-1))*Sim1f(it4,it3,it2,11)
                        end if
                    end if
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F12.5,F12.5)') (19+it2)*1d0, dum2
        end do

        !Single female welfare by age after 65


        do it4=1,Tret

            dum2=0d0
            dum3=0d0

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,T,10)<0.5) then
                        if(exp1f(it2,it,T,2)<3) then
                            dum2=dum2+(SimR1f(it2,it,it4,19)+(beta**(45+(it4-1)))*SimR1f(it2,it,it4,11))*WeightRet(it4)
                            dum3=dum3+1d0*WeightRet(it4)
                        end if
                    end if
                end do
            end do



            dum2=dum2/dum3

            write (1,'(F12.5,F12.5)') (64+it4)*1d0, dum2

        end do

        close(1)


        !Married Female Welfare

        open(1,file=trim(results_folder)//trim(tax_folder)//'welfarepathmarried_female.txt')

        do it2=1,T
            dum2=0.0d0
            dum3=0.0d0
            do it4=1,nsim2
                do it3=1,nsim
                    if(Sim1f(it4,it3,it2,10)>0.5d0) then
                        dum3=dum3+1d0
                        dum2=dum2+Sim1f(it4,it3,it2,16)+(beta**(it2-1))*Sim1f(it4,it3,it2,11)
                    end if
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F12.5,F12.5)') (19+it2)*1d0, dum2
        end do

        !Married female welfare by age after 65


        do it4=1,Tret

            dum2=0d0
            dum3=0d0

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,T,10)>0.5) then
                        dum2=dum2+(SimR1f(it2,it,it4,19)+(beta**(45+(it4-1)))*SimR1f(it2,it,it4,11))*WeightRet(it4)
                        dum3=dum3+1d0*WeightRet(it4)
                    end if
                end do
            end do



            dum2=dum2/dum3

            write (1,'(F12.5,F12.5)') (64+it4)*1d0, dum2

        end do

        close(1)


        !Married female welfare ability 1

        open(1,file=trim(results_folder)//trim(tax_folder)//'welfarepathmarried_female_a1.txt')


        do it2=1,T
            dum2=0.0d0
            dum3=0.0d0
            do it4=1,nsim2
                do it3=1,nsim
                    if(Sim1f(it4,it3,it2,10)>0.5d0) then
                        if(exp1f(it4,it3,it2,2)<2) then
                            dum3=dum3+1d0
                            dum2=dum2+Sim1f(it4,it3,it2,16)+(beta**(it2-1))*Sim1f(it4,it3,it2,11)
                        end if
                    end if
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F12.5,F12.5)') (19+it2)*1d0, dum2
        end do

        !Married female welfare by age after 65


        do it4=1,Tret

            dum2=0d0
            dum3=0d0

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,T,10)>0.5) then
                        if(exp1f(it2,it,T,2)<2) then
                            dum2=dum2+(SimR1f(it2,it,it4,19)+(beta**(45+(it4-1)))*SimR1f(it2,it,it4,11))*WeightRet(it4)
                            dum3=dum3+1d0*WeightRet(it4)
                        end if
                    end if
                end do
            end do



            dum2=dum2/dum3

            write (1,'(F12.5,F12.5)') (64+it4)*1d0, dum2

        end do

        close(1)

        !Married female welfare ability 1 and 2

        open(1,file=trim(results_folder)//trim(tax_folder)//'welfarepathmarried_female_a12.txt')


        do it2=1,T
            dum2=0.0d0
            dum3=0.0d0
            do it4=1,nsim2
                do it3=1,nsim
                    if(Sim1f(it4,it3,it2,10)>0.5d0) then
                        if(exp1f(it4,it3,it2,2)<3) then
                            dum3=dum3+1d0
                            dum2=dum2+Sim1f(it4,it3,it2,16)+(beta**(it2-1))*Sim1f(it4,it3,it2,11)
                        end if
                    end if
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F12.5,F12.5)') (19+it2)*1d0, dum2
        end do

        !Married female welfare by age after 65


        do it4=1,Tret

            dum2=0d0
            dum3=0d0

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,T,10)>0.5) then
                        if(exp1f(it2,it,T,2)<3) then
                            dum2=dum2+(SimR1f(it2,it,it4,19)+(beta**(45+(it4-1)))*SimR1f(it2,it,it4,11))*WeightRet(it4)
                            dum3=dum3+1d0*WeightRet(it4)
                        end if
                    end if
                end do
            end do



            dum2=dum2/dum3

            write (1,'(F12.5,F12.5)') (64+it4)*1d0, dum2

        end do

        close(1)

        !Single Male Welfare

        open(1,file=trim(results_folder)//trim(tax_folder)//'welfarepathsingle_male.txt')


        do it2=1,T
            dum2=0.0d0
            dum3=0.0d0
            do it4=1,nsim2
                do it3=1,nsim
                    if(Sim1m(it4,it3,it2,10)<0.5d0) then
                        dum3=dum3+1d0
                        dum2=dum2+Sim1m(it4,it3,it2,16)+(beta**(it2-1))*Sim1m(it4,it3,it2,11)
                    end if
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F12.5,F12.5)') (19+it2)*1d0, dum2
        end do

        !Single male welfare by age after 65


        do it4=1,Tret

            dum2=0d0
            dum3=0d0

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,T,10)<0.5) then
                        dum2=dum2+(SimR1m(it2,it,it4,19)+(beta**(45+(it4-1)))*SimR1m(it2,it,it4,11))*WeightRet(it4)
                        dum3=dum3+1d0*WeightRet(it4)
                    end if
                end do
            end do



            dum2=dum2/dum3

            write (1,'(F12.5,F12.5)') (64+it4)*1d0, dum2

        end do

        close(1)


        !Single male welfare ability 1

        open(1,file=trim(results_folder)//trim(tax_folder)//'welfarepathsingle_male_a1.txt')


        do it2=1,T
            dum2=0.0d0
            dum3=0.0d0
            do it4=1,nsim2
                do it3=1,nsim
                    if(Sim1m(it4,it3,it2,10)<0.5d0) then
                        if(exp1m(it4,it3,it2,2)<2) then
                            dum3=dum3+1d0
                            dum2=dum2+Sim1m(it4,it3,it2,16)+(beta**(it2-1))*Sim1m(it4,it3,it2,11)
                        end if
                    end if
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F12.5,F12.5)') (19+it2)*1d0, dum2
        end do

        !Single male welfare by age after 65


        do it4=1,Tret

            dum2=0d0
            dum3=0d0

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,T,10)<0.5) then
                        if(exp1m(it2,it,T,2)<2) then
                            dum2=dum2+(SimR1m(it2,it,it4,19)+(beta**(45+(it4-1)))*SimR1m(it2,it,it4,11))*WeightRet(it4)
                            dum3=dum3+1d0*WeightRet(it4)
                        end if
                    end if
                end do
            end do



            dum2=dum2/dum3

            write (1,'(F12.5,F12.5)') (64+it4)*1d0, dum2

        end do

        close(1)

        !Single male welfare ability 1 and 2

        open(1,file=trim(results_folder)//trim(tax_folder)//'welfarepathsingle_male_a12.txt')


        do it2=1,T
            dum2=0.0d0
            dum3=0.0d0
            do it4=1,nsim2
                do it3=1,nsim
                    if(Sim1m(it4,it3,it2,10)<0.5d0) then
                        if(exp1m(it4,it3,it2,2)<3) then
                            dum3=dum3+1d0
                            dum2=dum2+Sim1m(it4,it3,it2,16)+(beta**(it2-1))*Sim1m(it4,it3,it2,11)
                        end if
                    end if
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F12.5,F12.5)') (19+it2)*1d0, dum2
        end do

        !Single male welfare by age after 65


        do it4=1,Tret

            dum2=0d0
            dum3=0d0

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,T,10)<0.5) then
                        if(exp1m(it2,it,T,2)<3) then
                            dum2=dum2+(SimR1m(it2,it,it4,19)+(beta**(45+(it4-1)))*SimR1m(it2,it,it4,11))*WeightRet(it4)
                            dum3=dum3+1d0*WeightRet(it4)
                        end if
                    end if
                end do
            end do



            dum2=dum2/dum3

            write (1,'(F12.5,F12.5)') (64+it4)*1d0, dum2

        end do

        close(1)



        open(1,file=trim(results_folder)//trim(tax_folder)//'welfarepathmarried_male.txt')

        do it2=1,T
            dum2=0.0d0
            dum3=0.0d0
            do it4=1,nsim2
                do it3=1,nsim
                    if(Sim1m(it4,it3,it2,10)>0.5d0) then
                        dum3=dum3+1d0
                        dum2=dum2+Sim1m(it4,it3,it2,16)+(beta**(it2-1))*Sim1m(it4,it3,it2,11)
                    end if
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F12.5,F12.5)') (19+it2)*1d0, dum2
        end do

        !Married male welfare by age after 65


        do it4=1,Tret

            dum2=0d0
            dum3=0d0

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,T,10)>0.5) then
                        dum2=dum2+(SimR1m(it2,it,it4,19)+(beta**(45+(it4-1)))*SimR1m(it2,it,it4,11))*WeightRet(it4)
                        dum3=dum3+1d0*WeightRet(it4)
                    end if
                end do
            end do



            dum2=dum2/dum3

            write (1,'(F12.5,F12.5)') (64+it4)*1d0, dum2

        end do


        close(1)

        !Married male welfare ability 1

        open(1,file=trim(results_folder)//trim(tax_folder)//'welfarepathmarried_male_a1.txt')


        do it2=1,T
            dum2=0.0d0
            dum3=0.0d0
            do it4=1,nsim2
                do it3=1,nsim
                    if(Sim1m(it4,it3,it2,10)>0.5d0) then
                        if(exp1m(it4,it3,it2,2)<2) then
                            dum3=dum3+1d0
                            dum2=dum2+Sim1m(it4,it3,it2,16)+(beta**(it2-1))*Sim1m(it4,it3,it2,11)
                        end if
                    end if
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F12.5,F12.5)') (19+it2)*1d0, dum2
        end do

        !Married male welfare by age after 65


        do it4=1,Tret

            dum2=0d0
            dum3=0d0

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,T,10)>0.5) then
                        if(exp1m(it2,it,T,2)<2) then
                            dum2=dum2+(SimR1m(it2,it,it4,19)+(beta**(45+(it4-1)))*SimR1m(it2,it,it4,11))*WeightRet(it4)
                            dum3=dum3+1d0*WeightRet(it4)
                        end if
                    end if
                end do
            end do



            dum2=dum2/dum3

            write (1,'(F12.5,F12.5)') (64+it4)*1d0, dum2

        end do

        close(1)

        !Married male welfare ability 1 and 2

        open(1,file=trim(results_folder)//trim(tax_folder)//'welfarepathmarried_male_a12.txt')


        do it2=1,T
            dum2=0.0d0
            dum3=0.0d0
            do it4=1,nsim2
                do it3=1,nsim
                    if(Sim1m(it4,it3,it2,10)>0.5d0) then
                        if(exp1m(it4,it3,it2,2)<3) then
                            dum3=dum3+1d0
                            dum2=dum2+Sim1m(it4,it3,it2,16)+(beta**(it2-1))*Sim1m(it4,it3,it2,11)
                        end if
                    end if
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F12.5,F12.5)') (19+it2)*1d0, dum2
        end do

        !Married male welfare by age after 65


        do it4=1,Tret

            dum2=0d0
            dum3=0d0

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,T,10)>0.5) then
                        if(exp1m(it2,it,T,2)<3) then
                            dum2=dum2+(SimR1m(it2,it,it4,19)+(beta**(45+(it4-1)))*SimR1m(it2,it,it4,11))*WeightRet(it4)
                            dum3=dum3+1d0*WeightRet(it4)
                        end if
                    end if
                end do
            end do



            dum2=dum2/dum3

            write (1,'(F12.5,F12.5)') (64+it4)*1d0, dum2

        end do

        close(1)



        open(1,file=trim(results_folder)//trim(tax_folder)//'kpath.txt')

        do it2=1,T
            dum2=0d0
            dum3=0d0
            do it4=1,nsim2
                do it3=1,10000
                    dum3=dum3+2d0
                    dum2=dum2+Sim1m(it4,it3,it2,1)
                    if(Sim1f(it4,it3,it2,10)<0.5d0) then
                        dum2=dum2+Sim1f(it4,it3,it2,1)
                    end if
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F8.3,F8.3)') it2*1d0+19d0, dum2
        end do

        do it2=1,Tret
            dum2=0.0
            dum3=0.0
            do it4=1,nsim2
                do it3=1,10000
                    dum3=dum3+2d0
                    dum2=dum2+SimR1m(it4,it3,it2,1)
                    if(Sim1f(it4,it3,T,10)<0.5d0) then
                        dum2=dum2+SimR1f(it4,it3,it2,1)
                    end if
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F8.3,F8.3)') it2*1d0+64d0, dum2
        end do

        close(1)

        open(1,file=trim(results_folder)//trim(tax_folder)//'earningspathm.txt')

        do it2=1,T
            dum2=0d0
            dum3=0d0
            do it4=1,nsim2
                do it3=1,10000
                    dum3=dum3+1d0
                    dum2=dum2+Sim1m(it4,it3,it2,5)
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do

        close(1)

        open(1,file=trim(results_folder)//trim(tax_folder)//'earningspathm2.txt')

        do it2=1,T
            dum2=0d0
            dum3=0d0
            do it4=1,nsim2
                do it3=1,10000
                    if(Sim1m(it4,it3,it2,4)>0.001d0) then
                        dum3=dum3+1d0
                        dum2=dum2+Sim1m(it4,it3,it2,5)
                    end if
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do

        close(1)

        open(1,file=trim(results_folder)//trim(tax_folder)//'earningspathf.txt')

        do it2=1,T
            dum2=0d0
            dum3=0d0
            do it4=1,nsim2
                do it3=1,10000
                    dum3=dum3+1d0
                    dum2=dum2+Sim1f(it4,it3,it2,5)
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do

        close(1)

        open(1,file=trim(results_folder)//trim(tax_folder)//'cpathm.txt')

        do it2=1,T
            dum2=0.0
            do it4=1,nsim2
                do it3=1,10000
                    dum2=dum2+Sim1m(it4,it3,it2,2)
                end do
            end do
            dum2=dum2/160000d0
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do

        do it2=1,Tret
            dum2=0.0
            do it4=1,nsim2
                do it3=1,10000
                    dum2=dum2+SimR1m(it4,it3,it2,2)
                end do
            end do
            dum2=dum2/160000d0
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do

        close(1)

        open(1,file=trim(results_folder)//trim(tax_folder)//'cpathf.txt')

        do it2=1,T
            dum2=0.0
            do it4=1,nsim2
                do it3=1,10000
                    dum2=dum2+Sim1f(it4,it3,it2,2)
                end do
            end do
            dum2=dum2/160000d0
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do

        do it2=1,Tret
            dum2=0.0
            do it4=1,nsim2
                do it3=1,10000
                    dum2=dum2+SimR1f(it4,it3,it2,2)
                end do
            end do
            dum2=dum2/160000d0
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do

        open(1,file=trim(results_folder)//trim(tax_folder)//'kpathm.txt')

        do it2=1,T
            dum2=0.0
            do it4=1,nsim2
                do it3=1,10000
                    dum2=dum2+Sim1m(it4,it3,it2,1)
                end do
            end do
            dum2=dum2/160000d0
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do

        do it2=1,Tret
            dum2=0.0
            do it4=1,nsim2
                do it3=1,10000
                    dum2=dum2+SimR1m(it4,it3,it2,1)
                end do
            end do
            dum2=dum2/160000d0
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do

        close(1)

        open(1,file=trim(results_folder)//trim(tax_folder)//'kpathf.txt')

        do it2=1,T
            dum2=0.0
            do it4=1,nsim2
                do it3=1,10000
                    dum2=dum2+Sim1f(it4,it3,it2,1)
                end do
            end do
            dum2=dum2/160000d0
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do

        do it2=1,Tret
            dum2=0.0
            do it4=1,nsim2
                do it3=1,10000
                    dum2=dum2+SimR1f(it4,it3,it2,1)
                end do
            end do
            dum2=dum2/160000d0
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do

        close(1)

        open(1,file=trim(results_folder)//trim(tax_folder)//'npathm.txt')

        do it2=1,T
            dum2=0.0
            do it4=1,nsim2
                do it3=1,10000
                    dum2=dum2+Sim1m(it4,it3,it2,4)
                end do
            end do
            dum2=dum2/160000d0
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do

        close(1)

        open(1,file=trim(results_folder)//trim(tax_folder)//'npathf.txt')

        do it2=1,T
            dum2=0.0
            do it4=1,nsim2
                do it3=1,10000
                    dum2=dum2+Sim1f(it4,it3,it2,4)
                end do
            end do
            dum2=dum2/160000d0
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do

        close(1)

        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppathsingle.txt')


        do it2=1,T
            dum2=0.0d0
            dum3=0.0d0
            do it4=1,nsim2
                do it3=1,10000
                    if(Sim1f(it4,it3,it2,10)<0.5d0) then
                        dum3=dum3+1d0
                        if(Sim1f(it4,it3,it2,4)>0.001d0) then
                            dum2=dum2+1d0
                        end if
                    end if
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F8.5,F8.5)') it2*1d0, dum2
        end do

        !Single female labor force participation by age after 65


        do it4=1,Tret

            dum2=0d0
            dum3=0d0

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,T,10)<0.5) then
                        if(SimR1f(it2,it,it4,5)>1d-3) then
                            dum2=dum2+(1d0)*WeightRet(it4)
                        end if
                        dum3=dum3+1d0*WeightRet(it4)
                    end if
                end do
            end do



            dum2=dum2/dum3

            write (1,'(F8.5,F8.5)') (64+it4)*1d0, dum2

        end do

        close(1)

        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppathmarried.txt')

        do it2=1,T
            dum2=0.0d0
            dum3=0.0d0
            do it4=1,nsim2
                do it3=1,10000
                    if(Sim1f(it4,it3,it2,10)>0.5d0) then
                        dum3=dum3+1d0
                        if(Sim1f(it4,it3,it2,4)>0.001d0) then
                            dum2=dum2+1d0
                        end if
                    end if
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F8.5,F8.5)') it2*1d0, dum2
        end do

        !Married female labor force participation by age after 65

        do it4=1,Tret

            dum2=0d0
            dum3=0d0    

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,T,10)>0.5) then
                        if(SimR1f(it2,it,it4,5)>1d-3) then
                            dum2=dum2+(1d0)*WeightRet(it4)
                        end if
                        dum3=dum3+1d0*WeightRet(it4)
                    end if
                end do
            end do

            dum2=dum2/dum3

            write (1,'(F8.5,F8.5)') (64+it4)*1d0, dum2

        end do

        close(1)



        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppathsingle_male.txt')

        do i=1,T

            dum2=0.0d0
            dum3=0.0d0    

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,i,10)<0.5) then
                        if(Sim1m(it2,it,i,4)>1d-3) then
                            dum2=dum2+(1d0)*WeightActive(i)
                        end if
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end do
            end do

            dum2=dum2/dum3
            write (1,'(F8.5,F8.5)') i*1d0, dum2

        end do

        !Single male labor force participation by age after 65

        do i=1,Tret

            dum2=0d0
            dum3=0d0

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,T,10)<0.5) then
                        if(SimR1m(it2,it,i,5)>1d-3) then
                            dum2=dum2+(1d0)*WeightRet(i)
                        end if
                        dum3=dum3+1d0*WeightRet(i)
                    end if
                end do
            end do

            dum2=dum2/dum3
            write (1,'(F8.5,F8.5)') (64+i)*1d0, dum2

        end do


        close(1)

        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppathmarried_male.txt')

        do i=1,T

            dum2=0d0
            dum3=0d0

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,i,10)>0.5) then
                        if(Sim1m(it2,it,i,4)>1d-3) then
                            dum2=dum2+(1d0)*WeightActive(i)
                        end if
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end do
            end do

            dum2=dum2/dum3

            write (1,'(F8.5,F8.5)') i*1d0, dum2

        end do


        !Married male labor force participation by age after 65


        do i=1,Tret

            dum2=0d0
            dum3=0d0

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,T,10)>0.5) then
                        if(SimR1m(it2,it,i,5)>1d-3) then
                            dum2=dum2+(1d0)*WeightRet(i)
                        end if
                        dum3=dum3+1d0*WeightRet(i)
                    end if
                end do
            end do

            dum2=dum2/dum3

            write (1,'(F8.5,F8.5)') (64+i)*1d0, dum2

        end do


        close(1)


        !!!!!!!!!!!!!!!LFP by Fixed Cost !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !Single Females

        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppathsinglefc.txt')

        do it6=1,nfc

            do it2=1,T


                dum2=0.0d0
                dum3=0.0d0
                do it4=1,nsim2
                    do it3=1,10000
                        if(exp1f(it4,it3,it2,6)==it6) then

                            if(Sim1f(it4,it3,it2,10)<0.5d0) then
                                dum3=dum3+1d0
                                if(Sim1f(it4,it3,it2,4)>0.001d0) then
                                    dum2=dum2+1d0
                                end if
                            end if

                        end if
                    end do
                end do
                dum2=dum2/dum3

                write (1,'(F12.5,F12.5,F12.5,F12.5)') it6*1d0, (it2+19)*1d0, Fc(2,it6), dum2

            end do

            do it4=1,Tret

                dum2=0d0
                dum3=0d0

                do it2=1,nsim2
                    do it=1,nsim
                        if(exp1f(it2,it,T,6)==it6) then
                            if(Sim1f(it2,it,T,10)<0.5) then
                                if(SimR1f(it2,it,it4,5)>1d-3) then
                                    dum2=dum2+(1d0)*WeightRet(it4)
                                end if
                                dum3=dum3+1d0*WeightRet(it4)
                            end if
                        end if
                    end do
                end do



                dum2=dum2/dum3

                write (1,'(F12.5,F12.5,F12.5,F12.5)') it6*1d0, (64+it4)*1d0, Fc(2,it6), dum2

            end do

        end do

        close(1)




        !Married Females

        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppathmarriedfc.txt')

        do it6=1,nfc

            do it2=1,T


                dum2=0.0d0
                dum3=0.0d0
                do it4=1,nsim2
                    do it3=1,10000
                        if(exp1f(it4,it3,it2,6)==it6) then

                            if(Sim1f(it4,it3,it2,10)>0.5d0) then
                                dum3=dum3+1d0
                                if(Sim1f(it4,it3,it2,4)>0.001d0) then
                                    dum2=dum2+1d0
                                end if
                            end if

                        end if
                    end do
                end do
                dum2=dum2/dum3

                write (1,'(F12.5,F12.5,F12.5,F12.5)') it6*1d0, (it2+19)*1d0, Fc(1,it6), dum2

            end do

            do it4=1,Tret

                dum2=0d0
                dum3=0d0

                do it2=1,nsim2
                    do it=1,nsim
                        if(exp1f(it2,it,T,6)==it6) then
                            if(Sim1f(it2,it,T,10)>0.5) then
                                if(SimR1f(it2,it,it4,5)>1d-3) then
                                    dum2=dum2+(1d0)*WeightRet(it4)
                                end if
                                dum3=dum3+1d0*WeightRet(it4)
                            end if
                        end if
                    end do
                end do



                dum2=dum2/dum3

                write (1,'(F12.5,F12.5,F12.5,F12.5)') it6*1d0, (64+it4)*1d0, Fc(1,it6), dum2

            end do

        end do

        close(1)




        !Single Males

        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppathsinglemalefc.txt')

        do it6=1,nfcm

            do it2=1,T


                dum2=0.0d0
                dum3=0.0d0
                do it4=1,nsim2
                    do it3=1,10000
                        if(exp1m(it4,it3,it2,6)==it6) then

                            if(Sim1m(it4,it3,it2,10)<0.5d0) then
                                dum3=dum3+1d0
                                if(Sim1m(it4,it3,it2,4)>0.001d0) then
                                    dum2=dum2+1d0
                                end if
                            end if

                        end if
                    end do
                end do
                dum2=dum2/dum3

                write (1,'(F12.5,F12.5,F12.5,F12.5)') it6*1d0, (it2+19)*1d0, Fcm(2,it6), dum2

            end do

            do it4=1,Tret

                dum2=0d0
                dum3=0d0

                do it2=1,nsim2
                    do it=1,nsim
                        if(exp1m(it2,it,T,6)==it6) then
                            if(Sim1m(it2,it,T,10)<0.5) then
                                if(SimR1m(it2,it,it4,5)>1d-3) then
                                    dum2=dum2+(1d0)*WeightRet(it4)
                                end if
                                dum3=dum3+1d0*WeightRet(it4)
                            end if
                        end if
                    end do
                end do



                dum2=dum2/dum3

                write (1,'(F12.5,F12.5,F12.5,F12.5)') it6*1d0, (64+it4)*1d0, Fcm(2,it6), dum2

            end do

        end do

        close(1)




        !Married Males

        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppathmarriedmalefc.txt')

        do it6=1,nfcm

            do it2=1,T


                dum2=0.0d0
                dum3=0.0d0
                do it4=1,nsim2
                    do it3=1,10000
                        if(exp1m(it4,it3,it2,6)==it6) then

                            if(Sim1m(it4,it3,it2,10)>0.5d0) then
                                dum3=dum3+1d0
                                if(Sim1m(it4,it3,it2,4)>0.001d0) then
                                    dum2=dum2+1d0
                                end if
                            end if

                        end if
                    end do
                end do
                dum2=dum2/dum3

                write (1,'(F12.5,F12.5,F12.5,F12.5)') it6*1d0, (it2+19)*1d0, Fcm(1,it6), dum2

            end do

            do it4=1,Tret

                dum2=0d0
                dum3=0d0

                do it2=1,nsim2
                    do it=1,nsim
                        if(exp1m(it2,it,T,6)==it6) then
                            if(Sim1m(it2,it,T,10)>0.5) then
                                if(SimR1m(it2,it,it4,5)>1d-3) then
                                    dum2=dum2+(1d0)*WeightRet(it4)
                                end if
                                dum3=dum3+1d0*WeightRet(it4)
                            end if
                        end if
                    end do
                end do



                dum2=dum2/dum3

                write (1,'(F12.5,F12.5,F12.5,F12.5)') it6*1d0, (64+it4)*1d0, Fcm(1,it6), dum2

            end do

        end do

        close(1)

        !!!!!!!!!!!!!!!!!LFP by a,u,f !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Single Females

        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppathsingle_auf.txt')

        do it8=1,nfc

            do it7=1,na

                do it6=1,nu

                    do it2=1,T


                        dum2=0.0d0
                        dum3=0.0d0
                        dum4=0.0d0
                        dum5=0.0d0

                        do it4=1,nsim2
                            do it3=1,10000

                                if(exp1f(it4,it3,it2,6)==it8) then

                                    if(exp1f(it4,it3,it2,2)==it7) then

                                        if(exp1f(it4,it3,it2,3)==it6) then

                                            if(Sim1f(it4,it3,it2,10)<0.5d0) then
                                                dum3=dum3+1d0
                                                dum4=dum4+Sim1f(it4,it3,it2,1)
                                                dum5=dum5+Sim1f(it4,it3,it2,2)
                                                if(Sim1f(it4,it3,it2,4)>0.001d0) then
                                                    dum2=dum2+1d0
                                                end if
                                            end if

                                        end if

                                    end if

                                end if

                            end do
                        end do
                        dum2=dum2/dum3
                        dum4=dum4/dum3
                        dum5=dum5/dum3

                        write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5)') it8*1d0,it7*1d0, it6*1d0, (it2+19)*1d0, Fc(2,it8), A(2,it7), U(2,it6), dum2, dum4, dum5,0d0

                    end do

                    do it4=1,Tret

                        dum2=0d0
                        dum3=0d0
                        dum4=0d0
                        dum5=0d0
                        dum6=0d0

                        do it2=1,nsim2
                            do it=1,nsim
                                if(exp1f(it2,it,T,6)==it8) then
                                    if(exp1f(it2,it,T,2)==it7) then
                                        if(exp1f(it2,it,T,3)==it6) then
                                            if(Sim1f(it2,it,T,10)<0.5) then
                                                if(SimR1f(it2,it,it4,5)>1d-3) then
                                                    dum2=dum2+(1d0)
                                                end if
                                                dum3=dum3+1d0
                                                dum4=dum4+SimR1f(it2,it,it4,1)
                                                dum5=dum5+SimR1f(it2,it,it4,2)
                                                dum6=dum6+SimR1f(it2,it,it4,13)
                                            end if
                                        end if
                                    end if
                                end if
                            end do
                        end do



                        dum2=dum2/dum3
                        dum4=dum4/dum3
                        dum5=dum5/dum3
                        dum6=dum6/dum3

                        write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5)') it8*1d0, it7*1d0, it6*1d0, (64+it4)*1d0, Fc(2,it8), A(2,it7), U(2,it6), dum2, dum4, dum5, dum6

                    end do

                end do

            end do

        end do

        close(1)


        !Single Males

        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppathsinglemale_auf.txt')

        do it8=1,nfcm

            do it7=1,na

                do it6=1,nu

                    do it2=1,T


                        dum2=0.0d0
                        dum3=0.0d0
                        dum4=0.0d0
                        dum5=0.0d0

                        do it4=1,nsim2
                            do it3=1,10000

                                if(exp1m(it4,it3,it2,6)==it8) then

                                    if(exp1m(it4,it3,it2,2)==it7) then

                                        if(exp1m(it4,it3,it2,3)==it6) then

                                            if(Sim1m(it4,it3,it2,10)<0.5d0) then
                                                dum3=dum3+1d0
                                                dum4=dum4+Sim1m(it4,it3,it2,1)
                                                dum5=dum5+Sim1m(it4,it3,it2,2)

                                                if(Sim1m(it4,it3,it2,4)>0.001d0) then
                                                    dum2=dum2+1d0
                                                end if
                                            end if

                                        end if

                                    end if

                                end if

                            end do
                        end do
                        dum2=dum2/dum3
                        dum4=dum4/dum3
                        dum5=dum5/dum3

                        write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5)') it8*1d0, it7*1d0, it6*1d0, (it2+19)*1d0, Fcm(2,it8), A(1,it7), U(1,it6), dum2, dum4, dum5, 0d0

                    end do

                    do it4=1,Tret

                        dum2=0d0
                        dum3=0d0
                        dum4=0d0
                        dum5=0d0
                        dum6=0d0

                        do it2=1,nsim2
                            do it=1,nsim
                                if(exp1m(it2,it,T,6)==it8) then
                                    if(exp1m(it2,it,T,2)==it7) then
                                        if(exp1m(it2,it,T,3)==it6) then
                                            if(Sim1m(it2,it,T,10)<0.5) then
                                                if(SimR1m(it2,it,it4,5)>1d-3) then
                                                    dum2=dum2+(1d0)
                                                end if
                                                dum3=dum3+1d0
                                                dum4=dum4+SimR1m(it2,it,it4,1)
                                                dum5=dum5+SimR1m(it2,it,it4,2)
                                                dum6=dum6+SimR1m(it2,it,it4,13)
                                            end if
                                        end if
                                    end if
                                end if
                            end do
                        end do



                        dum2=dum2/dum3
                        dum4=dum4/dum3
                        dum5=dum5/dum3
                        dum6=dum6/dum3

                        write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5)') it8*1d0, it7*1d0, it6*1d0, (64+it4)*1d0, Fcm(2,it8), A(1,it7), U(1,it6), dum2, dum4, dum5, dum6

                    end do

                end do

            end do

        end do

        close(1)

        !Married Females

        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppathmarried_auf.txt')

        do it8=1,nfc

            do it7=1,na

                do it6=1,nu

                    do it2=1,T


                        dum2=0.0d0
                        dum3=0.0d0
                        dum4=0.0d0
                        dum5=0.0d0

                        do it4=1,nsim2
                            do it3=1,10000

                                if(exp1f(it4,it3,it2,6)==it8) then

                                    if(exp1f(it4,it3,it2,2)==it7) then

                                        if(exp1f(it4,it3,it2,3)==it6) then

                                            if(Sim1f(it4,it3,it2,10)>0.5d0) then
                                                dum3=dum3+1d0
                                                dum4=dum4+Sim1f(it4,it3,it2,1)
                                                dum5=dum5+Sim1f(it4,it3,it2,2)
                                                if(Sim1f(it4,it3,it2,4)>0.001d0) then
                                                    dum2=dum2+1d0
                                                end if
                                            end if

                                        end if

                                    end if

                                end if

                            end do
                        end do
                        dum2=dum2/dum3
                        dum4=dum4/dum3
                        dum5=dum5/dum3

                        write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5)') it8*1d0,it7*1d0, it6*1d0, (it2+19)*1d0, Fc(2,it8), A(2,it7), U(2,it6), dum2, dum4, dum5,0d0

                    end do

                    do it4=1,Tret

                        dum2=0d0
                        dum3=0d0
                        dum4=0d0
                        dum5=0d0
                        dum6=0d0

                        do it2=1,nsim2
                            do it=1,nsim
                                if(exp1f(it2,it,T,6)==it8) then
                                    if(exp1f(it2,it,T,2)==it7) then
                                        if(exp1f(it2,it,T,3)==it6) then
                                            if(Sim1f(it2,it,T,10)>0.5) then
                                                if(SimR1f(it2,it,it4,5)>1d-3) then
                                                    dum2=dum2+(1d0)
                                                end if
                                                dum3=dum3+1d0
                                                dum4=dum4+SimR1f(it2,it,it4,1)
                                                dum5=dum5+SimR1f(it2,it,it4,2)
                                                dum6=dum6+SimR1f(it2,it,it4,13)
                                            end if
                                        end if
                                    end if
                                end if
                            end do
                        end do



                        dum2=dum2/dum3
                        dum4=dum4/dum3
                        dum5=dum5/dum3
                        dum6=dum6/dum3

                        write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5)') it8*1d0, it7*1d0, it6*1d0, (64+it4)*1d0, Fc(2,it8), A(2,it7), U(2,it6), dum2, dum4, dum5, dum6

                    end do

                end do

            end do

        end do

        close(1)


        !Married Males

        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppathmarriedmale_auf.txt')

        do it8=1,nfcm

            do it7=1,na

                do it6=1,nu

                    do it2=1,T


                        dum2=0.0d0
                        dum3=0.0d0
                        dum4=0.0d0
                        dum5=0.0d0

                        do it4=1,nsim2
                            do it3=1,10000

                                if(exp1m(it4,it3,it2,6)==it8) then

                                    if(exp1m(it4,it3,it2,2)==it7) then

                                        if(exp1m(it4,it3,it2,3)==it6) then

                                            if(Sim1m(it4,it3,it2,10)>0.5d0) then
                                                dum3=dum3+1d0
                                                dum4=dum4+Sim1m(it4,it3,it2,1)
                                                dum5=dum5+Sim1m(it4,it3,it2,2)

                                                if(Sim1m(it4,it3,it2,4)>0.001d0) then
                                                    dum2=dum2+1d0
                                                end if
                                            end if

                                        end if

                                    end if

                                end if

                            end do
                        end do
                        dum2=dum2/dum3
                        dum4=dum4/dum3
                        dum5=dum5/dum3

                        write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5)') it8*1d0, it7*1d0, it6*1d0, (it2+19)*1d0, Fcm(2,it8), A(1,it7), U(1,it6), dum2, dum4, dum5, 0d0

                    end do

                    do it4=1,Tret

                        dum2=0d0
                        dum3=0d0
                        dum4=0d0
                        dum5=0d0
                        dum6=0d0

                        do it2=1,nsim2
                            do it=1,nsim
                                if(exp1m(it2,it,T,6)==it8) then
                                    if(exp1m(it2,it,T,2)==it7) then
                                        if(exp1m(it2,it,T,3)==it6) then
                                            if(Sim1m(it2,it,T,10)>0.5) then
                                                if(SimR1m(it2,it,it4,5)>1d-3) then
                                                    dum2=dum2+(1d0)
                                                end if
                                                dum3=dum3+1d0
                                                dum4=dum4+SimR1m(it2,it,it4,1)
                                                dum5=dum5+SimR1m(it2,it,it4,2)
                                                dum6=dum6+SimR1m(it2,it,it4,13)
                                            end if
                                        end if
                                    end if
                                end if
                            end do
                        end do



                        dum2=dum2/dum3
                        dum4=dum4/dum3
                        dum5=dum5/dum3
                        dum6=dum6/dum3

                        write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5)') it8*1d0, it7*1d0, it6*1d0, (64+it4)*1d0, Fcm(2,it8), A(1,it7), U(1,it6), dum2, dum4, dum5, dum6

                    end do

                end do

            end do

        end do

        close(1)



        !!!!!!!!!!!!!!!!!LFP by a !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Single Females

        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppathsingle_a.txt')


        do it7=1,na

            do it2=1,T


                dum2=0.0d0
                dum3=0.0d0
                dum4=0.0d0
                dum5=0.0d0

                do it4=1,nsim2
                    do it3=1,10000

                        if(exp1f(it4,it3,it2,2)==it7) then

                            if(Sim1f(it4,it3,it2,10)<0.5d0) then
                                dum3=dum3+1d0
                                dum4=dum4+Sim1f(it4,it3,it2,1)
                                dum5=dum5+Sim1f(it4,it3,it2,2)
                                if(Sim1f(it4,it3,it2,4)>0.001d0) then
                                    dum2=dum2+1d0
                                end if
                            end if

                        end if

                    end do
                end do
                dum2=dum2/dum3
                dum4=dum4/dum3
                dum5=dum5/dum3

                write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5)') it8*1d0,it7*1d0, it6*1d0, (it2+19)*1d0, A(2,it7), dum2, dum4, dum5,0d0

            end do

            do it4=1,Tret

                dum2=0d0
                dum3=0d0
                dum4=0d0
                dum5=0d0
                dum6=0d0

                do it2=1,nsim2
                    do it=1,nsim
                        if(exp1f(it2,it,T,2)==it7) then
                            if(Sim1f(it2,it,T,10)<0.5) then
                                if(SimR1f(it2,it,it4,5)>1d-3) then
                                    dum2=dum2+(1d0)
                                end if
                                dum3=dum3+1d0
                                dum4=dum4+SimR1f(it2,it,it4,1)
                                dum5=dum5+SimR1f(it2,it,it4,2)
                                dum6=dum6+SimR1f(it2,it,it4,13)
                            end if
                        end if
                    end do
                end do


                dum2=dum2/dum3
                dum4=dum4/dum3
                dum5=dum5/dum3
                dum6=dum6/dum3

                write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5)') it8*1d0, it7*1d0, it6*1d0, (64+it4)*1d0, A(2,it7), dum2, dum4, dum5, dum6

            end do

        end do


        close(1)


        !Single Males

        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppathsinglemale_a.txt')

        do it7=1,na


            do it2=1,T


                dum2=0.0d0
                dum3=0.0d0
                dum4=0.0d0
                dum5=0.0d0

                do it4=1,nsim2
                    do it3=1,10000


                        if(exp1m(it4,it3,it2,2)==it7) then


                            if(Sim1m(it4,it3,it2,10)<0.5d0) then
                                dum3=dum3+1d0
                                dum4=dum4+Sim1m(it4,it3,it2,1)
                                dum5=dum5+Sim1m(it4,it3,it2,2)

                                if(Sim1m(it4,it3,it2,4)>0.001d0) then
                                    dum2=dum2+1d0
                                end if
                            end if

                        end if

                    end do
                end do
                dum2=dum2/dum3
                dum4=dum4/dum3
                dum5=dum5/dum3

                write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5)') it8*1d0, it7*1d0, it6*1d0, (it2+19)*1d0, A(1,it7), dum2, dum4, dum5, 0d0

            end do

            do it4=1,Tret

                dum2=0d0
                dum3=0d0
                dum4=0d0
                dum5=0d0
                dum6=0d0

                do it2=1,nsim2
                    do it=1,nsim
                        if(exp1m(it2,it,T,2)==it7) then
                            if(Sim1m(it2,it,T,10)<0.5) then
                                if(SimR1m(it2,it,it4,5)>1d-3) then
                                    dum2=dum2+(1d0)
                                end if
                                dum3=dum3+1d0
                                dum4=dum4+SimR1m(it2,it,it4,1)
                                dum5=dum5+SimR1m(it2,it,it4,2)
                                dum6=dum6+SimR1m(it2,it,it4,13)
                            end if
                        end if
                    end do
                end do



                dum2=dum2/dum3
                dum4=dum4/dum3
                dum5=dum5/dum3
                dum6=dum6/dum3

                write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5)') it8*1d0, it7*1d0, it6*1d0, (64+it4)*1d0, A(1,it7), dum2, dum4, dum5, dum6

            end do

        end do

        close(1)

        !Married Females

        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppathmarried_a.txt')

        do it7=1,na

            do it2=1,T


                dum2=0.0d0
                dum3=0.0d0
                dum4=0.0d0
                dum5=0.0d0

                do it4=1,nsim2
                    do it3=1,10000

                        if(exp1f(it4,it3,it2,2)==it7) then


                            if(Sim1f(it4,it3,it2,10)>0.5d0) then
                                dum3=dum3+1d0
                                dum4=dum4+Sim1f(it4,it3,it2,1)
                                dum5=dum5+Sim1f(it4,it3,it2,2)
                                if(Sim1f(it4,it3,it2,4)>0.001d0) then
                                    dum2=dum2+1d0
                                end if
                            end if

                        end if


                    end do
                end do
                dum2=dum2/dum3
                dum4=dum4/dum3
                dum5=dum5/dum3

                write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5)') it8*1d0,it7*1d0, it6*1d0, (it2+19)*1d0, A(2,it7), dum2, dum4, dum5,0d0

            end do

            do it4=1,Tret

                dum2=0d0
                dum3=0d0
                dum4=0d0
                dum5=0d0
                dum6=0d0

                do it2=1,nsim2
                    do it=1,nsim
                        if(exp1f(it2,it,T,2)==it7) then
                            if(Sim1f(it2,it,T,10)>0.5) then
                                if(SimR1f(it2,it,it4,5)>1d-3) then
                                    dum2=dum2+(1d0)
                                end if
                                dum3=dum3+1d0
                                dum4=dum4+SimR1f(it2,it,it4,1)
                                dum5=dum5+SimR1f(it2,it,it4,2)
                                dum6=dum6+SimR1f(it2,it,it4,13)
                            end if
                        end if
                    end do
                end do



                dum2=dum2/dum3
                dum4=dum4/dum3
                dum5=dum5/dum3
                dum6=dum6/dum3

                write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5)') it8*1d0, it7*1d0, it6*1d0, (64+it4)*1d0, A(2,it7), dum2, dum4, dum5, dum6

            end do

        end do


        close(1)


        !Married Males

        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppathmarriedmale_a.txt')

        do it7=1,na

            do it2=1,T


                dum2=0.0d0
                dum3=0.0d0
                dum4=0.0d0
                dum5=0.0d0

                do it4=1,nsim2
                    do it3=1,10000

                        if(exp1m(it4,it3,it2,2)==it7) then


                            if(Sim1m(it4,it3,it2,10)>0.5d0) then
                                dum3=dum3+1d0
                                dum4=dum4+Sim1m(it4,it3,it2,1)
                                dum5=dum5+Sim1m(it4,it3,it2,2)

                                if(Sim1m(it4,it3,it2,4)>0.001d0) then
                                    dum2=dum2+1d0
                                end if
                            end if

                        end if

                    end do
                end do
                dum2=dum2/dum3
                dum4=dum4/dum3
                dum5=dum5/dum3

                write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5)') it8*1d0, it7*1d0, it6*1d0, (it2+19)*1d0, A(1,it7), dum2, dum4, dum5, 0d0

            end do

            do it4=1,Tret

                dum2=0d0
                dum3=0d0
                dum4=0d0
                dum5=0d0
                dum6=0d0

                do it2=1,nsim2
                    do it=1,nsim
                        if(exp1m(it2,it,T,2)==it7) then
                            if(Sim1m(it2,it,T,10)>0.5) then
                                if(SimR1m(it2,it,it4,5)>1d-3) then
                                    dum2=dum2+(1d0)
                                end if
                                dum3=dum3+1d0
                                dum4=dum4+SimR1m(it2,it,it4,1)
                                dum5=dum5+SimR1m(it2,it,it4,2)
                                dum6=dum6+SimR1m(it2,it,it4,13)
                            end if
                        end if
                    end do
                end do


                dum2=dum2/dum3
                dum4=dum4/dum3
                dum5=dum5/dum3
                dum6=dum6/dum3

                write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5)') it8*1d0, it7*1d0, it6*1d0, (64+it4)*1d0, A(1,it7), dum2, dum4, dum5, dum6

            end do

        end do


        close(1)

        !!!!!!!!!!!!!!!!!LFP by a (not by marital status) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Single Females

        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppath_females_a.txt')


        do it7=1,na

            do it2=1,T


                dum2=0.0d0
                dum3=0.0d0
                dum4=0.0d0
                dum5=0.0d0

                do it4=1,nsim2
                    do it3=1,10000

                        if(exp1f(it4,it3,it2,2)==it7) then

                            dum3=dum3+1d0
                            dum4=dum4+Sim1f(it4,it3,it2,1)
                            dum5=dum5+Sim1f(it4,it3,it2,2)
                            if(Sim1f(it4,it3,it2,4)>0.001d0) then
                                dum2=dum2+1d0
                            end if
                        end if

                    end do
                end do
                dum2=dum2/dum3
                dum4=dum4/dum3
                dum5=dum5/dum3

                write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5)') it8*1d0,it7*1d0, it6*1d0, (it2+19)*1d0, A(2,it7), dum2, dum4, dum5,0d0

            end do

            do it4=1,Tret

                dum2=0d0
                dum3=0d0
                dum4=0d0
                dum5=0d0
                dum6=0d0

                do it2=1,nsim2
                    do it=1,nsim
                        if(exp1f(it2,it,T,2)==it7) then
                            if(SimR1f(it2,it,it4,5)>1d-3) then
                                dum2=dum2+(1d0)
                            end if
                            dum3=dum3+1d0
                            dum4=dum4+SimR1f(it2,it,it4,1)
                            dum5=dum5+SimR1f(it2,it,it4,2)
                            dum6=dum6+SimR1f(it2,it,it4,13)
                        end if
                    end do
                end do


                dum2=dum2/dum3
                dum4=dum4/dum3
                dum5=dum5/dum3
                dum6=dum6/dum3

                write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5)') it8*1d0, it7*1d0, it6*1d0, (64+it4)*1d0, A(2,it7), dum2, dum4, dum5, dum6

            end do

        end do


        close(1)


        !Single Males

        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppath_males_a.txt')

        do it7=1,na


            do it2=1,T


                dum2=0.0d0
                dum3=0.0d0
                dum4=0.0d0
                dum5=0.0d0

                do it4=1,nsim2
                    do it3=1,10000


                        if(exp1m(it4,it3,it2,2)==it7) then

                            dum3=dum3+1d0
                            dum4=dum4+Sim1m(it4,it3,it2,1)
                            dum5=dum5+Sim1m(it4,it3,it2,2)

                            if(Sim1m(it4,it3,it2,4)>0.001d0) then
                                dum2=dum2+1d0
                            end if
                        end if

                    end do
                end do
                dum2=dum2/dum3
                dum4=dum4/dum3
                dum5=dum5/dum3

                write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5)') it8*1d0, it7*1d0, it6*1d0, (it2+19)*1d0, A(1,it7), dum2, dum4, dum5, 0d0

            end do

            do it4=1,Tret

                dum2=0d0
                dum3=0d0
                dum4=0d0
                dum5=0d0
                dum6=0d0

                do it2=1,nsim2
                    do it=1,nsim
                        if(exp1m(it2,it,T,2)==it7) then
                            if(SimR1m(it2,it,it4,5)>1d-3) then
                                dum2=dum2+(1d0)
                            end if
                            dum3=dum3+1d0
                            dum4=dum4+SimR1m(it2,it,it4,1)
                            dum5=dum5+SimR1m(it2,it,it4,2)
                            dum6=dum6+SimR1m(it2,it,it4,13)
                        end if
                    end do
                end do



                dum2=dum2/dum3
                dum4=dum4/dum3
                dum5=dum5/dum3
                dum6=dum6/dum3

                write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5,F12.5)') it8*1d0, it7*1d0, it6*1d0, (64+it4)*1d0, A(1,it7), dum2, dum4, dum5, dum6

            end do

        end do

        close(1)


        !Earners in couples

        open(1,file=trim(results_folder)//trim(tax_folder)//'earners_in_couples.txt')

        open(2,file=trim(results_folder)//trim(tax_folder)//'ability_couple_f_works.txt')

        do it2=1,T


            dum2=0.0d0
            dum3=0.0d0
            dum4=0.0d0
            dum5=0.0d0
            dum7=0d0
            a_couple=0d0

            do it4=1,nsim2
                do it3=1,nsim

                    if(Sim1m(it4,it3,it2,10)>0.5) then

                        it7=exp1m(it4,it3,it2,4)
                        dum3=dum3+1d0

                        if((Sim1m(it4,it3,it2,4)>0.001).AND.(Sim1f(it4,it7,it2,4)>0.001d0)) then
                            dum2=dum2+1d0
                        elseif((Sim1m(it4,it3,it2,4)>0.001).AND.(Sim1f(it4,it7,it2,4)<0.001d0)) then
                            dum4=dum4+1d0
                        elseif((Sim1m(it4,it3,it2,4)<0.001).AND.(Sim1f(it4,it7,it2,4)>0.001d0)) then
                            dum5=dum5+1d0
                            it8=exp1m(it4,it3,it2,2)
                            it9=exp1f(it4,it7,it2,2)
                            a_couple(it9,it8)=a_couple(it9,it8)+1d0
                        else
                            dum7=dum7+1d0
                        end if
                    end if

                end do
            end do

            a_couple=a_couple/dum5
            dum2=dum2/dum3
            dum4=dum4/dum3
            dum5=dum5/dum3
            dum7=dum7/dum3

            write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5)') (it2+19)*1d0, dum2, dum4, dum5, dum7
            write (2,'(26F12.5)') (it2+19)*1d0, a_couple(1,1), a_couple(1,2), a_couple(1,3), a_couple(1,4), a_couple(1,5), a_couple(2,1), a_couple(2,2), a_couple(2,3), a_couple(2,4), a_couple(2,5), a_couple(3,1), a_couple(3,2), a_couple(3,3), a_couple(3,4), a_couple(3,5), a_couple(4,1), a_couple(4,2), a_couple(4,3), a_couple(4,4), a_couple(4,5), a_couple(5,1), a_couple(5,2), a_couple(5,3), a_couple(5,4), a_couple(5,5)


        end do

        do it2=1,Tret

            dum2=0.0d0
            dum3=0.0d0
            dum4=0.0d0
            dum5=0.0d0
            dum7=0d0
            a_couple=0d0

            do it4=1,nsim2
                do it3=1,nsim

                    if(Sim1m(it4,it3,T,10)>0.5) then

                        it7=exp1m(it4,it3,T,4)
                        dum3=dum3+1d0

                        if((SimR1m(it4,it3,it2,5)>0.001).AND.(SimR1f(it4,it7,it2,5)>0.001d0)) then
                            dum2=dum2+1d0
                        elseif((SimR1m(it4,it3,it2,5)>0.001).AND.(SimR1f(it4,it7,it2,5)<0.001d0)) then
                            dum4=dum4+1d0
                        elseif((SimR1m(it4,it3,it2,5)<0.001).AND.(SimR1f(it4,it7,it2,5)>0.001d0)) then
                            dum5=dum5+1d0
                            it8=exp1m(it4,it3,T,2)
                            it9=exp1f(it4,it7,T,2)
                            a_couple(it9,it8)=a_couple(it9,it8)+1d0
                        else
                            dum7=dum7+1d0
                        end if
                    end if

                end do
            end do

            a_couple=a_couple/dum5
            dum2=dum2/dum3
            dum4=dum4/dum3
            dum5=dum5/dum3
            dum7=dum7/dum3

            write (1,'(F12.5,F12.5,F12.5,F12.5,F12.5)') (it2+64)*1d0, dum2, dum4, dum5, dum7
            write (2,'(26F12.5)') (it2+64)*1d0, a_couple(1,1), a_couple(1,2), a_couple(1,3), a_couple(1,4), a_couple(1,5), a_couple(2,1), a_couple(2,2), a_couple(2,3), a_couple(2,4), a_couple(2,5), a_couple(3,1), a_couple(3,2), a_couple(3,3), a_couple(3,4), a_couple(3,5), a_couple(4,1), a_couple(4,2), a_couple(4,3), a_couple(4,4), a_couple(4,5), a_couple(5,1), a_couple(5,2), a_couple(5,3), a_couple(5,4), a_couple(5,5)

        end do


        close(1)

        close(2)   

        !Persistence of Employment

        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppathsingle2.txt')


        do it2=4,T
            dum2=0.0d0
            dum3=0.0d0
            dum4=0.0d0
            dum5=0.0d0

            do it4=1,nsim2
                do it3=1,10000
                    if(Sim1f(it4,it3,it2,10)<0.5d0) then
                        if(Sim1f(it4,it3,it2,4)>0.001d0) then
                            dum3=dum3+1d0
                            if(Sim1f(it4,it3,it2-1,4)>0.001d0) then    
                                dum2=dum2+1d0
                            end if
                            if((Sim1f(it4,it3,it2-1,4)>0.001d0).AND.(Sim1f(it4,it3,it2-2,4)>0.001d0)) then    
                                dum4=dum4+1d0
                            end if
                            if((Sim1f(it4,it3,it2-1,4)>0.001d0).AND.(Sim1f(it4,it3,it2-2,4)>0.001d0).AND.(Sim1f(it4,it3,it2-3,4)>0.001d0)) then    
                                dum5=dum5+1d0
                            end if
                        end if
                    end if
                end do
            end do
            dum2=dum2/dum3
            dum4=dum4/dum3
            dum5=dum5/dum3
            write (1,'(F8.5,F8.5,F8.5,F8.5)') (19+it2)*1d0, dum2, dum4, dum5
        end do

        !Single female labor force participation by age after 65

        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0

        it4=1

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1f(it2,it,T,10)<0.5) then
                    if(SimR1f(it2,it,it4,5)>1d-3) then
                        dum3=dum3+1d0*WeightRet(it4)
                        if(Sim1f(it2,it,T,4)>1d-3) then
                            dum2=dum2+(1d0)*WeightRet(it4)
                        end if
                        if((Sim1f(it2,it,T,4)>1d-3).AND.(Sim1f(it2,it,T-1,4)>1d-3)) then
                            dum4=dum4+(1d0)*WeightRet(it4)
                        end if
                        if((Sim1f(it2,it,T,4)>1d-3).AND.(Sim1f(it2,it,T-1,4)>1d-3).AND.(Sim1f(it2,it,T-2,4)>1d-3)) then
                            dum5=dum5+(1d0)*WeightRet(it4)
                        end if   
                    end if
                end if
            end do
        end do



        dum2=dum2/dum3
        dum4=dum4/dum3
        dum5=dum5/dum3

        write (1,'(F8.5,F8.5,F8.5,F8.5)') (64+it4)*1d0, dum2, dum4, dum5

        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0

        it4=2

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1f(it2,it,T,10)<0.5) then
                    if(SimR1f(it2,it,it4,5)>1d-3) then
                        dum3=dum3+1d0*WeightRet(it4)
                        if(SimR1f(it2,it,it4-1,5)>1d-3) then
                            dum2=dum2+(1d0)*WeightRet(it4)
                        end if
                        if((SimR1f(it2,it,it4-1,5)>1d-3).AND.(Sim1f(it2,it,T,4)>1d-3)) then
                            dum4=dum4+(1d0)*WeightRet(it4)
                        end if
                        if((SimR1f(it2,it,it4-1,5)>1d-3).AND.(Sim1f(it2,it,T,4)>1d-3).AND.(Sim1f(it2,it,T-1,4)>1d-3)) then
                            dum5=dum5+(1d0)*WeightRet(it4)
                        end if   
                    end if
                end if
            end do
        end do



        dum2=dum2/dum3
        dum4=dum4/dum3
        dum5=dum5/dum3

        write (1,'(F8.5,F8.5,F8.5,F8.5)') (64+it4)*1d0, dum2, dum4, dum5


        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0

        it4=3

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1f(it2,it,T,10)<0.5) then
                    if(SimR1f(it2,it,it4,5)>1d-3) then
                        dum3=dum3+1d0*WeightRet(it4)
                        if(SimR1f(it2,it,it4-1,5)>1d-3) then
                            dum2=dum2+(1d0)*WeightRet(it4)
                        end if
                        if((SimR1f(it2,it,it4-1,5)>1d-3).AND.(SimR1f(it2,it,it4-2,5)>1d-3)) then
                            dum4=dum4+(1d0)*WeightRet(it4)
                        end if
                        if((SimR1f(it2,it,it4-1,5)>1d-3).AND.(SimR1f(it2,it,it4-2,5)>1d-3).AND.(Sim1f(it2,it,T,4)>1d-3)) then
                            dum5=dum5+(1d0)*WeightRet(it4)
                        end if   
                    end if
                end if
            end do
        end do



        dum2=dum2/dum3
        dum4=dum4/dum3
        dum5=dum5/dum3

        write (1,'(F8.5,F8.5,F8.5,F8.5)') (64+it4)*1d0, dum2, dum4, dum5



        do it4=4,Tret

            dum2=0d0
            dum3=0d0
            dum4=0d0
            dum5=0d0

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,T,10)<0.5) then
                        if(SimR1f(it2,it,it4,5)>1d-3) then
                            dum3=dum3+1d0*WeightRet(it4)
                            if(SimR1f(it2,it,it4-1,5)>1d-3) then
                                dum2=dum2+(1d0)*WeightRet(it4)
                            end if
                            if((SimR1f(it2,it,it4-1,5)>1d-3).AND.(SimR1f(it2,it,it4-2,5)>1d-3)) then
                                dum4=dum4+(1d0)*WeightRet(it4)
                            end if
                            if((SimR1f(it2,it,it4-1,5)>1d-3).AND.(SimR1f(it2,it,it4-2,5)>1d-3).AND.(SimR1f(it2,it,it4-3,5)>1d-3)) then
                                dum5=dum5+(1d0)*WeightRet(it4)
                            end if   
                        end if
                    end if
                end do
            end do



            dum2=dum2/dum3
            dum4=dum4/dum3
            dum5=dum5/dum3

            write (1,'(F8.5,F8.5,F8.5,F8.5)') (64+it4)*1d0, dum2, dum4, dum5

        end do



        close(1)

        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppathmarried2.txt')


        do it2=4,T
            dum2=0.0d0
            dum3=0.0d0
            dum4=0.0d0
            dum5=0.0d0

            do it4=1,nsim2
                do it3=1,10000
                    if(Sim1f(it4,it3,it2,10)>0.5d0) then
                        if(Sim1f(it4,it3,it2,4)>0.001d0) then
                            dum3=dum3+1d0
                            if(Sim1f(it4,it3,it2-1,4)>0.001d0) then    
                                dum2=dum2+1d0
                            end if
                            if((Sim1f(it4,it3,it2-1,4)>0.001d0).AND.(Sim1f(it4,it3,it2-2,4)>0.001d0)) then    
                                dum4=dum4+1d0
                            end if
                            if((Sim1f(it4,it3,it2-1,4)>0.001d0).AND.(Sim1f(it4,it3,it2-2,4)>0.001d0).AND.(Sim1f(it4,it3,it2-3,4)>0.001d0)) then    
                                dum5=dum5+1d0
                            end if
                        end if
                    end if
                end do
            end do
            dum2=dum2/dum3
            dum4=dum4/dum3
            dum5=dum5/dum3
            write (1,'(F8.5,F8.5,F8.5,F8.5)') (19+it2)*1d0, dum2, dum4, dum5
        end do

        !Married female labor force participation by age after 65

        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0

        it4=1

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1f(it2,it,T,10)>0.5) then
                    if(SimR1f(it2,it,it4,5)>1d-3) then
                        dum3=dum3+1d0*WeightRet(it4)
                        if(Sim1f(it2,it,T,4)>1d-3) then
                            dum2=dum2+(1d0)*WeightRet(it4)
                        end if
                        if((Sim1f(it2,it,T,4)>1d-3).AND.(Sim1f(it2,it,T-1,4)>1d-3)) then
                            dum4=dum4+(1d0)*WeightRet(it4)
                        end if
                        if((Sim1f(it2,it,T,4)>1d-3).AND.(Sim1f(it2,it,T-1,4)>1d-3).AND.(Sim1f(it2,it,T-2,4)>1d-3)) then
                            dum5=dum5+(1d0)*WeightRet(it4)
                        end if   
                    end if
                end if
            end do
        end do



        dum2=dum2/dum3
        dum4=dum4/dum3
        dum5=dum5/dum3

        write (1,'(F8.5,F8.5,F8.5,F8.5)') (64+it4)*1d0, dum2, dum4, dum5

        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0

        it4=2

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1f(it2,it,T,10)>0.5) then
                    if(SimR1f(it2,it,it4,5)>1d-3) then
                        dum3=dum3+1d0*WeightRet(it4)
                        if(SimR1f(it2,it,it4-1,5)>1d-3) then
                            dum2=dum2+(1d0)*WeightRet(it4)
                        end if
                        if((SimR1f(it2,it,it4-1,5)>1d-3).AND.(Sim1f(it2,it,T,4)>1d-3)) then
                            dum4=dum4+(1d0)*WeightRet(it4)
                        end if
                        if((SimR1f(it2,it,it4-1,5)>1d-3).AND.(Sim1f(it2,it,T,4)>1d-3).AND.(Sim1f(it2,it,T-1,4)>1d-3)) then
                            dum5=dum5+(1d0)*WeightRet(it4)
                        end if   
                    end if
                end if
            end do
        end do



        dum2=dum2/dum3
        dum4=dum4/dum3
        dum5=dum5/dum3

        write (1,'(F8.5,F8.5,F8.5,F8.5)') (64+it4)*1d0, dum2, dum4, dum5


        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0

        it4=3

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1f(it2,it,T,10)>0.5) then
                    if(SimR1f(it2,it,it4,5)>1d-3) then
                        dum3=dum3+1d0*WeightRet(it4)
                        if(SimR1f(it2,it,it4-1,5)>1d-3) then
                            dum2=dum2+(1d0)*WeightRet(it4)
                        end if
                        if((SimR1f(it2,it,it4-1,5)>1d-3).AND.(SimR1f(it2,it,it4-2,5)>1d-3)) then
                            dum4=dum4+(1d0)*WeightRet(it4)
                        end if
                        if((SimR1f(it2,it,it4-1,5)>1d-3).AND.(SimR1f(it2,it,it4-2,5)>1d-3).AND.(Sim1f(it2,it,T,4)>1d-3)) then
                            dum5=dum5+(1d0)*WeightRet(it4)
                        end if   
                    end if
                end if
            end do
        end do



        dum2=dum2/dum3
        dum4=dum4/dum3
        dum5=dum5/dum3

        write (1,'(F8.5,F8.5,F8.5,F8.5)') (64+it4)*1d0, dum2, dum4, dum5



        do it4=4,Tret

            dum2=0d0
            dum3=0d0
            dum4=0d0
            dum5=0d0

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,T,10)>0.5) then
                        if(SimR1f(it2,it,it4,5)>1d-3) then
                            dum3=dum3+1d0*WeightRet(it4)
                            if(SimR1f(it2,it,it4-1,5)>1d-3) then
                                dum2=dum2+(1d0)*WeightRet(it4)
                            end if
                            if((SimR1f(it2,it,it4-1,5)>1d-3).AND.(SimR1f(it2,it,it4-2,5)>1d-3)) then
                                dum4=dum4+(1d0)*WeightRet(it4)
                            end if
                            if((SimR1f(it2,it,it4-1,5)>1d-3).AND.(SimR1f(it2,it,it4-2,5)>1d-3).AND.(SimR1f(it2,it,it4-3,5)>1d-3)) then
                                dum5=dum5+(1d0)*WeightRet(it4)
                            end if   
                        end if
                    end if
                end do
            end do



            dum2=dum2/dum3
            dum4=dum4/dum3
            dum5=dum5/dum3

            write (1,'(F8.5,F8.5,F8.5,F8.5)') (64+it4)*1d0, dum2, dum4, dum5

        end do



        close(1) 






        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppathsingle_male2.txt')

        do it2=4,T
            dum2=0.0d0
            dum3=0.0d0
            dum4=0.0d0
            dum5=0.0d0

            do it4=1,nsim2
                do it3=1,10000
                    if(Sim1m(it4,it3,it2,10)<0.5d0) then
                        if(Sim1m(it4,it3,it2,4)>0.001d0) then
                            dum3=dum3+1d0
                            if(Sim1m(it4,it3,it2-1,4)>0.001d0) then    
                                dum2=dum2+1d0
                            end if
                            if((Sim1m(it4,it3,it2-1,4)>0.001d0).AND.(Sim1m(it4,it3,it2-2,4)>0.001d0)) then    
                                dum4=dum4+1d0
                            end if
                            if((Sim1m(it4,it3,it2-1,4)>0.001d0).AND.(Sim1m(it4,it3,it2-2,4)>0.001d0).AND.(Sim1m(it4,it3,it2-3,4)>0.001d0)) then    
                                dum5=dum5+1d0
                            end if
                        end if
                    end if
                end do
            end do
            dum2=dum2/dum3
            dum4=dum4/dum3
            dum5=dum5/dum3
            write (1,'(F8.5,F8.5,F8.5,F8.5)') (19+it2)*1d0, dum2, dum4, dum5
        end do

        !Single male labor force participation by age after 65

        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0

        it4=1

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1m(it2,it,T,10)<0.5) then
                    if(SimR1m(it2,it,it4,5)>1d-3) then
                        dum3=dum3+1d0*WeightRet(it4)
                        if(Sim1m(it2,it,T,4)>1d-3) then
                            dum2=dum2+(1d0)*WeightRet(it4)
                        end if
                        if((Sim1m(it2,it,T,4)>1d-3).AND.(Sim1m(it2,it,T-1,4)>1d-3)) then
                            dum4=dum4+(1d0)*WeightRet(it4)
                        end if
                        if((Sim1m(it2,it,T,4)>1d-3).AND.(Sim1m(it2,it,T-1,4)>1d-3).AND.(Sim1m(it2,it,T-2,4)>1d-3)) then
                            dum5=dum5+(1d0)*WeightRet(it4)
                        end if   
                    end if
                end if
            end do
        end do



        dum2=dum2/dum3
        dum4=dum4/dum3
        dum5=dum5/dum3

        write (1,'(F8.5,F8.5,F8.5,F8.5)') (64+it4)*1d0, dum2, dum4, dum5

        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0

        it4=2

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1m(it2,it,T,10)<0.5) then
                    if(SimR1m(it2,it,it4,5)>1d-3) then
                        dum3=dum3+1d0*WeightRet(it4)
                        if(SimR1m(it2,it,it4-1,5)>1d-3) then
                            dum2=dum2+(1d0)*WeightRet(it4)
                        end if
                        if((SimR1m(it2,it,it4-1,5)>1d-3).AND.(Sim1m(it2,it,T,4)>1d-3)) then
                            dum4=dum4+(1d0)*WeightRet(it4)
                        end if
                        if((SimR1m(it2,it,it4-1,5)>1d-3).AND.(Sim1m(it2,it,T,4)>1d-3).AND.(Sim1m(it2,it,T-1,4)>1d-3)) then
                            dum5=dum5+(1d0)*WeightRet(it4)
                        end if   
                    end if
                end if
            end do
        end do



        dum2=dum2/dum3
        dum4=dum4/dum3
        dum5=dum5/dum3

        write (1,'(F8.5,F8.5,F8.5,F8.5)') (64+it4)*1d0, dum2, dum4, dum5


        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0

        it4=3

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1m(it2,it,T,10)<0.5) then
                    if(SimR1m(it2,it,it4,5)>1d-3) then
                        dum3=dum3+1d0*WeightRet(it4)
                        if(SimR1m(it2,it,it4-1,5)>1d-3) then
                            dum2=dum2+(1d0)*WeightRet(it4)
                        end if
                        if((SimR1m(it2,it,it4-1,5)>1d-3).AND.(SimR1m(it2,it,it4-2,5)>1d-3)) then
                            dum4=dum4+(1d0)*WeightRet(it4)
                        end if
                        if((SimR1m(it2,it,it4-1,5)>1d-3).AND.(SimR1m(it2,it,it4-2,5)>1d-3).AND.(Sim1m(it2,it,T,4)>1d-3)) then
                            dum5=dum5+(1d0)*WeightRet(it4)
                        end if   
                    end if
                end if
            end do
        end do



        dum2=dum2/dum3
        dum4=dum4/dum3
        dum5=dum5/dum3

        write (1,'(F8.5,F8.5,F8.5,F8.5)') (64+it4)*1d0, dum2, dum4, dum5



        do it4=4,Tret

            dum2=0d0
            dum3=0d0
            dum4=0d0
            dum5=0d0

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,T,10)<0.5) then
                        if(SimR1m(it2,it,it4,5)>1d-3) then
                            dum3=dum3+1d0*WeightRet(it4)
                            if(SimR1m(it2,it,it4-1,5)>1d-3) then
                                dum2=dum2+(1d0)*WeightRet(it4)
                            end if
                            if((SimR1m(it2,it,it4-1,5)>1d-3).AND.(SimR1m(it2,it,it4-2,5)>1d-3)) then
                                dum4=dum4+(1d0)*WeightRet(it4)
                            end if
                            if((SimR1m(it2,it,it4-1,5)>1d-3).AND.(SimR1m(it2,it,it4-2,5)>1d-3).AND.(SimR1m(it2,it,it4-3,5)>1d-3)) then
                                dum5=dum5+(1d0)*WeightRet(it4)
                            end if   
                        end if
                    end if
                end do
            end do



            dum2=dum2/dum3
            dum4=dum4/dum3
            dum5=dum5/dum3

            write (1,'(F8.5,F8.5,F8.5,F8.5)') (64+it4)*1d0, dum2, dum4, dum5

        end do


        close(1)

        open(1,file=trim(results_folder)//trim(tax_folder)//'lfppathmarried_male2.txt')

        do it2=4,T
            dum2=0.0d0
            dum3=0.0d0
            dum4=0.0d0
            dum5=0.0d0

            do it4=1,nsim2
                do it3=1,10000
                    if(Sim1m(it4,it3,it2,10)>0.5d0) then
                        if(Sim1m(it4,it3,it2,4)>0.001d0) then
                            dum3=dum3+1d0
                            if(Sim1m(it4,it3,it2-1,4)>0.001d0) then    
                                dum2=dum2+1d0
                            end if
                            if((Sim1m(it4,it3,it2-1,4)>0.001d0).AND.(Sim1m(it4,it3,it2-2,4)>0.001d0)) then    
                                dum4=dum4+1d0
                            end if
                            if((Sim1m(it4,it3,it2-1,4)>0.001d0).AND.(Sim1m(it4,it3,it2-2,4)>0.001d0).AND.(Sim1m(it4,it3,it2-3,4)>0.001d0)) then    
                                dum5=dum5+1d0
                            end if
                        end if
                    end if
                end do
            end do
            dum2=dum2/dum3
            dum4=dum4/dum3
            dum5=dum5/dum3
            write (1,'(F8.5,F8.5,F8.5,F8.5)') (19+it2)*1d0, dum2, dum4, dum5
        end do

        !Married male labor force participation by age after 65

        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0

        it4=1

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1m(it2,it,T,10)>0.5) then
                    if(SimR1m(it2,it,it4,5)>1d-3) then
                        dum3=dum3+1d0*WeightRet(it4)
                        if(Sim1m(it2,it,T,4)>1d-3) then
                            dum2=dum2+(1d0)*WeightRet(it4)
                        end if
                        if((Sim1m(it2,it,T,4)>1d-3).AND.(Sim1m(it2,it,T-1,4)>1d-3)) then
                            dum4=dum4+(1d0)*WeightRet(it4)
                        end if
                        if((Sim1m(it2,it,T,4)>1d-3).AND.(Sim1m(it2,it,T-1,4)>1d-3).AND.(Sim1m(it2,it,T-2,4)>1d-3)) then
                            dum5=dum5+(1d0)*WeightRet(it4)
                        end if   
                    end if
                end if
            end do
        end do



        dum2=dum2/dum3
        dum4=dum4/dum3
        dum5=dum5/dum3

        write (1,'(F8.5,F8.5,F8.5,F8.5)') (64+it4)*1d0, dum2, dum4, dum5

        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0

        it4=2

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1m(it2,it,T,10)>0.5) then
                    if(SimR1m(it2,it,it4,5)>1d-3) then
                        dum3=dum3+1d0*WeightRet(it4)
                        if(SimR1m(it2,it,it4-1,5)>1d-3) then
                            dum2=dum2+(1d0)*WeightRet(it4)
                        end if
                        if((SimR1m(it2,it,it4-1,5)>1d-3).AND.(Sim1m(it2,it,T,4)>1d-3)) then
                            dum4=dum4+(1d0)*WeightRet(it4)
                        end if
                        if((SimR1m(it2,it,it4-1,5)>1d-3).AND.(Sim1m(it2,it,T,4)>1d-3).AND.(Sim1m(it2,it,T-1,4)>1d-3)) then
                            dum5=dum5+(1d0)*WeightRet(it4)
                        end if   
                    end if
                end if
            end do
        end do



        dum2=dum2/dum3
        dum4=dum4/dum3
        dum5=dum5/dum3

        write (1,'(F8.5,F8.5,F8.5,F8.5)') (64+it4)*1d0, dum2, dum4, dum5


        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0

        it4=3

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1m(it2,it,T,10)>0.5) then
                    if(SimR1m(it2,it,it4,5)>1d-3) then
                        dum3=dum3+1d0*WeightRet(it4)
                        if(SimR1m(it2,it,it4-1,5)>1d-3) then
                            dum2=dum2+(1d0)*WeightRet(it4)
                        end if
                        if((SimR1m(it2,it,it4-1,5)>1d-3).AND.(SimR1m(it2,it,it4-2,5)>1d-3)) then
                            dum4=dum4+(1d0)*WeightRet(it4)
                        end if
                        if((SimR1m(it2,it,it4-1,5)>1d-3).AND.(SimR1m(it2,it,it4-2,5)>1d-3).AND.(Sim1m(it2,it,T,4)>1d-3)) then
                            dum5=dum5+(1d0)*WeightRet(it4)
                        end if   
                    end if
                end if
            end do
        end do



        dum2=dum2/dum3
        dum4=dum4/dum3
        dum5=dum5/dum3

        write (1,'(F8.5,F8.5,F8.5,F8.5)') (64+it4)*1d0, dum2, dum4, dum5



        do it4=4,Tret

            dum2=0d0
            dum3=0d0
            dum4=0d0
            dum5=0d0

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,T,10)>0.5) then
                        if(SimR1m(it2,it,it4,5)>1d-3) then
                            dum3=dum3+1d0*WeightRet(it4)
                            if(SimR1m(it2,it,it4-1,5)>1d-3) then
                                dum2=dum2+(1d0)*WeightRet(it4)
                            end if
                            if((SimR1m(it2,it,it4-1,5)>1d-3).AND.(SimR1m(it2,it,it4-2,5)>1d-3)) then
                                dum4=dum4+(1d0)*WeightRet(it4)
                            end if
                            if((SimR1m(it2,it,it4-1,5)>1d-3).AND.(SimR1m(it2,it,it4-2,5)>1d-3).AND.(SimR1m(it2,it,it4-3,5)>1d-3)) then
                                dum5=dum5+(1d0)*WeightRet(it4)
                            end if   
                        end if
                    end if
                end do
            end do



            dum2=dum2/dum3
            dum4=dum4/dum3
            dum5=dum5/dum3

            write (1,'(F8.5,F8.5,F8.5,F8.5)') (64+it4)*1d0, dum2, dum4, dum5

        end do



        close(1)

        !open(1, file='Probm.txt')
        !do it2=1,T
        !    write (1,'(F12.4,F12.4)') (it2+19)*1d0, probm(it2)
        !end do
        !close(1)
        !
        !open(1, file='Probd.txt')
        !do it2=1,T
        !    write (1,'(F12.4,F12.4)') (it2+19)*1d0, probd(it2)
        !end do
        !close(1)
        
        if (gen_earnings_distr == 1) then
            !Variables are age, ID number, earnings
            !open(1,file=trim(results_folder)//trim(tax_folder)//'lfppathmarried_male2.txt')
            open(1, file=trim(results_folder)//trim(tax_folder)//'Simulated_earnings_single_male.txt')

            do it2=1,T
                do it3=1,nsim2
                    do it4=1,nsim
                        if(Sim1m(it4,it3,it2,10)<0.5d0) then
                            it7=exp1m(it3,it4,it2,2)
                            !if(it7==1) then
                                write (1,'(F12.4,F12.4,F12.4)') (it2+19)*1d0, nsim*(it3-1)*1d0+it4*1d0, Sim1m(it3,it4,it2,5)
                            !end if
                        end if
                    end do
                end do
            end do
    
            close(1)
    
            !Variables are age, ID number, earnings
    
            open(1, file=trim(results_folder)//trim(tax_folder)//'Simulated_earnings_single_female.txt')

            do it2=1,T
                do it3=1,nsim2
                    do it4=1,nsim
                        if(Sim1f(it4,it3,it2,10)<0.5d0) then
                            it7=exp1f(it3,it4,it2,2)
                            !if(it7==1) then
                                write (1,'(F12.4,F12.4,F12.4)') (it2+19)*1d0, nsim*(it3-1)*1d0+it4*1d0, Sim1f(it3,it4,it2,5)
                            !end if
                        end if
                    end do
                end do
            end do
    
            close(1)
    
            !Variables are age, male ID number, household earnings
    
            open(1, file=trim(results_folder)//trim(tax_folder)//'Simulated_earnings_married.txt')

            do it2=1,T
                do it3=1,nsim2
                    do it4=1,nsim
                        if(Sim1m(it4,it3,it2,10)>0.5d0) then
                            it7=exp1m(it3,it4,it2,4)                
                            if((exp1m(it3,it4,it2,2)==1).AND.(exp1f(it3,it7,it2,2)==1)) then
                                write (1,'(F12.4,F12.4,F12.4)') (it2+19)*1d0, nsim*(it3-1)*1d0+it4*1d0, Sim1m(it3,it4,it2,6)
                            end if
                        end if
                    end do
                end do
            end do
    
            close(1)    
            
        end if


    end subroutine generate_paths

end program Laffer
