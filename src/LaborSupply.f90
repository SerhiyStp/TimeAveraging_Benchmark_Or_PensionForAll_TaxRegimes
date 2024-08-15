module LaborSupply
    
    implicit none
    
    real(8) :: chi1, chi2, cmod, eta1, eta2, wage1, wage2
    
contains

    subroutine Solve_Hours(tax_regime)
        use Model_Parameters, only: nc, Deduct_cutoff, Deduct_cutoff_Mar, yhat, yhat_mar
        integer, intent(in) :: tax_regime
        integer :: ik
        integer :: tax_regime_nit
        
        if (tax_regime == 1) then
            !$OMP PARALLEL PRIVATE(ik)
            !$OMP DO SCHEDULE(DYNAMIC)
            do ik=1,nc
                call lsupply_baseline(ik)
            end do
            !$OMP END DO    
            !$OMP END PARALLEL  
            
        else if (tax_regime == 2) then
            if (Deduct_cutoff <= yhat .and. Deduct_cutoff_Mar <= yhat_mar) then
                tax_regime_nit = 1
            else if (Deduct_cutoff > yhat .and. Deduct_cutoff_Mar > yhat_mar) then
                tax_regime_nit = 2
            else
                print *, 'ERROR with tax regime'
                stop
            end if 
            if (tax_regime_nit == 1) then
                !$OMP PARALLEL PRIVATE(ik)
                !$OMP DO SCHEDULE(DYNAMIC)
                do ik=1,nc
                    call lsupply_nit_regime1(ik)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL
            else
                !$OMP PARALLEL PRIVATE(ik)
                !$OMP DO SCHEDULE(DYNAMIC)
                do ik=1,nc
                    call lsupply_nit_regime2(ik)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL
            end if            
        else
            print *, 'WARNING: wrong tax regime'
        end if
        
    end subroutine Solve_Hours

    
    subroutine lsupply_baseline(ik)
        !This subroutine computes optimal policies durin active worklife.
        use Model_Parameters
        use PolicyFunctions
        use glob0
        use hours
        use root_module

        implicit none

        integer, INTENT(IN) :: ik
        integer :: ium, iuf

        real(8) :: test_h_lo, h_cutoff, test_h_cutoff1, test_h_cutoff2, test_h_hi1, test_h_hi2, h_lo, h_hi, h_sol
        real(8) :: h_a, h_b, test_a, test_b
        integer :: i_gend
        real(8) :: test_hm_lo, test_hf_lo, test_hm_hi, test_hf_hi, hf_tmp
        real(8) :: hm_lo, hf_lo
        real(8) :: hm_hi, hf_hi
        real(8) :: hm_d, hf_d
        real(8) :: hm_sol, hf_sol
        real(8) :: test_h_d   

        real(8) :: hsol, fsol_test
        integer :: ic


        !CALL ERSET(IERSVR, IPACT, ISACT)

        dum=c_grid(ik)
        cmod = dum

        ! Married couples:
        do ium=1,nw
            do iuf=1,nw

                wagem=wage_grid(ium)
                wagef=wage_grid(iuf)

                hm_lo = 1d-9
                hf_lo = hours_ratio(hm_lo, wagef, wagem, chif, chim, etaf, etam)
                wage1 = wagem
                chi1 = chim
                eta1 = etam
                test_hm_lo = foc_hrs_mar_below_deduct(hm_lo)
                if (test_hm_lo > 0d0) then
                    ! Corner solution with hm = 0, hf = 0
                    hm_sol = hm_lo
                    hf_sol = hf_lo
                else

                    ! Check if there is hm < 1 or hf < 1 such that wm*hm + wf*hf <= Deduct
                    hm_hi = 1d0
                    hf_hi = hours_ratio(hm_hi, wagef, wagem, chif, chim, etaf, etam)


                    if (hm_hi > hf_hi) then
                        i_gend = 2
                        wage1 = wagef
                        wage2 = wagem
                        chi1 = chif
                        chi2 = chim
                        eta1 = etaf
                        eta2 = etam

                        h_a = 1d-9
                        h_b = 1d0
                        test_a = mar_hours_deduct(h_a)
                        if (test_a > 0d0) then
                            hf_d = h_a
                            hm_d = min(1.0d0, hours_ratio(hf_d,wagem,wagef,chim,chif,etam,etaf))
                        else
                            test_b = mar_hours_deduct(h_b)
                            if (test_a*test_b < 0d0) then
                                ! call zbren(mar_hours_deduct, h_a, h_b) 

                                call zeroin(mar_hours_deduct, h_a, h_b, tol=1d-8, xzero=hsol, fzero=fsol_test, iflag=ic)
                                ! if (abs(h_b-xsol_test) > 1d-5) then
                                !     write(file_id,*) 'warning'
                                ! end if

                                hf_d = hsol !h_b
                                hm_d = min(1.0d0, hours_ratio(hf_d,wagem,wagef,chim,chif,etam,etaf))
                            else
                                hm_d = 100d0
                                hf_d = 100d0
                            end if
                        end if
                    else if (hf_hi > hm_hi) then
                        i_gend = 1
                        wage1 = wagem
                        wage2 = wagef
                        chi1 = chim
                        chi2 = chif
                        eta1 = etam
                        eta2 = etaf

                        h_a = 1d-9
                        h_b = 1d0
                        test_a = mar_hours_deduct(h_a)
                        if (test_a > 0d0) then
                            hm_d = h_a
                            hf_d = min(1.0d0, hours_ratio(hm_d,wagef,wagem,chif,chim,etaf,etam))
                        else
                            test_b = mar_hours_deduct(h_b)
                            if (test_a*test_b < 0d0) then
                                ! call zbren(mar_hours_deduct, h_a, h_b) 

                                call zeroin(mar_hours_deduct, h_a, h_b, tol=1d-8, xzero=hsol, fzero=fsol_test, iflag=ic)
                                ! if (abs(h_b-xsol_test) > 1d-5) then
                                !     write(file_id,*) 'warning'
                                ! end if                            

                                hm_d = hsol !h_b
                                hf_d = min(1.0d0, hours_ratio(hm_d,wagef,wagem,chif,chim,etaf,etam))
                            else
                                hm_d = 100d0
                                hf_d = 100d0
                            end if   
                        end if

                    else

                        print *, 'WARNING: hm_hi == hf_hi'
                        hf_hi = hours_ratio(hm_hi, wagef, wagem, chif, chim, etaf, etam)

                    end if

                    if (hm_d < 100d0 .and. hf_d < 100d0) then
                        ! Look both below and above Deduct 
                        if (i_gend == 1) then
                            wage1 = wagem
                            wage2 = wagef
                            ! chi1 = chim
                            ! chi2 = chif
                            ! eta1 = etam
                            ! eta2 = etaf
                            test_h_cutoff1 = foc_hrs_mar_below_deduct(hm_d)    
                            test_h_cutoff2 = foc_hrs_mar_above_deduct(hm_d + 1d-6)
                            h_cutoff = hm_d
                        else
                            wage1 = wagef
                            wage2 = wagem
                            ! chi1 = chif
                            ! chi2 = chim
                            ! eta1 = etaf
                            ! eta2 = etam
                            test_h_cutoff1 = foc_hrs_mar_below_deduct(hf_d) 
                            test_h_cutoff2 = foc_hrs_mar_above_deduct(hf_d + 1d-6)
                            h_cutoff = hf_d
                        end if

                        if (test_h_cutoff1 > 0d0) then
                            ! h is between 0 and h_cutoff, find using bisection    
                            h_a = 1d-9
                            h_b = h_cutoff
                            test_a = foc_hrs_mar_below_deduct(h_a)
                            test_b = foc_hrs_mar_below_deduct(h_b)
                            if (test_a*test_b < 0d0) then
                                ! call zbren(foc_hrs_mar_below_deduct, h_a, h_b) 

                                call zeroin(foc_hrs_mar_below_deduct, h_a, h_b, tol=1d-8, xzero=hsol, fzero=fsol_test, iflag=ic)
                                ! if (abs(h_b-xsol_test) > 1d-5) then
                                !     write(file_id,*) 'warning'
                                ! end if                            
                                h_b = hsol  
                            else
                                print *, 'warning: f(a) and f(b) have the same sign'
                            end if
                            if (i_gend == 1) then
                                hm_sol = h_b
                                hf_sol = min(1.0d0, hours_ratio(hm_sol, wagef, wagem, chif, chim, etaf, etam))
                            else
                                hf_sol = h_b
                                hm_sol = min(1.0d0, hours_ratio(hf_sol, wagem, wagef, chim, chif, etam, etaf))
                            end if
                        else if (test_h_cutoff2 > 0d0) then
                            ! h = h_d
                            if (i_gend == 1) then
                                hm_sol = hm_d
                                hf_sol = min(1.0d0, hours_ratio(hm_sol, wagef, wagem, chif, chim, etaf, etam))
                            else
                                hf_sol = hf_d
                                hm_sol = min(1.0d0, hours_ratio(hf_sol, wagem, wagef, chim, chif, etam, etaf))
                            end if
                        else
                            ! h_d < h <= 1
                            if (i_gend == 1) then
                                wage1 = wagem
                                wage2 = wagef
                                ! chi1 = chim
                                ! chi2 = chif
                                ! eta1 = etam
                                ! eta2 = etaf
                                hm_hi = 1d0
                                test_hm_hi = foc_hrs_mar_above_deduct(hm_hi)
                                if (test_hm_hi < 0d0) then
                                    hm_sol = 1d0
                                    hf_sol = 1d0
                                else
                                    h_a = hm_d + 1d-6
                                    h_b = 1d0
                                    test_a = foc_hrs_mar_above_deduct(h_a)
                                    test_b = foc_hrs_mar_above_deduct(h_b)  
                                    if (test_a*test_b < 0d0) then
                                        ! call zbren(foc_hrs_mar_above_deduct, h_a, h_b) 

                                        call zeroin(foc_hrs_mar_above_deduct, h_a, h_b, tol=1d-8, xzero=hsol, fzero=fsol_test, iflag=ic)
                                        ! if (abs(h_b-xsol_test) > 1d-5) then
                                        !     write(file_id,*) 'warning'
                                        ! end if
                                        h_b = hsol
                                    else
                                        print *, 'warning: f(a) and f(b) have the same sign'
                                    end if
                                    hm_sol = h_b
                                    hf_sol = min(1.0d0, hours_ratio(hm_sol, wagef, wagem, chif, chim, etaf, etam))
                                end if
                            else
                                wage1 = wagef
                                wage2 = wagem
                                ! chi1 = chif
                                ! chi2 = chim
                                ! eta1 = etaf
                                ! eta2 = etam
                                hf_hi = 1d0
                                test_hm_hi = foc_hrs_mar_above_deduct(hf_hi)
                                if (test_hm_hi < 0d0) then
                                    hm_sol = 1d0
                                    hf_sol = 1d0
                                else
                                    h_a = hf_d + 1d-6
                                    h_b = 1d0
                                    test_a = foc_hrs_mar_above_deduct(h_a)
                                    test_b = foc_hrs_mar_above_deduct(h_b)  
                                    if (test_a*test_b < 0d0) then
                                        ! call zbren(foc_hrs_mar_above_deduct, h_a, h_b) 

                                        call zeroin(foc_hrs_mar_above_deduct, h_a, h_b, tol=1d-8, xzero=hsol, fzero=fsol_test, iflag=ic)
                                        ! if (abs(h_b-xsol_test) > 1d-5) then
                                        !     write(file_id,*) 'warning'
                                        ! end if
                                        h_b = hsol
                                    else
                                        print *, 'warning: f(a) and f(b) have the same sign'
                                    end if
                                    hf_sol = h_b
                                    hm_sol = min(1.0d0, hours_ratio(hf_sol, wagem, wagef, chim, chif, etam, etaf))
                                end if

                            end if
                        end if

                    else 
                        ! Look only below Deduct
                        if (i_gend == 1) then
                            wage1 = wagem
                            wage2 = wagef
                            ! chi1 = chim
                            ! chi2 = chif
                            ! eta1 = etam
                            ! eta2 = etaf
                            hm_hi = 1d0
                            test_hm_hi = foc_hrs_mar_below_deduct(hm_hi)
                            if (test_hm_hi < 0d0) then
                                hm_sol = 1d0
                                hf_sol = 1d0
                            else
                                h_a = 1d-9
                                h_b = 1d0
                                test_a = foc_hrs_mar_below_deduct(h_a)
                                test_b = foc_hrs_mar_below_deduct(h_b)  
                                if (test_a*test_b < 0d0) then
                                    ! call zbren(foc_hrs_mar_below_deduct, h_a, h_b) 

                                    call zeroin(foc_hrs_mar_below_deduct, h_a, h_b, tol=1d-8, xzero=hsol, fzero=fsol_test, iflag=ic)
                                    ! if (abs(h_b-xsol_test) > 1d-5) then
                                    !     write(file_id,*) 'warning'
                                    ! end if
                                    h_b = hsol
                                else
                                    print *, 'warning: f(a) and f(b) have the same sign'
                                end if
                                hm_sol = h_b
                                hf_sol = min(1.0d0, hours_ratio(hm_sol, wagef, wagem, chif, chim, etaf, etam))
                            end if
                        else
                            wage1 = wagef
                            wage2 = wagem
                            ! chi1 = chif
                            ! chi2 = chim
                            ! eta1 = etaf
                            ! eta2 = etam
                            hf_hi = 1d0
                            test_hm_hi = foc_hrs_mar_below_deduct(hf_hi)
                            if (test_hm_hi < 0d0) then
                                hm_sol = 1d0
                                hf_sol = 1d0
                            else
                                h_a = 1d-9
                                h_b = 1d0
                                test_a = foc_hrs_mar_below_deduct(h_a)
                                test_b = foc_hrs_mar_below_deduct(h_b)  
                                if (test_a*test_b < 0d0) then
                                    ! call zbren(foc_hrs_mar_below_deduct, h_a, h_b) 

                                    call zeroin(foc_hrs_mar_below_deduct, h_a, h_b, tol=1d-8, xzero=hsol, fzero=fsol_test, iflag=ic)
                                    ! if (abs(h_b-xsol_test) > 1d-5) then
                                    !     write(file_id,*) 'warning'
                                    ! end if
                                    h_b = hsol
                                else
                                    print *, 'warning: f(a) and f(b) have the same sign'
                                end if
                                hf_sol = h_b
                                hm_sol = min(1.0d0, hours_ratio(hf_sol, wagem, wagef, chim, chif, etam, etaf))
                            end if

                        end if

                    end if

                end if        

                laborm(ik,ium,iuf)=hm_sol
                laborf(ik,ium,iuf)=hf_sol            


            end do
        end do

        ! Only one spouse works:
        do i_gend = 1, 2

            if (i_gend == 1) then
                chi1 = chim
                eta1 = etam
            else
                chi1 = chif
                eta1 = etaf
            end if

            do ium=1,nw
                wagem=wage_grid(ium)
                wage1 = wagem


                h_lo = 1d-9
                test_h_lo = foc_hrs_mar_below_deduct(h_lo)
                h_cutoff = Deduct_Cutoff/wage1

                h_hi = 1.0d0
                test_h_hi1 = foc_hrs_mar_below_deduct(h_hi)

                if (h_cutoff < h_hi) then
                    test_h_cutoff2 = foc_hrs_mar_above_deduct_1w(h_cutoff + 1d-5)
                    test_h_hi2 = foc_hrs_mar_above_deduct_1w(h_hi)
                end if

                if (test_h_lo > 0d0) then
                    ! Uh > w*Uc, so  h = 0d0
                    h_sol = h_lo
                else
                    if (h_cutoff < h_hi) then 
                        if (test_h_cutoff1 > 0d0) then
                            ! h is between 0 and h_cutoff, find using bisection
                            h_a = h_lo
                            h_b = h_cutoff
                            test_a = foc_hrs_mar_below_deduct(h_a)
                            test_b = foc_hrs_mar_below_deduct(h_b)  
                            if (test_a*test_b < 0d0) then
                                ! call zbren(foc_hrs_mar_below_deduct, h_a, h_b)

                                call zeroin(foc_hrs_mar_below_deduct, h_a, h_b, tol=1d-8, xzero=hsol, fzero=fsol_test, iflag=ic)
                                ! if (abs(h_b-xsol_test) > 1d-5) then
                                !     write(file_id,*) 'warning'
                                ! end if
                                h_b = hsol
                            else
                                print *, 'warning: f(a) and f(b) have the same sign'
                            end if
                            h_sol = h_b
                        else if (test_h_cutoff2 > 0d0) then
                            ! h = h_cutoff
                            h_sol = h_cutoff 
                        else if (test_h_hi2 > 0d0) then
                            ! h is between h_cutoff and 1, find using bisection
                            h_a = h_cutoff + 1d-5
                            h_b = h_hi
                            test_a = foc_hrs_mar_above_deduct_1w(h_a)
                            test_b = foc_hrs_mar_above_deduct_1w(h_b)
                            if (test_a*test_b < 0d0) then
                                ! call zbren(foc_hrs_mar_above_deduct_1w, h_a, h_b)

                                call zeroin(foc_hrs_mar_above_deduct_1w, h_a, h_b, tol=1d-8, xzero=hsol, fzero=fsol_test, iflag=ic)
                                ! if (abs(h_b-xsol_test) > 1d-5) then
                                !     write(file_id,*) 'warning'
                                ! end if
                                h_b = hsol
                            else
                                print *, 'warning: f(a) and f(b) have the same sign'
                            end if
                            h_sol = h_b
                        else
                            h_sol = h_hi  
                        end if

                    else
                        if (test_h_hi1 > 0d0) then
                            ! h is between 0 and h_hi, find using bisection
                            h_a = h_lo
                            h_b = h_hi
                            test_a = foc_hrs_mar_below_deduct(h_a)
                            test_b = foc_hrs_mar_below_deduct(h_b)  
                            if (test_a*test_b < 0d0) then
                                ! call zbren(foc_hrs_mar_below_deduct, h_a, h_b)

                                call zeroin(foc_hrs_mar_below_deduct, h_a, h_b, tol=1d-8, xzero=hsol, fzero=fsol_test, iflag=ic)
                                ! if (abs(h_b-xsol_test) > 1d-5) then
                                !     write(file_id,*) 'warning'
                                ! end if
                                h_b = hsol
                            else
                                print *, 'warning: f(a) and f(b) have the same sign'
                            end if
                            h_sol = h_b
                        else
                            h_sol = h_hi
                        end if

                    end if
                end if

                if (i_gend == 1) then
                    labormwork(ik,ium)=h_sol
                else
                    laborfwork(ik,ium)=h_sol
                end if

            end do
        end do


        ! Singles:
        do i_gend = 1, 2

            if (i_gend == 1) then
                chi1 = chims
                eta1 = etam
            else
                chi1 = chifs
                eta1 = etaf
            end if

            do ium=1,nw

                wagem=wage_grid(ium)
                wage1 = wagem



                h_lo = 0.0001d0
                test_h_lo = foc_hrs_single_below_deduct(h_lo)
                h_cutoff = Deduct_Cutoff/wage1
                test_h_cutoff1 = foc_hrs_single_below_deduct(h_cutoff)

                h_hi = 1.0d0
                test_h_hi1 = foc_hrs_single_below_deduct(h_hi)
                if (h_cutoff < h_hi) then
                    test_h_cutoff2 = foc_hrs_single_above_deduct(h_cutoff + 1d-5)
                    test_h_hi2 = foc_hrs_single_above_deduct(h_hi)
                end if

                if (test_h_lo > 0d0) then
                    ! Uh > w*Uc, so  h = 0d0
                    h_sol = h_lo
                else
                    if (h_cutoff < h_hi) then 
                        if (test_h_cutoff1 > 0d0) then
                            ! h is between 0 and h_cutoff, find using bisection
                            h_a = h_lo
                            h_b = h_cutoff
                            ! call zbren(foc_hrs_single_below_deduct, h_a, h_b)

                            call zeroin(foc_hrs_mar_below_deduct, h_a, h_b, tol=1d-8, xzero=hsol, fzero=fsol_test, iflag=ic)
                            ! if (abs(h_b-xsol_test) > 1d-5) then
                            !     write(file_id,*) 'warning'
                            ! end if
                            h_sol = hsol !h_b
                        else if (test_h_cutoff2 > 0d0) then
                            ! h = h_cutoff
                            h_sol = h_cutoff 
                        else if (test_h_hi2 > 0d0) then
                            ! h is between h_cutoff and 1, find using bisection
                            h_a = h_cutoff + 1d-5
                            h_b = h_hi
                            ! call zbren(foc_hrs_single_above_deduct, h_a, h_b)

                            call zeroin(foc_hrs_single_above_deduct, h_a, h_b, tol=1d-8, xzero=hsol, fzero=fsol_test, iflag=ic)
                            ! if (abs(h_b-xsol_test) > 1d-5) then
                            !     write(file_id,*) 'warning'
                            ! end if
                            h_sol = hsol !h_b
                        else
                            h_sol = h_hi  
                        end if

                    else
                        if (test_h_hi1 > 0d0) then
                            ! h is between 0 and h_hi, find using bisection
                            h_a = h_lo
                            h_b = h_hi
                            ! call zbren(foc_hrs_single_below_deduct, h_a, h_b)

                            call zeroin(foc_hrs_single_below_deduct, h_a, h_b, tol=1d-8, xzero=hsol, fzero=fsol_test, iflag=ic)
                            ! if (abs(h_b-xsol_test) > 1d-5) then
                            !     write(file_id,*) 'warning'
                            ! end if
                            h_sol = hsol !h_b
                        else
                            h_sol = h_hi
                        end if

                    end if
                end if   

                if (i_gend == 1) then
                    laborsinglem(ik,ium)=h_sol
                else
                    laborsinglef(ik,ium)=h_sol
                end if

            end do

        end do

    end subroutine lsupply_baseline

    subroutine lsupply_nit_regime1(ik)
        ! for the case where d <= yhat
        use Model_Parameters
        use PolicyFunctions
        use glob0
        use root_module
        
        implicit none
        
        integer, INTENT(IN) :: ik
        integer :: ium, iuf
        
        real(8) :: test_h_lo, h_cutoff, test_h_cutoff1, test_h_cutoff2, test_h_hi1, test_h_hi2, h_lo, h_hi, h_sol
        real(8) :: h_a, h_b, test_a, test_b
        integer :: i_gend
        real(8) :: test_hm_lo, test_hf_lo, test_hm_hi, test_hf_hi, hf_tmp
        real(8) :: hm_lo, hf_lo
        real(8) :: hm_hi, hf_hi
        real(8) :: hm_d, hf_d
        real(8) :: hm_sol, hf_sol
        real(8) :: test_h_d   
        
        real(8) :: hsol, fsol_test
        integer :: ic
        real(8) :: hhat, htmp1, htmp2, htmp3
        real(8) :: htmp1m, htmp2m, htmp3m, htmp1f, htmp2f, htmp3f, ytot1, ytot2
        
        
        !CALL ERSET(IERSVR, IPACT, ISACT)
        
        dum=c_grid(ik)
        cmod = dum
        h_hi = 1.0d0
        
        ! Married couples:
        do ium=1,nw
            do iuf=1,nw
                
                wagem=wage_grid(ium)
                wagef=wage_grid(iuf)
                
                wage1 = wagem
                wage2 = wagef
                chi1 = chim
                chi2 = chif
                eta1 = etam
                eta2 = etaf           
                
                htmp1m = ( (1d0-t_const)*wage1*(1d0-s_nit-t_employee) / chi1 / cmod / (1d0+tc) )**(1d0/eta1)
                htmp2m = ( (1d0-t_const)*wage1*(1d0-tau_L-s_nit-t_employee) / chi1 / cmod / (1d0+tc) )**(1d0/eta1) 
                htmp3m = ( (1d0-t_const)*wage1*(1d0-tau_L-t_employee) / chi1 / cmod / (1d0+tc) )**(1d0/eta1)
                
                
                htmp1f = ( (1d0-t_const)*wage2*(1d0-s_nit-t_employee) / chi2 / cmod / (1d0+tc) )**(1d0/eta2)
                htmp2f = ( (1d0-t_const)*wage2*(1d0-tau_L-s_nit-t_employee) / chi2 / cmod / (1d0+tc) )**(1d0/eta2)
                htmp3f = ( (1d0-t_const)*wage2*(1d0-tau_L-t_employee) / chi2 / cmod / (1d0+tc) )**(1d0/eta2)
                
                ytot1 = wage1*htmp1m + wage2*htmp1f
                ytot2 = wage1*htmp2m + wage2*htmp2f
                
                !if (ytot  <= yhat) then
                if (ytot1  <= Deduct_cutoff_Mar) then
                    hm_sol = htmp1m
                    hf_sol = htmp1f
                else if (ytot2 <= yhat_mar) then
                    hm_sol = htmp2m
                    hf_sol = htmp2f                
                else
                    hm_sol = min(h_hi, htmp3m)
                    hf_sol = min(h_hi, htmp3f)
                end if     
                
                laborm(ik,ium,iuf)=hm_sol
                laborf(ik,ium,iuf)=hf_sol            
                
                
            end do
        end do
        
        ! Only one spouse works:
        do i_gend = 1, 2
            
            if (i_gend == 1) then
                chi1 = chim
                eta1 = etam
            else
                chi1 = chif
                eta1 = etaf
            end if
            
            do ium=1,nw
                wagem=wage_grid(ium)
                wage1 = wagem         
                
                htmp1 = ( (1d0-t_const)*wage1*(1d0-s_nit-t_employee) / chi1 / cmod / (1d0+tc) )**(1d0/eta1)
                htmp2 = ( (1d0-t_const)*wage1*(1d0-tau_L-s_nit-t_employee) / chi1 / cmod / (1d0+tc) )**(1d0/eta1)    
                htmp3 = ( (1d0-t_const)*wage1*(1d0-tau_L-t_employee) / chi1 / cmod / (1d0+tc) )**(1d0/eta1)
                
                ytot1 = htmp1*wage1
                ytot2 = htmp2*wage1
                
                if (ytot1 <= Deduct_cutoff_Mar) then
                    h_sol = htmp1
                else if (ytot2  <= yhat_mar) then
                    h_sol = htmp2
                else
                    h_sol = min(h_hi, htmp3)
                end if            
                
                if (i_gend == 1) then
                    labormwork(ik,ium)=h_sol
                else
                    laborfwork(ik,ium)=h_sol
                end if
                
            end do
        end do
        
        
        ! Singles:
        do i_gend = 1, 2
            
            if (i_gend == 1) then
                chi1 = chims
                eta1 = etam
            else
                chi1 = chifs
                eta1 = etaf
            end if
            
            do ium=1,nw
                
                wagem=wage_grid(ium)
                wage1 = wagem          
                
                
                htmp1 = ( (1d0-t_const)*wage1*(1d0-s_nit-t_employee) / chi1 / cmod / (1d0+tc) )**(1d0/eta1)
                htmp2 = ( (1d0-t_const)*wage1*(1d0-tau_L-s_nit-t_employee) / chi1 / cmod / (1d0+tc) )**(1d0/eta1)   
                htmp3 = ( (1d0-t_const)*wage1*(1d0-tau_L-t_employee) / chi1 / cmod / (1d0+tc) )**(1d0/eta1) 
                ytot1 = wage1*htmp1
                ytot2 = wage1*htmp2
                
                if (ytot1 <= Deduct_cutoff) then
                    h_sol = htmp1
                else if (ytot2 <= yhat) then
                    h_sol = htmp2    
                else
                    h_sol = min(h_hi, htmp3)
                end if   
                
                if (i_gend == 1) then
                    laborsinglem(ik,ium)=h_sol
                else
                    laborsinglef(ik,ium)=h_sol
                end if
                
            end do
            
        end do
        
    end subroutine lsupply_nit_regime1

    subroutine lsupply_nit_regime2(ik)
        ! for the case where d > yhat
        ! for the case where d <= yhat
        use Model_Parameters
        use PolicyFunctions
        use glob0
        use root_module
        
        implicit none
        
        integer, INTENT(IN) :: ik
        integer :: ium, iuf
        
        real(8) :: test_h_lo, h_cutoff, test_h_cutoff1, test_h_cutoff2, test_h_hi1, test_h_hi2, h_lo, h_hi, h_sol
        real(8) :: h_a, h_b, test_a, test_b
        integer :: i_gend
        real(8) :: test_hm_lo, test_hf_lo, test_hm_hi, test_hf_hi, hf_tmp
        real(8) :: hm_lo, hf_lo
        real(8) :: hm_hi, hf_hi
        real(8) :: hm_d, hf_d
        real(8) :: hm_sol, hf_sol
        real(8) :: test_h_d   
        
        real(8) :: hsol, fsol_test
        integer :: ic
        real(8) :: hhat, htmp1, htmp2, htmp3
        real(8) :: htmp1m, htmp2m, htmp3m, htmp1f, htmp2f, htmp3f, ytot1, ytot2
        
        dum=c_grid(ik)
        cmod = dum
        h_hi = 1.0d0
        
        ! Married couples:
        do ium=1,nw
            do iuf=1,nw
                
                wagem=wage_grid(ium)
                wagef=wage_grid(iuf)
                
                wage1 = wagem
                wage2 = wagef
                chi1 = chim
                chi2 = chif
                eta1 = etam
                eta2 = etaf           
                
                htmp1m = ( (1d0-t_const)*wage1*(1d0-s_nit-t_employee) / chi1 / cmod / (1d0+tc) )**(1d0/eta1)
                htmp2m = ( (1d0-t_const)*wage1*(1d0-t_employee) / chi1 / cmod / (1d0+tc) )**(1d0/eta1) 
                htmp3m = ( (1d0-t_const)*wage1*(1d0-tau_L-t_employee) / chi1 / cmod / (1d0+tc) )**(1d0/eta1)
                
                
                htmp1f = ( (1d0-t_const)*wage2*(1d0-s_nit-t_employee) / chi2 / cmod / (1d0+tc) )**(1d0/eta2)
                htmp2f = ( (1d0-t_const)*wage2*(1d0-t_employee) / chi2 / cmod / (1d0+tc) )**(1d0/eta2)
                htmp3f = ( (1d0-t_const)*wage2*(1d0-tau_L-t_employee) / chi2 / cmod / (1d0+tc) )**(1d0/eta2)
                
                ytot1 = wage1*htmp1m + wage2*htmp1f
                ytot2 = wage1*htmp2m + wage2*htmp2f
                
                if (ytot1  <= Deduct_cutoff_Mar) then
                    hm_sol = htmp1m
                    hf_sol = htmp1f
                else if (ytot2 <= yhat_mar) then
                    hm_sol = htmp2m
                    hf_sol = htmp2f                
                else
                    hm_sol = min(h_hi, htmp3m)
                    hf_sol = min(h_hi, htmp3f)
                end if     
                
                laborm(ik,ium,iuf)=hm_sol
                laborf(ik,ium,iuf)=hf_sol            
                
                
            end do
        end do
        
        ! Only one spouse works:
        do i_gend = 1, 2
            
            if (i_gend == 1) then
                chi1 = chim
                eta1 = etam
            else
                chi1 = chif
                eta1 = etaf
            end if
            
            do ium=1,nw
                wagem=wage_grid(ium)
                wage1 = wagem         
                
                htmp1 = ( (1d0-t_const)*wage1*(1d0-s_nit-t_employee) / chi1 / cmod / (1d0+tc) )**(1d0/eta1)
                htmp2 = ( (1d0-t_const)*wage1*(1d0-t_employee) / chi1 / cmod / (1d0+tc) )**(1d0/eta1)    
                htmp3 = ( (1d0-t_const)*wage1*(1d0-tau_L-t_employee) / chi1 / cmod / (1d0+tc) )**(1d0/eta1)
                
                ytot1 = htmp1*wage1
                ytot2 = htmp2*wage1
                
                if (ytot1 <= Deduct_cutoff) then
                    h_sol = htmp1
                else if (ytot2  <= yhat) then
                    h_sol = htmp2
                else
                    h_sol = min(h_hi, htmp3)
                end if            
                
                if (i_gend == 1) then
                    labormwork(ik,ium)=h_sol
                else
                    laborfwork(ik,ium)=h_sol
                end if
                
            end do
        end do
        
        
        ! Singles:
        do i_gend = 1, 2
            
            if (i_gend == 1) then
                chi1 = chims
                eta1 = etam
            else
                chi1 = chifs
                eta1 = etaf
            end if
            
            do ium=1,nw
                
                wagem=wage_grid(ium)
                wage1 = wagem
                
                htmp1 = ( (1d0-t_const)*wage1*(1d0-s_nit-t_employee) / chi1 / cmod / (1d0+tc) )**(1d0/eta1)
                htmp2 = ( (1d0-t_const)*wage1*(1d0-t_employee) / chi1 / cmod / (1d0+tc) )**(1d0/eta1)   
                htmp3 = ( (1d0-t_const)*wage1*(1d0-tau_L-t_employee) / chi1 / cmod / (1d0+tc) )**(1d0/eta1) 
                ytot1 = wage1*htmp1
                ytot2 = wage1*htmp2
                
                if (ytot1 <= Deduct_cutoff) then
                    h_sol = htmp1
                else if (ytot2 <= yhat) then
                    h_sol = htmp2    
                else
                    h_sol = min(h_hi, htmp3)
                end if   
                
                if (i_gend == 1) then
                    laborsinglem(ik,ium)=h_sol
                else
                    laborsinglef(ik,ium)=h_sol
                end if
                
            end do
            
        end do
        
    end subroutine lsupply_nit_regime2
    
end module LaborSupply
