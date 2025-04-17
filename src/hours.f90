module hours
    
    implicit none
    
    real(8) :: lab_inc_spouse
    real(8) :: wage1, wage2, chi1, chi2, eta1, eta2
    real(8) :: cmod
    !$OMP THREADPRIVATE(lab_inc_spouse)
    !$OMP THREADPRIVATE(wage1, wage2, chi1, chi2, eta1, eta2, cmod)
    
contains

    function mar_hours_deduct(h1)
        USE Model_Parameters
        !USE glob0
        implicit none
        real(8), intent(in) :: h1
        real(8) :: h2
        real(8) :: mar_hours_deduct
        
        h2 = min(1.0d0, hours_ratio(h1, wage2, wage1, chi2, chi1, eta2, eta1))
        mar_hours_deduct = h1*wage1 + h2*wage2 - Deduct_Cutoff
    end function mar_hours_deduct

    function hours_ratio(h2,w1,w2,chi1,chi2,eta1,eta2)
        implicit none
        real(8) :: h2, w1, w2, chi1, chi2, eta1, eta2
        real(8) :: hours_ratio
        
        ! hf = ( chim*wf / chif*wm )^(1/etaf) * hm^(etam/etaf)
        ! 1 = f, 2 = m
        ! hf = ( chi2*w1 / chi1*w2 )^(1/eta1) * hm^(eta2/eta1)
        
        hours_ratio = (chi2*w1/chi1/w2)**(1d0/eta1)*(h2**(eta2/eta1))
    end function hours_ratio
    
    function foc_hrs_mar_below_deduct(h)
        USE Model_Parameters
        !USE glob0
        implicit none
        real(8), intent(in) :: h
        !real(8) :: C
        real(8) :: foc_hrs_mar_below_deduct
        
        foc_hrs_mar_below_deduct = chi1*(h**eta1)*cmod*(1d0+tc) - (1d0-t_const)*wage1*(1.0d0-t_employee)
    
    end function foc_hrs_mar_below_deduct
    
    function foc_hrs_mar_above_deduct(h)
        USE Model_Parameters
        !USE glob0
        implicit none
        real(8), intent(in) :: h
        real(8) :: h2, tot_lab_inc, lab_inc_spouse
        !real(8) :: C
        real(8) :: foc_hrs_mar_above_deduct
        
        h2 = min(1.0d0, hours_ratio(h, wage2, wage1, chi2, chi1, eta2, eta1))
        lab_inc_spouse = wage2*h2
        tot_lab_inc = lab_inc_spouse + wage1*h
        
        foc_hrs_mar_above_deduct = chi1*(h**eta1)*cmod*(1d0+tc) - (1d0-t_const)*wage1*(theta(1)*(1d0-theta(2))*(tot_lab_inc-Deduct_Cutoff)**(-theta(2))-t_employee)
    
    end function foc_hrs_mar_above_deduct    
    
    function foc_hrs_mar_above_deduct_1w(h)
        USE Model_Parameters
        !USE glob0
        implicit none
        real(8), intent(in) :: h
        real(8) :: tot_lab_inc
        !real(8) :: C
        real(8) :: foc_hrs_mar_above_deduct_1w
        
        !h2 = min(1.0d0, hours_ratio(h, wage2, wage1, chi2, chi1, eta2, eta1))
        !lab_inc_spouse = wage2*h2
        tot_lab_inc = wage1*h
        
        foc_hrs_mar_above_deduct_1w = chi1*(h**eta1)*cmod*(1d0+tc) - (1d0-t_const)*wage1*(theta(1)*(1d0-theta(2))*(tot_lab_inc-Deduct_Cutoff)**(-theta(2))-t_employee)
    
    end function foc_hrs_mar_above_deduct_1w        

    function foc_hrs_single_below_deduct(h)
        USE Model_Parameters
        !USE glob0
    
        implicit none
    
        real(8), intent(in) :: h
        real(8) :: C
        real(8) :: foc_hrs_single_below_deduct
    
        !C = dum
    
        foc_hrs_single_below_deduct = chi1*(h**eta1)*cmod*(1d0+tc) - (1d0-t_const)*wage1*(1.0d0-t_employee)
    
    end function foc_hrs_single_below_deduct

    function foc_hrs_single_above_deduct(h)
        USE Model_Parameters
        !USE glob0
    
        implicit none
    
        real(8), intent(in) :: h
        real(8) :: C
        real(8) :: foc_hrs_single_above_deduct
    
        !C = dum
    
        foc_hrs_single_above_deduct = chi1*(h**eta1)*cmod*(1d0+tc) - (1d0-t_const)*wage1*(thetas(1)*(1d0-thetas(2))*(wage1*h-Deduct_Cutoff)**(-thetas(2)) - t_employee)
    
    end function foc_hrs_single_above_deduct

    
end module hours