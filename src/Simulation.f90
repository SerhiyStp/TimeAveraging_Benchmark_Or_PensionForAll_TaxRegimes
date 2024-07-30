module Simulations_mod
    
    implicit none
    
    integer, parameter :: nsim   = 12000 !12000 !10 !10000         ! Number of households used for simulation
    integer, parameter :: nsim2   = 80 !80         ! Number of simulations    
    
    real(8), dimension (:,:,:,:), allocatable :: Sim1m,Sim1f,exp2m,exp2f ! Array to hold simulation results
    integer, dimension (:,:,:,:), allocatable :: exp1m,exp1f,expR1f,expR1m ! Array to hold simulation results
    real(8), dimension (:,:,:,:), allocatable :: SimR1m, SimR1f ! Array to hold simulation results for retired
    real(8), dimension (:,:,:), allocatable :: Random3m, Random3f, marstatm, marstatf, partshock, partshock2
    real(8), dimension (:,:), allocatable :: Random1m, Random2m, Random1f, Random2f, marstatm_init, marstatf_init    
    
    integer, parameter :: SIM_MARSTAT = 10
    integer, parameter :: SIM_PARTID = 4    
    
contains

    subroutine InitSimulation()
        use Model_Parameters, only: prob_u, T, Tret, prob_a, prob_fc, nfc, trans_u
        
        integer :: ik
        integer :: it2, iu2, i, um, uf
        real(8) :: dum5
        integer :: rand_seed_size
        integer, allocatable :: seed(:)        
        
        allocate(exp1m(nsim2,nsim,T+1,6))
        allocate(exp1f(nsim2,nsim,T+1,6))
        allocate(exp2m(nsim2,nsim,T+1,6))
        allocate(exp2f(nsim2,nsim,T+1,6))
        allocate(expR1m(nsim2,nsim,Tret+1,4))
        allocate(expR1f(nsim2,nsim,Tret+1,4))
        allocate(Random3m(nsim2,nsim,T+Tret))
        allocate(Random3f(nsim2,nsim,T+Tret))   
        allocate(partshock(nsim2,nsim,1))
        allocate(partshock2(nsim2,nsim,1))
        allocate(Random1m(nsim2,nsim))
        allocate(Random1f(nsim2,nsim))
        allocate(Random2m(nsim2,nsim))
        allocate(Random2f(nsim2,nsim))

        call random_seed(SIZE=rand_seed_size)
        allocate(seed(rand_seed_size))
        seed = 123456789
        call random_seed(PUT=seed(1:rand_seed_size))

        
        call random_number(random1m)
        call random_number(random2m)
        call random_number(random3m)
        call random_number(random1f)
        call random_number(random2f)
        call random_number(random3f)
        call random_number(partshock)
        call random_number(partshock2)
        
       
        !$OMP PARALLEL PRIVATE(ik, it2, iu2, dum5, i, um, uf)
        !$OMP DO SCHEDULE(DYNAMIC)        
        do ik = 1, nsim2
            exp2m(ik,:,:,1)=0d0
            exp2f(ik,:,:,1)=0d0 
            !Initialzing the idiosyncratic shocks to wages
            do it2=1,nsim
                if(Random1m(ik,it2)<Prob_u(1,1)) then
                    exp1m(ik,it2,1,3)=1
                elseif((Prob_u(1,1)<Random1m(ik,it2)).AND.(Random1m(ik,it2)<(Prob_u(1,1)+Prob_u(1,2)))) then
                    exp1m(ik,it2,1,3)=2
                elseif(((Prob_u(1,1)+Prob_u(1,2))<Random1m(ik,it2)).AND.(Random1m(ik,it2)<(Prob_u(1,1)+Prob_u(1,2)+Prob_u(1,3)))) then
                    exp1m(ik,it2,1,3)=3
                elseif(((Prob_u(1,1)+Prob_u(1,2)+Prob_u(1,3))<Random1m(ik,it2)).AND.(Random1m(ik,it2)<(Prob_u(1,1)+Prob_u(1,2)+Prob_u(1,3)+Prob_u(1,4)))) then
                    exp1m(ik,it2,1,3)=4
                else
                    exp1m(ik,it2,1,3)=5
                end if
            end do 
        
            do it2=1,nsim
                if(Random1f(ik,it2)<Prob_u(2,1)) then
                    exp1f(ik,it2,1,3)=1
                elseif((Prob_u(2,1)<Random1f(ik,it2)).AND.(Random1f(ik,it2)<(Prob_u(2,1)+Prob_u(2,2)))) then
                    exp1f(ik,it2,1,3)=2
                elseif(((Prob_u(2,1)+Prob_u(2,2))<Random1f(ik,it2)).AND.(Random1f(ik,it2)<(Prob_u(2,1)+Prob_u(2,2)+Prob_u(2,3)))) then
                    exp1f(ik,it2,1,3)=3
                elseif(((Prob_u(2,1)+Prob_u(2,2)+Prob_u(2,3))<Random1f(ik,it2)).AND.(Random1f(ik,it2)<(Prob_u(2,1)+Prob_u(2,2)+Prob_u(2,3)+Prob_u(2,4)))) then
                    exp1f(ik,it2,1,3)=4
                else
                    exp1f(ik,it2,1,3)=5
                end if
            end do

            !Initialzing the distribution of abilit2ies
            do it2=1,nsim
                if(Random2m(ik,it2)<Prob_a(1,1)) then
                    exp1m(ik,it2,:,2)=1
                elseif((Prob_a(1,1)<Random2m(ik,it2)).AND.(Random2m(ik,it2)<(Prob_a(1,1)+Prob_a(1,2)))) then
                    exp1m(ik,it2,:,2)=2
                elseif(((Prob_a(1,1)+Prob_a(1,2))<Random2m(ik,it2)).AND.(Random2m(ik,it2)<(Prob_a(1,1)+Prob_a(1,2)+Prob_a(1,3)))) then
                    exp1m(ik,it2,:,2)=3
                elseif(((Prob_a(1,1)+Prob_a(1,2)+Prob_a(1,3))<Random2m(ik,it2)).AND.(Random2m(ik,it2)<(Prob_a(1,1)+Prob_a(1,2)+Prob_a(1,3)+Prob_a(1,4)))) then
                    exp1m(ik,it2,:,2)=4
                else
                    exp1m(ik,it2,:,2)=5
                end if
            end do

            do it2=1,nsim
                if(Random2f(ik,it2)<Prob_a(2,1)) then
                    exp1f(ik,it2,:,2)=1
                elseif((Prob_a(2,1)<Random2f(ik,it2)).AND.(Random2f(ik,it2)<(Prob_a(2,1)+Prob_a(2,2)))) then
                    exp1f(ik,it2,:,2)=2
                elseif(((Prob_a(2,1)+Prob_a(2,2))<Random2f(ik,it2)).AND.(Random2f(ik,it2)<(Prob_a(2,1)+Prob_a(2,2)+Prob_a(2,3)))) then
                    exp1f(ik,it2,:,2)=3
                elseif(((Prob_a(2,1)+Prob_a(2,2)+Prob_a(2,3))<Random2f(ik,it2)).AND.(Random2f(ik,it2)<(Prob_a(2,1)+Prob_a(2,2)+Prob_a(2,3)+Prob_a(2,4)))) then
                    exp1f(ik,it2,:,2)=4
                else
                    exp1f(ik,it2,:,2)=5
                end if
            end do 
            
            !Initialzing the participation costs 
            do it2=1,nsim
                iu2=1
                dum5=Prob_fc(1,iu2)
                do while((dum5<Partshock(ik,it2,1)).AND.(iu2<nfc))
                    iu2=iu2+1
                    dum5=dum5+Prob_fc(1,iu2)
                end do
                exp1f(ik,it2,:,6)=iu2
            end do            
            
            do it2=1,nsim
                exp1m(ik,it2,:,6)=1
            end do            
            
            !Uppdating the idiosyncratic wage shock
            do i = 1, T
                do it2=1,nsim
                    um=exp1m(ik,it2,i,3)
                    if(Random3m(ik,it2,i)<trans_u(1,um,1)) then
                        exp1m(ik,it2,i+1,3)=1
                    elseif((trans_u(1,um,1)<Random3m(ik,it2,i)).AND.(Random3m(ik,it2,i)<(trans_u(1,um,1)+trans_u(1,um,2)))) then
                        exp1m(ik,it2,i+1,3)=2
                    elseif(((trans_u(1,um,1)+trans_u(1,um,2))<Random3m(ik,it2,i)).AND.(Random3m(ik,it2,i)<(trans_u(1,um,1)+trans_u(1,um,2)+trans_u(1,um,3)))) then
                        exp1m(ik,it2,i+1,3)=3
                    elseif(((trans_u(1,um,1)+trans_u(1,um,2)+trans_u(1,um,3))<Random3m(ik,it2,i)).AND.(Random3m(ik,it2,i)<(trans_u(1,um,1)+trans_u(1,um,2)+trans_u(1,um,3)+trans_u(1,um,4)))) then
                        exp1m(ik,it2,i+1,3)=4
                    else
                        exp1m(ik,it2,i+1,3)=5
                    end if
                end do
                do it2=1,nsim
                    uf=exp1f(ik,it2,i,3)
                    if(Random3f(ik,it2,i)<trans_u(2,uf,1)) then
                        exp1f(ik,it2,i+1,3)=1
                    elseif((trans_u(2,uf,1)<Random3f(ik,it2,i)).AND.(Random3f(ik,it2,i)<(trans_u(2,uf,1)+trans_u(2,uf,2)))) then
                        exp1f(ik,it2,i+1,3)=2
                    elseif(((trans_u(2,uf,1)+trans_u(2,uf,2))<Random3f(ik,it2,i)).AND.(Random3f(ik,it2,i)<(trans_u(2,uf,1)+trans_u(2,uf,2)+trans_u(2,uf,3)))) then
                        exp1f(ik,it2,i+1,3)=3
                    elseif(((trans_u(2,uf,1)+trans_u(2,uf,2)+trans_u(2,uf,3))<Random3f(ik,it2,i)).AND.(Random3f(ik,it2,i)<(trans_u(2,uf,1)+trans_u(2,uf,2)+trans_u(2,uf,3)+trans_u(2,uf,4)))) then
                        exp1f(ik,it2,i+1,3)=4
                    else
                        exp1f(ik,it2,i+1,3)=5
                    end if
                end do
            end do
            
            expR1m(ik,:,1,2)=exp1m(ik,:,T+1,3)
            expR1f(ik,:,1,2)=exp1f(ik,:,T+1,3)
            
            do i=1,Tret
                do it2=1,nsim
                    um=expR1m(ik,it2,i,2)
                    if(Random3m(ik,it2,T+i)<trans_u(1,um,1)) then
                        expR1m(ik,it2,i+1,2)=1
                    elseif((trans_u(1,um,1)<Random3m(ik,it2,T+i)).AND.(Random3m(ik,it2,T+i)<(trans_u(1,um,1)+trans_u(1,um,2)))) then
                        expR1m(ik,it2,i+1,2)=2
                    elseif(((trans_u(1,um,1)+trans_u(1,um,2))<Random3m(ik,it2,T+i)).AND.(Random3m(ik,it2,T+i)<(trans_u(1,um,1)+trans_u(1,um,2)+trans_u(1,um,3)))) then
                        expR1m(ik,it2,i+1,2)=3
                    elseif(((trans_u(1,um,1)+trans_u(1,um,2)+trans_u(1,um,3))<Random3m(ik,it2,T+i)).AND.(Random3m(ik,it2,T+i)<(trans_u(1,um,1)+trans_u(1,um,2)+trans_u(1,um,3)+trans_u(1,um,4)))) then
                        expR1m(ik,it2,i+1,2)=4
                    else
                        expR1m(ik,it2,i+1,2)=5
                    end if
                end do
            
                do it2=1,nsim
                    uf=expR1f(ik,it2,i,2)
                    if(Random3f(ik,it2,T+i)<trans_u(2,uf,1)) then
                        expR1f(ik,it2,i+1,2)=1
                    elseif((trans_u(2,uf,1)<Random3f(ik,it2,T+i)).AND.(Random3f(ik,it2,T+i)<(trans_u(2,uf,1)+trans_u(2,uf,2)))) then
                        expR1f(ik,it2,i+1,2)=2
                    elseif(((trans_u(2,uf,1)+trans_u(2,uf,2))<Random3f(ik,it2,T+i)).AND.(Random3f(ik,it2,T+i)<(trans_u(2,uf,1)+trans_u(2,uf,2)+trans_u(2,uf,3)))) then
                        expR1f(ik,it2,i+1,2)=3
                    elseif(((trans_u(2,uf,1)+trans_u(2,uf,2)+trans_u(2,uf,3))<Random3f(ik,it2,T+i)).AND.(Random3f(ik,it2,T+i)<(trans_u(2,uf,1)+trans_u(2,uf,2)+trans_u(2,uf,3)+trans_u(2,uf,4)))) then
                        expR1f(ik,it2,i+1,2)=4
                    else
                        expR1f(ik,it2,i+1,2)=5
                    end if
                end do
            end do
        end do
        !$OMP END DO    
        !$OMP END PARALLEL 
        
        deallocate(Random1m)
        deallocate(Random1f)
        deallocate(Random2m)
        deallocate(Random2f)
        deallocate(Random3m)
        deallocate(Random3f)
        
        allocate(Sim1m(nsim2,nsim,T+1,17))
        allocate(Sim1f(nsim2,nsim,T+1,17))
        allocate(SimR1m(nsim2,nsim,Tret+1,19))
        allocate(SimR1f(nsim2,nsim,Tret+1,19))

        allocate(marstatm(nsim2,nsim,T))
        allocate(marstatf(nsim2,nsim,T))
        allocate(marstatm_init(nsim2,nsim))
        allocate(marstatf_init(nsim2,nsim))        
        
        call random_number(marstatm)
        call random_number(marstatf)
        call random_number(marstatm_init)
        call random_number(marstatf_init)        
        
        !$OMP PARALLEL PRIVATE(ik)
        !$OMP DO SCHEDULE(DYNAMIC)         
        do ik = 1, nsim2
            Sim1m(ik,:,:,:)=0d0
            Sim1f(ik,:,:,:)=0d0            
            exp1m(ik,:,:,4)=0
            exp1f(ik,:,:,4)=0            
            call Marriages(ik)            
        end do
        !$OMP END DO    
        !$OMP END PARALLEL    
        
        deallocate(marstatm)
        deallocate(marstatf)
        deallocate(marstatm_init)
        deallocate(marstatf_init)           
        
    end subroutine InitSimulation
    
    subroutine Marriages(ik)
        use Model_Parameters, only: a, match, probd, probm, T
        use SortH
        integer, INTENT(IN) :: ik
        real(8), allocatable :: mixm(:,:), mixf(:,:), mixdum(:,:)
        integer :: i, it2
        integer :: count_sing, count_mar
        integer :: i_m, i_f
        integer :: i_husb, i_wife
        
        ! Initial period
        allocate(mixm(nsim,3), mixf(nsim,3), mixdum(nsim,3))
        mixm=0d0
        mixf=0d0
        do it2=1,nsim
            !Single men
            mixm(it2,1) = marstatm_init(ik, it2)
            mixm(it2,2) = marstatm_init(ik, it2) + match*A(1,exp1m(ik,it2,1,2)) 
            mixm(it2,3) = it2*1d0
            mixf(it2,1) = marstatf_init(ik, it2)
            mixf(it2,2) = marstatf_init(ik, it2) + match*A(2,exp1f(ik,it2,1,2))
            mixf(it2,3) = it2*1d0
        end do
        
        count_mar = int(0.1259d0*nsim)
        !Sorting single men and women by marriage shock       
        call SortMultiArrayDescending(mixm, mixdum, 1)
        mixm(1:count_mar,:) = mixdum(1:count_mar,:)
        call SortMultiArrayDescending(mixm(1:count_mar,:), mixdum(1:count_mar,:), 2)
        mixm(1:count_mar,:) = mixdum(1:count_mar,:)
        call SortMultiArrayDescending(mixf, mixdum, 1)
        mixf(1:count_mar,:) = mixdum(1:count_mar,:) 
        call SortMultiArrayDescending(mixf(1:count_mar,:), mixdum(1:count_mar,:), 2)
        mixf(1:count_mar,:) = mixdum (1:count_mar,:) 
        
        do it2 = 1, count_mar
            i_husb = int(mixm(it2, 3)+0.01)
            i_wife = int(mixf(it2, 3)+0.01)
            exp1m(ik,i_husb,1,SIM_PARTID) = i_wife
            exp1f(ik,i_wife,1,SIM_PARTID) = i_husb
            Sim1m(ik,i_husb,1,SIM_MARSTAT) = 1d0
            Sim1f(ik,i_wife,1,SIM_MARSTAT) = 1d0
        end do    
        

        
        ! All other periods
        do i=1, T
            mixm=0d0
            mixf=0d0
            ! Married
            do it2=1,nsim
                if (Sim1m(ik,it2,i,SIM_MARSTAT)>0.5d0) then
                    i_husb = it2
                    i_wife = exp1m(ik,i_husb,i,SIM_PARTID)
                    if (marstatm(ik,it2,i)>Probd(i)) then
                        !staying married
                        exp1m(ik,i_husb,i+1,SIM_PARTID)=i_wife
                        exp1f(ik,i_wife,i+1,SIM_PARTID)=i_husb
                        Sim1m(ik,i_husb,i+1,SIM_MARSTAT)=1d0
                        Sim1f(ik,i_wife,i+1,SIM_MARSTAT)=1d0
                    else
                        !getting divorced
                        exp1m(ik,i_husb,i+1,SIM_PARTID)=0
                        exp1f(ik,i_wife,i+1,SIM_PARTID)=0
                        Sim1m(ik,i_husb,i+1,SIM_MARSTAT)=0d0
                        Sim1f(ik,i_wife,i+1,SIM_MARSTAT)=0d0
                        !Sim1m(ik,i_husb,i+1,1)=Sim1m(ik,it2,i+1,1)/2d0
                        !Sim1f(ik,i_wife,i+1,1)=Sim1f(ik,it2,i+1,1)/2d0  
                        !exp1f(ik,i_wife,i+1,5)=0
                    end if
                end if                
            
            end do      
            
            ! Single
            count_sing = count( Sim1m(ik,:,i,10)<0.5d0 )
            count_mar=int(count_sing*Probm(i))
            i_m = 1
            i_f = 1
            do it2=1, nsim
                !Single men
                if (Sim1m(ik,it2,i,SIM_MARSTAT)<0.5d0) then
                    mixm(i_m, 1) = marstatm(ik,it2,i)
                    mixm(i_m, 2) = marstatm(ik,it2,i) + match*A(1,exp1m(ik,it2,i,2))
                    mixm(i_m, 3) = it2*1d0
                    i_m = i_m + 1
                end if
                !Single women
                if (Sim1f(ik,it2,i,SIM_MARSTAT)<0.5d0) then
                    mixf(i_f, 1) = marstatf(ik,it2,i)
                    mixf(i_f, 2) = marstatf(ik,it2,i) + match*A(2,exp1f(ik,it2,i,2))
                    mixf(i_f, 3) = it2*1d0
                    i_f = i_f + 1
                end if
            end do
            !Sorting single men and women by marriage shock Mn
            call SortMultiArrayDescending(mixm(1:count_sing,:), mixdum(1:count_sing,:), 1)
            mixm(1:count_mar,:) = mixdum(1:count_mar,:)
            call SortMultiArrayDescending(mixm(1:count_mar,:), mixdum(1:count_mar,:), 2)
            mixm(1:count_mar,:) = mixdum(1:count_mar,:)            
            call SortMultiArrayDescending(mixf(1:count_sing,:), mixdum(1:count_sing,:), 1)
            mixf(1:count_mar,:) = mixdum(1:count_mar,:)
            call SortMultiArrayDescending(mixf(1:count_mar,:), mixdum(1:count_mar,:), 2)
            mixf(1:count_mar,:) = mixdum(1:count_mar,:)            

            do it2 = 1, count_mar
                i_husb = int(mixm(it2,3)+0.01)
                i_wife = int(mixf(it2,3)+0.01)
                exp1m(ik,i_husb,i+1,SIM_PARTID) = i_wife
                exp1f(ik,i_wife,i+1,SIM_PARTID) = i_husb
                Sim1m(ik,i_husb,i+1,SIM_MARSTAT)=1d0
                Sim1f(ik,i_wife,i+1,SIM_MARSTAT)=1d0
                !Sim1m(ik,i_husb,i+1,1)=Sim1m(ik,i_husb,i+1,1)+Sim1f(ik,i_wife,i+1,1)
                !Sim1f(ik,i_wife,i+1,1)=Sim1m(ik,i_husb,i+1,1)
            end do               
            
        end do
        
        deallocate(mixm, mixf, mixdum)
        
    end subroutine Marriages

    subroutine simulation(ik)
        !This subroutine simulates the lifecycle for nsim households
        use Model_Parameters
        use PolicyFunctions
        use Utilities
        use PolicyFunctions_obj, only: pol_v_mar_lfp, pol_v_sing_lfp, pol_vs_ret
        use GlobParams
        use SortH
        use, intrinsic :: IEEE_ARITHMETIC    

        implicit none

        integer, INTENT(IN) :: ik
        integer :: i,iam,ium,iaf,iuf,um,uf,it2,it3,it4,it5,j,count2,count3,ifc,ifcm,iu2,ir,irm,iu3
        real(8) :: d1,d2,ix,ixm,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum15,y,r_ret,pnt1(2),pnt12(2),pnt2(3),exp_grid_dum(nexp),exp_grid_dum2(nexp),INTERP2D(nk,nexp),INTERP3D(nk,nexp,nexp)
        real(8) :: P1,P2,P3,P4 !, mixm(nsim,3), mixf(nsim,3),mixdum(nsim,3),mixdum2(nsim,3)
        real(8) :: test1, test2
        real(8) :: tmp
        real(8) :: LFP_M, LFP_F, LFP_S
        integer :: LFP_M_I, LFP_F_I, LFP_S_I    
        real(8) :: dum3_m, dum3_f
        integer :: count_sing, count_mar
        integer :: i_m, i_f
        integer :: i_husb, i_wife
        

        !dum2=0.0

        !Sim1m(ik,:,:,:)=0d0
        !Sim1f(ik,:,:,:)=0d0
        !
        !exp1m(ik,:,:,4)=0
        !exp1f(ik,:,:,4)=0
        
        !call Marriages(ik)

        !Starting to simulate age 25-64, using piecewise linear interpolation of the policy functions

        do i=1,T

            exp_grid_dum=exp_grid(:,i)

            do it2=1,nsim

                !print *, it2, nsim

                if(Sim1m(ik,it2,i,10)>0.5d0) then

                    it3=exp1m(ik,it2,i,4)
                    dum2=Sim1m(ik,it2,i,1)
                    iam=exp1m(ik,it2,i,2)
                    ium=exp1m(ik,it2,i,3)
                    ixm=exp2m(ik,it2,i,1)
                    ifcm=exp1m(ik,it2,i,6)
                    iaf=exp1f(ik,it3,i,2)
                    iuf=exp1f(ik,it3,i,3)
                    ix=exp2f(ik,it3,i,1)
                    ifc=exp1f(ik,it3,i,6)
                    pnt2 = (/dum2, ix, ixm/)
                    
                    !Total years of marriage
                    Sim1m(ik,it2,i+1,17)=Sim1m(ik,it2,i,17)+1d0
                    Sim1f(ik,it3,i+1,17)=Sim1f(ik,it3,i,17)+1d0
                    
                    
                    INTERP3D=lfpm(:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                    LFP_M = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                    if (LFP_M >= 0.5d0) then
                        LFP_M_I = LFP_1
                    else
                        LFP_M_I = LFP_0
                    end if
                    INTERP3D=lfpf(:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                    LFP_F = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)     
                    if (LFP_F >= 0.5d0) then
                        LFP_F_I = LFP_1
                    else
                        LFP_F_I = LFP_0
                    end if                  
                    

                    !Next period's capital
                    INTERP3D=k_lfp(:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm,LFP_M_I,LFP_F_I)
                    Sim1m(ik,it2,i+1,1) = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                    Sim1f(ik,it3,i+1,1) = Sim1m(ik,it2,i+1,1)

                    !This period's consumption
                    INTERP3D=c_lfp(:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm,LFP_M_I,LFP_F_I)
                    dum4 = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                    Sim1m(ik,it2,i,2) = dum4
                    Sim1f(ik,it3,i,2) = dum4
                    
                    if (dum4 <= 0d0) then
                        print *, 'warning: c=0'
                    end if                

                    ! Male wage, work hours and earnings
                    Sim1m(ik,it2,i,3) = (1d0/(1d0+exp(kappa*(i-agestart))))*wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)
                    if (LFP_M_I == LFP_0) then
                        dum5 = 0d0
                    else
                        INTERP3D=nm_lfp(:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm,LFP_M_I,LFP_F_I)
                        dum5=trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                        dum5=max(dum5,0d0)
                        if(dum5<(minhours/2d0)) then
                            dum5=0d0
                        else
                            dum5=max(minhours,dum5)
                        end if
                    end if                                
                    Sim1m(ik,it2,i,4) = dum5
                    Sim1m(ik,it2,i,5) = dum5*(1d0/(1d0+exp(kappa*(i-agestart))))*wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)

                    ! Female wage, work hours and earnings
                    Sim1f(ik,it3,i,3) = (1d0/(1d0+exp(kappa*(i-agestart))))*wage(2,a(2,iaf),dble(ix),u(2,iuf))/(1d0+t_employer)
                    if (LFP_F_I == LFP_0) then
                        dum6 = 0d0
                    else
                        INTERP3D=nf_lfp(:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm,LFP_M_I,LFP_F_I)
                        dum6=trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                        dum6=max(dum6,0d0)
                        if(dum6<(minhours/2d0)) then
                            dum6=0d0
                        else
                            dum6=max(minhours,dum6)
                        end if
                    end if                                
                    Sim1f(ik,it3,i,4) = dum6
                    Sim1f(ik,it3,i,5) = dum6*(1d0/(1d0+exp(kappa*(i-agestart))))*wage(2,a(2,iaf),dble(ix),u(2,iuf))/(1d0+t_employer)

                    !Household income and taxes
                    dum3=dum5*(1d0/(1d0+exp(kappa*(i-agestart))))*wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)+dum6*(1d0/(1d0+exp(kappa*(i-agestart))))*wage(2,a(2,iaf),dble(ix),u(2,iuf))/(1d0+t_employer)
                    Sim1m(ik,it2,i,6) = dum3
                    Sim1f(ik,it3,i,6) = dum3
                    
                    tmp = labor_tax_married(dum3)
                    Sim1m(ik,it2,i,7) = tmp
                    Sim1f(ik,it3,i,7) = tmp


                    Sim1m(ik,it2,i,8)= dum4*tc
                    Sim1f(ik,it3,i,8)= dum4*tc
                    Sim1m(ik,it2,i,9)=dum3*t_employee+t_employer*dum3
                    Sim1f(ik,it3,i,9)=dum3*t_employee+t_employer*dum3


                    !Individual marginal labor income tax rate including social security
                    if(Sim1m(ik,it2,i,6)>Deduct_cutoff_Mar) then
                        Sim1m(ik,it2,i,12)=1d0-(1d0-(t_employee+t_employer+(1d0-(1d0-theta(2))*theta(1)*Sim1m(ik,it2,i,6)**(-theta(2))))/(1d0*(1d0+t_employer)))*(1d0-t_const)
                        Sim1f(ik,it3,i,12)=1d0-(1d0-(t_employee+t_employer+(1d0-(1d0-theta(2))*theta(1)*Sim1f(ik,it3,i,6)**(-theta(2))))/(1d0*(1d0+t_employer)))*(1d0-t_const)
                    else if(Sim1m(ik,it2,i,6)>0d0) then
                        Sim1m(ik,it2,i,12)=1d0-(1d0-(t_employee+t_employer)/(1d0*(1d0+t_employer)))*(1d0-t_const)
                        Sim1f(ik,it2,i,12)=1d0-(1d0-(t_employee+t_employer)/(1d0*(1d0+t_employer)))*(1d0-t_const)
                    else
                        Sim1m(ik,it2,i,12)= 0d0
                        Sim1f(ik,it3,i,12)= 0d0
                    end if

                    if(Sim1m(ik,it2,i,6)>0d0) then
                        Sim1m(ik,it2,i,14)=1d0-(1d0-(Sim1m(ik,it2,i,9)+Sim1m(ik,it2,i,7))/(Sim1m(ik,it2,i,6)*(1d0+t_employer)))*(1d0-t_const)
                        Sim1f(ik,it3,i,14)=1d0-(1d0-(Sim1f(ik,it3,i,9)+Sim1f(ik,it3,i,7))/(Sim1f(ik,it3,i,6)*(1d0+t_employer)))*(1d0-t_const)
                    else
                        Sim1m(ik,it2,i,14)= 0d0
                        Sim1f(ik,it3,i,14)= 0d0
                    end if


                    if(Sim1f(ik,it3,i,4)>1d-3) then
                        exp2f(ik,it3,i+1,1)=exp2f(ik,it3,i,1)+1d0
                    else
                        exp2f(ik,it3,i+1,1)=exp2f(ik,it3,i,1)*(1d0-deltaexp)
                    end if

                    if(Sim1m(ik,it2,i,4)>1d-3) then
                        exp2m(ik,it2,i+1,1)=exp2m(ik,it2,i,1)+1d0
                    else
                        exp2m(ik,it2,i+1,1)=exp2m(ik,it2,i,1)*(1d0-deltaexp)
                    end if

                    !Social Welfare
                    dum3_m = pol_v_mar_lfp(iam, ium, iaf, iuf, i, ifc, ifcm, LFP_M_I, LFP_F_I, MEN)%eval(pnt2)
                    Sim1m(ik,it2,i,11) = dum3_m
                    dum3_f = pol_v_mar_lfp(iam, ium, iaf, iuf, i, ifc, ifcm, LFP_M_I, LFP_F_I, WOMEN)%eval(pnt2)
                    Sim1f(ik,it3,i,11) = dum3_f                   

                    !Ex-post social welfare
                    if((dum5>0.001d0).AND.(dum6>0.001d0)) then
                        Sim1m(ik,it2,i+1,16) = Sim1m(ik,it2,i,16)+ (beta**(i-1d0))*(Uc(dum4)+Ul(dum5,dum6)-fc(1,ifc)-fcm(1,ifcm))
                        Sim1f(ik,it3,i+1,16) = Sim1f(ik,it3,i,16)+ (beta**(i-1d0))*(Uc(dum4)+Ul(dum5,dum6)-fc(1,ifc)-fcm(1,ifcm))
                    elseif((dum5>0.001d0).AND.(dum6<0.001d0)) then
                        Sim1m(ik,it2,i+1,16) = Sim1m(ik,it2,i,16)+ (beta**(i-1d0))*(Uc(dum4)+Ul(dum5,dum6)-fcm(1,ifcm))
                        Sim1f(ik,it3,i+1,16) = Sim1f(ik,it3,i,16)+ (beta**(i-1d0))*(Uc(dum4)+Ul(dum5,dum6)-fcm(1,ifcm))
                    elseif((dum5<0.001d0).AND.(dum6>0.001d0)) then
                        Sim1m(ik,it2,i+1,16) = Sim1m(ik,it2,i,16)+ (beta**(i-1d0))*(Uc(dum4)+Ul(dum5,dum6)-fc(1,ifc))
                        Sim1f(ik,it3,i+1,16) = Sim1f(ik,it3,i,16)+ (beta**(i-1d0))*(Uc(dum4)+Ul(dum5,dum6)-fc(1,ifc))
                    else
                        Sim1m(ik,it2,i+1,16) = Sim1m(ik,it2,i,16)+ (beta**(i-1d0))*(Uc(dum4))
                        Sim1f(ik,it3,i+1,16) = Sim1f(ik,it3,i,16)+ (beta**(i-1d0))*(Uc(dum4))
                    end if
                    
                    if (.not. ieee_is_finite(Sim1m(ik,it2,i+1,16))) then
                        print *, 'warning: infinite!'
                    end if
                    if (.not. ieee_is_finite(Sim1f(ik,it2,i+1,16))) then
                        print *, 'warning: infinite!'
                    end if                  
                        
                        
                    !Extra tax revenue
                    if(Sim1m(ik,it2,i,6)>0d0) then
                        dum3=Sim1m(ik,it2,i,6)
                        tmp = t_const*(after_tax_labor_inc_married(dum3) - t_employee*dum3)
                        Sim1m(ik,it2,i,13)=tmp
                        Sim1f(ik,it3,i,13)=tmp
                    else
                        Sim1m(ik,it2,i,13)= 0d0
                        Sim1f(ik,it3,i,13)= 0d0
                    end if

                else

                    dum2=Sim1m(ik,it2,i,1)
                    iam=exp1m(ik,it2,i,2)
                    ium=exp1m(ik,it2,i,3)
                    ixm=exp2m(ik,it2,i,1)
                    ifcm=exp1m(ik,it2,i,6)
                    j=1
                    pnt1 = (/dum2, ixm/)

                    !Total years of marriage
                    Sim1m(ik,it2,i+1,17)=Sim1m(ik,it2,i,17)
                    
                    INTERP2D=lfps(j,:,:,iam,ium,i,ifcm)
                    LFP_S=bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    if (LFP_S >= 0.5d0) then
                        LFP_S_I = LFP_1
                    else
                        LFP_S_I = LFP_0
                    end if                 
                    
                    !Next period's capital
                    INTERP2D=ks_lfp(j,:,:,iam,ium,i,ifcm,LFP_S_I)
                    Sim1m(ik,it2,i+1,1)=bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)

                    !This period's consumption
                    INTERP2D=cs_lfp(j,:,:,iam,ium,i,ifcm,LFP_S_I)
                    dum4 = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    Sim1m(ik,it2,i,2) = dum4

                    ! Male wage, work hours and earnings
                    Sim1m(ik,it2,i,3) = (1d0/(1d0+exp(kappa*(i-agestart))))*wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)
                    if (LFP_S_I == LFP_0) then
                        dum5 = 0d0
                    else                
                        INTERP2D=ns(j,:,:,iam,ium,i,ifcm)
                        dum5=bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                        dum5=max(dum5,0d0)
                        if(dum5<(minhours/2d0)) then
                            dum5=0d0
                        else
                            dum5=max(minhours,dum5)
                        end if
                    end if                 
                    Sim1m(ik,it2,i,4) = dum5
                    Sim1m(ik,it2,i,5) = dum5*(1d0/(1d0+exp(kappa*(i-agestart))))*wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)

                    !Household income and taxes

                    dum3=dum5*(1d0/(1d0+exp(kappa*(i-agestart))))*wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)
                    Sim1m(ik,it2,i,6) = dum3
                    if(dum3>0d0) then
                        Sim1m(ik,it2,i,7)= labor_tax_single(dum3)
                    else
                        Sim1m(ik,it2,i,7)=0d0
                    end if

                    Sim1m(ik,it2,i,8)= dum4*tc
                    Sim1m(ik,it2,i,9)=dum3*t_employee+t_employer*dum3

                    !Individual marginal labor income tax rate including social security

                    if(Sim1m(ik,it2,i,6)>Deduct_cutoff) then
                        Sim1m(ik,it2,i,12)=1d0-(1d0-(t_employee+t_employer+(1d0-(1d0-thetas(2))*thetas(1)*Sim1m(ik,it2,i,6)**(-thetas(2))))/(1d0*(1d0+t_employer)))*(1d0-t_const)
                    else if(Sim1m(ik,it2,i,6)>0d0) then
                        Sim1m(ik,it2,i,12)=1d0-(1d0-(t_employee+t_employer)/(1d0*(1d0+t_employer)))*(1d0-t_const)
                    else
                        Sim1m(ik,it2,i,12)= 0d0
                    end if

                    if(Sim1m(ik,it2,i,6)>Deduct_cutoff) then
                        Sim1m(ik,it2,i,14)=1d0-(1d0-(Sim1m(ik,it2,i,9)+Sim1m(ik,it2,i,7))/(Sim1m(ik,it2,i,6)*(1d0+t_employer)))*(1d0-t_const)
                    else
                        Sim1m(ik,it2,i,14)= 0d0
                    end if


                    if(Sim1m(ik,it2,i,4)>1d-3) then
                        exp2m(ik,it2,i+1,1)=exp2m(ik,it2,i,1)+1d0
                    else
                        exp2m(ik,it2,i+1,1)=exp2m(ik,it2,i,1)*(1d0-deltaexp)
                    end if

                    !Social welfare
                    dum3 = pol_v_sing_lfp(j, iam, ium, i, ifcm, LFP_S_I)%eval(pnt1)
                    Sim1m(ik,it2,i,11) = dum3
                    
                    !Ex-post social welfare
                    if(dum5>0.001d0) then
                        Sim1m(ik,it2,i+1,16) = Sim1m(ik,it2,i,16) + (beta**(i-1d0))*(Uc(dum4)-chims*(dum5**(1d0+etam))/(1d0+etam)-fcm(2,ifcm))
                    else
                        Sim1m(ik,it2,i+1,16) = Sim1m(ik,it2,i,16) + (beta**(i-1d0))*(Uc(dum4))
                    end if

                    !Extra tax revenue
                    if(Sim1m(ik,it2,i,6)>0d0) then
                        dum3=Sim1m(ik,it2,i,6)
                        Sim1m(ik,it2,i,13)= t_const*( after_tax_labor_inc_married(dum3) - t_employee*dum3 ) 
                    else
                        Sim1m(ik,it2,i,13)= 0d0
                    end if

                end if

                !Single women

                if(Sim1f(ik,it2,i,10)<0.5d0) then

                    dum2=Sim1f(ik,it2,i,1)
                    iam=exp1f(ik,it2,i,2)
                    ium=exp1f(ik,it2,i,3)
                    ix=exp2f(ik,it2,i,1)
                    ifc=exp1f(ik,it2,i,6)
                    j=2
                    pnt1 = (/dum2, ix/)

                    !Total years of marriage
                    Sim1f(ik,it2,i+1,17)=Sim1f(ik,it2,i,17)
                    
                    INTERP2D=lfps(j,:,:,iam,ium,i,ifc)
                    LFP_S=bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    if (LFP_S >= 0.5d0) then
                        LFP_S_I = LFP_1
                    else
                        LFP_S_I = LFP_0
                    end if                  
                    
                    !Next period's capital
                    INTERP2D=ks_lfp(j,:,:,iam,ium,i,ifc,LFP_S_I)  
                    Sim1f(ik,it2,i+1,1)=bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)

                    !This period's consumption
                    INTERP2D=cs_lfp(j,:,:,iam,ium,i,ifc,LFP_S_I) 
                    dum4=bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    Sim1f(ik,it2,i,2) = dum4    

                    ! Female wage, work hours and earnings
                    Sim1f(ik,it2,i,3) = (1d0/(1d0+exp(kappa*(i-agestart))))*wage(2,a(2,iam),dble(ix),u(2,ium))/(1d0+t_employer)
                    if (LFP_S_I == LFP_0) then
                        dum5 = 0d0
                    else      
                        INTERP2D=ns(j,:,:,iam,ium,i,ifc)
                        dum5=bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                        dum5=max(dum5,0d0)
                        if(dum5<(minhours/2d0)) then
                            dum5=0d0
                        else
                            dum5=max(minhours,dum5)
                        end if
                    end if                                
                    Sim1f(ik,it2,i,4) = dum5
                    Sim1f(ik,it2,i,5) = dum5*(1d0/(1d0+exp(kappa*(i-agestart))))*wage(2,a(2,iam),dble(ix),u(2,ium))/(1d0+t_employer)

                    !Household income and taxes
                    dum3=dum5*(1d0/(1d0+exp(kappa*(i-agestart))))*wage(2,a(2,iam),dble(ix),u(2,ium))/(1d0+t_employer)
                    Sim1f(ik,it2,i,6) = dum3
                    if(dum3>0d0) then
                        Sim1f(ik,it2,i,7)= labor_tax_single(dum3)
                    else
                        Sim1f(ik,it2,i,7)=0d0
                    end if

                    Sim1f(ik,it2,i,8)= dum4*tc
                    Sim1f(ik,it2,i,9)=dum3*t_employee+t_employer*dum3

                    !Individual marginal labor income tax rate including social security
                    if(Sim1f(ik,it2,i,6)>Deduct_cutoff) then
                        Sim1f(ik,it2,i,12)=1d0-(1d0-(t_employee+t_employer+(1d0-(1d0-thetas(2))*thetas(1)*Sim1f(ik,it2,i,6)**(-thetas(2))))/(1d0*(1d0+t_employer)))*(1d0-t_const)
                    else if (Sim1f(ik,it2,i,6)>0d0) then
                        Sim1f(ik,it2,i,12)= 1d0-(1d0-(t_employee+t_employer)/(1d0*(1d0+t_employer)))*(1d0-t_const)
                    end if

                    if(Sim1f(ik,it2,i,6)>0d0) then
                        Sim1f(ik,it2,i,14)=1d0-(1d0-(Sim1f(ik,it2,i,9)+Sim1f(ik,it2,i,7))/(Sim1f(ik,it2,i,6)*(1d0+t_employer)))*(1d0-t_const)
                    else
                        Sim1f(ik,it2,i,14)= 0d0
                    end if


                    if(Sim1f(ik,it2,i,4)>1d-3) then
                        y=Sim1f(ik,it2,i,6)
                        exp2f(ik,it2,i+1,1)=exp2f(ik,it2,i,1)+1d0
                        Sim1f(ik,it2,i,15)=Sim1f(ik,it2,i+1,1)-((dum2 + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+(1d0-t_const)*(after_tax_labor_inc_single(y)-y*tSS_employee(y))-Sim1f(ik,it2,i,2)*(1d0+tc))/(1d0+mu)
                    else
                        exp2f(ik,it2,i+1,1)=exp2f(ik,it2,i,1)*(1d0-deltaexp)
                        Sim1f(ik,it2,i,15)=Sim1f(ik,it2,i+1,1)-((dum2+Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+Unemp_benefit-Sim1f(ik,it2,i,2)*(1d0+tc))/(1d0+mu)
                    end if

                    !Social Welfare
                    dum3 = pol_v_sing_lfp(j, iam, ium, i, ifc, LFP_S_I)%eval(pnt1)
                    Sim1f(ik,it2,i,11) = dum3
                    
                    !Ex-post social welfare
                    if(dum5>0.001d0) then
                        Sim1f(ik,it2,i+1,16) = Sim1f(ik,it2,i,16) + (beta**(i-1d0))*(Uc(dum4)-chifs*(dum5**(1d0+etaf))/(1d0+etaf)-fc(2,ifc))
                    else
                        Sim1f(ik,it2,i+1,16) = Sim1f(ik,it2,i,16) + (beta**(i-1d0))*(Uc(dum4))
                    endif

                    !Extra tax revenue
                    if(Sim1f(ik,it2,i,6)>0d0) then
                        dum3=Sim1f(ik,it2,i,6)
                        Sim1f(ik,it2,i,13)=t_const*(after_tax_labor_inc_married(dum3)- t_employee*dum3)
                    else
                        Sim1f(ik,it2,i,13)= 0d0
                    end if

                end if

            end do

            !mixm=0d0
            !mixf=0d0
        
            do it2=1,nsim
                if (Sim1m(ik,it2,i,SIM_MARSTAT) > 0.5d0) then
                    if (Sim1m(ik,it2,i+1,SIM_MARSTAT) < 0.5d0) then
                        i_husb = it2    
                        i_wife = exp1m(ik,i_husb,i,SIM_PARTID)
                        Sim1m(ik,i_husb,i+1,1)=Sim1m(ik,it2,i+1,1)/2d0
                        Sim1f(ik,i_wife,i+1,1)=Sim1f(ik,it2,i+1,1)/2d0
                        exp1f(ik,i_wife,i+1,5)=0
                    end if
                end if              
            end do            
            
            
            do it2=1, nsim
                !Single men
                if (Sim1m(ik,it2,i,SIM_MARSTAT)<0.5d0) then
                    if (Sim1m(ik,it2,i+1,SIM_MARSTAT)>0.5d0) then
                        i_husb = it2
                        i_wife = exp1m(ik,i_husb,i+1,SIM_PARTID) 
                        Sim1m(ik,i_husb,i+1,1) = Sim1m(ik,i_husb,i+1,1) + Sim1f(ik,i_wife,i+1,1)
                        Sim1f(ik,i_wife,i+1,1) = Sim1m(ik,i_husb,i+1,1)
                    end if
                end if
            end do
                              
            dum5=maxval(Sim1m(ik,:,i+1,1))

        end do    


        !Beginning the simulation for retired households

        SimR1m(ik,:,1,1)=Sim1m(ik,:,T+1,1)
        SimR1m(ik,:,1,12)=exp2m(ik,:,T+1,1)
        SimR1m(ik,:,1,19)=Sim1m(ik,:,T+1,16)
        expR1m(ik,:,1,2)=exp1m(ik,:,T+1,3)
        expR1m(ik,:,1,3)=2
        SimR1f(ik,:,1,1)=Sim1f(ik,:,T+1,1)
        SimR1f(ik,:,1,12)=exp2f(ik,:,T+1,1)
        SimR1f(ik,:,1,19)=Sim1f(ik,:,T+1,16)
        expR1f(ik,:,1,2)=exp1f(ik,:,T+1,3)
        expR1f(ik,:,1,3)=2

        do i=1,Tret
            expR1m(ik,:,i,4)=exp1m(ik,:,T+1,6)
            expR1m(ik,:,i,1)=exp1m(ik,:,T+1,2)
            expR1f(ik,:,i,4)=exp1f(ik,:,T+1,6)
            expR1f(ik,:,i,1)=exp1f(ik,:,T+1,2)
        end do


        do i=1,Tret

            r_ret=((1d0+r)/OmegaRet2(i))-1d0
            exp_grid_dum=exp_grid(:,T+i)
            if(i<Tret) then
                exp_grid_dum2=exp_grid(:,T+i+1)
            end if

            do it2=1,nsim

                if(Sim1m(ik,it2,T,10)>0.5) then

                    it3=exp1m(ik,it2,T,4)
                    dum2=SimR1m(ik,it2,i,1)
                    ixm=SimR1m(ik,it2,i,12)
                    iam=expR1m(ik,it2,i,1)
                    ium=expR1m(ik,it2,i,2)
                    irm=expR1m(ik,it2,i,3)
                    ifcm=expR1m(ik,it2,i,4)
                    ix=SimR1f(ik,it3,i,12)
                    iaf=expR1f(ik,it3,i,1)
                    iuf=expR1f(ik,it3,i,2)
                    ir=expR1f(ik,it3,i,3)
                    ifc=expR1f(ik,it3,i,4)
                    pnt2 = (/dum2, ix, ixm/)

                    
                    !Retirement decision
                    
                    if((irm==2).AND.(ir==2)) then
                        
                        INTERP3D=V_ret2(2,2,:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                        P4=trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                        INTERP3D=V_ret2(2,1,:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                        P3=trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                        
                        If(P4<P3) then
                            P4=P3
                            irm=1
                            ir=2
                        end if
                        
                        INTERP3D=V_ret2(1,2,:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                        P2=trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                            
                        If(P4<P2) then
                            P4=P2
                            irm=2
                            ir=1
                        end if
                        
                        INTERP3D=V_ret2(1,1,:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                        P1=trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                            
                        If(P4<P1) then
                            P4=P1
                            irm=1
                            ir=1
                        end if
                        
                    elseif((irm==1).AND.(ir==2)) then
                        INTERP3D=V_ret2(2,1,:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                        P4=trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                        INTERP3D=V_ret2(1,1,:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                        P1=trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                            
                        If(P4<P1) then
                            P4=P1
                            irm=1
                            ir=1
                        end if
                    
                    elseif((irm==2).AND.(ir==1)) then
                        INTERP3D=V_ret2(1,2,:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                        P4=trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                        INTERP3D=V_ret2(1,1,:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                        P1=trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                        
                        If(P4<P1) then
                            P4=P1
                            irm=1
                            ir=1
                        end if
                        
                    else
                        ir=1
                        irm=1
                    end if
                    
                    expR1m(ik,it2,i,3)=irm
                    expR1m(ik,it2,i+1,3)=irm
                    expR1f(ik,it3,i,3)=ir
                    expR1f(ik,it3,i+1,3)=ir
                    
                    !Next period's capital
                    INTERP3D=k_ret(ir,irm,:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                    SimR1m(ik,it2,i+1,1) = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                    SimR1f(ik,it3,i+1,1) = SimR1m(ik,it2,i+1,1)

                    !Next periods's consumption
                    INTERP3D=c_ret(ir,irm,:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                    dum4 = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                    SimR1m(ik,it2,i,2) = dum4
                    SimR1f(ik,it3,i,2) = dum4
                    SimR1m(ik,it2,i,3)= dum4*tc
                    SimR1f(ik,it3,i,3)= dum4*tc

                    ! Male wage, work hours and earnings

                    INTERP3D=nm_ret(ir,irm,:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                    SimR1m(ik,it2,i,4) = (1d0/(1d0+exp(kappa*(44+i-(agestart-1)))))*wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)
                    dum5=trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                    dum5=max(dum5,0d0)
                    if(dum5<(minhours/2d0)) then
                        dum5=0d0
                    else
                        dum5=max(minhours,dum5)
                    end if
                    SimR1m(ik,it2,i,5) = dum5
                    SimR1m(ik,it2,i,6) = dum5*(1d0/(1d0+exp(kappa*(44+i-(agestart-1)))))*wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)

                    ! Female wage, work hours and earnings

                    INTERP3D=nf_ret(ir,irm,:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                    SimR1f(ik,it3,i,4) = (1d0/(1d0+exp(kappa*(44+i-(agestart-1)))))*wage(2,a(2,iaf),dble(ix),u(2,iuf))/(1d0+t_employer)
                    dum6=trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                    dum6=max(dum6,0d0)
                    if(dum6<(minhours/2d0)) then
                        dum6=0d0
                    else
                        dum6=max(minhours,dum6)
                    end if
                    SimR1f(ik,it3,i,5) = dum6
                    SimR1f(ik,it3,i,6) = dum6*(1d0/(1d0+exp(kappa*(44+i-(agestart-1)))))*wage(2,a(2,iaf),dble(ix),u(2,iuf))/(1d0+t_employer)

                    !Household income and taxes

                    dum3=dum5*(1d0/(1d0+exp(kappa*(44+i-(agestart-1)))))*wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)+dum6*(1d0/(1d0+exp(kappa*(44+i-(agestart-1)))))*wage(2,a(2,iaf),dble(ix),u(2,iuf))/(1d0+t_employer)
                    SimR1f(ik,it3,i,7) = dum3
                    SimR1m(ik,it2,i,7) = dum3
                    if(dum3>0d0) then
                        tmp =  labor_tax_married(dum3)
                        SimR1f(ik,it3,i,8) = tmp
                        SimR1m(ik,it2,i,8) = tmp
                    else
                        SimR1f(ik,it3,i,8) = 0d0
                        SimR1m(ik,it2,i,8) = 0d0
                    end if

                    SimR1f(ik,it3,i,9)=dum3*t_employee+t_employer*dum3
                    SimR1m(ik,it2,i,9)=dum3*t_employee+t_employer*dum3

                    if(SimR1m(ik,it2,i,7)>Deduct_cutoff_Mar) then
                        SimR1m(ik,it2,i,15)=1d0-(1d0-(t_employee+t_employer+(1d0-(1d0-theta(2))*theta(1)*SimR1m(ik,it2,i,7)**(-theta(2))))/(1d0*(1d0+t_employer)))*(1d0-t_const)
                        SimR1f(ik,it3,i,15)=1d0-(1d0-(t_employee+t_employer+(1d0-(1d0-theta(2))*theta(1)*SimR1f(ik,it3,i,7)**(-theta(2))))/(1d0*(1d0+t_employer)))*(1d0-t_const)
                    else if(SimR1m(ik,it2,i,7)>0d0) then
                        SimR1m(ik,it2,i,15)=1d0-(1d0-(t_employee+t_employer)/(1d0*(1d0+t_employer)))*(1d0-t_const)
                        SimR1f(ik,it3,i,15)=1d0-(1d0-(t_employee+t_employer)/(1d0*(1d0+t_employer)))*(1d0-t_const)
                    else
                        SimR1m(ik,it2,i,15)= 0d0
                        SimR1f(ik,it3,i,15)= 0d0
                    end if

                    if(SimR1m(ik,it2,i,7)>0d0) then
                        SimR1m(ik,it2,i,17)=1d0-(1d0-(SimR1m(ik,it2,i,9)+SimR1m(ik,it2,i,8))/(SimR1m(ik,it2,i,7)*(1d0+t_employer)))*(1d0-t_const)
                        SimR1f(ik,it3,i,17)=1d0-(1d0-(SimR1f(ik,it3,i,9)+SimR1f(ik,it3,i,8))/(SimR1f(ik,it3,i,7)*(1d0+t_employer)))*(1d0-t_const)
                    else
                        SimR1m(ik,it2,i,17)= 0d0
                        SimR1f(ik,it3,i,17)= 0d0
                    end if

                    
                    !Experience
                    if(SimR1f(ik,it3,i,5)>1d-3) then
                        SimR1f(ik,it3,i+1,12)=SimR1f(ik,it3,i,12)+1d0
                    else
                        SimR1f(ik,it3,i+1,12)=SimR1f(ik,it3,i,12)*(1d0-deltaexp)
                    end if

                    if(SimR1m(ik,it2,i,5)>1d-3) then
                        SimR1m(ik,it2,i+1,12)=SimR1m(ik,it2,i,12)+1d0
                    else
                        SimR1m(ik,it2,i+1,12)=SimR1m(ik,it2,i,12)*(1d0-deltaexp)
                    end if
                    
                    
                    !Social Security

                    if(irm==2) then
                        SimR1m(ik,it2,i,14)=0d0
                    else
                        !Pension depends on expected wage conditional on ability and experience
                        SimR1m(ik,it2,i,14)=(psi0+psi1*av_earnings(1,1,iam)*min(1d0,dble(ixm)/35d0))*(1d0-t_const)
                    end if

                    if(ir==2) then
                        SimR1f(ik,it3,i,14)=0d0
                    else
                        SimR1f(ik,it3,i,14)=(psi0+psi1*av_earnings(2,1,iaf)*min(1d0,dble(ix)/35d0))*(1d0-t_const)
                    end if

                    !Budget
                    if(SimR1m(ik,it2,i,7)>0.001d0) then
                        SimR1m(ik,it2,i,18)=SimR1m(ik,it2,i+1,1)-((SimR1m(ik,it2,i,1)+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk)) + SimR1f(ik,it3,i,14)+ SimR1m(ik,it2,i,14) +lumpsum*0.5d0+ (1d0-t_const)*(after_tax_labor_inc_married(SimR1m(ik,it2,i,7))-SimR1m(ik,it2,i,7)*tSS_employee(SimR1m(ik,it2,i,7)))-SimR1m(ik,it2,i,2)*(1d0+tc))/(1d0+mu)
                    else
                        SimR1m(ik,it2,i,18)=SimR1m(ik,it2,i+1,1)-((SimR1m(ik,it2,i,1)+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk)) + SimR1f(ik,it3,i,14)+ SimR1m(ik,it2,i,14) +lumpsum*0.5d0-SimR1m(ik,it2,i,2)*(1d0+tc))/(1d0+mu)
                    end if
                    
                    !Social welfare
                    INTERP3D=v_ret(ir,irm,:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                    dum3 = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                    SimR1m(ik,it2,i,11) = dum3
                    SimR1f(ik,it3,i,11) = dum3

                    !Ex-post social welfare
                    if((dum5>0.001d0).AND.(dum6>0.001d0)) then
                        SimR1m(ik,it2,i+1,19) = SimR1m(ik,it2,i,19)+ (beta**(45+i-1d0))*(Uc(dum4)+Ul(dum5,dum6)-fc(1,ifc)-fcm(1,ifcm))
                        SimR1f(ik,it3,i+1,19) = SimR1f(ik,it3,i,19)+ (beta**(45+i-1d0))*(Uc(dum4)+Ul(dum5,dum6)-fc(1,ifc)-fcm(1,ifcm))
                    elseif((dum5>0.001d0).AND.(dum6<0.001d0)) then
                        SimR1m(ik,it2,i+1,19) = SimR1m(ik,it2,i,19)+ (beta**(45+i-1d0))*(Uc(dum4)+Ul(dum5,dum6)-fcm(1,ifcm))
                        SimR1f(ik,it3,i+1,19) = SimR1f(ik,it3,i,19)+ (beta**(45+i-1d0))*(Uc(dum4)+Ul(dum5,dum6)-fcm(1,ifcm))
                    elseif((dum5<0.001d0).AND.(dum6>0.001d0)) then
                        SimR1m(ik,it2,i+1,19) = SimR1m(ik,it2,i,19)+ (beta**(i-1d0))*(Uc(dum4)+Ul(dum5,dum6)-fc(1,ifc))
                        SimR1f(ik,it3,i+1,19) = SimR1f(ik,it3,i,19)+ (beta**(i-1d0))*(Uc(dum4)+Ul(dum5,dum6)-fc(1,ifc))
                    else
                        SimR1m(ik,it2,i+1,19) = SimR1m(ik,it2,i,19)+ (beta**(i-1d0))*(Uc(dum4))
                        SimR1f(ik,it3,i+1,19) = SimR1f(ik,it3,i,19)+ (beta**(i-1d0))*(Uc(dum4))
                    end if
                    
                    !Extra tax revenue
                    if(SimR1m(ik,it2,i,7)>0d0) then
                        dum3=SimR1m(ik,it2,i,7)
                        !SimR1m(ik,it2,i,16)=t_const*(1d0-(t_employee+tax_labor(dum3)))*dum3
                        !SimR1f(ik,it3,i,16)=t_const*(1d0-(t_employee+tax_labor(dum3)))*dum3
                        tmp=t_const*(after_tax_labor_inc_married(dum3)-dum3*t_employee)*dum3
                        SimR1m(ik,it2,i,16)=tmp
                        SimR1f(ik,it3,i,16)=tmp
                    else
                        SimR1m(ik,it2,i,16)= 0d0
                        SimR1f(ik,it3,i,16)= 0d0
                    end if

                else

                    dum2=SimR1m(ik,it2,i,1)
                    ixm=SimR1m(ik,it2,i,12)
                    iam=expR1m(ik,it2,i,1)
                    ium=expR1m(ik,it2,i,2)
                    irm=expR1m(ik,it2,i,3)
                    ifcm=expR1m(ik,it2,i,4)
                    j=1
                    pnt1 = (/dum2, ixm/)
                    
                    !Retirement decision
                    
                    if(irm==2) then
                        
                        INTERP2D=Vs_ret2(2,j,:,:,iam,ium,i,ifcm)
                        P2=bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                        INTERP2D=Vs_ret2(1,j,:,:,iam,ium,i,ifcm)
                        P1=bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                        
                        if(P1<P2) then
                            expR1m(ik,it2,i+1,3)=2
                        else
                            irm=1
                            expR1m(ik,it2,i,3)=1
                            expR1m(ik,it2,i+1,3)=1
                        end if
                        
                    else
                        expR1m(ik,it2,i+1,3)=1
                    end if

                    !Next period's capital
                    INTERP2D=ks_ret(irm,j,:,:,iam,ium,i,ifcm)
                    SimR1m(ik,it2,i+1,1) = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)

                    !This periods's consumption
                    INTERP2D=cs_ret(irm,j,:,:,iam,ium,i,ifcm)
                    dum4 = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    SimR1m(ik,it2,i,2) = dum4
                    SimR1m(ik,it2,i,3)= dum4*tc


                    !>Single male wage, work hours and earnings

                    INTERP2D=ns_ret(irm,j,:,:,iam,ium,i,ifcm)
                    SimR1m(ik,it2,i,4) = (1d0/(1d0+exp(kappa*(44+i-(agestart-1)))))*wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)
                    dum5=bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    dum5=max(dum5,0d0)
                    if(dum5<(minhours/2d0)) then
                        dum5=0d0
                    else
                        dum5=max(minhours,dum5)
                    end if
                    SimR1m(ik,it2,i,5) = dum5
                    SimR1m(ik,it2,i,6) = dum5*(1d0/(1d0+exp(kappa*(44+i-(agestart-1)))))*wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)

                    !Household income and taxes

                    dum3=dum5*(1d0/(1d0+exp(kappa*(44+i-(agestart-1)))))*wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)
                    SimR1m(ik,it2,i,7) = dum3
                    if(dum3>0d0) then
                        !SimR1m(ik,it2,i,8) = tax_labors(dum3)*dum3
                        !test1 = labor_tax_single(dum3)
                        SimR1m(ik,it2,i,8) = labor_tax_single(dum3)
                    else
                        SimR1m(ik,it2,i,8) = 0d0
                    end if

                    SimR1m(ik,it2,i,9)=dum3*t_employee+t_employer*dum3

                    !Individual marginal labor income tax rate including social security
                    if(SimR1m(ik,it2,i,7)>Deduct_Cutoff) then
                        SimR1m(ik,it2,i,15)=1d0-(1d0-(t_employee+t_employer+(1d0-(1d0-thetas(2))*thetas(1)*SimR1m(ik,it2,i,7)**(-thetas(2))))/(1d0*(1d0+t_employer)))*(1d0-t_const)
                    else if(SimR1m(ik,it2,i,7)>0d0) then
                        SimR1m(ik,it2,i,15)=1d0-(1d0-(t_employee+t_employer)/(1d0*(1d0+t_employer)))*(1d0-t_const)
                    else
                        SimR1m(ik,it2,i,15)= 0d0
                    end if

                    if(SimR1m(ik,it2,i,7)>0d0) then
                        SimR1m(ik,it2,i,17)=1d0-(1d0-(SimR1m(ik,it2,i,9)+SimR1m(ik,it2,i,8))/(SimR1m(ik,it2,i,7)*(1d0+t_employer)))*(1d0-t_const)
                    else
                        SimR1m(ik,it2,i,17)= 0d0
                    end if

                    !Experience

                    if(SimR1m(ik,it2,i,5)>1d-3) then
                        SimR1m(ik,it2,i+1,12)=SimR1m(ik,it2,i,12)+1d0
                    else
                        SimR1m(ik,it2,i+1,12)=SimR1m(ik,it2,i,12)*(1d0-deltaexp)
                    end if


                    if(irm==2) then
                        SimR1m(ik,it2,i,14)=0d0
                        if(SimR1m(ik,it2,i,7)>0.0001d0) then
                            y=SimR1m(ik,it2,i,7)
                            SimR1m(ik,it2,i,18)=SimR1m(ik,it2,i+1,1)-((dum2+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk))+lumpsum*0.5d0+(1d0-t_const)*(after_tax_labor_inc_single(y)-y*tSS_employee(y))-SimR1m(ik,it2,i,2)*(1d0+tc))/(1d0+mu)
                        else
                            SimR1m(ik,it2,i,18)=SimR1m(ik,it2,i+1,1)-((SimR1m(ik,it2,i,1)+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk)) + SimR1m(ik,it2,i,14)+lumpsum*0.5d0+(SimR1m(ik,it2,i,7)*(1d0-t_employee)-SimR1m(ik,it2,i,8))-SimR1m(ik,it2,i,2)*(1d0+tc))/(1d0+mu)
                        end if
                    else
                        SimR1m(ik,it2,i,14)=(psi0+psi1*av_earnings(1,2,iam)*min(1d0,dble(ixm)/35d0))*(1d0-t_const)
                        SimR1m(ik,it2,i,18)=SimR1m(ik,it2,i+1,1)-((SimR1m(ik,it2,i,1)+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk)) + SimR1m(ik,it2,i,14)+lumpsum*0.5d0+(SimR1m(ik,it2,i,7)*(1d0-t_employee)-SimR1m(ik,it2,i,8))-SimR1m(ik,it2,i,2)*(1d0+tc))/(1d0+mu)
                    end if

                    !Social welfare
                    INTERP2D=vs_ret(irm,j,:,:,iam,ium,i,ifcm)
                    dum3 = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    SimR1m(ik,it2,i,11) = dum3

                    
                    !Ex-post social welfare
                    if(dum5>0.001d0) then
                        SimR1m(ik,it2,i+1,19) = SimR1m(ik,it2,i,19) + (beta**(i-1d0))*(Uc(dum4)-chims*(dum5**(1d0+etam))/(1d0+etam)-fcm(2,ifcm))
                        !dum15=(beta**(i-1d0))*(Uc(dum4)-chims*(dum5**(1d0+etam))/(1d0+etam)-fcm(2,ifcm))
                    else
                        SimR1m(ik,it2,i+1,19) = SimR1m(ik,it2,i,19) + (beta**(i-1d0))*(Uc(dum4))
                        !dum15=(beta**(i-1d0))*(Uc(dum4))
                    end if
                    
                    !Euler error
                    if(i>2) then
                        dum8=SimR1m(ik,it2,i-1,2)
                        SimR1m(ik,it2,i,13) = dUc(dum8)-beta*OmegaRet(i)*((1d0+r*(1d0-tk))/(1d0+mu))*dUc(dum4)
                    end if

                    !Extra tax revenue
                    if(SimR1m(ik,it2,i,7)>0d0) then
                        dum3=SimR1m(ik,it2,i,7)
                        !SimR1m(ik,it2,i,16)=t_const*(1d0-(t_employee+tax_labor(dum3)))*dum3
                        !test1 = t_const*(after_tax_labor_inc_married(dum3)-dum3*t_employee)
                        SimR1m(ik,it2,i,16)= t_const*(after_tax_labor_inc_married(dum3)-dum3*t_employee)
                    else
                        SimR1m(ik,it2,i,16)= 0d0
                    end if

                end if


                !Single women
                if(Sim1f(ik,it2,T,10)<0.5) then

                    dum2=SimR1f(ik,it2,i,1)
                    ix=SimR1f(ik,it2,i,12)
                    iaf=expR1f(ik,it2,i,1)
                    iuf=expR1f(ik,it2,i,2)
                    ir=expR1f(ik,it2,i,3)
                    ifc=expR1f(ik,it2,i,4)
                    j=2
                    pnt1 = (/dum2, ix/)

                    !Retirement decision
                    
                    if(ir==2) then
                        
                        INTERP2D=Vs_ret2(2,j,:,:,iaf,iuf,i,ifc)
                        P2=bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                        INTERP2D=Vs_ret2(1,j,:,:,iaf,iuf,i,ifc)
                        P1=bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                        
                        if(P1<P2) then
                            expR1f(ik,it2,i+1,3)=2
                        else
                            ir=1
                            expR1f(ik,it2,i,3)=1
                            expR1f(ik,it2,i+1,3)=1
                        end if
                        
                    else
                        expR1f(ik,it2,i+1,3)=1
                    end if
                    
                    !if(i==1) then
                    !    dum2=k_grid(2)
                    !    ix=exp_grid_dum(3)
                    !    iaf=5
                    !    iuf=4
                    !    ir=2
                    !    pnt1 = (/dum2, ix/)
                    !end if
                    
                    !Next period's capital
                    INTERP2D=ks_ret(ir,j,:,:,iaf,iuf,i,ifc)
                    SimR1f(ik,it2,i+1,1) = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)

                    !Next periods's consumption
                    INTERP2D=cs_ret(ir,j,:,:,iaf,iuf,i,ifc)
                    dum4 = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    SimR1f(ik,it2,i,2) = dum4
                    SimR1f(ik,it2,i,3)= dum4*tc

                    ! Female wage, work hours and earnings

                    INTERP2D=ns_ret(ir,j,:,:,iaf,iuf,i,ifc)
                    SimR1f(ik,it2,i,4) = (1d0/(1d0+exp(kappa*(44+i-(agestart-1)))))*wage(2,a(2,iaf),dble(ix),u(2,iuf))/(1d0+t_employer)
                    dum5=bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    dum5=max(dum5,0d0)
                    if(dum5<(minhours/2d0)) then
                        dum5=0d0
                    else
                        dum5=max(minhours,dum5)
                    end if
                    SimR1f(ik,it2,i,5) = dum5
                    SimR1f(ik,it2,i,6) = dum5*(1d0/(1d0+exp(kappa*(44+i-(agestart-1)))))*wage(2,a(2,iaf),dble(ix),u(2,iuf))/(1d0+t_employer)
                    
                    !Household income and taxes

                    dum3=dum5*(1d0/(1d0+exp(kappa*(44+i-(agestart-1)))))*wage(2,a(2,iaf),dble(ix),u(2,iuf))/(1d0+t_employer)
                    SimR1f(ik,it2,i,7) = dum3
                    if(dum3>0d0) then
                        !SimR1f(ik,it2,i,8) = tax_labors(dum3)*dum3
                        !test1 = labor_tax_single(dum3)
                        SimR1f(ik,it2,i,8) = labor_tax_single(dum3)
                    else
                        SimR1f(ik,it2,i,8) = 0d0
                    end if

                    SimR1f(ik,it2,i,9)=dum3*t_employee+t_employer*dum3

                    !Individual marginal labor income tax rate including social security

                    if(SimR1f(ik,it2,i,7)>Deduct_Cutoff) then
                        SimR1f(ik,it2,i,15)=1d0-(1d0-(t_employee+t_employer+(1d0-(1d0-thetas(2))*thetas(1)*SimR1f(ik,it2,i,7)**(-thetas(2))))/(1d0*(1d0+t_employer)))*(1d0-t_const)
                    else if(SimR1f(ik,it2,i,7)>0d0) then
                        SimR1f(ik,it2,i,15)=1d0-(1d0-(t_employee+t_employer)/(1d0*(1d0+t_employer)))*(1d0-t_const)
                    else
                        SimR1f(ik,it2,i,15)= 0d0
                    end if

                    if(SimR1f(ik,it2,i,7)>0d0) then
                        SimR1f(ik,it2,i,17)=1d0-(1d0-(SimR1f(ik,it2,i,9)+SimR1f(ik,it2,i,8))/(SimR1f(ik,it2,i,7)*(1d0+t_employer)))*(1d0-t_const)
                    else
                        SimR1f(ik,it2,i,17)= 0d0
                    end if

                    !Experience
                    if(SimR1f(ik,it2,i,5)>1d-3) then
                        SimR1f(ik,it2,i+1,12)=SimR1f(ik,it2,i,12)+1d0
                    else
                        SimR1f(ik,it2,i+1,12)=SimR1f(ik,it2,i,12)*(1d0-deltaexp)
                    end if
                    
                    !Social security and budget constraint

                    
                    if(ir==2) then
                        !expR1f(ik,it2,i+1,3)=2
                        SimR1f(ik,it2,i,14)=0d0
                        !Budget Constraint
                        if(SimR1f(ik,it2,i,7)>0.0001d0) then
                            y=SimR1f(ik,it2,i,6)
                            !SimR1f(ik,it2,i,18)=SimR1f(ik,it2,i+1,1)-((dum2+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk))+lumpsum*0.5d0+y*(1d0-t_const)*(1d0-tax_labors(y)-tSS_employee(y))-SimR1f(ik,it2,i,2)*(1d0+tc))/(1d0+mu)
                            !test1=SimR1f(ik,it2,i+1,1)-((dum2+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk))+lumpsum*0.5d0+(1d0-t_const)*(after_tax_labor_inc_single(y)-y*tSS_employee(y))-SimR1f(ik,it2,i,2)*(1d0+tc))/(1d0+mu)
                            SimR1f(ik,it2,i,18)=SimR1f(ik,it2,i+1,1)-((dum2+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk))+lumpsum*0.5d0+(1d0-t_const)*(after_tax_labor_inc_single(y)-y*tSS_employee(y))-SimR1f(ik,it2,i,2)*(1d0+tc))/(1d0+mu)
                        else
                            SimR1f(ik,it2,i,18)=SimR1f(ik,it2,i+1,1)-((dum2+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk)) + SimR1f(ik,it2,i,14)+lumpsum*0.5d0-SimR1f(ik,it2,i,2)*(1d0+tc))/(1d0+mu)
                        end if
                            
                            
                    else
                        !expR1f(ik,it2,i+1,3)=1
                        SimR1f(ik,it2,i,14)=(psi0+psi1*av_earnings(2,2,iaf)*min(1d0,dble(ix)/35d0))*(1d0-t_const)
                        !Budget Constraint
                        SimR1f(ik,it2,i,18)=SimR1f(ik,it2,i+1,1)-((dum2+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk)) + SimR1f(ik,it2,i,14)+lumpsum*0.5d0-SimR1f(ik,it2,i,2)*(1d0+tc))/(1d0+mu)
                    end if

                    !Social welfare
                    INTERP2D=vs_ret(ir,j,:,:,iaf,iuf,i,ifc)
                    dum3 = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    SimR1f(ik,it2,i,11) = dum3
                    
                    !Ex-post social welfare
                    
                    if(dum5>0.001d0) then
                        SimR1f(ik,it2,i+1,19) = SimR1f(ik,it2,i,19) + (beta**(i-1d0))*(Uc(dum4)-chifs*(dum5**(1d0+etaf))/(1d0+etaf)-fc(2,ifc))
                    else
                        SimR1f(ik,it2,i+1,19) = SimR1f(ik,it2,i,19) + (beta**(i-1d0))*(Uc(dum4))
                    endif

                    !Euler error
                    if(i>2) then
                        dum8=SimR1f(ik,it2,i-1,2)
                        SimR1f(ik,it2,i,13) = dUc(dum8)-beta*OmegaRet(i)*((1d0+r*(1d0-tk))/(1d0+mu))*dUc(dum4)
                    end if

                    !Extra tax revenue
                    if(SimR1f(ik,it2,i,7)>0d0) then
                        dum3=SimR1f(ik,it2,i,7)
                        SimR1f(ik,it2,i,16)=t_const*(after_tax_labor_inc_married(dum3)-dum3*t_employee)
                    else
                        SimR1f(ik,it2,i,16)= 0d0
                    end if

                end if



            end do

            !do it2=1,nsim
            !    um=expR1m(ik,it2,i,2)
            !    if(Random3m(ik,it2,T+i)<trans_u(1,um,1)) then
            !        expR1m(ik,it2,i+1,2)=1
            !    elseif((trans_u(1,um,1)<Random3m(ik,it2,T+i)).AND.(Random3m(ik,it2,T+i)<(trans_u(1,um,1)+trans_u(1,um,2)))) then
            !        expR1m(ik,it2,i+1,2)=2
            !    elseif(((trans_u(1,um,1)+trans_u(1,um,2))<Random3m(ik,it2,T+i)).AND.(Random3m(ik,it2,T+i)<(trans_u(1,um,1)+trans_u(1,um,2)+trans_u(1,um,3)))) then
            !        expR1m(ik,it2,i+1,2)=3
            !    elseif(((trans_u(1,um,1)+trans_u(1,um,2)+trans_u(1,um,3))<Random3m(ik,it2,T+i)).AND.(Random3m(ik,it2,T+i)<(trans_u(1,um,1)+trans_u(1,um,2)+trans_u(1,um,3)+trans_u(1,um,4)))) then
            !        expR1m(ik,it2,i+1,2)=4
            !    else
            !        expR1m(ik,it2,i+1,2)=5
            !    end if
            !end do
            !
            !do it2=1,nsim
            !    uf=expR1f(ik,it2,i,2)
            !    if(Random3f(ik,it2,T+i)<trans_u(2,uf,1)) then
            !        expR1f(ik,it2,i+1,2)=1
            !    elseif((trans_u(2,uf,1)<Random3f(ik,it2,T+i)).AND.(Random3f(ik,it2,T+i)<(trans_u(2,uf,1)+trans_u(2,uf,2)))) then
            !        expR1f(ik,it2,i+1,2)=2
            !    elseif(((trans_u(2,uf,1)+trans_u(2,uf,2))<Random3f(ik,it2,T+i)).AND.(Random3f(ik,it2,T+i)<(trans_u(2,uf,1)+trans_u(2,uf,2)+trans_u(2,uf,3)))) then
            !        expR1f(ik,it2,i+1,2)=3
            !    elseif(((trans_u(2,uf,1)+trans_u(2,uf,2)+trans_u(2,uf,3))<Random3f(ik,it2,T+i)).AND.(Random3f(ik,it2,T+i)<(trans_u(2,uf,1)+trans_u(2,uf,2)+trans_u(2,uf,3)+trans_u(2,uf,4)))) then
            !        expR1f(ik,it2,i+1,2)=4
            !    else
            !        expR1f(ik,it2,i+1,2)=5
            !    end if
            !end do

        end do

    end subroutine Simulation

end module Simulations_mod
