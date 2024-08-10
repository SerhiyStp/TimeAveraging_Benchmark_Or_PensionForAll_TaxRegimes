module Summarize_Simulation
    
    implicit none
    
contains
    
    subroutine Statistics(file_id)

        !This subroutine computes aggregate statistics from the simulation

        use Model_Parameters
        use PolicyFunctions
        use Utilities
        use Simulations_mod

        implicit none

        integer :: i,country,ia,ia2,iu,ix,um,it2,ik,ifc, NVAR=2,it3,it4,ICOPT=2,ik2
        real(8) :: dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,dum10,dum11,dum12,dum13,dum14,dum15,dum16,dum17,SST,SSE,COV(2,2)
        real(8) :: dum18,dum19,dum20,dum21,dum22,dum23,dum24,dum25,r_ret,ss_tax,ss_expense,new_revenue,population_mass
        real(8), allocatable :: spousewage(:,:), spousewage2(:,:)
        real(8) :: test_corr
        real(8) :: dum2a(na,3)
        real(8) :: share_single(T, na+1)        
        integer :: file_id

        !Computing the weight of each generation

        WeightActive(1)=1d0
        do i=2,T
            WeightActive(i)=WeightActive(i-1)*OmegaActive(i-1)
        end do

        WeightRet(1)=WeightActive(T)*OmegaActive(T)
        do i=2,Tret
            WeightRet(i)=WeightRet(i-1)*OmegaRet(i-1)
        end do

        dum2=0d0
        dum3=0d0

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    if(SimR1m(it2,it,i,5)>1d-3) then
                        dum2=dum2+SimR1m(it2,it,i,6)*WeightRet(i)
                        dum3=dum3+1d0*WeightRet(i)
                    end if
                end do
            end do

        end do


        Print *,'Average male earnings above 65 is',dum2/dum3


        dum2=0d0
        dum3=0d0

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    if(SimR1f(it2,it,i,5)>1d-3) then
                        dum2=dum2+SimR1f(it2,it,i,6)*WeightRet(i)
                        dum3=dum3+1d0*WeightRet(i)
                    end if
                end do
            end do

        end do

        Print *,'Average female earnings above 65 is',dum2/dum3


        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0

        !Population Labor Supply

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+Sim1m(it2,it,i,4)*WeightActive(i)
                    dum2=dum2+Sim1f(it2,it,i,4)*WeightActive(i)
                    dum3=dum3+2d0*WeightActive(i)
                end do
            end do

        end do

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+SimR1m(it2,it,i,5)*WeightRet(i)
                    dum2=dum2+SimR1f(it2,it,i,5)*WeightRet(i)
                    dum3=dum3+2d0*WeightRet(i)
                end do
            end do

        end do

        population_mass=dum3
        dum2=dum2/dum3

        Print *,'Population Labor Supply is',dum2

        dum4=dum2

        dum2=0d0
        dum3=0d0

        !Population Employment Rate

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,i,4)>1d-3) then
                        dum2=dum2+1d0*WeightActive(i)
                    end if
                    if(Sim1f(it2,it,i,4)>1d-3) then
                        dum2=dum2+1d0*WeightActive(i)
                    end if
                    dum3=dum3+2d0*WeightActive(i)
                end do
            end do

        end do

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    if(SimR1m(it2,it,i,5)>1d-3) then
                        dum2=dum2+1d0*WeightRet(i)
                    end if
                    if(SimR1f(it2,it,i,5)>1d-3) then
                        dum2=dum2+1d0*WeightRet(i)
                    end if
                    dum3=dum3+2d0*WeightRet(i)
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Population Employment rate is',dum2

        Print *,'Population Intensive margin is',dum4/dum2


        dum2=0d0
        dum3=0d0

        !Employment Rate 65+


        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    if(SimR1m(it2,it,i,5)>1d-3) then
                        dum2=dum2+1d0*WeightRet(i)
                    end if
                    if(SimR1f(it2,it,i,5)>1d-3) then
                        dum2=dum2+1d0*WeightRet(i)
                    end if
                    dum3=dum3+2d0*WeightRet(i)
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Employment rate 65+ is',dum2  

        dum2=0d0
        dum3=0d0

        !Employment Rate Below 65

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,i,4)>1d-3) then
                        dum2=dum2+1d0*WeightActive(i)
                    end if
                    if(Sim1f(it2,it,i,4)>1d-3) then
                        dum2=dum2+1d0*WeightActive(i)
                    end if
                    dum3=dum3+2d0*WeightActive(i)
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Employment rate below 65 is',dum2

        !Labor Supply before 65

        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0


        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+Sim1m(it2,it,i,4)*WeightActive(i)
                    dum2=dum2+Sim1f(it2,it,i,4)*WeightActive(i)
                    dum3=dum3+2d0*WeightActive(i)
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Labor supply below 65 is',dum2

        dum2=0d0
        dum3=0d0

        !Labor supply of males 65 and older

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+SimR1m(it2,it,i,5)*WeightRet(i)
                    dum3=dum3+WeightRet(i)
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Labor supply of males 65+ is',dum2

        dum2=0d0
        dum3=0d0

        !Male labor force participation after 65

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    if(SimR1m(it2,it,i,5)>1d-3) then
                        dum2=dum2+(1d0)*WeightRet(i)
                    end if
                    dum3=dum3+1d0*WeightRet(i)
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'LFP of males 65+ is',dum2

        !Single male labor force participation at age 65

        dum2=0d0
        dum3=0d0

        i=1

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

        !dum10=dum10+((dum2-0.2)/0.2)**2

        Print *,'LFP of single males at age 65 is',dum2

        !Single male labor force participation after 65

        dum2=0d0
        dum3=0d0

        do i=1,Tret

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

        end do

        dum2=dum2/dum3

        Print *,'LFP of single males 65+ is',dum2

        !Single male labor supply after 65

        dum2=0d0
        dum3=0d0

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,T,10)<0.5) then
                        dum2=dum2+SimR1m(it2,it,i,5)*WeightRet(i)
                        dum3=dum3+1d0*WeightRet(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Labor supply of single males 65+ is',dum2

        !Single male labor force participation by age after 65


        !do i=1,Tret
        !
        !dum2=0d0
        !dum3=0d0
        !    
        !do it2=1,nsim2
        !do it=1,nsim
        !    if(Sim1m(it2,it,T,10)<0.5) then
        !        if(SimR1m(it2,it,i,5)>1d-3) then
        !        dum2=dum2+(1d0)*WeightRet(i)
        !    end if
        !    dum3=dum3+1d0*WeightRet(i)
        !    end if
        !end do
        !end do
        !
        !
        !
        !dum2=dum2/dum3
        !
        !Print *,'LFP of single males at age', 64+i, dum2
        !
        !end do

        !Married male labor force participation at age 65

        dum2=0d0
        dum3=0d0

        i=1

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

        !dum10=dum10+((dum2-0.2)/0.2)**2

        Print *,'LFP of married males at age 65 is',dum2

        !Married male labor force participation after 65

        dum2=0d0
        dum3=0d0

        do i=1,Tret

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

        end do

        dum2=dum2/dum3

        Print *,'LFP of married males 65+ is',dum2

        !Married male labor supply after 65

        dum2=0d0
        dum3=0d0

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,T,10)>0.5) then
                        dum2=dum2+SimR1m(it2,it,i,5)*WeightRet(i)
                        dum3=dum3+1d0*WeightRet(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Labor supply of married males 65+ is',dum2

        !!Married male labor force participation by age after 65
        !
        !do i=1,Tret
        !
        !dum2=0d0
        !dum3=0d0    
        !    
        !do it2=1,nsim2
        !do it=1,nsim
        !    if(Sim1m(it2,it,T,10)>0.5) then
        !        if(SimR1m(it2,it,i,5)>1d-3) then
        !        dum2=dum2+(1d0)*WeightRet(i)
        !    end if
        !    dum3=dum3+1d0*WeightRet(i)
        !    end if
        !end do
        !end do
        !
        !dum2=dum2/dum3
        !
        !Print *,'LFP of married males at age', 64+i, dum2
        !
        !end do


        !!Male retirement age
        !
        !dum2=0d0
        !dum3=0d0
        !
        !do i=1,Tret
        !
        !    do it2=1,nsim2
        !        do it=1,nsim
        !            if((expR1m(it2,it,i,3)==2).AND.(expR1m(it2,it,i+1,3)==1)) then
        !                dum2=dum2+(1d0)*i*WeightRet(i)
        !                dum3=dum3+1d0*WeightRet(i)
        !            end if
        !        end do
        !    end do
        !
        !end do
        !
        !dum2=dum2/dum3
        !
        !Print *,'Male retirement age is',dum2+64d0

        !!Single male retirement age
        !
        !dum2=0d0
        !dum3=0d0
        !
        !do i=1,Tret
        !
        !    do it2=1,nsim2
        !        do it=1,nsim
        !            if(Sim1m(it2,it,T,10)<0.5) then
        !                if((expR1m(it2,it,i,3)==2).AND.(expR1m(it2,it,i+1,3)==1)) then
        !                    dum2=dum2+(1d0)*i*WeightRet(i)
        !                    dum3=dum3+1d0*WeightRet(i)
        !                end if
        !            end if
        !        end do
        !    end do
        !
        !end do
        !
        !dum2=dum2/dum3
        !
        !Print *,'Single male retirement age is',dum2+65d0

        !!Single male retirement by age
        !
        !do i=1,Tret
        !
        !dum2=0d0
        !dum3=0d0
        !    
        !do it2=1,nsim2
        !do it=1,nsim
        !    if(Sim1m(it2,it,T,10)<0.5) then
        !        dum3=dum3+1d0*WeightRet(i)
        !        if((expR1m(it2,it,i,3)==2).AND.(SimR1m(it2,it,i,5)<1d-3)) then
        !            dum2=dum2+(1d0)*WeightRet(i)
        !        end if
        !    end if
        !end do
        !end do
        !
        !dum2=dum2/dum3
        !
        !Print *,'Single male retirement at age',i+64,dum2
        !
        !end do

        !!Married male retirement age
        !
        !dum2=0d0
        !dum3=0d0
        !
        !do i=1,Tret
        !
        !    do it2=1,nsim2
        !        do it=1,nsim
        !            if(Sim1m(it2,it,T,10)>0.5) then
        !                if((expR1m(it2,it,i,3)==2).AND.(expR1m(it2,it,i+1,3)==1)) then
        !                    dum2=dum2+(1d0)*i*WeightRet(i)
        !                    dum3=dum3+1d0*WeightRet(i)
        !                end if
        !            end if
        !        end do
        !    end do
        !
        !end do
        !
        !dum2=dum2/dum3
        !
        !Print *,'Married male retirement age is',dum2+65d0

        !Married male retirement by age

        !do i=1,Tret
        !
        !dum2=0d0
        !dum3=0d0
        !    
        !do it2=1,nsim2
        !do it=1,nsim
        !    if(Sim1m(it2,it,T,10)>0.5) then
        !        dum3=dum3+1d0*WeightRet(i)
        !        if((expR1m(it2,it,i,3)==2).AND.(SimR1m(it2,it,i,5)<1d-3)) then
        !            dum2=dum2+(1d0)*WeightRet(i)
        !        end if
        !    end if
        !end do
        !end do
        !
        !dum2=dum2/dum3
        !
        !Print *,'Single male retirement at age',i+64,dum2
        !
        !end do



        !Labor supply of females 65 and older

        dum2=0d0
        dum3=0d0

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+SimR1f(it2,it,i,5)*WeightRet(i)
                    dum3=dum3+WeightRet(i)
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Labor supply of females 65+ is',dum2

        dum2=0d0
        dum3=0d0

        !Female labor force participation after 65

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    if(SimR1f(it2,it,i,5)>1d-3) then
                        dum2=dum2+(1d0)*WeightRet(i)
                    end if
                    dum3=dum3+1d0*WeightRet(i)
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'LFP of females 65+ is',dum2

        !Single female labor force participation after 65

        dum2=0d0
        dum3=0d0

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,T,10)<0.5) then
                        if(SimR1f(it2,it,i,5)>1d-3) then
                            dum2=dum2+(1d0)*WeightRet(i)
                        end if
                        dum3=dum3+1d0*WeightRet(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'LFP of single females 65+ is',dum2

        !Single female labor supply after 65

        dum2=0d0
        dum3=0d0

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,T,10)<0.5) then
                        dum2=dum2+SimR1f(it2,it,i,5)*WeightRet(i)
                        dum3=dum3+1d0*WeightRet(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Labor supply of single females 65+ is',dum2

        !Single female labor force participation by age after 65


        !do i=1,Tret
        !
        !dum2=0d0
        !dum3=0d0
        !    
        !do it2=1,nsim2
        !do it=1,nsim
        !    if(Sim1f(it2,it,T,10)<0.5) then
        !        if(SimR1f(it2,it,i,5)>1d-3) then
        !        dum2=dum2+(1d0)*WeightRet(i)
        !    end if
        !    dum3=dum3+1d0*WeightRet(i)
        !    end if
        !end do
        !end do
        !
        !
        !
        !dum2=dum2/dum3
        !
        !Print *,'LFP of single females at age', 64+i, dum2
        !
        !end do

        !Married female labor force participation after 65

        dum2=0d0
        dum3=0d0

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,T,10)>0.5) then
                        if(SimR1f(it2,it,i,5)>1d-3) then
                            dum2=dum2+(1d0)*WeightRet(i)
                        end if
                        dum3=dum3+1d0*WeightRet(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'LFP of married females 65+ is',dum2

        !Married female labor supply after 65

        dum2=0d0
        dum3=0d0

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,T,10)>0.5) then
                        dum2=dum2+SimR1f(it2,it,i,5)*WeightRet(i)
                        dum3=dum3+1d0*WeightRet(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Labor supply of married females 65+ is',dum2

        !!Married female labor force participation by age after 65
        !
        !do i=1,Tret
        !
        !dum2=0d0
        !dum3=0d0    
        !    
        !do it2=1,nsim2
        !do it=1,nsim
        !    if(Sim1f(it2,it,T,10)>0.5) then
        !        if(SimR1f(it2,it,i,5)>1d-3) then
        !        dum2=dum2+(1d0)*WeightRet(i)
        !    end if
        !    dum3=dum3+1d0*WeightRet(i)
        !    end if
        !end do
        !end do
        !
        !dum2=dum2/dum3
        !
        !Print *,'LFP of married females at age', 64+i, dum2
        !
        !end do


        !!Female retirement age
        !
        !dum2=0d0
        !dum3=0d0
        !
        !do i=1,Tret
        !
        !    do it2=1,nsim2
        !        do it=1,nsim
        !            if((expR1f(it2,it,i,3)==2).AND.(expR1f(it2,it,i+1,3)==1)) then
        !                dum2=dum2+(1d0)*i*WeightRet(i)
        !                dum3=dum3+1d0*WeightRet(i)
        !            end if
        !        end do
        !    end do
        !
        !end do
        !
        !dum2=dum2/dum3
        !
        !Print *,'Female retirement age is',dum2+65d0

        !!Single female retirement age
        !
        !dum2=0d0
        !dum3=0d0
        !
        !do i=1,Tret
        !
        !    do it2=1,nsim2
        !        do it=1,nsim
        !            if(Sim1f(it2,it,T,10)<0.5) then
        !                if((expR1f(it2,it,i,3)==2).AND.(expR1f(it2,it,i+1,3)==1)) then
        !                    dum2=dum2+(1d0)*i*WeightRet(i)
        !                    dum3=dum3+1d0*WeightRet(i)
        !                end if
        !            end if
        !        end do
        !    end do
        !
        !end do
        !
        !dum2=dum2/dum3
        !
        !Print *,'Single female retirement age is',dum2+65d0

        !!Single female retirement by age
        !
        !do i=1,Tret
        !
        !dum2=0d0
        !dum3=0d0
        !    
        !do it2=1,nsim2
        !do it=1,nsim
        !    if(Sim1f(it2,it,T,10)<0.5) then
        !        dum3=dum3+1d0*WeightRet(i)
        !        if((expR1f(it2,it,i,3)==2).AND.(SimR1f(it2,it,i,5)<1d-3)) then
        !            dum2=dum2+(1d0)*WeightRet(i)
        !        end if
        !    end if
        !end do
        !end do
        !
        !dum2=dum2/dum3
        !
        !Print *,'Single female retirement at age',i+64,dum2
        !
        !end do

        !!Married female retirement age
        !
        !dum2=0d0
        !dum3=0d0
        !
        !do i=1,Tret
        !
        !    do it2=1,nsim2
        !        do it=1,nsim
        !            if(Sim1f(it2,it,T,10)>0.5) then
        !                if((expR1f(it2,it,i,3)==2).AND.(expR1f(it2,it,i+1,3)==1)) then
        !                    dum2=dum2+(1d0)*i*WeightRet(i)
        !                    dum3=dum3+1d0*WeightRet(i)
        !                end if
        !            end if
        !        end do
        !    end do
        !
        !end do
        !
        !dum2=dum2/dum3
        !
        !Print *,'Married female retirement age is',dum2+65d0

        !Married female retirement by age

        !do i=1,Tret
        !
        !dum2=0d0
        !dum3=0d0
        !    
        !do it2=1,nsim2
        !do it=1,nsim
        !    if(Sim1f(it2,it,T,10)>0.5) then
        !        dum3=dum3+1d0*WeightRet(i)
        !        if((expR1f(it2,it,i,3)==2).AND.(SimR1f(it2,it,i,5)<1d-3)) then
        !            dum2=dum2+(1d0)*WeightRet(i)
        !        end if
        !    end if
        !end do
        !end do
        !
        !dum2=dum2/dum3
        !
        !Print *,'Single female retirement at age',i+64,dum2
        !
        !end do

        !Welfare of single women after 65

        !dum2=0d0
        !dum3=0d0
        !dum4=0d0
        !dum5=0d0
        !
        !!do i=1,Tret
        !
        !i=1
        !
        !do it2=1,nsim2
        !do it=1,nsim
        !    if(Sim1f(it2,it,T,10)<0.5) then
        !        dum2=dum2+SimR1f(it2,it,i,11)*WeightRet(i)
        !        dum3=dum3+1d0*WeightRet(i)
        !    end if
        !end do
        !end do
        !
        !!end do
        !
        !dum2=dum2/dum3
        !
        !Print *,'Average welfare of single women +65 is',dum2

        !variance of welfare

        !do i=1,Tret

        !do it2=1,nsim2
        !do it=1,nsim
        !    if(Sim1f(it2,it,T,10)<0.5) then
        !        dum4=dum4+((SimR1f(it2,it,i,11)-dum2)**(2d0))*WeightRet(i)
        !    end if
        !end do
        !end do

        !end do

        !dum4=dum4/dum3
        !
        !Print *,'Variance of welfare single women +65 is',dum4
        !
        !Print *,'Coefficient of variation welfare single women +65 is',SQRT(dum4)/dum2

        !Consumption of married women after 65

        !dum2=0d0
        !dum3=0d0
        !dum4=0d0
        !dum5=0d0

        !do i=1,Tret

        !do it2=1,nsim2
        !do it=1,nsim
        !    if(Sim1f(it2,it,T,10)>0.5) then
        !        dum2=dum2+SimR1f(it2,it,i,11)*WeightRet(i)
        !        dum3=dum3+1d0*WeightRet(i)
        !    end if
        !end do
        !end do

        !end do

        !dum2=dum2/dum3
        !
        !Print *,'Average welfare of narried women +65 is',dum2

        !variance of consumption

        !do i=1,Tret

        !do it2=1,nsim2
        !do it=1,nsim
        !    if(Sim1f(it2,it,T,10)>0.5) then
        !        dum4=dum4+((SimR1f(it2,it,i,11)-dum2)**(2d0))*WeightRet(i)
        !    end if
        !end do
        !end do

        !end do

        !dum4=dum4/dum3
        !
        !Print *,'Variance of welfare married women +65 is',dum4
        !
        !Print *,'Coefficient of variation of welfare married women +65 is',SQRT(dum4)/dum2


        dum2=0d0
        dum3=0d0

        !Male Labor Supply

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+Sim1m(it2,it,i,4)*WeightActive(i)
                    dum3=dum3+1d0*WeightActive(i)
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Male labor supply is',dum2

        !Single Male Labor Supply

        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,i,10)<0.5) then
                        dum2=dum2+Sim1m(it2,it,i,4)*WeightActive(i)
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Single Male labor supply is',dum2

        dum10=((dum2-0.281917216)/0.281917216)**2

        !Single male labor force participation

        dum2=0d0
        dum3=0d0

        do i=1,T

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

        end do

        dum2=dum2/dum3

        dum10=dum10+((dum2-0.799)/0.799)**2

        Print *,'Single male labor force participation is',dum2

        !YVAR=sqrt(-1.0)
        !XVARS=sqrt(-1.0)

        !do i=2,T

        !do it2=1,nsim2
        !do it=1,nsim
        !if(Sim1m(it2,it,i,10)<0.5) then
        !if(Sim1m(it2,it,i-1,4)>1d-3) then
        !XVARs((i-1)*it*it2+(it2-1)*it+it,1)=1d0
        !else
        !XVARs((i-1)*it*it2+(it2-1)*it+it,1)=0d0
        !end if    
        !if(Sim1m(it2,it,i,4)>1d-3) then
        !YVAR((i-1)*it*it2+(it2-1)*it+it)=1d0
        !else
        !YVAR((i-1)*it*it2+(it2-1)*it+it)=0d0
        !end if      
        !end if
        !end do
        !end do

        !end do

        !CALL RLSE (YVAR, XVARS, BREG, SST=SST, SSE=SSE)
        !
        !dum2=1d0-SSE/SST

        !dum10=dum10+((dum2-0.4081)/0.4081)**2

        !Print *,'BREG is',BREG

        !Print *,'Single male LFP R2 is',dum2


        !Married Male Labor Supply

        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,i,10)>0.5) then
                        dum2=dum2+Sim1m(it2,it,i,4)*WeightActive(i)
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Married Male labor supply is',dum2

        dum10=dum10+((dum2-0.359934432)/0.359934432)**2


        !Married male labor force participation

        dum2=0d0
        dum3=0d0

        do i=1,T

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

        end do

        dum2=dum2/dum3

        dum10=dum10+((dum2-0.879)/0.879)**2

        Print *,'Married male labor force participation is',dum2

        !YVAR=sqrt(-1.0)
        !XVARS=sqrt(-1.0)

        !do i=2,T

        !do it2=1,nsim2
        !do it=1,nsim
        !if(Sim1m(it2,it,i,10)>0.5) then
        !if(Sim1m(it2,it,i-1,4)>1d-3) then
        !XVARs((i-1)*it*it2+(it2-1)*it+it,1)=1d0
        !else
        !XVARs((i-1)*it*it2+(it2-1)*it+it,1)=0d0
        !end if    
        !if(Sim1m(it2,it,i,4)>1d-3) then
        !YVAR((i-1)*it*it2+(it2-1)*it+it)=1d0
        !else
        !YVAR((i-1)*it*it2+(it2-1)*it+it)=0d0
        !end if      
        !end if
        !end do
        !end do

        !end do

        !CALL RLSE (YVAR, XVARS, BREG, SST=SST, SSE=SSE)
        !
        !dum2=1d0-SSE/SST

        !dum10=dum10+((dum2-0.4573)/0.4573)**2

        !Print *,'BREG is',BREG

        !Print *,'Married male LFP R2 is',dum2


        !Female labor supply

        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+Sim1f(it2,it,i,4)*WeightActive(i)
                    dum3=dum3+1d0*WeightActive(i)
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Female labor supply is',dum2

        !Single female labor supply

        dum2=0d0
        dum3=0d0
        dum4=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,i,10)<0.5) then
                        dum2=dum2+Sim1f(it2,it,i,4)*WeightActive(i)
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Single female labor supply is',dum2

        dum10=dum10+((dum2-0.25116337)/0.25116337)**2

        !Variance of single female hours

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,i,10)<0.5) then
                        dum4=dum4+((Sim1f(it2,it,i,4)-dum2)**2)*WeightActive(i)
                    end if
                end do
            end do

        end do

        dum4=dum4/dum3

        Print *,'Stdev single female labor supply is',SQRT(dum4)


        !!Variance of single female hours
        !do i=10,10
        !do it2=1,nsim2
        !do it=1,nsim
        !if(Sim1f(it2,it,i,10)<0.5) then
        !dum4=dum4+((Sim1f(it2,it,i,4)-dum2)**2)*WeightActive(i)
        !end if
        !end do
        !end do
        !end do

        !dum4=dum4/dum3

        !Print *,'Stdev single female labor supply at age 30 is',SQRT(dum4)

        !Married female labor supply

        dum2=0d0
        dum3=0d0
        dum4=0d0

        do i=1,T
            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,i,10)>0.5) then
                        dum2=dum2+Sim1f(it2,it,i,4)*WeightActive(i)
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end do
            end do
        end do

        dum2=dum2/dum3

        Print *,'Married female labor supply is',dum2

        dum10=dum10+((dum2-0.224398901)/0.224398901)**2

        !Variance of married female hours

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,i,10)>0.5) then
                        dum4=dum4+((Sim1f(it2,it,i,4)-dum2)**2)*WeightActive(i)
                    end if
                end do
            end do

        end do

        dum4=dum4/dum3

        Print *,'Stdev married female labor supply is',SQRT(dum4)


        !Female labor force participation

        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,i,4)>1d-3) then
                        dum2=dum2+(1d0)*WeightActive(i)
                    end if
                    dum3=dum3+1d0*WeightActive(i)
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Female labor force participation is',dum2

        !Single female labor force participation

        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,i,10)<0.5) then
                        if(Sim1f(it2,it,i,4)>1d-3) then
                            dum2=dum2+(1d0)*WeightActive(i)
                        end if
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Single female labor force participation is',dum2

        dum10=dum10+((dum2-0.7600114)/0.7600114)**2

        !YVAR=sqrt(-1.0)
        !XVARS=sqrt(-1.0)

        !do i=2,T

        !do it2=1,nsim2
        !do it=1,nsim
        !if(Sim1f(it2,it,i,10)<0.5) then
        !if(Sim1f(it2,it,i-1,4)>1d-3) then
        !XVARs((i-1)*it*it2+(it2-1)*it+it,1)=1d0
        !else
        !XVARs((i-1)*it*it2+(it2-1)*it+it,1)=0d0
        !end if    
        !if(Sim1f(it2,it,i,4)>1d-3) then
        !YVAR((i-1)*it*it2+(it2-1)*it+it)=1d0
        !else
        !YVAR((i-1)*it*it2+(it2-1)*it+it)=0d0
        !end if      
        !end if
        !end do
        !end do

        !end do

        !CALL RLSE (YVAR, XVARS, BREG, SST=SST, SSE=SSE)
        !
        !dum2=1d0-SSE/SST
        !
        !
        !Print *,'Single female LFP R2 is',dum2


        !Married female labor force participation

        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,i,10)>0.5) then
                        if(Sim1f(it2,it,i,4)>1d-3) then
                            dum2=dum2+(1d0)*WeightActive(i)
                        end if
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Married female labor force participation is',dum2

        dum10=dum10+((dum2-0.6755489)/0.6755489)**2

        !YVAR=sqrt(-1.0)
        !XVARS=sqrt(-1.0)

        !do i=2,T

        !do it2=1,nsim2
        !do it=1,nsim
        !if(Sim1f(it2,it,i,10)>0.5) then
        !if(Sim1f(it2,it,i-1,4)>1d-3) then
        !XVARS((i-1)*it*it2+(it2-1)*it+it,1)=1d0
        !else
        !XVARS((i-1)*it*it2+(it2-1)*it+it,1)=0d0
        !end if      
        !if(Sim1f(it2,it,i,4)>1d-3) then
        !YVAR((i-1)*it*it2+(it2-1)*it+it)=1d0
        !else
        !YVAR((i-1)*it*it2+(it2-1)*it+it)=0d0
        !end if      
        !end if
        !end do
        !end do

        !end do

        !CALL RLSE (YVAR, XVARS, BREG, SST=SST, SSE=SSE)
        !
        !dum2=1d0-SSE/SST

        !dum10=dum10+((dum2-0.5532)/0.5532)**2

        !Print *,'BREG is',BREG

        !Print *,'Married female LFP R2 is',dum2

        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,i,4)>1d-3) then
                        dum2=dum2+Sim1f(it2,it,i,4)*WeightActive(i)
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Female intensive margin is',dum2

        !Variance of log male earnings

        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,i,5)>0d0) then
                        dum2=dum2+log(Sim1m(it2,it,i,5))*WeightActive(i)
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        do i=1,T

            do it2=1,nsim2    
                do it=1,nsim
                    if(Sim1m(it2,it,i,5)>0d0) then
                        dum4=dum4+((log(Sim1m(it2,it,i,5))-dum2)**2)*WeightActive(i)
                    end if
                end do
            end do

        end do

        dum4=dum4/dum3

        Print *,'Variance of log male earnings is',dum4

        dum2=0d0
        dum3=0d0
        dum4=0d0

        !Variance of log female earnings

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,i,5)>1d-3) then
                        dum2=dum2+log(Sim1f(it2,it,i,5))*WeightActive(i)
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,i,5)>1d-3) then
                        dum4=dum4+((log(Sim1f(it2,it,i,5))-dum2)**2)*WeightActive(i)
                    end if
                end do
            end do

        end do

        dum4=dum4/dum3

        Print *,'Variance of log female earnings is',dum4

        !Variance of log male wage

        dum2=0d0
        dum3=0d0
        dum4=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+log(Sim1m(it2,it,i,3))*WeightActive(i)
                    dum3=dum3+1d0*WeightActive(i)
                end do
            end do

        end do

        dum2=dum2/dum3

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    dum4=dum4+((log(Sim1m(it2,it,i,3))-dum2)**2)*WeightActive(i)
                end do
            end do

        end do

        dum4=dum4/dum3

        Print *,'Variance of log male wage is',dum4

        dum2=0d0
        dum3=0d0
        dum4=0d0

        !Variance of log female wage

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+log(Sim1f(it2,it,i,3))*WeightActive(i)
                    dum3=dum3+1d0*WeightActive(i)
                end do
            end do

        end do

        dum2=dum2/dum3

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    dum4=dum4+((log(Sim1f(it2,it,i,3))-dum2)**2)*WeightActive(i)
                end do
            end do

        end do

        dum4=dum4/dum3

        Print *,'Variance of log female wage is',dum4

        !Male experience

        dum2=0d0
        dum3=0d0
        dum4=0d0

        do it2=1,nsim2
            do it=1,nsim
                dum2=dum2+1d0*(exp2m(it2,it,T+1,1)-1)
                dum3=dum3+1d0
            end do
        end do

        dum2=dum2/dum3

        Print *,'Male experience at age 65 is',dum2

        !Female experience

        dum2=0d0
        dum3=0d0
        dum4=0d0

        do it2=1,nsim2
            do it=1,nsim
                dum2=dum2+1d0*(exp2f(it2,it,T+1,1)-1)
                dum3=dum3+1d0
            end do
        end do

        dum2=dum2/dum3

        Print *,'Female experience at age 65 is',dum2

        !Male experience at age 99

        dum2=0d0
        dum3=0d0
        dum4=0d0

        do it2=1,nsim2
            do it=1,nsim
                dum2=dum2+1d0*(SimR1m(it2,it,Tret-1,12)-1)
                dum3=dum3+1d0
            end do
        end do

        dum2=dum2/dum3

        Print *,'Male experience at age 99 is',dum2    

        !Female experience at age 99

        dum2=0d0
        dum3=0d0
        dum4=0d0

        do it2=1,nsim2
            do it=1,nsim
                dum2=dum2+1d0*(SimR1f(it2,it,Tret-1,12)-1)
                dum3=dum3+1d0
            end do
        end do

        dum2=dum2/dum3

        Print *,'Female experience at age 99 is',dum2    


        !Male wage

        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,i,4)>1d-3) then
                        dum2=dum2+Sim1m(it2,it,i,3)*WeightActive(i)
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Average male wage is',dum2

        dum7=dum2

        !Female wage

        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,i,4)>1d-3) then
                        dum2=dum2+Sim1f(it2,it,i,3)*WeightActive(i)
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Average female wage is',dum2

        Print *,'Male_wage/Female_wage is',dum7/dum2


        !Wages at 0.5AE and 2AE

        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim

                    if((Sim1m(it2,it,i,5)>0.4d0).AND.(Sim1m(it2,it,i,5)<0.6d0)) then
                        dum2=dum2+Sim1m(it2,it,i,3)*WeightActive(i)
                        dum3=dum3+1d0*WeightActive(i)
                    end if

                    if((Sim1f(it2,it,i,5)>0.4d0).AND.(Sim1f(it2,it,i,5)<0.6d0)) then
                        dum2=dum2+Sim1f(it2,it,i,3)*WeightActive(i)
                        dum3=dum3+1d0*WeightActive(i)
                    end if

                    if((Sim1m(it2,it,i,5)>1.9d0).AND.(Sim1m(it2,it,i,5)<2.1d0)) then
                        dum4=dum4+Sim1m(it2,it,i,3)*WeightActive(i)
                        dum5=dum5+1d0*WeightActive(i)
                    end if

                    if((Sim1f(it2,it,i,5)>1.9d0).AND.(Sim1f(it2,it,i,5)<2.1d0)) then
                        dum4=dum4+Sim1f(it2,it,i,3)*WeightActive(i)
                        dum5=dum5+1d0*WeightActive(i)
                    end if

                end do
            end do

        end do

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim

                    if((SimR1m(it2,it,i,6)>0.4d0).AND.(SimR1m(it2,it,i,6)<0.6d0)) then
                        dum2=dum2+SimR1m(it2,it,i,4)*WeightRet(i)
                        dum3=dum3+1d0*WeightRet(i)
                    end if

                    if((SimR1f(it2,it,i,6)>0.4d0).AND.(SimR1f(it2,it,i,6)<0.6d0)) then
                        dum2=dum2+SimR1f(it2,it,i,4)*WeightRet(i)
                        dum3=dum3+1d0*WeightRet(i)
                    end if

                    if((SimR1m(it2,it,i,6)>1.9d0).AND.(SimR1m(it2,it,i,6)<2.1d0)) then
                        dum4=dum4+SimR1m(it2,it,i,4)*WeightRet(i)
                        dum5=dum5+1d0*WeightRet(i)
                    end if

                    if((SimR1f(it2,it,i,6)>1.9d0).AND.(SimR1f(it2,it,i,6)<2.1d0)) then
                        dum4=dum4+SimR1f(it2,it,i,4)*WeightRet(i)
                        dum5=dum5+1d0*WeightRet(i)
                    end if

                end do
            end do

        end do

        Print *,'Wage rate of people making 0.5AE is',dum2/dum3
        Print *,'Wage rate of people making 2AE is',dum4/dum5

        !Male earnings

        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,i,4)>1d-3) then
                        dum2=dum2+Sim1m(it2,it,i,5)*WeightActive(i)
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Average married male earnings is',dum2

        dum7=dum2

        !Female earnings

        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,i,4)>1d-3) then
                        dum2=dum2+Sim1f(it2,it,i,5)*WeightActive(i)
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        dum7=dum7/dum2

        Print *,'Average married female earnings is',dum2

        Print *,'Male_earnings/Female_earnings is',dum7

        dum10=dum10+((dum7-1.5689894)/1.5689894)**2

        !Correlation in spousal ability
        it4=0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,i,10)>0.5) then
                        !if((Sim1m(it2,it,i,4)>0.001).AND.(Sim1f(it3,it,i,4)>0.001d0)) then
                        it4=it4+1
                        !end if
                    end if
                end do
            end do

        end do

        allocate(spousewage(it4,2))

        it4=0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,i,10)>0.5) then

                        it3=exp1m(it2,it,i,4)
                        !if((Sim1m(it2,it,i,4)>0.001).AND.(Sim1f(it3,it,i,4)>0.001d0)) then
                        it4=it4+1
                        Spousewage(it4,1)=A(1,exp1m(it2,it,i,2))
                        Spousewage(it4,2)=A(2,exp1f(it2,it3,i,2))
                        !end if
                    end if
                end do
            end do

        end do

        !CALL D_CORVC(NVAR, Spousewage, COV, ICOPT=ICOPT)
        test_corr = corr_xy(Spousewage(:,1), Spousewage(:,2), it4)
        !Print *,'Correlation of spousal ability is',COV(1,2)
        Print *, 'Correlation of spousal ability is', test_corr

        dum10=dum10+((test_corr-0.646)/0.646)**2

        it4=0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,i,10)>0.5) then

                        it3=exp1m(it2,it,i,4)
                        if((Sim1m(it2,it,i,4)>0.001).AND.(Sim1f(it2,it3,i,4)>0.001d0)) then
                            it4=it4+1
                        end if
                    end if
                end do
            end do

        end do

        allocate(spousewage2(it4,2))

        it4=0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,i,10)>0.5) then

                        it3=exp1m(it2,it,i,4)
                        if((Sim1m(it2,it,i,4)>0.001).AND.(Sim1f(it2,it3,i,4)>0.001d0)) then
                            it4=it4+1
                            Spousewage2(it4,1)=Sim1m(it2,it,i,3)
                            Spousewage2(it4,2)=Sim1f(it2,it3,i,3)
                        end if
                    end if
                end do
            end do

        end do

        !CALL D_CORVC(NVAR, Spousewage2, COV, ICOPT=ICOPT)
        test_corr = corr_xy(Spousewage2(:,1), Spousewage2(:,2), it4)

        !Print *,'Correlation of spousal wages is',COV(1,2)
        Print *,'Correlation of spousal wages is', test_corr


        !dum10=dum10+((test_corr-0.4070)/0.4070)**2

        !Total labor income taxes and tax revenue

        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0
        dum7=0d0
        dum15=0d0
        dum16=0d0
        dum17=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    dum5=dum5+2d0*WeightActive(i)
                    dum2=dum2+Sim1m(it2,it,i,6)*(1d0+t_employer)*WeightActive(i)
                    dum7=dum7+Sim1m(it2,it,i,7)*WeightActive(i)
                    dum15=dum15+Sim1m(it2,it,i,9)*WeightActive(i)
                    dum4=dum4+Sim1m(it2,it,i,13)*WeightActive(i)+Sim1m(it2,it,i,9)*WeightActive(i)+Sim1m(it2,it,i,7)*WeightActive(i)
                    dum3=dum3+Sim1m(it2,it,i,13)*WeightActive(i)+Sim1m(it2,it,i,9)*WeightActive(i)+Sim1m(it2,it,i,7)*WeightActive(i)+Sim1m(it2,it,i,8)*WeightActive(i)+Sim1m(it2,it,i,1)*WeightActive(i)*r*tk
                    if(Sim1f(it2,it,i,10)<0.5d0) then
                        dum2=dum2+Sim1f(it2,it,i,6)*(1d0+t_employer)*WeightActive(i)
                        dum7=dum7+Sim1f(it2,it,i,7)*WeightActive(i)
                        dum15=dum15+Sim1f(it2,it,i,9)*WeightActive(i)
                        dum4=dum4+Sim1f(it2,it,i,13)*WeightActive(i)+Sim1f(it2,it,i,9)*WeightActive(i)+Sim1f(it2,it,i,7)*WeightActive(i)
                        dum3=dum3+Sim1f(it2,it,i,13)*WeightActive(i)+Sim1f(it2,it,i,9)*WeightActive(i)+Sim1f(it2,it,i,7)*WeightActive(i)+Sim1f(it2,it,i,8)*WeightActive(i)+Sim1f(it2,it,i,1)*WeightActive(i)*r*tk
                    end if
                end do
            end do

        end do

        do i=1,Tret

            r_ret=((1d0+r)/OmegaRet2(i))-1d0

            do it2=1,nsim2
                do it=1,nsim
                    dum5=dum5+2d0*WeightRet(i)
                    dum3=dum3+SimR1m(it2,it,i,16)*WeightRet(i)+SimR1m(it2,it,i,3)*WeightRet(i)+SimR1m(it2,it,i,8)*WeightRet(i)+SimR1m(it2,it,i,9)*WeightRet(i)+SimR1m(it2,it,i,1)*WeightRet(i)*r_ret*tk
                    dum2=dum2+SimR1m(it2,it,i,7)*(1d0+t_employer)*WeightRet(i)
                    dum4=dum4+SimR1m(it2,it,i,16)*WeightRet(i)+SimR1m(it2,it,i,8)*WeightRet(i)+SimR1m(it2,it,i,9)*WeightRet(i)
                    dum15=dum15+SimR1m(it2,it,i,9)*WeightRet(i)
                    if(Sim1f(it2,it,T,10)<0.5d0) then
                        dum3=dum3+SimR1f(it2,it,i,16)*WeightRet(i)+SimR1f(it2,it,i,3)*WeightRet(i)+SimR1f(it2,it,i,8)*WeightRet(i)+SimR1f(it2,it,i,9)*WeightRet(i)+SimR1f(it2,it,i,1)*WeightRet(i)*r_ret*tk
                        dum2=dum2+SimR1f(it2,it,i,7)*(1d0+t_employer)*WeightRet(i)
                        dum4=dum4+SimR1f(it2,it,i,16)*WeightRet(i)+SimR1f(it2,it,i,8)*WeightRet(i)+SimR1f(it2,it,i,9)*WeightRet(i)
                        dum15=dum15+SimR1f(it2,it,i,9)*WeightRet(i)
                    end if
                end do
            end do

        end do

        Print *,'Labor income tax rate including TSS is',dum4/dum2

        Print *,'Tax revenue per capita including TSS is',dum3/dum5

        !Print *,'Labor Income tax per capita is',dum7/dum5

        !Print *,'Social security tax per capita is',dum15/dum5

        !Individual tax rates

        dum2=max(maxval(Sim1m(:,:,:,12)),maxval(Sim1f(:,:,:,12)))

        dum2=max(dum2,maxval(SimR1m(:,:,:,15)),maxval(SimR1f(:,:,:,15)))

        Print *,'Max individual marginal tax rate is',dum2


        dum2=max(maxval(Sim1m(:,:,:,14)),maxval(Sim1f(:,:,:,14)))

        dum2=max(dum2,maxval(SimR1m(:,:,:,17)),maxval(SimR1f(:,:,:,17)))

        Print *,'Max individual average tax rate is',dum2



        !Social security

        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0
        dum7=0d0
        dum15=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    dum3=dum3+2d0*WeightActive(i)
                    dum4=dum4+Sim1m(it2,it,i,9)*WeightActive(i)

                    if(Sim1f(it2,it,i,10)<0.5d0) then
                        dum4=dum4+Sim1f(it2,it,i,9)*WeightActive(i)
                    end if

                end do
            end do

        end do

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    dum3=dum3+2d0*WeightRet(i)
                    dum4=dum4+SimR1m(it2,it,i,9)*WeightRet(i)

                    if(Sim1f(it2,it,i,10)<0.5d0) then
                        dum4=dum4+SimR1f(it2,it,i,9)*WeightRet(i)
                    end if

                end do
            end do

        end do

        ss_tax=dum4/dum3

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim

                    if(expR1m(it2,it,i,3)==1) then
                        dum2=dum2+1d0*WeightRet(i)
                    end if

                    dum5=dum5+SimR1m(it2,it,i,14)*WeightRet(i)

                    if(expR1f(it2,it,i,3)==1) then
                        dum2=dum2+1d0*WeightRet(i)
                    end if

                    dum5=dum5+SimR1f(it2,it,i,14)*WeightRet(i)

                end do
            end do

        end do

        Print *,'SS tax per retiree is',dum4/dum2
        Print *,'Average pension is',dum5/dum2

        dum5=dum5/dum3

        ss_expense=dum5

        Print *,'Social Security expenses per capita is',dum5
        !Print *,'Pension',Psi_pension/2d0

        !epsilon=Psi_pension-dum4*2d0

        epsilon=0d0
        !Psi_pension=Psi_pension-0.1d0*(Psi_pension-dum4*2d0)

        !Savings

        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+Sim1m(it2,it,i,1)*WeightActive(i)
                    dum3=dum3+2d0*WeightActive(i)

                    if(Sim1f(it2,it,i,10)<0.5d0) then
                        dum2=dum2+Sim1f(it2,it,i,1)*WeightActive(i)
                    end if

                end do
            end do

        end do

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+SimR1m(it2,it,i,1)*WeightRet(i)
                    dum3=dum3+2d0*WeightRet(i)

                    if(Sim1f(it2,it,T,10)<0.5d0) then
                        dum2=dum2+SimR1f(it2,it,i,1)*WeightRet(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Savings per capita is',dum2
        dum6=dum2

        ! Assets for redistribution

        !dum2=0d0
        !dum3=0d0
        !
        !do i=1,T
        !
        !do it2=1,nsim2
        !do it=1,nsim
        !    dum2=dum2+Sim1m(it2,it,i+1,1)*WeightActive(i)*(1d0-OmegaActive(i))
        !    dum3=dum3+1d0*WeightActive(i)
        !    
        !    if(Sim1f(it2,it,i,10)<0.5d0) then
        !        dum2=dum2+Sim1f(it2,it,i+1,1)*WeightActive(i)*(1d0-OmegaActive(i))
        !    end if
        !    
        !end do
        !end do
        !
        !end do
        !
        !do i=1,Tret
        !
        !do it2=1,nsim2
        !do it=1,nsim
        !    dum2=dum2+SimR1m(it2,it,i+1,1)*WeightRet(i)*(1d0-OmegaRet(i))
        !    dum3=dum3+1d0*WeightRet(i)
        !    
        !    if(Sim1f(it2,it,T,10)<0.5d0) then
        !        dum2=dum2+SimR1f(it2,it,i+1,1)*WeightRet(i)*(1d0-OmegaRet(i))
        !    end if
        !    
        !end do
        !end do
        !
        !end do
        !
        !dum2=dum2/dum3

        !Print *,'Assets per capita to be redistributed',dum2/2d0
        !Print *,'Gamma_redistr',Gamma_redistr/2d0

        !epsilon2=Gamma_redistr-dum2
        !
        !Gamma_redistr=Gamma_redistr-0.5d0*(Gamma_redistr-dum2)

        ! Capital tax

        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2 
                do it=1,nsim
                    dum2=dum2+Sim1m(it2,it,i,1)*WeightActive(i)*r*tk
                    dum3=dum3+2d0*WeightActive(i)

                    if(Sim1f(it2,it,i,10)<0.5d0) then
                        dum2=dum2+Sim1f(it2,it,i,1)*WeightActive(i)*r*tk
                    end if

                end do
            end do

        end do

        do i=1,Tret

            r_ret=((1d0+r)/OmegaRet2(i))-1d0    

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+SimR1m(it2,it,i,1)*WeightRet(i)*r_ret*tk
                    dum3=dum3+2d0*WeightRet(i)

                    if(Sim1f(it2,it,T,10)<0.5d0) then
                        dum2=dum2+SimR1f(it2,it,i,1)*WeightRet(i)*r_ret*tk
                    end if

                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Capital tax per capita is',dum2

        dum7=dum2

        !Consumption

        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+Sim1m(it2,it,i,2)*WeightActive(i)
                    dum3=dum3+2d0*WeightActive(i)

                    if(Sim1f(it2,it,i,10)<0.5d0) then
                        dum2=dum2+Sim1f(it2,it,i,2)*WeightActive(i)
                    end if

                end do
            end do

        end do

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+SimR1m(it2,it,i,2)*WeightRet(i)
                    dum3=dum3+2d0*WeightRet(i)

                    if(Sim1f(it2,it,T,10)<0.5d0) then
                        dum2=dum2+SimR1f(it2,it,i,2)*WeightRet(i)
                    end if

                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Consumption per capita is',dum2

        dum2=0d0
        dum3=0d0

        !Consumption tax

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+Sim1m(it2,it,i,8)*WeightActive(i)
                    dum3=dum3+2d0*WeightActive(i)

                    if(Sim1f(it2,it,i,10)<0.5d0) then
                        dum2=dum2+Sim1f(it2,it,i,8)*WeightActive(i)
                    end if

                end do
            end do

        end do

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+SimR1m(it2,it,i,3)*WeightRet(i)
                    dum3=dum3+2d0*WeightRet(i)

                    if(Sim1f(it2,it,T,10)<0.5d0) then
                        dum2=dum2+SimR1f(it2,it,i,3)*WeightRet(i)
                    end if

                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Consumption tax per capita is',dum2

        !Labor income tax

        dum4=0d0
        dum5=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    dum4=dum4+Sim1m(it2,it,i,7)*WeightActive(i)
                    dum5=dum5+2d0*WeightActive(i)

                    if(Sim1f(it2,it,i,10)<0.5d0) then
                        dum4=dum4+Sim1f(it2,it,i,7)*WeightActive(i)
                    end if

                end do
            end do

        end do


        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    dum5=dum5+2d0*WeightRet(i)
                    dum4=dum4+SimR1m(it2,it,i,8)*WeightRet(i)
                    if(Sim1f(it2,it,T,10)<0.5d0) then
                        dum4=dum4+SimR1f(it2,it,i,8)*WeightRet(i)
                    end if
                end do
            end do

        end do

        dum4=dum4/dum5

        Print *,'Labor income tax per capita is',dum4

        Print *,'Tax revenue per capita is',dum4+dum2+dum7

        dum5=dum4+dum2

        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+Sim1m(it2,it,i,13)*WeightActive(i)
                    dum3=dum3+2d0*WeightActive(i)

                    if(Sim1f(it2,it,i,10)<0.5d0) then
                        dum2=dum2+Sim1f(it2,it,i,13)*WeightActive(i)
                    end if

                end do
            end do

        end do


        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    dum3=dum3+2d0*WeightRet(i)
                    dum2=dum2+SimR1m(it2,it,i,16)*WeightRet(i)
                    if(Sim1f(it2,it,T,10)<0.5d0) then
                        dum2=dum2+SimR1f(it2,it,i,16)*WeightRet(i)
                    end if
                end do
            end do

        end do

        new_revenue=dum2/dum3

        Print *,'New Revenue is',new_revenue
        Print *,'t_const is',t_const

        dum2=0d0
        dum3=0d0
        dum15=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    dum3=dum3+2d0*WeightActive(i)
                    if(Sim1f(it2,it,i,4)<0.001) then
                        dum15=dum15+Unemp_benefit*WeightActive(i)
                    end if
                end do
            end do

        end do

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    dum3=dum3+2d0*WeightRet(i)
                    if(SimR1f(it2,it,i,5)<1d-3) then
                        dum15=dum15+Unemp_benefit*WeightRet(i)
                    end if
                end do
            end do

        end do

        dum15=dum15/dum3

        !GDP per capita

        dum9=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    dum9=dum9+(Sim1m(it2,it,i,6)*(1d0+t_employer)/w)*WeightActive(i)

                    if(Sim1f(it2,it,i,10)<0.5d0) then
                        dum9=dum9+(Sim1f(it2,it,i,6)*(1d0+t_employer)/w)*WeightActive(i)
                    end if

                end do
            end do

        end do



        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    dum9=dum9+SimR1m(it2,it,i,7)*((1d0+t_employer)/w)*WeightRet(i)

                    if(Sim1f(it2,it,T,10)<0.5d0) then
                        dum9=dum9+SimR1f(it2,it,i,7)*((1d0+t_employer)/w)*WeightRet(i)
                    end if

                end do
            end do

        end do


        !Print *,'Ltot is',dum9

        dum3=((ratio*dum9)**alpha)*(dum9**(1-alpha))/dum3

        Print *,'GDP per capita is',dum3

        lumpsumdum=0.02d0*dum3
        !
        Print *,'Lumpsumdum is',lumpsumdum/2d0
        Print *,'Lumpsum is',lumpsum/2d0

        !epsilon3=lumpsum-lumpsumdum

        !epsilon3=lumpsum-lumpsumdum

        epsilon3=0d0
        !lumpsum=lumpsum-0.1d0*(lumpsum-lumpsumdum)


        !Government Budget

        !lumpsumdum=(new_revenue+dum5+dum7+ss_tax)+mu*debttoGDP*dum3-(ss_expense+dum15+r*debttoGDP*dum3+(lumpsum/2d0))

        !lumpsumdum=(dum5+dum7+ss_tax)+mu*debttoGDP*dum3-(ss_expense+dum15+r*debttoGDP*dum3+lumpsum/2d0)

        lumpsumdum=(dum5+dum7+ss_tax)+mu*debttoGDP*dum3-(ss_expense+dum15+r*debttoGDP*dum3+lumpsum/2d0)
        
        Print *,'Net revenue is',lumpsumdum

        !epsilon3=lumpsumdum-0.119363895479142d0

        !lumpsumdum=new_revenue+0.00902d0
        !
        !lumpsumdum=lumpsumdum*2d0
        !
        !epsilon3=lumpsum-lumpsumdum
        !
        !Print *,'epsilon3 is',epsilon3

        !epsilon3=0d0
        !lumpsum=lumpsum-0.1d0*(lumpsum-lumpsumdum)


        !Labor income tax level

        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+Sim1m(it2,it,i,6)*WeightActive(i)
                    dum4=dum4+Sim1m(it2,it,i,7)*WeightActive(i)
                    dum3=dum3+2d0*WeightActive(i)

                    if(Sim1f(it2,it,i,10)<0.5d0) then
                        dum2=dum2+Sim1f(it2,it,i,6)*WeightActive(i)
                        dum4=dum4+Sim1f(it2,it,i,7)*WeightActive(i)
                    end if

                end do
            end do

        end do

        do i=1,Tret
        
            do it2=1,nsim2
                do it=1,nsim
                    dum3=dum3+2d0*WeightRet(i)
                    dum2=dum2+SimR1f(it2,it,i,7)*WeightRet(i)
                    dum4=dum4+SimR1f(it2,it,i,8)*WeightRet(i)
                end do
            end do
        
        end do


        Print *,'Average labor income tax rate is',dum4/dum2

        Print *,'Average individual earnings is',dum2/dum3
        
        !dum2=0d0
        !dum3=0d0
        !dum4=0d0
        !dum5=0d0
        !
        !do i=1,T
        !
        !    do it2=1,nsim2
        !        do it=1,nsim
        !        if(Sim1m(it2,it,i,10)>0.5d0) then
        !            dum2=dum2+Sim1m(it2,it,i,6)*WeightActive(i)
        !            dum3=dum3+1d0*WeightActive(i)
        !        end if
        !        end do
        !    end do
        !
        !    end do
        !    
        !    Print *,'Average married household earnings is',dum2/dum3
        !
        !STOP
        
        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,i,4)>0.001) then
                        dum2=dum2+Sim1m(it2,it,i,5)*WeightActive(i)
                        dum3=dum3+WeightActive(i)
                    end if
                    if(Sim1f(it2,it,i,4)>0.001) then
                        dum2=dum2+Sim1f(it2,it,i,5)*WeightActive(i)
                        dum3=dum3+WeightActive(i)
                    end if
                end do
            end do

        end do


        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    if(SimR1f(it2,it,i,5)>1d-3) then
                        dum2=dum2+SimR1f(it2,it,i,6)*WeightRet(i)
                        dum3=dum3+1d0*WeightRet(i)
                    end if
                    if(SimR1m(it2,it,i,5)>1d-3) then
                        dum2=dum2+SimR1m(it2,it,i,6)*WeightRet(i)
                        dum3=dum3+1d0*WeightRet(i)
                    end if
                end do
            end do

        end do


        Print *,'Average individual earnings for working people is',dum2/dum3

        Print *,'AE is',AE

        dum10=dum10+((dum2/dum3-1d0)/1d0)**2

        epsilon5=0d0

        !epsilon5=AE-dum2/dum3

        !AE=AE-0.2d0*(AE-dum2/dum3)

        !Unemp_benefit=0.201795*AE

        !Filling the average earnings matrix to be used for pensions


        !do ia=1,na
        !
        !    dum2=0d0
        !    dum3=0d0
        !    dum4=0d0
        !    dum5=0d0
        !
        !    do i=11,T
        !
        !        do it2=1,nsim2
        !            do it=1,nsim
        !
        !                if(exp1m(it2,it,i,2)==ia) then
        !                    if(Sim1m(it2,it,i,10)>0.5) then    
        !                        if(Sim1m(it2,it,i,4)>0.001) then
        !                            dum2=dum2+Sim1m(it2,it,i,5)*WeightActive(i)
        !                            dum3=dum3+WeightActive(i)
        !                        end if
        !                    end if
        !                end if
        !
        !                if(exp1f(it2,it,i,2)==ia) then
        !                    if(Sim1f(it2,it,i,10)>0.5) then
        !                        if(Sim1f(it2,it,i,4)>0.001) then
        !                            dum4=dum4+Sim1f(it2,it,i,5)*WeightActive(i)
        !                            dum5=dum5+WeightActive(i)
        !                        end if
        !                    end if
        !                end if
        !
        !
        !            end do
        !        end do
        !
        !    end do
        !
        !    if(dum3>0.01d0) then
        !        av_earnings(1,1,ia)=dum2/dum3
        !    else
        !        av_earnings(1,1,ia)=0d0
        !    end if
        !
        !
        !    if(dum5>0.01d0) then
        !        av_earnings(2,1,ia)=dum4/dum5
        !    else
        !        av_earnings(2,1,ia)=0d0
        !    end if
        !
        !end do
        !
        !do ia=1,na
        !
        !    dum2=0d0
        !    dum3=0d0
        !    dum4=0d0
        !    dum5=0d0
        !
        !    do i=11,T
        !
        !        do it2=1,nsim2
        !            do it=1,nsim
        !
        !                if(exp1m(it2,it,i,2)==ia) then
        !                    if(Sim1m(it2,it,i,10)<0.5) then    
        !                        if(Sim1m(it2,it,i,4)>0.001) then
        !                            dum2=dum2+Sim1m(it2,it,i,5)*WeightActive(i)
        !                            dum3=dum3+WeightActive(i)
        !                        end if
        !                    end if
        !                end if
        !
        !                if(exp1f(it2,it,i,2)==ia) then
        !                    if(Sim1f(it2,it,i,10)<0.5) then
        !                        if(Sim1f(it2,it,i,4)>0.001) then
        !                            dum4=dum4+Sim1f(it2,it,i,5)*WeightActive(i)
        !                            dum5=dum5+WeightActive(i)
        !                        end if
        !                    end if
        !                end if
        !
        !
        !            end do
        !        end do
        !
        !    end do
        !
        !    if(dum3>0.01d0) then
        !        av_earnings(1,2,ia)=dum2/dum3
        !    else
        !        av_earnings(1,2,ia)=0d0
        !    end if
        !
        !
        !    if(dum5>0.01d0) then
        !        av_earnings(2,2,ia)=dum4/dum5
        !    else
        !        av_earnings(2,2,ia)=0d0
        !    end if
        !
        !end do


        !Prices

        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+Sim1m(it2,it,i,1)*WeightActive(i)

                    if(Sim1f(it2,it,i,10)<0.5d0) then
                        dum2=dum2+Sim1f(it2,it,i,1)*WeightActive(i)
                    end if

                end do
            end do

        end do

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+(SimR1m(it2,it,i,1)/OmegaRet2(i))*WeightRet(i)

                    if(Sim1f(it2,it,T,10)<0.5d0) then
                        dum2=dum2+(SimR1f(it2,it,i,1)/OmegaRet2(i))*WeightRet(i)
                    end if

                end do
            end do

        end do

        dum6=dum2-debttoGDP*((ratio*dum9)**(alpha))*(dum9**(1d0-alpha))

        Print *,'capital stock is',dum6/population_mass

        ratiodum=dum6/dum9

        Print *,'Ratio between capital and labor is',ratiodum

        Print *,'Implied wage is',(1d0-alpha)*ratiodum**alpha
        Print *,'Implied interest is',alpha*ratiodum**(alpha-1d0)-delta

        Print *,'Wage is',w
        Print *,'Interest is',r

        Print *,'K/Y is',ratiodum**(1d0-alpha)

        dum10=dum10+(((ratiodum**(1d0-alpha))-2.6399)/2.6399)**2

        Print *,'FCN is',dum10

        ! beta=1.00251
        ! Fw=0.0125
        ! Chi=26.1d0
        ! gamma0=0.3288d0
        ! FCN=0.00000559

        !Print *,'I/Y is',delta*ratiodum**(1d0-alpha)

        dum9=max(maxval(Sim1m(:,:,:,1)),maxval(Sim1f(:,:,:,1)))

        dum9=max(dum9,maxval(SimR1m(:,:,:,1)),maxval(SimR1f(:,:,:,1)))

        Print *,'Max savings is',dum9

        dum2=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,i,10)>0.5d0) then
                        dum2=dum2+1d0*WeightActive(i)
                    end if

                end do
            end do

        end do

        Print *,'Fraction of Married females is',dum2/(1d0*T*nsim2*nsim)

        dum2=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,i,10)>0.5d0) then
                        dum2=dum2+1d0*WeightActive(i)
                    end if

                end do
            end do

        end do

        Print *,'Fraction of Married males is',dum2/(1d0*T*nsim2*nsim)


        dum2=0d0
        dum3=0d0

        i=min(T,11)

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1f(it2,it,i,10)>0.5d0) then
                    dum2=dum2+1d0*WeightActive(i)
                    if(Sim1f(it2,it,T,10)>0.5d0) then
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end if
            end do
        end do


        Print *,'Probability of being married at 65 conditional on married at 30 is',dum3/dum2

        dum2=0d0
        dum3=0d0

        i=min(T,11)

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1f(it2,it,i,10)>0.5d0) then
                    dum2=dum2+1d0*WeightActive(i)
                    if(Sim1f(it2,it,min(T,12),10)>0.5d0) then
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end if
            end do
        end do


        Print *,'Probability of being married at 31 conditional on married at 30 is',dum3/dum2

        dum2=0d0
        dum3=0d0

        i=min(T,11)

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1f(it2,it,i,10)>0.5d0) then
                    dum2=dum2+1d0*WeightActive(i)
                    if(Sim1f(it2,it,min(T,16),10)>0.5d0) then
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end if
            end do
        end do


        Print *,'Probability of being married at 35 conditional on married at 30 is',dum3/dum2

        dum2=0d0
        dum3=0d0

        i=min(T,11)

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1f(it2,it,i,10)>0.5d0) then
                    dum2=dum2+1d0*WeightActive(i)
                    if(Sim1f(it2,it,min(T,31),10)>0.5d0) then
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end if
            end do
        end do


        Print *,'Probability of being married at 50 conditional on married at 30 is',dum3/dum2

        dum2=0d0
        dum3=0d0

        i=min(T,11)

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1f(it2,it,i,10)>0.5d0) then
                    dum2=dum2+1d0*WeightActive(i)
                    if(Sim1f(it2,it,min(T,44),10)>0.5d0) then
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end if
            end do
        end do


        Print *,'Probability of being married at 63 conditional on married at 30 is',dum3/dum2

        dum2=0d0
        dum3=0d0

        i=min(T,11)

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1f(it2,it,i,10)>0.5d0) then
                    dum2=dum2+1d0*WeightActive(i)
                    if(Sim1f(it2,it,T,10)>0.5d0) then
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end if
            end do
        end do


        Print *,'Probability of being married at 64 conditional on married at 30 is',dum3/dum2

        dum2=0d0
        dum3=0d0

        i=min(T,11)

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1m(it2,it,i,10)>0.5d0) then
                    dum2=dum2+1d0*WeightActive(i)
                    if(Sim1m(it2,it,T,10)>0.5d0) then
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end if
            end do
        end do


        Print *,'Male Probability of being married at 64 conditional on married at 30 is',dum3/dum2


        dum2=0d0

        i=T

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1f(it2,it,i,10)>0.5d0) then
                    dum2=dum2+1d0*WeightActive(i)
                end if

            end do
        end do

        Print *,'Fraction of Married females at age 64 is',dum2/(1d0*nsim2*nsim)

        dum2=0d0

        i=min(T,31)

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1f(it2,it,i,10)>0.5d0) then
                    dum2=dum2+1d0*WeightActive(i)
                end if

            end do
        end do

        Print *,'Fraction of Married females at age 50 is',dum2/(1d0*nsim2*nsim)


        mpartnerdum=0d0
        fpartnerdum=0d0

        do ia2=1,na

            do i=1,T

                dum3=0d0    

                do it2=1,nsim2
                    do it=1,nsim

                        if((Sim1m(it2,it,i,10)<0.5).AND.(ia2==exp1m(it2,it,i,2))) then
                            dum2=Sim1m(it2,it,i,1)
                            ia=ia2
                            iu=exp1m(it2,it,i,3)
                            ifc=exp1m(it2,it,i,6)
                            dum5=exp2m(it2,it,i,1)

                            if(i==1) then
                                ix=1
                            end if

                            if(dum5<exp_grid(2,i)/2d0) then
                                ix=1
                            end if

                            do ik=2,nexp-1  
                                if(((exp_grid(ik,i)-(exp_grid(ik,i)-exp_grid(ik-1,i))/2d0)<dum5).AND.(dum5<(exp_grid(ik,i)+(exp_grid(ik+1,i)-exp_grid(ik,i))/2d0))) then
                                    ix=ik
                                end if
                            end do

                            if(dum5>(exp_grid(nexp,i)-(exp_grid(nexp,i)-exp_grid(nexp-1,i))/2d0)) then
                                ix=nexp
                            end if

                            if(dum2<k_grid(2)/2d0) then
                                mpartnerdum(1,ix,ia,iu,i,ifc)=mpartnerdum(1,ix,ia,iu,i,ifc)+1d0
                            end if

                            do ik=2,nk-1  
                                if(((k_grid(ik)-(k_grid(ik)-k_grid(ik-1))/2d0)<dum2).AND.(dum2<(k_grid(ik)+(k_grid(ik+1)-k_grid(ik))/2d0))) then
                                    mpartnerdum(ik,ix,ia,iu,i,ifc)=mpartnerdum(ik,ix,ia,iu,i,ifc)+1d0
                                end if
                            end do

                            if(dum2>(k_grid(nk)-(k_grid(nk)-k_grid(nk-1))/2d0)) then
                                mpartnerdum(nk,ix,ia,iu,i,ifc)=mpartnerdum(nk,ix,ia,iu,i,ifc)+1d0
                            end if

                            dum3=dum3+1d0
                        end if

                    end do
                end do

                mpartnerdum(:,:,ia2,:,i,:)=mpartnerdum(:,:,ia2,:,i,:)/dum3

            end do

        end do

        do ia2=1,na

            do i=1,T

                dum3=0d0     

                do it2=1,nsim2
                    do it=1,nsim

                        if((Sim1f(it2,it,i,10)<0.5).AND.(ia2==exp1f(it2,it,i,2))) then
                            dum2=Sim1f(it2,it,i,1)
                            ia=ia2
                            iu=exp1f(it2,it,i,3)
                            ifc=exp1f(it2,it,i,6)
                            dum5=exp2f(it2,it,i,1)
                            dum4=0d0

                            if(i==1) then
                                ix=1
                            end if

                            if(dum5<exp_grid(2,i)/2d0) then
                                ix=1
                            end if

                            do ik=2,nexp-1  
                                if(((exp_grid(ik,i)-(exp_grid(ik,i)-exp_grid(ik-1,i))/2d0)<dum5).AND.(dum5<(exp_grid(ik,i)+(exp_grid(ik+1,i)-exp_grid(ik,i))/2d0))) then
                                    ix=ik
                                end if
                            end do

                            if(dum5>(exp_grid(nexp,i)-(exp_grid(nexp,i)-exp_grid(nexp-1,i))/2d0)) then
                                ix=nexp
                            end if

                            if(dum2<k_grid(2)/2d0) then
                                fpartnerdum(1,ix,ia,iu,i,ifc)=fpartnerdum(1,ix,ia,iu,i,ifc)+1d0
                                dum4=1d0
                            end if

                            do ik=2,nk-1
                                if(((k_grid(ik)-(k_grid(ik)-k_grid(ik-1))/2d0)<dum2).AND.(dum2<(k_grid(ik)+(k_grid(ik+1)-k_grid(ik))/2d0))) then
                                    fpartnerdum(ik,ix,ia,iu,i,ifc)=fpartnerdum(ik,ix,ia,iu,i,ifc)+1d0
                                    dum4=1d0
                                end if
                            end do

                            if(dum2>(k_grid(nk)-(k_grid(nk)-k_grid(nk-1))/2d0)) then
                                fpartnerdum(nk,ix,ia,iu,i,ifc)=fpartnerdum(nk,ix,ia,iu,i,ifc)+1d0
                                dum4=1d0
                            end if

                            dum3=dum3+1d0
                        end if

                    end do
                end do

                fpartnerdum(:,:,ia2,:,i,:)=fpartnerdum(:,:,ia2,:,i,:)/dum3

            end do

        end do


        fpartnerdum2=abs(fpartnerdum-fpartner)

        mpartnerdum2=abs(mpartnerdum-mpartner)

        epsilon6=max(maxval(fpartnerdum2),maxval(mpartnerdum2))

        Print *,'max distance between distribution of singles is',epsilon6

        !if(iter>5) then
        epsilon6=0d0
        !end if

        fpartner=fpartnerdum

        mpartner=mpartnerdum

        dum2=0d0

        !do ifc=1,nfc
        !    do ik=1,nk
        !        do ia=1,na
        !            do ix=1,nexp
        !                do iu=1,nu
        !                    dum2=dum2+fpartner(ik,ix,ia,iu,30,ifc)
        !                end do
        !            end do
        !        end do
        !    end do
        !end do
        !
        !Print *,'sum fpartner is',dum2

        dum2=0d0

        !do ifc=1,nfcm
        !    do ik=1,nk
        !        do ia=1,na
        !            do ix=1,nexp
        !                do iu=1,nu
        !                    dum2=dum2+mpartner(ik,ix,ia,iu,30,ifc)
        !                end do
        !            end do
        !        end do
        !    end do
        !end do
        !
        !Print *,'sum mpartner is',dum2

        !Here we are computing the matrix of  probabilities for marrying someone of ability a' if you have probability a.

        ability_prob=0d0

        do ia=1,na

            dum3=0d0

            do i=1,T

                do it2=1,nsim2
                    do it=1,nsim

                        if((Sim1m(it2,it,i,10)>0.5).AND.(ia==exp1m(it2,it,i,2))) then

                            it4=exp1m(it2,it,i,4)
                            ia2=exp1f(it2,it4,i,2)

                            ability_prob(ia,ia2)=ability_prob(ia,ia2)+1d0
                            dum3=dum3+1d0

                        end if

                    end do
                end do

            end do

            ability_prob(ia,:)=ability_prob(ia,:)/dum3

        end do

        Print *,'ability_prob is',ability_prob(2,:)

        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0
        dum6=0d0
        dum7=0d0
        dum8=0d0
        dum9=0d0
        dum10=0d0
        dum11=0d0
        dum12=0d0
        dum13=0d0
        dum14=0d0
        dum15=0d0
        dum16=0d0
        dum17=0d0
        dum18=0d0
        dum19=0d0
        dum20=0d0
        dum21=0d0
        dum22=0d0
        dum23=0d0
        dum24=0d0
        dum25=0d0


        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    dum7=dum7+1d0*WeightActive(i)
                    if(exp1f(it2,it,i,6)==1) then
                        dum2=dum2+1d0*WeightActive(i)
                    elseif(exp1f(it2,it,i,6)==2) then
                        dum3=dum3+1d0*WeightActive(i)
                    elseif(exp1f(it2,it,i,6)==3) then
                        dum4=dum4+1d0*WeightActive(i)
                    elseif(exp1f(it2,it,i,6)==4) then
                        dum5=dum5+1d0*WeightActive(i)
                    elseif(exp1f(it2,it,i,6)==5) then
                        dum6=dum6+1d0*WeightActive(i)
                    elseif(exp1f(it2,it,i,6)==6) then
                        dum8=dum8+1d0*WeightActive(i)
                    elseif(exp1f(it2,it,i,6)==7) then
                        dum9=dum9+1d0*WeightActive(i)
                    elseif(exp1f(it2,it,i,6)==8) then
                        dum10=dum10+1d0*WeightActive(i)
                    elseif(exp1f(it2,it,i,6)==9) then
                        dum11=dum11+1d0*WeightActive(i)
                    elseif(exp1f(it2,it,i,6)==10) then
                        dum12=dum12+1d0*WeightActive(i)
                    elseif(exp1f(it2,it,i,6)==11) then
                        dum13=dum13+1d0*WeightActive(i)
                    elseif(exp1f(it2,it,i,6)==12) then
                        dum14=dum14+1d0*WeightActive(i)
                    elseif(exp1f(it2,it,i,6)==13) then
                        dum15=dum15+1d0*WeightActive(i)
                    elseif(exp1f(it2,it,i,6)==14) then
                        dum16=dum16+1d0*WeightActive(i)
                    elseif(exp1f(it2,it,i,6)==15) then
                        dum17=dum17+1d0*WeightActive(i)
                    elseif(exp1f(it2,it,i,6)==16) then
                        dum18=dum18+1d0*WeightActive(i)
                    elseif(exp1f(it2,it,i,6)==17) then
                        dum19=dum19+1d0*WeightActive(i)
                    elseif(exp1f(it2,it,i,6)==18) then
                        dum20=dum20+1d0*WeightActive(i)
                    elseif(exp1f(it2,it,i,6)==19) then
                        dum21=dum21+1d0*WeightActive(i)
                    elseif(exp1f(it2,it,i,6)==20) then
                        dum22=dum22+1d0*WeightActive(i)
                    elseif(exp1f(it2,it,i,6)==21) then
                        dum23=dum23+1d0*WeightActive(i)
                    elseif(exp1f(it2,it,i,6)==22) then
                        dum24=dum24+1d0*WeightActive(i)
                    else
                        dum25=dum25+1d0*WeightActive(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum7
        dum3=dum3/dum7
        dum4=dum4/dum7
        dum5=dum5/dum7
        dum6=dum6/dum7
        dum8=dum8/dum7
        dum9=dum9/dum7
        dum10=dum10/dum7
        dum11=dum11/dum7
        dum12=dum12/dum7
        dum13=dum13/dum7
        dum14=dum14/dum7
        dum15=dum15/dum7
        dum16=dum16/dum7
        dum17=dum17/dum7
        dum18=dum18/dum7
        dum19=dum19/dum7
        dum20=dum20/dum7
        dum21=dum21/dum7
        dum22=dum22/dum7
        dum23=dum23/dum7
        dum24=dum24/dum7
        dum25=dum25/dum7



        Print *,'Distribution of female fixed cost shock is',dum2,dum3,dum4,dum5,dum6,dum8,dum9,dum10,dum11,dum12,dum13,dum14,dum15,dum16,dum17,dum18,dum19,dum20,dum21,dum22,dum23,dum24,dum25

        ! Taxes by demographic group

        !Singles

        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0
        dum6=0d0
        dum7=0d0


        do i=1,T

            do it2=1,nsim2 
                do it=1,nsim

                    if(Sim1m(it2,it,i,10)<0.5d0) then
                        dum2=dum2+Sim1m(it2,it,i,1)*WeightActive(i)*r*tk
                        dum3=dum3+1d0*WeightActive(i)
                    end if

                    if(Sim1f(it2,it,i,10)<0.5d0) then
                        dum4=dum4+Sim1f(it2,it,i,1)*WeightActive(i)*r*tk
                        dum5=dum5+1d0*WeightActive(i)
                    end if

                end do
            end do

        end do



        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim

                    r_ret=((1d0+r)/OmegaRet2(i))-1d0

                    if(Sim1m(it2,it,T,10)<0.5d0) then
                        dum2=dum2+SimR1m(it2,it,i,1)*WeightRet(i)*r_ret*tk
                        dum3=dum3+1d0*WeightRet(i)
                    end if

                    if(Sim1f(it2,it,T,10)<0.5d0) then
                        dum4=dum4+SimR1f(it2,it,i,1)*WeightRet(i)*r_ret*tk
                        dum5=dum5+1d0*WeightRet(i)
                    end if

                end do
            end do

        end do

        dum2=dum2/dum3
        dum4=dum4/dum5

        Print *,'Capital tax for singles is',(dum2+dum4)/2d0
        Print *,'Capital tax for single men is',dum2
        Print *,'Capital tax for single women is',dum4

        dum6=dum2
        dum7=dum4


        !Consumption tax singles

        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim

                    if(Sim1m(it2,it,i,10)<0.5d0) then
                        dum2=dum2+Sim1m(it2,it,i,8)*WeightActive(i)
                        dum3=dum3+1d0*WeightActive(i)
                    end if

                    if(Sim1f(it2,it,i,10)<0.5d0) then
                        dum4=dum4+Sim1f(it2,it,i,8)*WeightActive(i)
                        dum5=dum5+1d0*WeightActive(i)
                    end if

                end do
            end do

        end do

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim

                    if(Sim1f(it2,it,T,10)<0.5d0) then
                        dum2=dum2+SimR1m(it2,it,i,3)*WeightRet(i)
                        dum3=dum3+1d0*WeightRet(i)
                    end if

                    if(Sim1f(it2,it,T,10)<0.5d0) then
                        dum4=dum4+SimR1f(it2,it,i,3)*WeightRet(i)
                        dum5=dum5+1d0*WeightRet(i)
                    end if

                end do
            end do

        end do

        dum2=dum2/dum3
        dum4=dum4/dum5

        Print *,'Consumption tax for singles is',(dum2+dum4)/2d0
        Print *,'Consumption tax for single men is',dum2
        Print *,'Consumption tax for single women is',dum4

        dum6=dum6+dum2
        dum7=dum7+dum4

        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0
        dum10=0d0
        dum11=0d0
        dum12=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim

                    if(Sim1m(it2,it,i,10)<0.5d0) then
                        dum2=dum2+Sim1m(it2,it,i,7)*WeightActive(i)
                        dum10=dum10+Sim1m(it2,it,i,6)*WeightActive(i)
                        dum3=dum3+1d0*WeightActive(i)
                    end if

                    if(Sim1f(it2,it,i,10)<0.5d0) then
                        dum4=dum4+Sim1f(it2,it,i,7)*WeightActive(i)
                        dum11=dum11+Sim1f(it2,it,i,6)*WeightActive(i)
                        dum5=dum5+1d0*WeightActive(i)
                    end if

                end do
            end do

        end do

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim

                    if(Sim1f(it2,it,T,10)<0.5d0) then
                        dum3=dum3+1d0*WeightRet(i)
                    end if

                    if(Sim1f(it2,it,T,10)<0.5d0) then
                        dum5=dum5+1d0*WeightRet(i)
                    end if

                end do
            end do

        end do

        dum2=dum2/dum3
        dum4=dum4/dum5

        dum10=dum10/dum3
        dum11=dum11/dum5

        Print *,'Labor income tax for singles is',(dum2+dum4)/2d0
        Print *,'Labor income tax for single men is',dum2
        Print *,'Labor income tax for single women is',dum4

        Print *,'Labor income tax rate for singles is',(dum2+dum4)/(dum10+dum11)
        Print *,'Labor income tax rate for single men is',dum2/dum10
        Print *,'Labor income tax rate for single women is',dum4/dum11

        Print *,'Total tax revenue for singles is',(dum2+dum4+dum6+dum7)/2d0
        Print *,'Total tax revenueor single men is',dum2+dum6
        Print *,'Total tax revenue for single women is',dum4+dum7

        !Married

        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0
        dum6=0d0
        dum7=0d0


        do i=1,T

            do it2=1,nsim2 
                do it=1,nsim

                    if(Sim1m(it2,it,i,10)>0.5d0) then
                        dum2=dum2+Sim1m(it2,it,i,1)*WeightActive(i)*r*tk
                        dum3=dum3+2d0*WeightActive(i)
                    end if

                end do
            end do

        end do



        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim

                    r_ret=((1d0+r)/OmegaRet2(i))-1d0

                    if(Sim1m(it2,it,T,10)>0.5d0) then
                        dum2=dum2+SimR1m(it2,it,i,1)*WeightRet(i)*r_ret*tk
                        dum3=dum3+2d0*WeightRet(i)
                    end if

                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Capital tax for married is',dum2

        dum6=dum2


        !Consumption tax married

        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim

                    if(Sim1m(it2,it,i,10)>0.5d0) then
                        dum2=dum2+Sim1m(it2,it,i,8)*WeightActive(i)
                        dum3=dum3+2d0*WeightActive(i)
                    end if

                end do
            end do

        end do


        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim

                    if(Sim1f(it2,it,T,10)>0.5d0) then
                        dum2=dum2+SimR1m(it2,it,i,3)*WeightRet(i)
                        dum3=dum3+2d0*WeightRet(i)
                    end if

                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Consumption tax for married is',dum2+dum4

        dum6=dum6+dum2

        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim

                    if(Sim1m(it2,it,i,10)>0.5d0) then
                        dum2=dum2+Sim1m(it2,it,i,7)*WeightActive(i)
                        dum12=dum12+Sim1m(it2,it,i,6)*WeightActive(i)
                        dum3=dum3+2d0*WeightActive(i)
                    end if

                end do
            end do

        end do

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim

                    if(Sim1f(it2,it,T,10)>0.5d0) then
                        dum3=dum3+2d0*WeightRet(i)
                    end if

                end do
            end do

        end do

        dum2=dum2/dum3
        dum12=dum12/dum3

        Print *,'Labor income tax for married is',dum2

        Print *,'Labor income tax rate for married is',dum2/dum12

        Print *,'Total tax revenue for married is',dum2+dum6

        !Single female labor income tax after 65

        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0


        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,T,10)<0.5) then
                        if(SimR1f(it2,it,i,5)>1d-3) then
                            dum2=dum2+SimR1f(it2,it,i,7)*WeightRet(i)
                            dum4=dum4+SimR1f(it2,it,i,8)*WeightRet(i)
                            dum3=dum3+1d0*WeightRet(i)
                        end if
                    end if
                end do
            end do

        end do

        dum5=dum4/dum2

        Print *,'Labor income tax rate of single females 65+ is',dum5

        !Married female labor income tax after 65

        dum2=0d0
        dum3=0d0
        dum4=0d0
        dum5=0d0


        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,T,10)>0.5) then
                        if(SimR1f(it2,it,i,5)>1d-3) then
                            dum2=dum2+SimR1f(it2,it,i,7)*WeightRet(i)
                            dum4=dum4+SimR1f(it2,it,i,8)*WeightRet(i)
                            dum3=dum3+1d0*WeightRet(i)
                        end if
                    end if
                end do
            end do

        end do

        dum5=dum4/dum2

        Print *,'Labor income tax rate of married females 65+ is',dum5

        !STOP

        !!Share of singles by age
        !do i=1,T
        !    dum2=0d0
        !    dum3=0d0
        !    do it2=1,nsim2
        !        do it=1,nsim
        !            if(Sim1m(it2,it,i,10)<0.5d0) then
        !                dum2=dum2+1d0*WeightActive(i)
        !            end if
        !            dum3=dum3+1d0*WeightActive(i)
        !        end do
        !    end do
        !    dum2=dum2/dum3
        !    share_single(i)=dum2
        !end do

        !Social Welfare

        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+Sim1m(it2,it,i,11)*WeightActive(i)
                    dum2=dum2+Sim1f(it2,it,i,11)*WeightActive(i)
                    dum3=dum3+2d0*WeightActive(i)

                end do
            end do

        end do

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    dum2=dum2+SimR1m(it2,it,i,11)*WeightRet(i)
                    dum2=dum2+SimR1f(it2,it,i,11)*WeightRet(i)
                    dum3=dum3+2d0*WeightRet(i)
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Social welfare of everyone is',dum2

        dum2=0d0
        dum3=0d0

        i=1

        do it2=1,nsim2
            do it=1,nsim
                dum2=dum2+Sim1m(it2,it,i,11)*WeightActive(i)
                dum2=dum2+Sim1f(it2,it,i,11)*WeightActive(i)
                dum3=dum3+2d0*WeightActive(i)

            end do
        end do

        dum2=dum2/dum3

        Print *,'Social welfare of 20-year olds is',dum2

    !Ex-post + ex-ante Social Welfare of 65-year old married males by ability

    do ia=1,na

        dum2=0d0
        dum3=0d0
        dum4=0d0

        i=T

        do it2=1,nsim2
            do it=1,nsim
                !if((Sim1m(it2,it,i,10)>0.5d0).AND.(Sim1m(it2,it,i,17)>34d0)) then
                if(Sim1m(it2,it,i,10)>0.5d0) then
                if(exp1m(it2,it,i,2)==ia) then
                    !dum2=dum2+Sim1m(it2,it,i+1,16)*WeightActive(i)
                    !dum4=dum4+(beta**i)*Sim1m(it2,it,i,11)*WeightActive(i)
                    dum2=dum2+Sim1m(it2,it,i,16)*WeightActive(i)
                    dum4=dum4+(beta**(i-1))*Sim1m(it2,it,i,11)*WeightActive(i)
                    dum3=dum3+1d0*WeightActive(i)
                end if
                end if
            end do
        end do

        dum2=dum2/dum3
        dum4=dum4/dum3

        Print *,'ability',ia,'Ex-post Social welfare of 65-year old married men is',dum2
        Print *,'ability',ia,'Ex-ante Social welfare of 65-year old married men is',dum4
        Print *,'ability',ia,'Ex-post + Ex-ante Social welfare of 65-year old married men is',dum2+dum4
            
    end do

    !Ex-post + ex-ante Social Welfare of 65-year old married females by ability

    do ia=1,na

        dum2=0d0
        dum3=0d0
        dum4=0d0

        i=T

        do it2=1,nsim2
            do it=1,nsim
                !if((Sim1m(it2,it,i,10)>0.5d0).AND.(Sim1m(it2,it,i,17)>34d0)) then
                if(Sim1f(it2,it,i,10)>0.5d0) then
                if(exp1f(it2,it,i,2)==ia) then
                    !dum2=dum2+Sim1f(it2,it,i+1,16)*WeightActive(i)
                    !dum4=dum4+(beta**i)*Sim1f(it2,it,i,11)*WeightActive(i)
                    dum2=dum2+Sim1f(it2,it,i,16)*WeightActive(i)
                    dum4=dum4+(beta**(i-1))*Sim1f(it2,it,i,11)*WeightActive(i)
                    dum3=dum3+1d0*WeightActive(i)
                end if
                end if
            end do
        end do

        dum2=dum2/dum3
        dum4=dum4/dum3

        Print *,'ability',ia,'Ex-post Social welfare of 65-year old married women is',dum2
        Print *,'ability',ia,'Ex-ante Social welfare of 65-year old married women is',dum4
        Print *,'ability',ia,'Ex-post + Ex-ante Social welfare of 65-year old married women is',dum2+dum4
            
    end do


    !Ex-post + ex-ante Social Welfare of 65-year old single males by ability

    do ia=1,na

        dum2=0d0
        dum3=0d0
        dum4=0d0

        i=T

        do it2=1,nsim2
            do it=1,nsim
                !if((Sim1m(it2,it,i,10)>0.5d0).AND.(Sim1m(it2,it,i,17)>34d0)) then
                if(Sim1m(it2,it,i,10)<0.5d0) then
                if(exp1m(it2,it,i,2)==ia) then
                    !dum2=dum2+Sim1m(it2,it,i+1,16)*WeightActive(i)
                    !dum4=dum4+(beta**i)*Sim1m(it2,it,i,11)*WeightActive(i)
                    dum2=dum2+Sim1m(it2,it,i,16)*WeightActive(i)
                    dum4=dum4+(beta**(i-1))*Sim1m(it2,it,i,11)*WeightActive(i)
                    dum3=dum3+1d0*WeightActive(i)
                end if
                end if
            end do
        end do

        dum2=dum2/dum3
        dum4=dum4/dum3

        Print *,'ability',ia,'Ex-post Social welfare of 65-year old single men is',dum2
        Print *,'ability',ia,'Ex-ante Social welfare of 65-year old single men is',dum4
        Print *,'ability',ia,'Ex-post + Ex-ante Social welfare of 65-year old single men is',dum2+dum4
            
    end do


    !Ex-post + ex-ante Social Welfare of 65-year old single females by ability

    do ia=1,na

        dum2=0d0
        dum3=0d0
        dum4=0d0

        i=T

        do it2=1,nsim2
            do it=1,nsim
                !if((Sim1m(it2,it,i,10)>0.5d0).AND.(Sim1m(it2,it,i,17)>34d0)) then
                if(Sim1f(it2,it,i,10)<0.5d0) then
                if(exp1f(it2,it,i,2)==ia) then
                    !dum2=dum2+Sim1f(it2,it,i+1,16)*WeightActive(i)
                    !dum4=dum4+(beta**i)*Sim1f(it2,it,i,11)*WeightActive(i)
                    dum2=dum2+Sim1f(it2,it,i,16)*WeightActive(i)
                    dum4=dum4+(beta**(i-1))*Sim1f(it2,it,i,11)*WeightActive(i)
                    dum3=dum3+1d0*WeightActive(i)
                end if
                end if
            end do
        end do

        dum2=dum2/dum3
        dum4=dum4/dum3

        Print *,'ability',ia,'Ex-post Social welfare of 65-year old single women is',dum2
        Print *,'ability',ia,'Ex-ante Social welfare of 65-year old single women is',dum4
        Print *,'ability',ia,'Ex-post + Ex-ante Social welfare of 65-year old single women is',dum2+dum4
            
    end do


        
    !Social Welfare of males by ability

    do ia=1,na

    dum2=0d0
    dum3=0d0

        i=1

        do it2=1,nsim2
            do it=1,nsim
                if(exp1m(it2,it,i,2)==ia) then
                    dum2=dum2+Sim1m(it2,it,i,11)*WeightActive(i)
                    dum3=dum3+WeightActive(i)
                end if
            end do
        end do

        dum2=dum2/dum3

        Print *,'ability',ia,'Social welfare of 20-year males',dum2
        
    end do


    !Social Welfare of females by ability
        
    do ia=1,na

    dum2=0d0
    dum3=0d0

        i=1

        do it2=1,nsim2
            do it=1,nsim
                if(exp1f(it2,it,i,2)==ia) then
                    dum2=dum2+Sim1f(it2,it,i,11)*WeightActive(i)
                    dum3=dum3+WeightActive(i)
                end if
            end do
        end do

        dum2=dum2/dum3

        Print *,'ability',ia,'Social welfare of 20-year females',dum2
        
    end do
        
        
        !Ex-post Social Welfare at 65

        dum2=0d0
        dum3=0d0

        i=T

        do it2=1,nsim2
            do it=1,nsim
                dum2=dum2+Sim1m(it2,it,i+1,16)*WeightActive(i)
                dum2=dum2+Sim1f(it2,it,i+1,16)*WeightActive(i)
                dum3=dum3+2d0*WeightActive(i)

            end do
        end do

        dum2=dum2/dum3

        Print *,'Ex-post Social welfare of everyone at 65 is',dum2


        !Social welfare of married and single at age 65

        !Social welfare Married men at 65

        dum2=0d0
        dum3=0d0
        dum4=0d0

        i=T

        do it2=1,nsim2
            do it=1,nsim
                !if((Sim1m(it2,it,i,10)>0.5d0).AND.(Sim1m(it2,it,i,17)>34d0)) then
                if(Sim1m(it2,it,i,10)>0.5d0) then
                    dum2=dum2+Sim1m(it2,it,i+1,16)*WeightActive(i)
                    dum4=dum4+(beta**i)*Sim1m(it2,it,i,11)*WeightActive(i)
                    dum3=dum3+1d0*WeightActive(i)
                end if
            end do
        end do

        dum2=dum2/dum3
        dum4=dum4/dum3

        Print *,'Ex-post Social welfare of 65-year old married men is',dum2
        Print *,'Ex-ante Social welfare of 65-year old married men is',dum4
        Print *,'Ex-post + Ex-ante Social welfare of 65-year old married men is',dum2+dum4
        Print *,'Number of 65-year old married men is',dum3


        !Social welfare single men at 65

        dum2=0d0
        dum3=0d0
        dum4=0d0

        i=T

        do it2=1,nsim2
            do it=1,nsim
                !if((Sim1m(it2,it,i,10)<0.5d0).AND.(Sim1m(it2,it,i,17)<10d0)) then
                if(Sim1m(it2,it,i,10)<0.5d0) then
                    dum2=dum2+Sim1m(it2,it,i+1,16)*WeightActive(i)
                    dum4=dum4+(beta**i)*Sim1m(it2,it,i,11)*WeightActive(i)
                    dum3=dum3+1d0*WeightActive(i)
                end if
            end do
        end do

        dum2=dum2/dum3
        dum4=dum4/dum3

        Print *,'Ex-post Social welfare of 65-year old single men is',dum2
        Print *,'Ex-ante Social welfare of 65-year old single men is',dum4
        Print *,'Ex-post + Ex-ante Social welfare of 65-year old single men is',dum2+dum4
        Print *,'Number of 65-year old single men is',dum3


        !Social welfare married women at 65


        dum2=0d0
        dum3=0d0
        dum4=0d0

        i=T

        do it2=1,nsim2
            do it=1,nsim
                !if((Sim1f(it2,it,i,10)>0.5d0).AND.(Sim1f(it2,it,i,17)>34d0)) then
                if(Sim1f(it2,it,i,10)>0.5d0) then
                    dum2=dum2+Sim1f(it2,it,i+1,16)*WeightActive(i)
                    dum4=dum4+(beta**i)*Sim1f(it2,it,i,11)*WeightActive(i)
                    dum3=dum3+1d0*WeightActive(i)
                end if

            end do
        end do

        dum2=dum2/dum3
        dum4=dum4/dum3

        Print *,'Ex-post Social welfare of 65-year old married women is',dum2
        Print *,'Ex-ante Social welfare of 65-year old married women is',dum4
        Print *,'Ex-post + Ex-ante Social welfare of 65-year old married women is',dum2+dum4
        Print *,'Number of 65-year old married women is',dum3

        !Social welfare single women at 65

        dum2=0d0
        dum3=0d0
        dum4=0d0

        i=T

        do it2=1,nsim2
            do it=1,nsim
                if(Sim1f(it2,it,i,10)<0.5d0) then
                    dum2=dum2+Sim1f(it2,it,i+1,16)*WeightActive(i)
                    dum4=dum4+(beta**i)*Sim1f(it2,it,i,11)*WeightActive(i)
                    dum3=dum3+1d0*WeightActive(i)
                end if
            end do
        end do

        dum2=dum2/dum3
        dum4=dum4/dum3

        Print *,'Ex-post Social welfare of 65-year old single women is',dum2
        Print *,'Ex-ante Social welfare of 65-year old single women is',dum4
        Print *,'Ex-post + Ex-ante Social welfare of 65-year old single women is',dum2+dum4
        Print *,'Number of 65-year old single women is',dum3

        !Budget constraint single women

        dum2=0d0
        dum3=0d0

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,T,10)<0.5d0) then
                        if(expR1f(it2,it,i,3)==1) then
                            dum2=dum2+abs(SimR1f(it2,it,i,18))*WeightRet(i)
                            dum3=dum3+1d0*WeightRet(i)
                        end if
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Error in budget constraint of retired single women is',dum2

        dum2=0d0
        dum3=0d0

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,T,10)<0.5d0) then
                        if((expR1f(it2,it,i,3)==2).AND.(SimR1f(it2,it,i,5)>minhours)) then
                            dum2=dum2+abs(SimR1f(it2,it,i,18))*WeightRet(i)
                            dum3=dum3+1d0*WeightRet(i)
                        end if
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Error in budget constraint of working single women is',dum2

        dum2=0d0
        dum3=0d0

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,T,10)<0.5d0) then
                        if((expR1f(it2,it,i,3)==2).AND.(SimR1f(it2,it,i,5)<minhours)) then
                            dum2=dum2+abs(SimR1f(it2,it,i,18))*WeightRet(i)
                            dum3=dum3+1d0*WeightRet(i)
                        end if
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Error in budget constraint of non-retired women working 0 hours',dum2

        !Single female budget before 65

        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,i,10)<0.5) then
                        if(Sim1f(it2,it,i,4)>1d-3) then
                            dum2=dum2+abs(Sim1f(it2,it,i,15))*WeightActive(i)
                        end if
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Error in budget constraint of working single females below 65 is',dum2

        dum2=0d0
        dum3=0d0

        do i=1,T

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1f(it2,it,i,10)<0.5) then
                        if(Sim1f(it2,it,i,4)<1d-3) then
                            dum2=dum2+abs(Sim1f(it2,it,i,15))*WeightActive(i)
                        end if
                        dum3=dum3+1d0*WeightActive(i)
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Error in budget constraint of unemployed single females below 65 is',dum2

        !Budget constraint single men

        dum2=0d0
        dum3=0d0

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,T,10)<0.5d0) then
                        if(expR1m(it2,it,i,3)==1) then
                            dum2=dum2+abs(SimR1m(it2,it,i,18))*WeightRet(i)
                            dum3=dum3+1d0*WeightRet(i)
                        end if
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Error in budget constraint of retired single men is',dum2

        dum2=0d0
        dum3=0d0

        do i=1,Tret-1

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,T,10)<0.5d0) then
                        if((expR1m(it2,it,i,3)==2).AND.(SimR1m(it2,it,i,5)>minhours)) then
                            dum2=dum2+abs(SimR1m(it2,it,i,18))*WeightRet(i)
                            dum3=dum3+1d0*WeightRet(i)
                        end if
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Error in budget constraint of working single men is',dum2

        dum2=0d0
        dum3=0d0

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,T,10)<0.5d0) then
                        if((expR1m(it2,it,i,3)==2).AND.(SimR1m(it2,it,i,5)<minhours)) then
                            dum2=dum2+abs(SimR1m(it2,it,i,18))*WeightRet(i)
                            dum3=dum3+1d0*WeightRet(i)
                        end if
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Error in budget constraint of non-retired men working 0 hours',dum2    

        !Budget constraint married

        dum2=0d0
        dum3=0d0

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim
                    if(Sim1m(it2,it,T,10)>0.5d0) then
                        if(SimR1m(it2,it,i,7)<0.001) then
                            dum2=dum2+abs(SimR1m(it2,it,i,18))*WeightRet(i)
                            dum3=dum3+1d0*WeightRet(i)
                        end if
                    end if
                end do
            end do

        end do

        dum2=dum2/dum3

        Print *,'Error in budget constraint of couples, both not working is',dum2

        dum2=0d0
        dum3=0d0

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim

                    if(Sim1m(it2,it,T,10)>0.5d0) then

                        if(SimR1m(it2,it,i,7)>0.001d0) then

                            dum2=dum2+abs(SimR1m(it2,it,i,18))*WeightRet(i)
                            if(SimR1m(it2,it,i,18)/=SimR1m(it2,it,i,18)) then
                                Print *,'SimR1m(it2,it,i,18) is',SimR1m(it2,it,i,18)
                                Print *,'it2 is',it2
                                Print *,'it2 is',it
                                Print *,'i is',i
                                PAUSE
                            end if
                            dum3=dum3+1d0*WeightRet(i)

                        end if

                    end if

                end do
            end do

        end do

        Print *,'dum2 is',dum2
        Print *,'dum3 is',dum3

        dum2=dum2/dum3

        Print *,'Error in budget constraint of couples with some labor income is',dum2

        dum2=0d0
        dum3=0d0    

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim

                    if(Sim1f(it2,it,T,10)>0.5d0) then

                        if(expR1f(it2,it,i,3)==2) then

                            dum2=dum2+1d0*WeightRet(i)


                        end if

                        dum3=dum3+1d0*WeightRet(i)

                    end if

                end do
            end do

        end do    

        !Print *,'dum2 is',dum2
        !    Print *,'dum3 is',dum3

        dum2=dum2/dum3

        Print *,'Fraction of non-retired married women above 65 is',dum2

        dum2=0d0
        dum3=0d0    

        do i=1,Tret

            do it2=1,nsim2
                do it=1,nsim

                    if(Sim1m(it2,it,T,10)>0.5d0) then

                        if(expR1m(it2,it,i,3)==2) then

                            dum2=dum2+1d0*WeightRet(i)


                        end if

                        if((SimR1m(it2,it,i,5)>0d0).AND.(expR1m(it2,it,i,3)==1)) then

                            Print *,'expR1m(it2,it,i,3) is',expR1m(it2,it,i,3)
                            Print *,'labor supply is',SimR1m(it2,it,i,5)
                            Print *,'it2 is',it2
                            Print *,'it is',it
                            Print *,'i is',i

                            PAUSE

                        end if

                        dum3=dum3+1d0*WeightRet(i)

                    end if

                end do
            end do

        end do    

        !Print *,'dum2 is',dum2
        !    Print *,'dum3 is',dum3

        dum2=dum2/dum3

        Print *,'Fraction of non-retired married men above 65 is',dum2

        open(21, file='share_single.txt')
    
        do i=1,T
            share_single(i,1) = count( Sim1m(:,:,i,10)<0.5d0 ) / dble(nsim2*nsim)
            do ia = 1, na
                !dum2a(ia,1)=dum2a(ia,1)/dum3a(ia)
                dum2 = count( Sim1m(:,:,i,10)<0.5d0 .and. exp1m(:,:,i,2) == ia )
                dum3 = count( exp1m(:,:,i,2) == ia )
                dum2a(ia,1) = dum2 / dum3
                share_single(i,1+ia) = dum2a(ia,1)
            end do
            write(21, '(i2,<na+1>f12.6)') i, share_single(i,:)
        end do
        close(21)  

    end subroutine Statistics

end module Summarize_Simulation
