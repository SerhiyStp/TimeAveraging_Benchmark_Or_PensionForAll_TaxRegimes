module SortH
    implicit none
    
contains
    
    subroutine SortMultiArrayDescending(Array, ArraySorted, ix)
        real(8), intent(in) :: Array(:,:)
        real(8), intent(inout) :: ArraySorted(:,:)
        real(8), allocatable :: ArrayTmp(:,:)
        integer, intent(in) :: ix
        integer :: n1, n2
        integer :: it2, it3, it4
        real(8) :: dum2
        
        n1 = size(Array, 1)
        n2 = size(Array, 2)
        allocate(ArrayTmp(n1,n2))
        
        ArraySorted = 0d0
        ArrayTmp = 0d0
        
        do it2=1, n1
            it3=0
            dum2=0d0
            do
                it3=it3+1
                if(Array(it2,ix)>ArraySorted(it3,ix)) then
                    do it4=it3, n1
                        ArrayTmp(it4,:)=ArraySorted(it4,:)
                    end do
                    ArraySorted(it3,:)=Array(it2,:)
                    do it4=(it3+1), n1
                        ArraySorted(it4,:)=ArrayTmp(it4-1,:)
                    end do
                    dum2=1d0
                end if
                if(it3==n1) then
                    ArraySorted(it3,:)=Array(it2,:)
                    dum2=1d0
                end if
                if(dum2>0.5d0) exit
            end do

        end do        
        deallocate(ArrayTmp)
    end subroutine SortMultiArrayDescending
    
end module SortH