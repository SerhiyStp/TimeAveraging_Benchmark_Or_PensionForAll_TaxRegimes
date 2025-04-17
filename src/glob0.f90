module glob0

    real(8) :: dum,kdum,itg,netinc,wagem,wagef,SOL(1), XGUESS(1), XLB(1), XUB(1), GRAD(1,1), RHS(1),SOL2(2), XGUESS2(2), XLB2(2), XUB2(2), GRAD2(1,2),OBJ,FNORM
    integer :: MAXFCN,ITMAX,ind,ind2,ind3
    !$OMP THREADPRIVATE(dum,kdum,itg,netinc,wagem,wagef,SOL,XGUESS,XLB,XUB,GRAD,RHS,SOL2, XGUESS2, XLB2, XUB2,GRAD2,OBJ,FNORM,MAXFCN,ind,ind2,ind3)

end module glob0