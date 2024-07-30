module PolicyFunctions

    use Model_Parameters
    !use bspline_oo_module
    !use pppack

    implicit none
    ! Policy functions for active period of life
    real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable :: v !,vdum
    real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable, target :: v_aux ! Pseudo-value fn for married proposed by Lars, that tracks separately continuation values for single men and women
    real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable, target :: ev_aux ! Pseudo-value fn for married proposed by Lars, that tracks separately continuation values for single men and women

    ! LFP regime specific functions for married
    real(8), dimension (:,:,:,:,:,:,:,:,:,:,:,:), allocatable :: c_lfp, k_lfp, nm_lfp, nf_lfp
    real(8), dimension (:,:,:,:,:,:,:,:,:,:,:,:,:), allocatable :: v_lfp
    real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable :: lfpm, lfpf  

    real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable :: ev
    real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable :: c !,cdum
    real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable :: k !,gkdum
    real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable :: nm,nf !,nmdum,nfdum
    real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable :: ev_spln_coefs
    !real(8), dimension (:,:,:,:,:,:,:,:,:,:), allocatable :: v_spln_coefs, v_spln_coefs_kdim
    
    ! LFP regime specific functions for single
    real(8), dimension (:,:,:,:,:,:,:,:), allocatable :: cs_lfp, ks_lfp, vs_lfp    
    real(8), dimension (:,:,:,:,:,:,:), allocatable :: lfps    

    real(8), dimension (:,:,:,:,:,:,:), allocatable :: vs
    real(8), dimension (:,:,:,:,:,:,:), allocatable :: evs
    real(8), dimension (:,:,:,:,:,:,:), allocatable :: evm !Expected value function in the case of marriage next period
    real(8), dimension (:,:,:,:,:,:,:), allocatable :: cs, edcs, Uprimes !, edcs_spln_coefs
    real(8), dimension (:,:,:,:,:,:,:), allocatable :: ks
    real(8), dimension (:,:,:,:,:,:,:), allocatable :: ns
    !real(8), dimension (:,:,:,:,:,:,:), allocatable :: evs_spln_coefs, vs_spln_coefs, vs_spln_coefs_kdim
    !real(8), dimension (:,:,:,:,:,:,:), allocatable :: evm_spln_coefs !Coefficients to interpolate expected value function in the case of marriage next period
    real(8), dimension (:,:,:,:,:,:), allocatable :: fpartner,fpartnerdum,fpartnerdum2
    real(8), dimension (:,:,:,:,:,:), allocatable :: mpartner,mpartnerdum,mpartnerdum2
    real(8), dimension (:,:), allocatable :: ability_prob
    real(8), dimension (:,:,:), allocatable :: av_earnings

    real(8), dimension (:), allocatable :: c_grid
    real(8), dimension (:), allocatable :: wage_grid
    real(8), dimension (:), allocatable :: k_grid
    real(8), dimension (:,:), allocatable, target :: exp_grid
    !real(8), dimension (:), allocatable :: K_KNOT
    !real(8), dimension (:,:), allocatable :: EXP_KNOT
    real(8), dimension (:,:,:), allocatable :: laborm,laborf
    real(8), dimension (:,:), allocatable :: labormwork,laborfwork
    real(8), dimension (:,:), allocatable :: laborsinglem,laborsinglef
    ! Policy functions for retired
    !real(8), dimension (:,:,:), allocatable :: ev_spln_coefs_ret,evs_spln_coefs_ret
    !real(8), dimension (:,:,:), allocatable :: p_ev_spln_coefs_ret,p_evs_spln_coefs_ret  ! for testing only
    real(8), dimension (:,:,:,:,:,:,:,:,:,:,:,:), allocatable :: lfpm_ret, lfpf_ret
    real(8), dimension (:,:,:,:,:,:,:,:,:,:,:,:,:,:), allocatable :: c_ret_lfp, k_ret_lfp, nf_ret_lfp, nm_ret_lfp
    real(8), dimension (:,:,:,:,:,:,:,:,:,:,:,:), allocatable :: c_ret, v_ret, v_ret2 !, ev_ret_spln_coefs
    real(8), dimension (:,:,:,:,:,:,:,:,:,:,:,:), allocatable, target :: ev_ret
    real(8), dimension (:,:,:,:,:,:,:,:,:,:,:,:), allocatable :: edc_ret_spln_coefs, Uprime_ret
    real(8), dimension (:,:,:,:,:,:,:,:,:,:,:,:), allocatable, target :: edc_ret
    real(8), dimension (:,:,:,:,:,:,:,:,:,:,:,:), allocatable :: k_ret, nf_ret, nm_ret, retf, retm
    real(8), dimension (:,:,:,:,:,:,:,:), allocatable :: evs_ret, evs_ret_spln_coefs, Eulers_ret
    !real(8), dimension (:,:,:,:,:,:,:,:), allocatable :: cs_ret_spln_coefs
    real(8), dimension (:,:,:,:,:,:,:,:), allocatable :: edcs_ret_spln_coefs, Uprimes_ret
    real(8), dimension (:,:,:,:,:,:,:,:), allocatable, target :: edcs_ret, vs_ret, vs_ret2, cs_ret 
    real(8), dimension (:,:,:,:,:,:,:,:), allocatable :: ks_ret, ns_ret, Incomes_ret, rets
    real(8), dimension (:,:,:,:,:,:,:,:), allocatable, target :: lfps_ret      
    real(8), dimension (:,:,:,:,:,:,:,:,:), allocatable :: ks_ret_lfp, ns_ret_lfp, cs_ret_lfp    
    real(8), dimension (:), allocatable :: break

    ! New splines:        
    integer, parameter :: kx = 3
    integer, parameter :: ky = 3
    integer, parameter :: kz = 3
    real(8) :: tx(nk+kx)
    real(8) :: ty(nexp+ky,T+Tret)
    real(8) :: tz(nexp+kz,T+Tret)
    type bs_coefs_2d
        real(8) :: coefs(nk, nexp)
    end type bs_coefs_2d
    type bs_coefs_3d
        real(8) :: coefs(nk, nexp, nexp)
    end type bs_coefs_3d

    type(bs_coefs_2d) :: evs_ret_bspl(2,2,na,nu,nfc)
    type(bs_coefs_2d) :: edcs_ret_bspl(2,2,na,nu,nfc)

    type(bs_coefs_2d) :: evs_bspl(2,na,nu,nfc)
    type(bs_coefs_2d), target :: vs_bspl(2,na,nu,nfc)
    type(bs_coefs_2d) :: evm_bspl(2,na,nu,nfc)

    type(bs_coefs_3d) :: ev_ret_bspl(2,2,na,nu,na,nu,nfc,nfcm)
    type(bs_coefs_3d) :: edc_ret_bspl(2,2,na,nu,na,nu,nfc,nfcm)

    type(bs_coefs_3d) :: ev_bspl(na,nu,na,nu,nfc,nfcm)
    type(bs_coefs_3d), target :: v_bspl(na,nu,na,nu,nfc,nfcm)

    !type(bspline_2d) :: evs_ret_oo(2,2,na,nu,nfc)
    !type(bspline_2d) :: edcs_ret_oo(2,2,na,nu,nfc)
    !type(bspline_3d) :: ev_ret_oo(2,2,na,nu,na,nu,nfc,nfcm)
    !type(bspline_3d) :: edc_ret_oo(2,2,na,nu,na,nu,nfc,nfcm)

    contains


end module PolicyFunctions
