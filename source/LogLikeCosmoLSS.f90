module LogLikeCosmoLSS_module
  use Interpolation
  implicit none
  integer, parameter :: mcp= KIND(1.d0)
  integer :: Curvature
  real(mcp) :: CPr
  Type :: myCMBparam
     real(mcp) :: H0
     real(mcp) :: h
     real(mcp) :: omdm
     real(mcp) :: omb
     real(mcp) :: omk
  end Type myCMBparam
  Type :: myMultipoleOverlaps
     real(mcp), dimension(:,:), allocatable ::  Convolution_matrix
     real(mcp) :: z_eff, k_min_theory, dk_theory, k_spacing_obs, k_min_obs
     INTEGER :: k_num_obs,size_convolution, k_num_conv
  end type myMultipoleOverlaps
  type this_type
 !    real(mcp), dimension(2,70) :: arraysjorig1,arraysjorig2,arraysjorig3,arraysjorig4
     real(mcp), dimension(2,70,4) :: arraysjfull
     real(mcp), dimension(2,29) :: arraysjlenslowz,arraysjlens2dfloz
     real(mcp), dimension(2,28) :: arraysjlenscmass,arraysjlens2dfhiz
     real(mcp), allocatable, dimension(:) :: xipm !4 tom-bins, 9 ang-bins, 2 for +/- gives 10*9*2 = 180
     real(mcp), allocatable, dimension(:,:) :: invcovxipm !, covxipm, invcovxipm,covxipminv !nzbins*(nzbins+1)*nangbins
!     real(mcp), dimension(9) :: thetacfhtini, thetaradcfhtini !nangbins
     real(mcp), dimension(58999) :: ellgentestarrini
     real(mcp), allocatable,dimension(:) :: maskelements
     real(mcp), dimension(9,58999) :: bes0arr,bes4arr,bes2arr
     integer :: size_cov,size_covmask,sizcov,sizcovpremask,klinesum,set_scenario,size_covallmask
     logical :: use_morell, use_rombint, use_conservative
     type(myMultipoleOverlaps) :: cmass_mp_overlap, lowz_mp_overlap, twodfloz_mp_overlap, twodfhiz_mp_overlap
     logical :: use_cmass_overlap, use_lowz_overlap, use_2dfloz_overlap, use_2dfhiz_overlap
     real(mcp), allocatable, dimension(:) :: exact_z
     integer :: num_z !, exact_z_index
  end type this_type

  Type(this_type) :: this
  type(myCMBparam) :: CMB
  type(SCubicSpline)           :: ConformalDistance, Hubble, GrowthRate, GrowthFunction
  type(TInterpGrid2D)          :: NonlinearPowerspectrum, LinearPowerspectrum

  contains

  subroutine Interpolation_Init_All(a, H, conf_dist, D, f, size_bg, kvec, zvec, pk, pk_nl, size_kvec, size_zvec)
    integer, intent(in) :: size_bg, size_kvec, size_zvec
    real(mcp), dimension(size_bg), intent(in) :: a, H, conf_dist, D, f
    real(mcp), dimension(size_kvec), intent(in) :: kvec
    real(mcp), dimension(size_zvec), intent(in) :: zvec
    real(mcp), dimension(size_kvec,size_zvec), intent(in) :: pk, pk_nl

    real(mcp) :: pk_test
    
    call ConformalDistance%SCubicSpline_Init(a,conf_dist,size_bg)
    call Hubble%SCubicSpline_Init(a,H,size_bg)
    call GrowthFunction%SCubicSpline_Init(a,D,size_bg)
    call GrowthRate%SCubicSpline_Init(a,f,size_bg)
    !write (*,*) 'Growth rate at z=0.56:',GrowthRate%Value(1d0/(1d0+0.56d0))
    !write (*,*) 'Hubble rate at z=0.56:',Hubble%Value(1d0/(1d0+0.56d0))
    call LinearPowerspectrum%Init(kvec,zvec,pk)
    call NonlinearPowerspectrum%Init(kvec,zvec,pk_nl)

    !pk_test = NonlinearPowerspectrum%Value(0.5d0, 0.2d0)
    !write (*,*) "pk_test:"
    !write (*,*) pk_test
    
  end subroutine Interpolation_Init_All

  subroutine Interpolation_Free_All()
    call SCubicSpline_Free(ConformalDistance)
    call SCubicSpline_Free(Hubble)
    call SCubicSpline_Free(GrowthFunction)
    call SCubicSpline_Free(GrowthRate)

    call TInterpGrid2D_Free(LinearPowerspectrum)
    call TInterpGrid2D_Free(NonlinearPowerspectrum)

  end subroutine Interpolation_Free_All

  subroutine set_mp_overlap(mp_overlap_ptr, ConMat, size_ConMat, intParams, realParams)
    integer, intent(in) :: size_ConMat
    real(mcp), dimension(size_ConMat,size_ConMat), intent(in) :: ConMat
    integer, intent(in), dimension(3) :: intParams
    real(mcp), intent(in), dimension(5) :: realParams
    type(myMultipoleOverlaps), intent(inout) :: mp_overlap_ptr
    !     real(mcp) :: z_eff, k_min_theory, dk_theory, k_spacing_obs, k_min_obs
    !       INTEGER :: k_num_obs,size_convolution, k_num_conv

    mp_overlap_ptr%z_eff = realParams(1)
    mp_overlap_ptr%k_min_theory = realParams(2)
    mp_overlap_ptr%dk_theory = realParams(3)
    mp_overlap_ptr%k_spacing_obs = realParams(4)
    mp_overlap_ptr%k_min_obs = realParams(5)

    mp_overlap_ptr%k_num_obs = intParams(1)
    mp_overlap_ptr%size_convolution = intParams(2)
    mp_overlap_ptr%k_num_conv = intParams(3)

    mp_overlap_ptr%Convolution_matrix = ConMat

  end subroutine set_mp_overlap

  
  function MPKPowerAt(k, z)
    !! Very important: k is in h/Mpc and Pk is in h/Mpc**3
    real(mcp) MPKPowerAt
    real(mcp), intent(in) :: k, z
    real(mcp) :: logk
    !I will assume interpolation is in log10(k) while the result is P(k), not log(P(k))
    logk = max(log10(k),LinearPowerspectrum%x(1))
    MPKPowerAt =LinearPowerspectrum%Value(logk,z)
  end function MPKPowerAt

  function NL_MPKPowerAt(k, z)
    real(mcp) NL_MPKPowerAt
    real(mcp), intent(in) :: k, z
    real(mcp) :: logk
    !I will assume interpolation is in log10(k) while the result is P(k), not log(P(k))
    logk = max(log10(k),LinearPowerspectrum%x(1))
    NL_MPKPowerAt =NonlinearPowerspectrum%Value(logk,z)
  end function NL_MPKPowerAt
    
  function ComovingRadialDistance(z)
    real(mcp) ComovingRadialDistance
    real(mcp), intent(in) :: z

    ComovingRadialDistance = ConformalDistance%Value(1d0/(1d0+z))

  end function ComovingRadialDistance

  function Hofz(z)
    !!non-comoving Hubble in MPC units, divide by MPC_in_sec to get in SI units
    real(mcp) Hofz, dtauda,a
    real(mcp), intent(in) :: z

    Hofz = Hubble%Value(1d0/(1d0+z))
    !write(*,*) z,'H(z)',Hubble%Value(1d0/(1d0+z))
    
  end function Hofz
  
  function AngularDiameterDistance(z)
    !This is the physical (non-comoving) angular diameter distance in Mpc
    real(mcp) AngularDiameterDistance
    real(mcp), intent(in) :: z

    AngularDiameterDistance = CPr/(1d0+z)*rofChi(ComovingRadialDistance(z) /CPr)

  end function AngularDiameterDistance

  function AngularDiameterDistance2(z1, z2) ! z1 < z2
    !From http://www.slac.stanford.edu/~amantz/work/fgas14/#cosmomc
    real(mcp) AngularDiameterDistance2
    real(mcp), intent(in) :: z1, z2

    AngularDiameterDistance2 = CPr/(1+z2)*rofChi(ComovingRadialDistance(z2)/CPr - ComovingRadialDistance(z1)/CPr)

  end function AngularDiameterDistance2

  function rofChi(Chi) !sinh(chi) for open, sin(chi) for closed.
    real(mcp) Chi,rofChi

    if (Curvature .eq. 1) then
        rofChi=sin(chi)
    else if (Curvature .eq. -1) then
        rofChi=sinh(chi)
    else
        rofChi=chi
    endif
  end function rofChi

  function f_K(x)
    real(mcp) :: f_K
    real(mcp), intent(in) :: x
    f_K = CPr*rofChi(x/CPr)
  end function f_K










  !----------------------------------------------------------------------------------------------------------
    ! Main likelihood calculation
    !----------------------------------------------------------------------------------------------------------
  function CosmoLSS_LnLike(DataParams)
    use omp_lib
    implicit none
    real(mcp), intent(in) :: DataParams(19)
!    class(CosmoLSSLikelihood) :: this
!    class(CMBParams) :: CMB
    !    class(TCosmoTheoryPredictions), target :: Theory
    real(mcp), parameter :: c = 2.99792458d8
    type my_type
        !type(TCosmoTheoryPredictions), allocatable :: theoryomp
        real(mcp) :: pisjomp,h0omp,homp,omdmomp,ombomp,omkomp,distlensflatomp,mumin2vomp,mumaxomp,tol_erroromp,ampiayesomp,redziayesomp,lumiayesomp,nellbinsomp,bbiaz0yesomp,bbiaz1yesomp,bbiaz2yesomp,bbiaz3yesomp,bbiazconstomp
        real(mcp), allocatable, dimension(:) :: ellarromp, aphotarromp,exact_zmodomp,lumarromp,exact_zmodomplenslowz,exact_zmodomplenscmass
        integer, allocatable, dimension(:) :: momp,m1arromp,m2arromp !exact_z_index,
        integer :: wchooseomp,ellgenomp,bnumomp,bnumomplens,uselensomp,kline,m1omp,m2omp,binwomp,wtrapmaxomp,wtrapmaxlensomp,setscenarioomp
        real(mcp) :: mumax1omp,mumax2omp,mumaxlensomp,mumaxwomp,mumaxcrossomplowz,mumin2vcrossomplowz,mumaxcrossompcmass,mumin2vcrossompcmass
        real(mcp), allocatable, dimension(:,:) :: weightarromp,weightpsarromp,weightnoarromp,weightarromplowz,weightarrompcmass,weightpsarromplowz,weightpsarrompcmass,weightnoarromplenslowz,weightnoarromplenscmass,weightarromp2dfloz,weightarromp2dfhiz,weightpsarromp2dfloz,weightpsarromp2dfhiz,weightnoarromplens2dfloz,weightnoarromplens2dfhiz
        real(mcp), allocatable, dimension(:) :: weightarromplenslowz,weightarromplenscmass,weightarromplens2dfloz,weightarromplens2dfhiz
        type(SCubicSpline),allocatable, dimension(:) :: psourcetypeomp
        type(SCubicSpline), allocatable :: plenstypeomplowz,plenstypeompcmass,plenstypeomp2dfloz,plenstypeomp2dfhiz
        type(SCubicSpline), allocatable :: gcubicomp
        Type(TCubicSpline), allocatable :: clinterptypeomp,clinterptypeiiomp,clinterptypegiomp,clinterptypecrosslowzomp,clinterptypecrosscmassomp,clinterptypecross2dflozomp,clinterptypecross2dfhizomp,clinterptypecrosslowzgiomp,clinterptypecrosscmassgiomp,clinterptypecross2dflozgiomp,clinterptypecross2dfhizgiomp !keep this as tcubicspline
    end type my_type
    real(mcp) rombint
    real(mcp) rombint_obj
    external rombint
    external rombint_obj
    type(my_type) obj
    real(mcp) :: CosmoLSS_LnLike,ampiayes,redziayes,lumiayes
    real(mcp) :: tol_error
    real(mcp) :: sigma_v1,sigma_v2,b1_bias_lowZ,b1_bias_highZ,bbiaz0yes,bbiaz1yes,bbiaz2yes,bbiaz3yes
    real(mcp) :: a1phot,a2phot,a3phot,a4phot
    real(mcp) :: mumax,mumin,mumin2v, pisj,likesj,chisqsj,chisqsjsellentin,start,finish,mumax1,mumax2,mumaxw,morellfac,ckmz,mumaxcrosslowz,mumin2vcrosslowz,mumaxcrosscmass,mumin2vcrosscmass,startrsd
    real(mcp), allocatable, dimension(:,:,:) :: nonlinspec
    INTEGER   :: zim,bbint,bbintigen,bbintmera,bbintmeraigen
    real(mcp), allocatable, dimension(:) :: xiplus,xiplusii,xiplusgi,xicrosslowz,xicrosscmass,xicrosslowzgi,xicrosscmassgi,xicross2dfloz,xicross2dfhiz,xicross2dflozgi,xicross2dfhizgi
    real(mcp), allocatable, dimension(:) :: ximinus,ximinusii,ximinusgi
    real(mcp), allocatable, dimension(:) :: xiplusminus, xiplusminusgg,xiplusminusii,xiplusminusgi,dataminustheory,dataminustheoryhihi,dataminustheoryhoho,finaltheory
    real(mcp), allocatable, dimension(:) :: xicrossarrlowz,xicrossarrlowzgi,xicrossarrcmass,xicrossarrcmassgi,xicrossarr2dfloz,xicrossarr2dflozgi,xicrossarr2dfhiz,xicrossarr2dfhizgi
    real(mcp), allocatable, dimension(:,:) :: weightarr,weightpsarr,weightnoarr,weightnoarrlenslowz,weightnoarrlenscmass,weightarrlowz,weightarrcmass,weightpsarrlowz,weightpsarrcmass,weightnoarrlens2dfloz,weightnoarrlens2dfhiz,weightarr2dfloz,weightarr2dfhiz,weightpsarr2dfloz,weightpsarr2dfhiz
    real(mcp), allocatable, dimension(:) :: weightarrlenslowz,weightarrlenscmass,weightarrlens2dfloz,weightarrlens2dfhiz
    real(mcp), allocatable, dimension(:) :: m1arr,m2arr,binwarr,exact_zmodlenslowz,exact_zmodlenscmass,exact_zmodlens2dfloz,exact_zmodlens2dfhiz
    real(mcp), allocatable, dimension(:) :: ellarr,trapzarr,trapzarrii,trapzarrgi,trapzarrcrosslowz,trapzarrcrosscmass,trapzarrcross2dfloz,trapzarrcross2dfhiz,trapzarrcrosslowzgi,trapzarrcrosscmassgi,trapzarrcross2dflozgi,trapzarrcross2dfhizgi
    real(mcp), allocatable, dimension(:) :: mumax1arr,mumax2arr, aphotarr,garr,cfhtlum,lumarr,mnomp
    real(mcp), allocatable, dimension(:,:) :: clarr,clarrii,clarrgi,clarrcrosslowz,clarrcrosscmass,clarrcross2dfloz,clarrcross2dfhiz,clarrcrosslowzgi,clarrcrosscmassgi,clarrcross2dflozgi,clarrcross2dfhizgi
    real(mcp), allocatable, dimension(:) :: clarrbig,ellgentestarr,clarrbigii,clarrbiggi,clarrbigcrosslowz,clarrbigcrosscmass,clarrbigcross2dfloz,clarrbigcross2dfhiz,clarrbigcrosslowzgi,clarrbigcrosscmassgi,clarrbigcross2dflozgi,clarrbigcross2dfhizgi
    integer :: hoj,bnum,gigi,how,howin
    integer :: loopee,yeye,ellgen,intnum,bin,biin,wchoose, ho, m1, m2, mint,ttt,wtrapmax,wtrapmaxlens,wttt,binw
    integer :: nangbins, nzbins, nellbins,  nell, nzp
    real(mcp) :: ellgentest
    logical :: use_cfht = .true.
    logical :: use_rombintz = .true. !if true use Romberg integration instead of trapezoidal for outer integral
!    logical :: use_rombintz = .false.
!    logical :: use_morellz = .true. !if true use 101 ell values, otherwise 31 ell values
     logical :: use_morellz = .false. !if true use 101 ell values, otherwise 31 ell values
    !!CMASS, LOWZ, 2dFLOZ, 2dFHIZ interpolation objects.
    Type(TCubicSpline) ::  P0_theory_spline_cmass_overlap,P2_theory_spline_cmass_overlap
    Type(TCubicSpline) ::  P0_theory_spline_lowz_overlap,P2_theory_spline_lowz_overlap
    Type(TCubicSpline) ::  P0_theory_spline_2dfloz_overlap,P2_theory_spline_2dfloz_overlap
    Type(TCubicSpline) ::  P0_theory_spline_2dfhiz_overlap,P2_theory_spline_2dfhiz_overlap
    !! CMASS overlap variables
    REAL(mcp), allocatable, dimension(:) :: P0P2_final_cmass_overlap,diff_cmass_overlap
    REAL(mcp) :: a_perp_cmass_overlap,a_par_cmass_overlap,k_val_overlap_cmass, z_cmass,sigma_v_cmass,b1_cmass,N_shot_cmass
    !! LOWZ overlap variables
    REAL(mcp), allocatable, dimension(:) :: P0P2_final_lowz_overlap,diff_lowz_overlap
    REAL(mcp) :: a_perp_lowz_overlap,a_par_lowz_overlap,k_val_overlap_lowz,z_lowz,sigv_lowz,b1_lowz,N_shot_lowz
    !! 2dfloz overlap variables
    REAL(mcp), allocatable, dimension(:) :: P0P2_final_2dfloz_overlap,diff_2dfloz_overlap
    REAL(mcp) :: a_perp_2dfloz_overlap,a_par_2dfloz_overlap,k_val_overlap_2dfloz,z_2dfloz,sigv_2dfloz,b1_2dfloz,N_shot_2dfloz
    !! 2dfhiz overlap variables
    REAL(mcp), allocatable, dimension(:) :: P0P2_final_2dfhiz_overlap,diff_2dfhiz_overlap
    REAL(mcp) :: a_perp_2dfhiz_overlap,a_par_2dfhiz_overlap,k_val_overlap_2dfhiz,z_2dfhiz,sigv_2dfhiz,b1_2dfhiz,N_shot_2dfhiz
    !    REAL(mcp) :: b1_2dfloz,b1_2dfhiz
    integer :: i
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    sigma_v_cmass = DataParams(1)
    b1_cmass = DataParams(2)
    N_shot_cmass = DataParams(3)
    sigv_lowz = DataParams(4)
    b1_lowz = DataParams(5)
    N_shot_lowz = DataParams(6)
    b1_2dfloz = DataParams(7)
    b1_2dfhiz = DataParams(8)
    sigv_2dfloz = DataParams(9)
    sigv_2dfhiz = DataParams(10)
    N_shot_2dfloz = DataParams(11)
    N_shot_2dfhiz = DataParams(12)
    ampiayes = DataParams(13)
    redziayes = DataParams(14)
    lumiayes = DataParams(15)
    a1phot = DataParams(16)
    a2phot = DataParams(17)
    a3phot = DataParams(18)
    a4phot = DataParams(19)



    CosmoLSS_LnLike = 0.0d0
    startrsd  = OMP_get_wtime()

    bbiaz0yes = b1_lowz
    bbiaz1yes = b1_cmass
    bbiaz2yes = b1_2dfloz
    bbiaz3yes = b1_2dfhiz

    obj%setscenarioomp = this%set_scenario

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(this%use_cmass_overlap) then
        allocate(P0P2_final_cmass_overlap(2*this%cmass_mp_overlap%k_num_obs)) !! Only want P0,P2 no measurement of P4
        allocate(diff_cmass_overlap(2*this%cmass_mp_overlap%k_num_obs))
    end if

    if(this%use_lowz_overlap) then
        allocate(P0P2_final_lowz_overlap(2*this%lowz_mp_overlap%k_num_obs)) !! Only want P0,P2 no measurement of P4
        allocate(diff_lowz_overlap(2*this%lowz_mp_overlap%k_num_obs))
    end if

    if(this%use_2dfloz_overlap) then
        allocate(P0P2_final_2dfloz_overlap(2*this%twodfloz_mp_overlap%k_num_obs)) !! Only want P0,P2 no measurement of P4
        allocate(diff_2dfloz_overlap(2*this%twodfloz_mp_overlap%k_num_obs))
    end if

    if(this%use_2dfhiz_overlap) then
        allocate(P0P2_final_2dfhiz_overlap(2*this%twodfhiz_mp_overlap%k_num_obs)) !! Only want P0,P2 no measurement of P4
        allocate(diff_2dfhiz_overlap(2*this%twodfhiz_mp_overlap%k_num_obs))
    end if


    !-----------------------------------------------------------------------------------------
    ! Compute theory
    !-----------------------------------------------------------------------------------------
    !write (*,*) 'Computing theory...'

    z_cmass = this%cmass_mp_overlap%z_eff
    z_lowz = this%lowz_mp_overlap%z_eff
    z_2dfloz = this%twodfloz_mp_overlap%z_eff
    z_2dfhiz = this%twodfhiz_mp_overlap%z_eff


    !write (*,*) 'Debug zone:'
    !write (*,*) CMB%H0, ComovingRadialDistance(0.5d0), AngularDiameterDistance(0.5d0), AngularDiameterDistance2(0.25d0, 1.75d0), MPKPowerAt(0.5d0,0.5d0), NL_MPKPowerAt(0.5d0,0.5d0)

    !write(*,*) 'use overlaps:',this%use_cmass_overlap,this%use_lowz_overlap,this%use_2dfloz_overlap,this%use_2dfhiz_overlap,this%twodfhiz_mp_overlap%k_num_obs,this%cmass_mp_overlap%k_num_obs,this%lowz_mp_overlap%k_num_obs,this%twodfloz_mp_overlap%k_num_obs
    !write (*,*) this%arraysjfull(1,:,1)
    !write (*,*) 'MPK test:', NL_MPKPowerAt(0.1d0*0.7d0,0.5d0), MPKPowerAt(0.1d0*0.7d0,0.5d0)
    
    !OPEN(UNIT=12, FILE="aoutput.txt", ACTION="write", STATUS="replace")
    !DO i=1,1000
    !   WRITE(12,*) 0.001*i, MPKPowerAt(0.001d0*i,0.0d0)
    !END DO
    !CLOSE(UNIT=12)
    !stop

    
    if(this%use_cmass_overlap) call Get_CMASS_overlap_TheoryMPK(P0P2_final_cmass_overlap)
    if(this%use_lowz_overlap) call Get_LOWZ_overlap_TheoryMPK(P0P2_final_lowz_overlap)
    if(this%use_2dfloz_overlap) call Get_2dfloz_overlap_TheoryMPK(P0P2_final_2dfloz_overlap)
    if(this%use_2dfhiz_overlap) call Get_2dfhiz_overlap_TheoryMPK(P0P2_final_2dfhiz_overlap)


    !---------------------------------------------------------------------------------------------------------
    ! Compute Delta = Obs - Theory data vector for CMASS overlap, LOWZ overlap, 2dfloz/2dfhiz overlap
    !----------------------------------------------------------------------------------------------------------

    if((obj%setscenarioomp == 3) .and. this%use_cmass_overlap) then
        diff_cmass_overlap(:) = P0P2_final_cmass_overlap(:)
    end if
    if((obj%setscenarioomp == 3) .and. this%use_lowz_overlap) then
        diff_lowz_overlap(:) = P0P2_final_lowz_overlap(:)
    end if
    if((obj%setscenarioomp == 3) .and. this%use_2dfloz_overlap) then
        diff_2dfloz_overlap(:) = P0P2_final_2dfloz_overlap(:)
    end if
    if((obj%setscenarioomp == 3) .and. this%use_2dfhiz_overlap) then
        diff_2dfhiz_overlap(:) = P0P2_final_2dfhiz_overlap(:)
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    use_morellz = this%use_morell
    use_rombintz = this%use_rombint
    nangbins = 9
    nzbins = 4
!    if(use_morellz == .false.) nellbins = 31
!    if(use_morellz == .true.) nellbins = 101
    if (use_morellz) then
       nellbins = 101
    else
       nellbins = 31
    end if
    obj%nellbinsomp = nellbins
    nell = 58999
    nzp = 70
    wtrapmax = 37
    wtrapmaxlens = 17 !assuming wtrapmaxlens < wtrapmax
    obj%wtrapmaxomp = wtrapmax
    obj%wtrapmaxlensomp = wtrapmaxlens

    !write (*,*) 'Allocating memory...'

!    allocate(obj%theoryomp)
    allocate(obj%ellarromp(nellbins))
!    allocate(obj%exact_z_index(this%num_z))
    allocate(obj%momp(2))
    allocate(obj%psourcetypeomp(nzbins))
    allocate(obj%plenstypeomplowz)
    allocate(obj%plenstypeompcmass)
    allocate(obj%plenstypeomp2dfloz)
    allocate(obj%plenstypeomp2dfhiz)
    allocate(obj%aphotarromp(nzbins))
    allocate(obj%lumarromp(nzbins))
    allocate(obj%exact_zmodomp(wtrapmax))
    allocate(obj%exact_zmodomplenslowz(wtrapmaxlens))
    allocate(obj%exact_zmodomplenscmass(wtrapmaxlens))
    allocate(obj%weightarromp(nzbins,wtrapmax-1))
    allocate(obj%weightarromplowz(nzbins,wtrapmaxlens-1)) !same as weightarromp just evaluated at lowz redshifts
    allocate(obj%weightarrompcmass(nzbins,wtrapmaxlens-1)) !same as weightarromp just evaluated at cmass redshifts
    allocate(obj%weightarromp2dfloz(nzbins,wtrapmaxlens-1)) !same as weightarromp just evaluated at 2dfloz redshifts
    allocate(obj%weightarromp2dfhiz(nzbins,wtrapmaxlens-1)) !same as weightarromp just evaluated at 2dfhiz redshifts
    allocate(obj%weightpsarromplowz(nzbins,wtrapmaxlens-1)) !same as weightarromp just evaluated at lowz redshifts
    allocate(obj%weightpsarrompcmass(nzbins,wtrapmaxlens-1)) !same as weightarromp just evaluated at cmass redshifts
    allocate(obj%weightpsarromp2dfloz(nzbins,wtrapmaxlens-1)) !same as weightarromp just evaluated at 2dfloz redshifts
    allocate(obj%weightpsarromp2dfhiz(nzbins,wtrapmaxlens-1)) !same as weightarromp just evaluated at 2dfhiz redshifts
    allocate(obj%weightpsarromp(nzbins,wtrapmax-1))
    allocate(obj%weightarromplenslowz(wtrapmaxlens-1))
    allocate(obj%weightarromplenscmass(wtrapmaxlens-1))
    allocate(obj%weightarromplens2dfloz(wtrapmaxlens-1))
    allocate(obj%weightarromplens2dfhiz(wtrapmaxlens-1))
    allocate(obj%m1arromp(nzbins*(nzbins+1)/2))
    allocate(obj%m2arromp(nzbins*(nzbins+1)/2))
    allocate(obj%weightnoarromp(nellbins,wtrapmax-1))
    allocate(obj%weightnoarromplenslowz(nellbins,wtrapmaxlens-1))
    allocate(obj%weightnoarromplenscmass(nellbins,wtrapmaxlens-1))
    allocate(obj%weightnoarromplens2dfloz(nellbins,wtrapmaxlens-1))
    allocate(obj%weightnoarromplens2dfhiz(nellbins,wtrapmaxlens-1))
    allocate(obj%clinterptypeomp)
    allocate(obj%clinterptypeiiomp)
    allocate(obj%clinterptypegiomp)
    allocate(obj%clinterptypecrosslowzomp)
    allocate(obj%clinterptypecrosscmassomp)
    allocate(obj%clinterptypecross2dflozomp)
    allocate(obj%clinterptypecross2dfhizomp)
    allocate(obj%clinterptypecrosslowzgiomp)
    allocate(obj%clinterptypecrosscmassgiomp)
    allocate(obj%clinterptypecross2dflozgiomp)
    allocate(obj%clinterptypecross2dfhizgiomp)
    allocate(obj%gcubicomp)

    allocate(clarrbig(nell))
    allocate(clarrbiggi(nell))
    allocate(clarrbigii(nell))
    allocate(clarrbigcrosslowz(nell))
    allocate(clarrbigcrosscmass(nell))
    allocate(clarrbigcross2dfloz(nell))
    allocate(clarrbigcross2dfhiz(nell))
    allocate(clarrbigcrosslowzgi(nell))
    allocate(clarrbigcrosscmassgi(nell))
    allocate(clarrbigcross2dflozgi(nell))
    allocate(clarrbigcross2dfhizgi(nell))
    allocate(ellgentestarr(nell))
    allocate(xiplusminus(nzbins*(nzbins+1)/2*nangbins*2))
    allocate(xiplusminusgg(nzbins*(nzbins+1)/2*nangbins*2))
    allocate(xiplusminusii(nzbins*(nzbins+1)/2*nangbins*2))
    allocate(xiplusminusgi(nzbins*(nzbins+1)/2*nangbins*2))
    allocate(dataminustheory(this%sizcov))
    allocate(dataminustheoryhihi(this%sizcov))
    allocate(dataminustheoryhoho(this%sizcov))
    allocate(finaltheory(this%sizcovpremask))
    allocate(xicrossarrlowz(nzbins*nangbins))
    allocate(xicrossarrlowzgi(nzbins*nangbins))
    allocate(xicrossarrcmass(nzbins*nangbins))
    allocate(xicrossarrcmassgi(nzbins*nangbins))
    allocate(xicrossarr2dfloz(nzbins*nangbins))
    allocate(xicrossarr2dflozgi(nzbins*nangbins))
    allocate(xicrossarr2dfhiz(nzbins*nangbins))
    allocate(xicrossarr2dfhizgi(nzbins*nangbins))
    allocate(clarr(nellbins,nzbins*(nzbins+1)/2))
    allocate(clarrii(nellbins,nzbins*(nzbins+1)/2))
    allocate(clarrgi(nellbins,nzbins*(nzbins+1)/2))
    allocate(clarrcrosslowz(nellbins,nzbins))
    allocate(clarrcrosscmass(nellbins,nzbins))
    allocate(clarrcross2dfloz(nellbins,nzbins))
    allocate(clarrcross2dfhiz(nellbins,nzbins))
    allocate(clarrcrosslowzgi(nellbins,nzbins))
    allocate(clarrcrosscmassgi(nellbins,nzbins))
    allocate(clarrcross2dflozgi(nellbins,nzbins))
    allocate(clarrcross2dfhizgi(nellbins,nzbins))
    allocate(ellarr(nellbins))
    allocate(xiplus(nangbins))
    allocate(xiplusii(nangbins))
    allocate(xiplusgi(nangbins))
    allocate(ximinus(nangbins))
    allocate(ximinusii(nangbins))
    allocate(ximinusgi(nangbins))
    allocate(xicrosslowz(nangbins))
    allocate(xicrosslowzgi(nangbins))
    allocate(xicrosscmass(nangbins))
    allocate(xicrosscmassgi(nangbins))
    allocate(xicross2dfloz(nangbins))
    allocate(xicross2dflozgi(nangbins))
    allocate(xicross2dfhiz(nangbins))
    allocate(xicross2dfhizgi(nangbins))
    allocate(mumax1arr(nzbins*(nzbins+1)/2))
    allocate(mumax2arr(nzbins*(nzbins+1)/2))
    allocate(aphotarr(nzbins))
    allocate(cfhtlum(nzbins))
    allocate(lumarr(nzbins))
    allocate(m1arr(nzbins*(nzbins+1)/2))
    allocate(m2arr(nzbins*(nzbins+1)/2))
    allocate(mnomp(2))
    allocate(weightarr(nzbins,wtrapmax-1))
    allocate(weightarrlowz(nzbins,wtrapmaxlens-1)) !same as weightarr just evaluated at lowz redshifts
    allocate(weightarrcmass(nzbins,wtrapmaxlens-1)) !same as weightarr just evaluated at cmass redshifts
    allocate(weightarr2dfloz(nzbins,wtrapmaxlens-1)) !same as weightarr just evaluated at 2dfloz redshifts
    allocate(weightarr2dfhiz(nzbins,wtrapmaxlens-1)) !same as weightarr just evaluated at 2dfhiz redshifts
    allocate(weightpsarrlowz(nzbins,wtrapmaxlens-1))
    allocate(weightpsarrcmass(nzbins,wtrapmaxlens-1))
    allocate(weightpsarr2dfloz(nzbins,wtrapmaxlens-1))
    allocate(weightpsarr2dfhiz(nzbins,wtrapmaxlens-1))
    allocate(weightpsarr(nzbins,wtrapmax-1))
    allocate(weightarrlenslowz(wtrapmaxlens-1))
    allocate(weightarrlenscmass(wtrapmaxlens-1))
    allocate(weightarrlens2dfloz(wtrapmaxlens-1))
    allocate(weightarrlens2dfhiz(wtrapmaxlens-1))
    allocate(binwarr(nzbins))
    allocate(trapzarr(wtrapmax-2))
    allocate(trapzarrii(wtrapmax-2))
    allocate(trapzarrgi(wtrapmax-2))
    allocate(trapzarrcrosslowz(wtrapmaxlens-2))
    allocate(trapzarrcrosscmass(wtrapmaxlens-2))
    allocate(trapzarrcross2dfloz(wtrapmaxlens-2))
    allocate(trapzarrcross2dfhiz(wtrapmaxlens-2))
    allocate(trapzarrcrosslowzgi(wtrapmaxlens-2))
    allocate(trapzarrcrosscmassgi(wtrapmaxlens-2))
    allocate(trapzarrcross2dflozgi(wtrapmaxlens-2))
    allocate(trapzarrcross2dfhizgi(wtrapmaxlens-2))
    allocate(weightnoarr(nellbins,wtrapmax-1))
    allocate(weightnoarrlenslowz(nellbins,wtrapmaxlens-1))
    allocate(weightnoarrlenscmass(nellbins,wtrapmaxlens-1))
    allocate(weightnoarrlens2dfloz(nellbins,wtrapmaxlens-1))
    allocate(weightnoarrlens2dfhiz(nellbins,wtrapmaxlens-1))
    allocate(garr(NonlinearPowerspectrum%ny))

    aphotarr(1) = a1phot
    aphotarr(2) = a2phot
    aphotarr(3) = a3phot
    aphotarr(4) = a4phot
    obj%aphotarromp(1) = a1phot
    obj%aphotarromp(2) = a2phot
    obj%aphotarromp(3) = a3phot
    obj%aphotarromp(4) = a4phot

    start  = OMP_get_wtime()
    cfhtlum = (/ 0.0174d0, 0.0694d0, 0.1518d0, 0.2182d0 /) !not used for KiDS
    if (use_cfht) lumarr = cfhtlum
    obj%lumarromp = lumarr

    obj%exact_zmodomp = this%exact_z
    tol_error = 0.005d0 !Numerical integration relative error, for Romberg
    pisj = 3.14159265359d0
    obj%pisjomp = 3.14159265359d0
    obj%h0omp = CMB%H0
    obj%homp = CMB%h
    obj%omdmomp = CMB%omdm
    obj%ombomp = CMB%omb
    obj%omkomp = CMB%omk
    !obj%theoryomp = Theory


    !!!trapezoid z-steps for lowz
    obj%exact_zmodomplenslowz(1) = 0.0d0
    obj%exact_zmodomplenslowz(2) = 0.15d0
    do zim=3,17
       obj%exact_zmodomplenslowz(zim) = 0.15d0 + (0.43d0-0.15d0)/15.0d0*(zim-2)
    end do
    obj%exact_zmodomplenslowz(17) = 0.43d0 !0.429999d0

    !!!trapezoid z-steps for cmass
    obj%exact_zmodomplenscmass(1) = 0.0d0
    obj%exact_zmodomplenscmass(2) = 0.43d0
    do zim=3,17
       obj%exact_zmodomplenscmass(zim) = 0.43 + (0.7d0-0.43d0)/15.0d0*(zim-2)
    end do
    obj%exact_zmodomplenscmass(17) = 0.7d0 !0.699999d0


    mumin=0.0d0
    mumin2v=0.025d0 !min integration redshift
    mumax=3.474999d0 !max integration redshift
    mumin2vcrosslowz=0.15 !0.150001d0 !min for lowz and 2dfloz
    mumaxcrosslowz=0.43d0 !0.4299999d0 !max for lowz and 2dfloz
    mumin2vcrosscmass=0.430001d0 !min for cmass and 2dfhiz
    mumaxcrosscmass=0.7d0 !0.6999999d0 !max for cmass and 2dfhiz

    !Generate array containing l-values where C(l) is evaluated
    !if(use_morellz == .false.) morellfac = 0.5
    if(use_morellz) then
       morellfac = 0.1
    else
       morellfac = 0.5
    end if
    do ellgen=1,nellbins
        if(ellgen < 10) ellarr(ellgen) = ellgen+1
        if(ellgen > 10 .or. ellgen == 10) ellarr(ellgen) = ellarr(ellgen-1)+morellfac*ellarr(ellgen-1)
    end do
    ellarr = int(ellarr)
    obj%ellarromp = ellarr
    obj%tol_erroromp = tol_error
    obj%mumin2vomp = mumin2v
    obj%mumaxomp = mumax
    obj%mumin2vcrossomplowz = mumin2vcrosslowz
    obj%mumaxcrossomplowz = mumaxcrosslowz
    obj%mumin2vcrossompcmass = mumin2vcrosscmass
    obj%mumaxcrossompcmass = mumaxcrosscmass


    obj%ampiayesomp = ampiayes
    obj%redziayesomp = redziayes
    obj%lumiayesomp = lumiayes
    obj%bbiaz0yesomp = bbiaz0yes
    obj%bbiaz1yesomp = bbiaz1yes
    obj%bbiaz2yesomp = bbiaz2yes
    obj%bbiaz3yesomp = bbiaz3yes
    !obj%exact_z_index = this%exact_z_index
    wchoose = 1
    obj%wchooseomp = 1

    !write(*,*) 'Initialising source redshift distributions...'
    !interpolation for source redshift distributions
    do ho=1,nzbins
        call obj%psourcetypeomp(ho)%Init(this%arraysjfull(1,:,ho),this%arraysjfull(2,:,ho),nzp)
    end do

    !interpolation for lens redshift distribution
    call obj%plenstypeomplowz%Init(this%arraysjlenslowz(1,:),this%arraysjlenslowz(2,:),29)
    call obj%plenstypeompcmass%Init(this%arraysjlenscmass(1,:),this%arraysjlenscmass(2,:),28)
    call obj%plenstypeomp2dfloz%Init(this%arraysjlens2dfloz(1,:),this%arraysjlens2dfloz(2,:),29)
    call obj%plenstypeomp2dfhiz%Init(this%arraysjlens2dfhiz(1,:),this%arraysjlens2dfhiz(2,:),28)

    obj%kline = 0

    !tomographic ordering
    m1arr = (/ 1, 1, 1, 1, 2, 2, 2, 3, 3, 4 /)
    m2arr = (/ 1, 2, 3, 4, 2, 3, 4, 3, 4, 4 /)

    !compute integration range accounting for photo-z error
    mint = 0
    do m1=1,nzbins
        do m2=1,nzbins
            if(m2 >= m1) then
                mint = mint + 1
                mumax1arr(mint) = mumax + aphotarr(m1arr(mint))
                mumax2arr(mint) = mumax + aphotarr(m2arr(mint))
            end if
        end do
    end do

    !write (*,*) 'For each redshift in NonlinearPowerspectrum, do some integral...'
    !requires the redshift in CosmologyTypes.f90 file to be set to z>3.5 in case of photo-z varying
    do gigi=1,NonlinearPowerspectrum%ny
       !!TT Check that NonlinearPowerspectrum%y(gigi) is maximum z.
       garr(gigi) = exp(rombint(sjgrowtha,1.0d0/(1.0d0+NonlinearPowerspectrum%y(gigi)),1.0d0,tol_error))
    end do
    !write (*,*) 'Initialise som interpolation routine'
    call obj%gcubicomp%Init(NonlinearPowerspectrum%y,garr,NonlinearPowerspectrum%ny)
    bbint = 0
    bbintigen = 0
    bbintmera = 0
    bbintmeraigen = 0

    !write (*,*) 'Outer integral, using Romberg:',use_rombintz

    if(use_rombintz) then !Romberg Integration for the outer integral
       
        !for nz tomographic bins, compute nz*(nz+1)/2 C(l)-combinations
        do bin = 1,nzbins*(nzbins+1)/2
            mint = 0
            do m1=1,nzbins
                do m2=1,nzbins
                    if(m2 >= m1) then 
                        mint = mint + 1
                        if(mint == bin) then
                            obj%m1omp=m1arr(mint)
                            obj%m2omp=m2arr(mint)
                            obj%momp(1)=m1arr(mint)
                            obj%momp(2)=m2arr(mint)
                        end if
                    end if
                end do
            end do
            mnomp = obj%momp
            bnum = bin
            obj%bnumomp = bin
            mumax1 = mumax1arr(bin)
            mumax2 = mumax2arr(bin)
            obj%mumax1omp = mumax1arr(bin)
            obj%mumax2omp = mumax2arr(bin)
            obj%mumaxlensomp = min(mumax1,mumax2)
    
            !!!##############COMPUTE LENSING POWER SPECTRA (GG, GI, II, Gg) AT DISTINCT L-VALUES. PARALLELIZATION EMPLOYED############
            !$OMP PARALLEL DO PRIVATE(ellgen) FIRSTPRIVATE(obj)
                do ellgen=1,nellbins
                    obj%ellgenomp = ellgen
                    clarr(ellgen,obj%bnumomp) = rombint_obj(obj,sjclsobjsf,1.0d0/(1.0d0+obj%mumaxlensomp),1.0d0/(1.0d0+obj%mumin2vomp),obj%tol_erroromp) !!good!!
                    clarrii(ellgen,obj%bnumomp) = rombint_obj(obj,sjclsiiobjsf,1.0d0/(1.0d0+obj%mumaxlensomp),1.0d0/(1.0d0+obj%mumin2vomp),obj%tol_erroromp) !good!!
                    clarrgi(ellgen,obj%bnumomp) = rombint_obj(obj,sjclsgiobjsf,1.0d0/(1.0d0+obj%mumaxlensomp),1.0d0/(1.0d0+obj%mumin2vomp),obj%tol_erroromp) !!good!!
                end do
            !$OMP END PARALLEL DO

            if((m1arr(bin) == m2arr(bin)) .and. (obj%setscenarioomp > 1)) then
                bbint = bbint+1
                obj%bnumomplens = bbint
                mumaxw = mumax1arr(bin)
                obj%mumaxwomp = mumaxw

                !$OMP PARALLEL DO PRIVATE(ellgen) FIRSTPRIVATE(obj)
                    do ellgen=1,nellbins
                        obj%ellgenomp = ellgen
                        obj%uselensomp = 0
                        obj%bbiazconstomp = obj%bbiaz0yesomp
                        clarrcrosslowz(ellgen,obj%bnumomplens) = rombint_obj(obj,sjclscrossobjsf,1.0d0/(1.0d0+obj%mumaxcrossomplowz),1.0d0/(1.0d0+obj%mumin2vcrossomplowz),obj%tol_erroromp)
                        clarrcrosslowzgi(ellgen,obj%bnumomplens) = rombint_obj(obj,sjclscrossgiobjsf,1.0d0/(1.0d0+obj%mumaxcrossomplowz),1.0d0/(1.0d0+obj%mumin2vcrossomplowz),obj%tol_erroromp)
                        !!!clarrcrosslowz(ellgen,obj%bnumomplens) = rombint_obj(obj,sjclscrossobj,obj%mumin2vcrossomplowz,obj%mumaxcrossomplowz,obj%tol_erroromp) !with z instead of sf
                        obj%uselensomp = 1
                        obj%bbiazconstomp = obj%bbiaz1yesomp
                        clarrcrosscmass(ellgen,obj%bnumomplens) = rombint_obj(obj,sjclscrossobjsf,1.0d0/(1.0d0+obj%mumaxcrossompcmass),1.0d0/(1.0d0+obj%mumin2vcrossompcmass),obj%tol_erroromp)
                        clarrcrosscmassgi(ellgen,obj%bnumomplens) = rombint_obj(obj,sjclscrossgiobjsf,1.0d0/(1.0d0+obj%mumaxcrossompcmass),1.0d0/(1.0d0+obj%mumin2vcrossompcmass),obj%tol_erroromp)
                        !!!clarrcrosscmass(ellgen,obj%bnumomplens) = rombint_obj(obj,sjclscrossobj,obj%mumin2vcrossompcmass,obj%mumaxcrossompcmass,obj%tol_erroromp) !with z instead of sf
                        obj%uselensomp = 2
                        obj%bbiazconstomp = obj%bbiaz2yesomp
                        clarrcross2dfloz(ellgen,obj%bnumomplens) = rombint_obj(obj,sjclscrossobjsf,1.0d0/(1.0d0+obj%mumaxcrossomplowz),1.0d0/(1.0d0+obj%mumin2vcrossomplowz),obj%tol_erroromp)
                        clarrcross2dflozgi(ellgen,obj%bnumomplens) = rombint_obj(obj,sjclscrossgiobjsf,1.0d0/(1.0d0+obj%mumaxcrossomplowz),1.0d0/(1.0d0+obj%mumin2vcrossomplowz),obj%tol_erroromp)
                        obj%uselensomp = 3
                        obj%bbiazconstomp = obj%bbiaz3yesomp
                        clarrcross2dfhiz(ellgen,obj%bnumomplens) = rombint_obj(obj,sjclscrossobjsf,1.0d0/(1.0d0+obj%mumaxcrossompcmass),1.0d0/(1.0d0+obj%mumin2vcrossompcmass),obj%tol_erroromp)
                        clarrcross2dfhizgi(ellgen,obj%bnumomplens) = rombint_obj(obj,sjclscrossgiobjsf,1.0d0/(1.0d0+obj%mumaxcrossompcmass),1.0d0/(1.0d0+obj%mumin2vcrossompcmass),obj%tol_erroromp)
                    end do
                !$OMP END PARALLEL DO
            end if
        end do
    end if 

    !write (*,*) 'Now trapezoidal...'
    
    !Trapezoid Integration for the outer integral!!! for the case use_rombintz=F !!! ------ do not use this option when varying photo-z params, in that case set use_rombintz = T
    if(.not. use_rombintz) then

        binwarr = (/ 1, 5, 8, 10 /) !4 tomographic bins, m1arr(i) = m2arr(i) for these indices --- need to change this array for different number of tomographic bins
        weightarr(:,:) = 0.0d0
        weightarrlowz(:,:) = 0.0d0
        weightarrcmass(:,:) = 0.0d0
        weightpsarrlowz(:,:) = 0.0d0
        weightpsarrcmass(:,:) = 0.0d0
        weightarrlenslowz(:) = 0.0d0
        weightarrlenscmass(:) = 0.0d0
        weightarr2dfloz(:,:) = 0.0d0
        weightarr2dfhiz(:,:) = 0.0d0
        weightpsarr2dfloz(:,:) = 0.0d0
        weightpsarr2dfhiz(:,:) = 0.0d0
        weightarrlens2dfloz(:) = 0.0d0
        weightarrlens2dfhiz(:) = 0.0d0
        weightpsarr(:,:) = 0.0d0
        weightnoarr(:,:) = 0.0d0
        weightnoarrlenslowz(:,:) = 0.0d0
        weightnoarrlenscmass(:,:) = 0.0d0
        weightnoarrlens2dfloz(:,:) = 0.0d0
        weightnoarrlens2dfhiz(:,:) = 0.0d0

        do binw = 1,nzbins
            mumaxw = mumax1arr(binwarr(binw))
            obj%mumaxwomp = mumaxw
            bnum = binwarr(binw)
            obj%bnumomp = bnum
            obj%m1omp=m1arr(binwarr(binw))
            obj%m2omp=m2arr(binwarr(binw))
            obj%momp(1)=m1arr(binwarr(binw))
            obj%momp(2)=m2arr(binwarr(binw))
            mnomp = obj%momp
            obj%binwomp = binw

            !$OMP PARALLEL DO PRIVATE(wttt) FIRSTPRIVATE(obj)
                do wttt=2,obj%wtrapmaxomp
                    weightarr(obj%binwomp,wttt-1) = sjclsobjonlyweight(obj,obj%exact_zmodomp(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                    weightpsarr(obj%binwomp,wttt-1) = sjclsiiobjonlyweight(obj,obj%exact_zmodomp(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                    if((wttt < (obj%wtrapmaxlensomp + 1)) .and. (obj%setscenarioomp > 1)) then
                        weightarrlowz(obj%binwomp,wttt-1) = sjclsobjonlyweight(obj,obj%exact_zmodomplenslowz(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                        weightarrcmass(obj%binwomp,wttt-1) = sjclsobjonlyweight(obj,obj%exact_zmodomplenscmass(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                        weightarr2dfloz(obj%binwomp,wttt-1) = weightarrlowz(obj%binwomp,wttt-1) !since same redshifts
                        weightarr2dfhiz(obj%binwomp,wttt-1) = weightarrcmass(obj%binwomp,wttt-1) !since same redshifts
                        weightpsarrlowz(obj%binwomp,wttt-1) = sjclsiiobjonlyweight(obj,obj%exact_zmodomplenslowz(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                        weightpsarrcmass(obj%binwomp,wttt-1) = sjclsiiobjonlyweight(obj,obj%exact_zmodomplenscmass(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                        weightpsarr2dfloz(obj%binwomp,wttt-1) = weightpsarrlowz(obj%binwomp,wttt-1) !since same redshifts
                        weightpsarr2dfhiz(obj%binwomp,wttt-1) = weightpsarrcmass(obj%binwomp,wttt-1) !since same redshifts
                    end if
                    if((obj%binwomp == 1) .and. (wttt < (obj%wtrapmaxlensomp+1)) .and. (obj%setscenarioomp > 1)) then
                        obj%uselensomp = 0
                        obj%bbiazconstomp = obj%bbiaz0yesomp
                        weightarrlenslowz(wttt-1) = sjclscrossobjonlyweight(obj,obj%exact_zmodomplenslowz(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                        obj%uselensomp = 1
                        obj%bbiazconstomp = obj%bbiaz1yesomp
                        weightarrlenscmass(wttt-1) = sjclscrossobjonlyweight(obj,obj%exact_zmodomplenscmass(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                        obj%uselensomp = 2
                        obj%bbiazconstomp = obj%bbiaz2yesomp
                        weightarrlens2dfloz(wttt-1) = sjclscrossobjonlyweight(obj,obj%exact_zmodomplenslowz(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                        obj%uselensomp = 3
                        obj%bbiazconstomp = obj%bbiaz3yesomp
                        weightarrlens2dfhiz(wttt-1) = sjclscrossobjonlyweight(obj,obj%exact_zmodomplenscmass(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                    end if
                end do
            !!!!OMP END PARALLEL DO
        end do

        do ellgen=1,nellbins
            obj%ellgenomp = ellgen
            !$OMP PARALLEL DO PRIVATE(wttt) FIRSTPRIVATE(obj)
                do wttt=2,obj%wtrapmaxomp
                    weightnoarr(obj%ellgenomp,wttt-1) = sjclsobjnoweight(obj,obj%exact_zmodomp(wttt))   !have to make sure mumaxw > this%exact_z(wttt)
                    if((wttt < (obj%wtrapmaxlensomp+1)) .and. (obj%setscenarioomp > 1)) then
                        weightnoarrlenslowz(obj%ellgenomp,wttt-1) = sjclscrossobjnoweight(obj,obj%exact_zmodomplenslowz(wttt))   !have to make sure mumaxw > this%exact_z(wttt)
                        weightnoarrlenscmass(obj%ellgenomp,wttt-1) = sjclscrossobjnoweight(obj,obj%exact_zmodomplenscmass(wttt))   !have to make sure mumaxw > this%exact_z(wttt)
                        weightnoarrlens2dfloz(obj%ellgenomp,wttt-1) = weightnoarrlenslowz(obj%ellgenomp,wttt-1)
                        weightnoarrlens2dfhiz(obj%ellgenomp,wttt-1) = weightnoarrlenscmass(obj%ellgenomp,wttt-1)
                    end if
                end do
            !!!!OMP END PARALLEL DO
        end do

        !for nz tomographic bins, compute nz*(nz+1)/2 C(l)-combinations
        do bin = 1,nzbins*(nzbins+1)/2
            if(m1arr(bin) == m2arr(bin)) then
                bbintigen = bbintigen+1
                obj%bnumomplens = bbintigen
            end if
            mint = 0
            do m1=1,nzbins
                do m2=1,nzbins
                    if(m2 >= m1) then 
                        mint = mint + 1
                        if(mint == bin) then
                            obj%m1omp=m1arr(mint)
                            obj%m2omp=m2arr(mint)
                            obj%momp(1)=m1arr(mint)
                            obj%momp(2)=m2arr(mint)
                        end if
                    end if
                end do
            end do
            mnomp = obj%momp
            bnum = bin
            obj%bnumomp = bin
            mumax1 = mumax1arr(bin)
            mumax2 = mumax2arr(bin)
            obj%mumax1omp = mumax1arr(bin)
            obj%mumax2omp = mumax2arr(bin)
            obj%mumaxlensomp = min(mumax1,mumax2)
            obj%weightarromp = weightarr
            obj%weightarromplowz = weightarrlowz
            obj%weightarrompcmass = weightarrcmass
            obj%weightpsarromplowz = weightpsarrlowz
            obj%weightpsarrompcmass = weightpsarrcmass
            obj%weightarromplenslowz = weightarrlenslowz
            obj%weightarromplenscmass = weightarrlenscmass
            obj%weightarromp2dfloz = weightarr2dfloz
            obj%weightarromp2dfhiz = weightarr2dfhiz
            obj%weightarromplens2dfloz = weightarrlens2dfloz
            obj%weightarromplens2dfhiz = weightarrlens2dfhiz
            obj%weightpsarromp2dfloz = weightpsarr2dfloz
            obj%weightpsarromp2dfhiz = weightpsarr2dfhiz
            obj%weightpsarromp = weightpsarr
            obj%weightnoarromp = weightnoarr
            obj%weightnoarromplenslowz = weightnoarrlenslowz
            obj%weightnoarromplenscmass = weightnoarrlenscmass
            obj%weightnoarromplens2dfloz = weightnoarrlens2dfloz
            obj%weightnoarromplens2dfhiz = weightnoarrlens2dfhiz
            obj%m1arromp = m1arr
            obj%m2arromp = m2arr

            !!!##############COMPUTE LENSING POWER SPECTRA (GG, GI, II, Gg) AT DISTINCT L-VALUES. PARALLELIZATION EMPLOYED############
            do ellgen=1,obj%nellbinsomp
                obj%ellgenomp = ellgen
                !$OMP PARALLEL DO PRIVATE(TTT) FIRSTPRIVATE(obj)
                    do ttt=2,36
                        trapzarr(ttt-1) = 0.5d0*(obj%weightnoarromp(obj%ellgenomp,ttt-1)*obj%weightarromp(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightarromp(obj%m2arromp(obj%bnumomp),ttt-1) + obj%weightnoarromp(obj%ellgenomp,ttt)*obj%weightarromp(obj%m1arromp(obj%bnumomp),ttt)*obj%weightarromp(obj%m2arromp(obj%bnumomp),ttt))*(obj%exact_zmodomp(ttt+1)-obj%exact_zmodomp(ttt))
                        trapzarrii(ttt-1) = 0.5d0*(obj%weightnoarromp(obj%ellgenomp,ttt-1)*obj%weightpsarromp(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightpsarromp(obj%m2arromp(obj%bnumomp),ttt-1) + obj%weightnoarromp(obj%ellgenomp,ttt)*obj%weightpsarromp(obj%m1arromp(obj%bnumomp),ttt)*obj%weightpsarromp(obj%m2arromp(obj%bnumomp),ttt))*(obj%exact_zmodomp(ttt+1)-obj%exact_zmodomp(ttt))
                        trapzarrgi(ttt-1) = 0.5d0*(obj%weightnoarromp(obj%ellgenomp,ttt-1)*(obj%weightarromp(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightpsarromp(obj%m2arromp(obj%bnumomp),ttt-1)+obj%weightarromp(obj%m2arromp(obj%bnumomp),ttt-1)*obj%weightpsarromp(obj%m1arromp(obj%bnumomp),ttt-1)) + obj%weightnoarromp(obj%ellgenomp,ttt)*(obj%weightarromp(obj%m1arromp(obj%bnumomp),ttt)*obj%weightpsarromp(obj%m2arromp(obj%bnumomp),ttt)+obj%weightarromp(obj%m2arromp(obj%bnumomp),ttt)*obj%weightpsarromp(obj%m1arromp(obj%bnumomp),ttt)))*(obj%exact_zmodomp(ttt+1)-obj%exact_zmodomp(ttt))
                    end do
                !$OMP END PARALLEL DO
                clarr(ellgen,obj%bnumomp)=sum(trapzarr)
                clarrii(ellgen,obj%bnumomp)=sum(trapzarrii)
                clarrgi(ellgen,obj%bnumomp)=sum(trapzarrgi)
                if((m1arr(bin) == m2arr(bin)) .and. (obj%setscenarioomp > 1)) then
                    !$OMP PARALLEL DO PRIVATE(ttt) FIRSTPRIVATE(obj)
                        do ttt=2,16
                            trapzarrcrosslowz(ttt-1) = 0.5d0*(obj%weightnoarromplenslowz(obj%ellgenomp,ttt-1)*obj%weightarromplowz(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightarromplenslowz(ttt-1) + obj%weightnoarromplenslowz(obj%ellgenomp,ttt)*obj%weightarromplowz(obj%m1arromp(obj%bnumomp),ttt)*obj%weightarromplenslowz(ttt))*(obj%exact_zmodomplenslowz(ttt+1)-obj%exact_zmodomplenslowz(ttt))
                            trapzarrcrosscmass(ttt-1) = 0.5d0*(obj%weightnoarromplenscmass(obj%ellgenomp,ttt-1)*obj%weightarrompcmass(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightarromplenscmass(ttt-1) + obj%weightnoarromplenscmass(obj%ellgenomp,ttt)*obj%weightarrompcmass(obj%m1arromp(obj%bnumomp),ttt)*obj%weightarromplenscmass(ttt))*(obj%exact_zmodomplenscmass(ttt+1)-obj%exact_zmodomplenscmass(ttt))
                            trapzarrcross2dfloz(ttt-1) = 0.5d0*(obj%weightnoarromplens2dfloz(obj%ellgenomp,ttt-1)*obj%weightarromp2dfloz(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightarromplens2dfloz(ttt-1) + obj%weightnoarromplens2dfloz(obj%ellgenomp,ttt)*obj%weightarromp2dfloz(obj%m1arromp(obj%bnumomp),ttt)*obj%weightarromplens2dfloz(ttt))*(obj%exact_zmodomplenslowz(ttt+1)-obj%exact_zmodomplenslowz(ttt))
                            trapzarrcross2dfhiz(ttt-1) = 0.5d0*(obj%weightnoarromplens2dfhiz(obj%ellgenomp,ttt-1)*obj%weightarromp2dfhiz(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightarromplens2dfhiz(ttt-1) + obj%weightnoarromplens2dfhiz(obj%ellgenomp,ttt)*obj%weightarromp2dfhiz(obj%m1arromp(obj%bnumomp),ttt)*obj%weightarromplens2dfhiz(ttt))*(obj%exact_zmodomplenscmass(ttt+1)-obj%exact_zmodomplenscmass(ttt))
                            trapzarrcrosslowzgi(ttt-1) = 0.5d0*(obj%weightnoarromplenslowz(obj%ellgenomp,ttt-1)*obj%weightpsarromplowz(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightarromplenslowz(ttt-1) + obj%weightnoarromplenslowz(obj%ellgenomp,ttt)*obj%weightpsarromplowz(obj%m1arromp(obj%bnumomp),ttt)*obj%weightarromplenslowz(ttt))*(obj%exact_zmodomplenslowz(ttt+1)-obj%exact_zmodomplenslowz(ttt))
                            trapzarrcrosscmassgi(ttt-1) = 0.5d0*(obj%weightnoarromplenscmass(obj%ellgenomp,ttt-1)*obj%weightpsarrompcmass(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightarromplenscmass(ttt-1) + obj%weightnoarromplenscmass(obj%ellgenomp,ttt)*obj%weightpsarrompcmass(obj%m1arromp(obj%bnumomp),ttt)*obj%weightarromplenscmass(ttt))*(obj%exact_zmodomplenscmass(ttt+1)-obj%exact_zmodomplenscmass(ttt))
                            trapzarrcross2dflozgi(ttt-1) = 0.5d0*(obj%weightnoarromplens2dfloz(obj%ellgenomp,ttt-1)*obj%weightpsarromp2dfloz(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightarromplens2dfloz(ttt-1) + obj%weightnoarromplens2dfloz(obj%ellgenomp,ttt)*obj%weightpsarromp2dfloz(obj%m1arromp(obj%bnumomp),ttt)*obj%weightarromplens2dfloz(ttt))*(obj%exact_zmodomplenslowz(ttt+1)-obj%exact_zmodomplenslowz(ttt))
                            trapzarrcross2dfhizgi(ttt-1) = 0.5d0*(obj%weightnoarromplens2dfhiz(obj%ellgenomp,ttt-1)*obj%weightpsarromp2dfhiz(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightarromplens2dfhiz(ttt-1) + obj%weightnoarromplens2dfhiz(obj%ellgenomp,ttt)*obj%weightpsarromp2dfhiz(obj%m1arromp(obj%bnumomp),ttt)*obj%weightarromplens2dfhiz(ttt))*(obj%exact_zmodomplenscmass(ttt+1)-obj%exact_zmodomplenscmass(ttt))
                        end do
                    !$OMP END PARALLEL DO
                    clarrcrosslowz(ellgen,obj%bnumomplens)=sum(trapzarrcrosslowz)
                    clarrcrosscmass(ellgen,obj%bnumomplens)=sum(trapzarrcrosscmass)
                    clarrcross2dfloz(ellgen,obj%bnumomplens)=sum(trapzarrcross2dfloz)
                    clarrcross2dfhiz(ellgen,obj%bnumomplens)=sum(trapzarrcross2dfhiz)
                    clarrcrosslowzgi(ellgen,obj%bnumomplens)=sum(trapzarrcrosslowzgi)
                    clarrcrosscmassgi(ellgen,obj%bnumomplens)=sum(trapzarrcrosscmassgi)
                    clarrcross2dflozgi(ellgen,obj%bnumomplens)=sum(trapzarrcross2dflozgi)
                    clarrcross2dfhizgi(ellgen,obj%bnumomplens)=sum(trapzarrcross2dfhizgi)
                end if
            end do
        end do
    end if
    !write (*, *) 'After Trapexoidal method, now initialise many interpolations...'
    
    !Interpolation for C(l)
    ellgentestarr = this%ellgentestarrini
    intnum=nellbins
    do biin = 1,nzbins*(nzbins+1)/2
        call obj%clinterptypeomp%Init(ellarr,clarr(:,biin),intnum)
        call obj%clinterptypeiiomp%Init(ellarr,clarrii(:,biin),intnum)
        call obj%clinterptypegiomp%Init(ellarr,clarrgi(:,biin),intnum)
        if((m1arr(biin) == m2arr(biin)) .and. (obj%setscenarioomp > 1)) then
            bbintmera = bbintmera+1
            call obj%clinterptypecrosslowzomp%Init(ellarr,clarrcrosslowz(:,bbintmera),intnum)
            call obj%clinterptypecrosscmassomp%Init(ellarr,clarrcrosscmass(:,bbintmera),intnum)
            call obj%clinterptypecross2dflozomp%Init(ellarr,clarrcross2dfloz(:,bbintmera),intnum)
            call obj%clinterptypecross2dfhizomp%Init(ellarr,clarrcross2dfhiz(:,bbintmera),intnum)
            call obj%clinterptypecrosslowzgiomp%Init(ellarr,clarrcrosslowzgi(:,bbintmera),intnum)
            call obj%clinterptypecrosscmassgiomp%Init(ellarr,clarrcrosscmassgi(:,bbintmera),intnum)
            call obj%clinterptypecross2dflozgiomp%Init(ellarr,clarrcross2dflozgi(:,bbintmera),intnum)
            call obj%clinterptypecross2dfhizgiomp%Init(ellarr,clarrcross2dfhizgi(:,bbintmera),intnum)
        end if
        ellgentest = 1.0d0
        !$OMP PARALLEL DO PRIVATE(ellgen) FIRSTPRIVATE(obj)
            do ellgen=1,58999
                clarrbig(ellgen) = obj%clinterptypeomp%Value(ellgen+1.0d0)
                clarrbigii(ellgen) = obj%clinterptypeiiomp%Value(ellgen+1.0d0)
                clarrbiggi(ellgen) = obj%clinterptypegiomp%Value(ellgen+1.0d0)
            end do
        !$OMP END PARALLEL DO
        if((m1arr(biin) == m2arr(biin)) .and. (obj%setscenarioomp > 1)) then
            ellgentest = 1.0d0
            !$OMP PARALLEL DO PRIVATE(ellgen) FIRSTPRIVATE(obj)
                do ellgen=1,58999
                    clarrbigcrosslowz(ellgen) = obj%clinterptypecrosslowzomp%Value(ellgen+1.0d0)
                    clarrbigcrosscmass(ellgen) = obj%clinterptypecrosscmassomp%Value(ellgen+1.0d0)
                    clarrbigcross2dfloz(ellgen) = obj%clinterptypecross2dflozomp%Value(ellgen+1.0d0)
                    clarrbigcross2dfhiz(ellgen) = obj%clinterptypecross2dfhizomp%Value(ellgen+1.0d0)
                    clarrbigcrosslowzgi(ellgen) = obj%clinterptypecrosslowzgiomp%Value(ellgen+1.0d0)
                    clarrbigcrosscmassgi(ellgen) = obj%clinterptypecrosscmassgiomp%Value(ellgen+1.0d0)
                    clarrbigcross2dflozgi(ellgen) = obj%clinterptypecross2dflozgiomp%Value(ellgen+1.0d0)
                    clarrbigcross2dfhizgi(ellgen) = obj%clinterptypecross2dfhizgiomp%Value(ellgen+1.0d0)
                end do
            !$OMP END PARALLEL DO
        end if

        xiplus(:) = 0.0d0
        ximinus(:) = 0.0d0
        xiplusii(:) = 0.0d0
        ximinusii(:) = 0.0d0
        xiplusgi(:) = 0.0d0
        ximinusgi(:) = 0.0d0
        if((m1arr(biin) == m2arr(biin)) .and. (obj%setscenarioomp > 1)) then
            xicrosslowz(:) = 0.0d0
            xicrosscmass(:) = 0.0d0
            xicross2dfloz(:) = 0.0d0
            xicross2dfhiz(:) = 0.0d0
             xicrosslowzgi(:) = 0.0d0
             xicrosscmassgi(:) = 0.0d0
             xicross2dflozgi(:) = 0.0d0
             xicross2dfhizgi(:) = 0.0d0
        end if
        hoj=0
        do yeye=1,nangbins
            do loopee=1,nell
                xiplus(yeye) = xiplus(yeye) + 1.0d0/pisj/2.0d0*clarrbig(loopee)*ellgentestarr(loopee)*this%bes0arr(yeye,loopee)
                ximinus(yeye) = ximinus(yeye) + 1.0d0/pisj/2.0d0*clarrbig(loopee)*ellgentestarr(loopee)*this%bes4arr(yeye,loopee)
                xiplusii(yeye) = xiplusii(yeye) + 1.0d0/pisj/2.0d0*clarrbigii(loopee)*ellgentestarr(loopee)*this%bes0arr(yeye,loopee)
                ximinusii(yeye) = ximinusii(yeye) + 1.0d0/pisj/2.0d0*clarrbigii(loopee)*ellgentestarr(loopee)*this%bes4arr(yeye,loopee)
                xiplusgi(yeye) = xiplusgi(yeye) + 1.0d0/pisj/2.0d0*clarrbiggi(loopee)*ellgentestarr(loopee)*this%bes0arr(yeye,loopee)
                ximinusgi(yeye) = ximinusgi(yeye) + 1.0d0/pisj/2.0d0*clarrbiggi(loopee)*ellgentestarr(loopee)*this%bes4arr(yeye,loopee)
                if((m1arr(biin) == m2arr(biin)) .and. (obj%setscenarioomp > 1)) then
                    xicrosslowz(yeye) = xicrosslowz(yeye) + 1.0d0/pisj/2.0d0*clarrbigcrosslowz(loopee)*ellgentestarr(loopee)*this%bes2arr(yeye,loopee)
                    xicrosscmass(yeye) = xicrosscmass(yeye) + 1.0d0/pisj/2.0d0*clarrbigcrosscmass(loopee)*ellgentestarr(loopee)*this%bes2arr(yeye,loopee)
                    xicross2dfloz(yeye) = xicross2dfloz(yeye) + 1.0d0/pisj/2.0d0*clarrbigcross2dfloz(loopee)*ellgentestarr(loopee)*this%bes2arr(yeye,loopee)
                    xicross2dfhiz(yeye) = xicross2dfhiz(yeye) + 1.0d0/pisj/2.0d0*clarrbigcross2dfhiz(loopee)*ellgentestarr(loopee)*this%bes2arr(yeye,loopee)
                    xicrosslowzgi(yeye) = xicrosslowzgi(yeye) + 1.0d0/pisj/2.0d0*clarrbigcrosslowzgi(loopee)*ellgentestarr(loopee)*this%bes2arr(yeye,loopee)
                    xicrosscmassgi(yeye) = xicrosscmassgi(yeye) + 1.0d0/pisj/2.0d0*clarrbigcrosscmassgi(loopee)*ellgentestarr(loopee)*this%bes2arr(yeye,loopee)
                    xicross2dflozgi(yeye) = xicross2dflozgi(yeye) + 1.0d0/pisj/2.0d0*clarrbigcross2dflozgi(loopee)*ellgentestarr(loopee)*this%bes2arr(yeye,loopee)
                    xicross2dfhizgi(yeye) = xicross2dfhizgi(yeye) + 1.0d0/pisj/2.0d0*clarrbigcross2dfhizgi(loopee)*ellgentestarr(loopee)*this%bes2arr(yeye,loopee)
                end if
            end do
        end do

        !old Catherine's scheme
        xiplusminus(1+(biin-1)*nangbins*2:nangbins+(biin-1)*nangbins*2) = xiplus + xiplusii + xiplusgi !full xi+ accounting for intrinsic alignments
        xiplusminus((nangbins+1)+(biin-1)*nangbins*2:nangbins*2+(biin-1)*nangbins*2) = ximinus + ximinusii + ximinusgi !full xi- accounting for intrinsic alignments
        xiplusminusgg(1+(biin-1)*nangbins*2:nangbins+(biin-1)*nangbins*2) = xiplus
        xiplusminusgg((nangbins+1)+(biin-1)*nangbins*2:nangbins*2+(biin-1)*nangbins*2) = ximinus
        xiplusminusii(1+(biin-1)*nangbins*2:nangbins+(biin-1)*nangbins*2) = xiplusii
        xiplusminusii((nangbins+1)+(biin-1)*nangbins*2:nangbins*2+(biin-1)*nangbins*2) = ximinusii
        xiplusminusgi(1+(biin-1)*nangbins*2:nangbins+(biin-1)*nangbins*2) = xiplusgi
        xiplusminusgi((nangbins+1)+(biin-1)*nangbins*2:nangbins*2+(biin-1)*nangbins*2) = ximinusgi
        if((m1arr(biin) == m2arr(biin)) .and. (obj%setscenarioomp > 1)) then
            bbintmeraigen = bbintmeraigen+1
             xicrossarrlowz(1+(bbintmeraigen-1)*nangbins:nangbins+(bbintmeraigen-1)*nangbins) = xicrosslowz + xicrosslowzgi
             xicrossarrcmass(1+(bbintmeraigen-1)*nangbins:nangbins+(bbintmeraigen-1)*nangbins) = xicrosscmass + xicrosscmassgi
             xicrossarr2dfloz(1+(bbintmeraigen-1)*nangbins:nangbins+(bbintmeraigen-1)*nangbins) = xicross2dfloz + xicross2dflozgi
             xicrossarr2dfhiz(1+(bbintmeraigen-1)*nangbins:nangbins+(bbintmeraigen-1)*nangbins) = xicross2dfhiz + xicross2dfhizgi
        end if

    end do

    this%klinesum = this%klinesum + obj%kline


    finaltheory(1:size(xiplusminus)) = xiplusminus
    if(obj%setscenarioomp > 1) then !ordering is 2dfloz, 2dfhiz, cmass, then lowz
        finaltheory(1+size(xiplusminus):size(xiplusminus)+size(xicrossarr2dfloz)) = xicrossarr2dfloz
        finaltheory(1+size(xiplusminus)+size(xicrossarr2dfloz):size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)) = xicrossarr2dfhiz
        finaltheory(1+size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz):size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass)) = xicrossarrcmass
        finaltheory(1+size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass):size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass)+size(xicrossarrlowz)) = xicrossarrlowz
    end if
    if(obj%setscenarioomp == 3) then
        finaltheory(1+size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass)+size(xicrossarrlowz):size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass)+size(xicrossarrlowz)+size(P0P2_final_2dfloz_overlap)) = P0P2_final_2dfloz_overlap
        finaltheory(1+size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass)+size(xicrossarrlowz)+size(P0P2_final_2dfloz_overlap):size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass)+size(xicrossarrlowz)+size(P0P2_final_2dfloz_overlap)+size(P0P2_final_2dfhiz_overlap)) = P0P2_final_2dfhiz_overlap
        finaltheory(1+size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass)+size(xicrossarrlowz)+size(P0P2_final_2dfloz_overlap)+size(P0P2_final_2dfhiz_overlap):size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass)+size(xicrossarrlowz)+size(P0P2_final_2dfloz_overlap)+size(P0P2_final_2dfhiz_overlap)+size(P0P2_final_cmass_overlap)) = P0P2_final_cmass_overlap
        finaltheory(1+size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass)+size(xicrossarrlowz)+size(P0P2_final_2dfloz_overlap)+size(P0P2_final_2dfhiz_overlap)+size(P0P2_final_cmass_overlap):size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass)+size(xicrossarrlowz)+size(P0P2_final_2dfloz_overlap)+size(P0P2_final_2dfhiz_overlap)+size(P0P2_final_cmass_overlap)+size(P0P2_final_lowz_overlap)) = P0P2_final_lowz_overlap
    end if

    howin = 1
    do how=1,this%sizcovpremask
        if(this%maskelements(how) == 1) then
            dataminustheory(howin) = this%xipm(howin) - finaltheory(how)
            dataminustheoryhoho(howin) = finaltheory(how)
            dataminustheoryhihi(howin) = this%xipm(howin)
            howin = howin + 1
        end if
    end do

    likesj = (Matrix_QuadForm(this%invcovxipm,dataminustheory))/2.0d0 !compute likelihood
    chisqsj = likesj*2.0
    chisqsjsellentin = 930.0d0*dlog(1.0d0 + likesj*2.0d0/(930.0d0-1.0d0)) !KiDS
    CosmoLSS_LnLike = 930.0d0/2.0d0*dlog(1.0d0 + likesj*2.0d0/(930.0d0-1.0d0)) !!KiDS: Sellentin-Heavens
!!!    CosmoLSS_LnLike = likesj !analytic or Hartlap

    finish  = OMP_get_wtime()
    !print *, 'joint chisq, 2*CosmoLSS_LnLike, finish-start, finish-startrsd', chisqsj, 2.0d0*CosmoLSS_LnLike, finish-start, finish-startrsd



    if(this%use_cmass_overlap) then
        DEALLOCATE(P0P2_final_cmass_overlap)
        DEALLOCATE(diff_cmass_overlap)
    end if
    if(this%use_lowz_overlap) then
        DEALLOCATE(P0P2_final_lowz_overlap)
        DEALLOCATE(diff_lowz_overlap)
    end if
    if(this%use_2dfloz_overlap) then
        DEALLOCATE(P0P2_final_2dfloz_overlap)
        DEALLOCATE(diff_2dfloz_overlap)
    end if
    if(this%use_2dfhiz_overlap) then
        DEALLOCATE(P0P2_final_2dfhiz_overlap)
        DEALLOCATE(diff_2dfhiz_overlap)
    end if
!!!!!!!!!!!!!!!!!!!
!    deallocate(obj%theoryomp)
    deallocate(obj%ellarromp)
!    deallocate(obj%exact_z_index)
    deallocate(obj%psourcetypeomp)
    deallocate(obj%plenstypeomplowz)
    deallocate(obj%plenstypeompcmass)
    deallocate(obj%plenstypeomp2dfloz)
    deallocate(obj%plenstypeomp2dfhiz)
    deallocate(obj%momp)
    deallocate(obj%aphotarromp)
    deallocate(obj%weightarromp)
    deallocate(obj%weightarromplowz)
    deallocate(obj%weightarrompcmass)
    deallocate(obj%weightarromplenslowz)
    deallocate(obj%weightarromplenscmass)
    deallocate(obj%weightarromp2dfloz)
    deallocate(obj%weightarromp2dfhiz)
    deallocate(obj%weightarromplens2dfloz)
    deallocate(obj%weightarromplens2dfhiz)
    deallocate(obj%weightpsarromp)
    deallocate(obj%weightnoarromp)
    deallocate(obj%weightnoarromplenslowz)
    deallocate(obj%weightnoarromplenscmass)
    deallocate(obj%weightnoarromplens2dfloz)
    deallocate(obj%weightnoarromplens2dfhiz)
    deallocate(obj%m1arromp)
    deallocate(obj%m2arromp)
    deallocate(obj%clinterptypeomp)
    deallocate(obj%clinterptypeiiomp)
    deallocate(obj%clinterptypegiomp)
    deallocate(obj%clinterptypecrosslowzomp)
    deallocate(obj%clinterptypecrosscmassomp)
    deallocate(obj%clinterptypecross2dflozomp)
    deallocate(obj%clinterptypecross2dfhizomp)
    deallocate(obj%exact_zmodomp)
    deallocate(obj%exact_zmodomplenslowz)
    deallocate(obj%exact_zmodomplenscmass)
    deallocate(obj%gcubicomp)
    deallocate(obj%lumarromp)

    deallocate(clarrbig)
    deallocate(clarrbiggi)
    deallocate(clarrbigii)
    deallocate(clarrbigcrosslowz)
    deallocate(clarrbigcrosscmass)
    deallocate(clarrbigcross2dfloz)
    deallocate(clarrbigcross2dfhiz)
    deallocate(ellgentestarr)
    deallocate(xiplusminus)
    deallocate(xiplusminusgg)
    deallocate(xiplusminusii)
    deallocate(xiplusminusgi)
    deallocate(dataminustheory)
    deallocate(dataminustheoryhihi)
    deallocate(dataminustheoryhoho)
    deallocate(finaltheory)
    deallocate(xicrossarrlowz)
    deallocate(xicrossarrlowzgi)
    deallocate(xicrossarrcmass)
    deallocate(xicrossarrcmassgi)
    deallocate(xicrossarr2dfloz)
    deallocate(xicrossarr2dflozgi)
    deallocate(xicrossarr2dfhiz)
    deallocate(xicrossarr2dfhizgi)
    deallocate(clarr)
    deallocate(clarrii)
    deallocate(clarrgi)
    deallocate(clarrcrosslowz)
    deallocate(clarrcrosscmass)
    deallocate(clarrcross2dfloz)
    deallocate(clarrcross2dfhiz)
    deallocate(ellarr)
    deallocate(xiplus)
    deallocate(xiplusii)
    deallocate(xiplusgi)
    deallocate(xicrosslowz)
    deallocate(xicrosslowzgi)
    deallocate(xicrosscmass)
    deallocate(xicrosscmassgi)
    deallocate(xicross2dfloz)
    deallocate(xicross2dflozgi)
    deallocate(xicross2dfhiz)
    deallocate(xicross2dfhizgi)
    deallocate(ximinus)
    deallocate(ximinusii)
    deallocate(ximinusgi)
    deallocate(mumax1arr)
    deallocate(mumax2arr)
    deallocate(aphotarr)
    deallocate(m1arr)
    deallocate(m2arr)
    deallocate(mnomp)
    deallocate(cfhtlum)
    deallocate(lumarr)
    deallocate(weightarr)
    deallocate(weightarrlowz)
    deallocate(weightarrcmass)
    deallocate(weightarrlenslowz)
    deallocate(weightarrlenscmass)
    deallocate(weightarr2dfloz)
    deallocate(weightarr2dfhiz)
    deallocate(weightarrlens2dfloz)
    deallocate(weightarrlens2dfhiz)
    deallocate(weightpsarr)
    deallocate(weightnoarr)
    deallocate(weightnoarrlenslowz)
    deallocate(weightnoarrlenscmass)
    deallocate(weightnoarrlens2dfloz)
    deallocate(weightnoarrlens2dfhiz)
    deallocate(binwarr)
    deallocate(trapzarr)
    deallocate(trapzarrii)
    deallocate(trapzarrgi)
    deallocate(trapzarrcrosslowz)
    deallocate(trapzarrcrosscmass)
    deallocate(trapzarrcross2dfloz)
    deallocate(trapzarrcross2dfhiz)
    deallocate(garr)

    deallocate(obj%weightpsarromplowz)
    deallocate(obj%weightpsarrompcmass)
    deallocate(obj%weightpsarromp2dfloz)
    deallocate(obj%weightpsarromp2dfhiz)
    deallocate(obj%clinterptypecrosslowzgiomp)
    deallocate(obj%clinterptypecrosscmassgiomp)
    deallocate(obj%clinterptypecross2dflozgiomp)
    deallocate(obj%clinterptypecross2dfhizgiomp)
    deallocate(clarrbigcrosslowzgi)
    deallocate(clarrbigcrosscmassgi)
    deallocate(clarrbigcross2dflozgi)
    deallocate(clarrbigcross2dfhizgi)
    deallocate(clarrcrosslowzgi)
    deallocate(clarrcrosscmassgi)
    deallocate(clarrcross2dflozgi)
    deallocate(clarrcross2dfhizgi)
    deallocate(weightpsarrlowz)
    deallocate(weightpsarrcmass)
    deallocate(weightpsarr2dfloz)
    deallocate(weightpsarr2dfhiz)
    deallocate(trapzarrcrosslowzgi)
    deallocate(trapzarrcrosscmassgi)
    deallocate(trapzarrcross2dflozgi)
    deallocate(trapzarrcross2dfhizgi)

  contains
!----------------------------------------------------------------------------------------------------------
! Theory calculation CMASS - overlap region
!----------------------------------------------------------------------------------------------------------

    subroutine Get_CMASS_overlap_TheoryMPK(P0P2_final_cmass_overlap)
    REAL(mcp), allocatable, dimension(:) :: P0P2_final_cmass_overlap
    REAL(mcp), allocatable, dimension(:)   :: P0_multipole_theory,P2_multipole_theory,P4_multipole_theory
    REAL(mcp), allocatable, dimension(:)   :: P0P2P4_Conv
    REAL(mcp), allocatable, dimension(:)   :: theory_vec,k_values_obs,k_values_conv
    REAL(mcp) :: mu_max,mu_min,mu_value,d_mu,nu_total,tol_error
    REAL(mcp) :: qnum,wnum
    REAL(mcp) :: scale_multipole_L0,scale_multipole_L2,scale_multipole_L4,AP_scaling
    REAL(mcp), parameter :: H_fid_no_h=95.49, DA_fid_no_h=1345.9,h_fiducial=0.7    !! Fiducial parameters at z= 0.57 units [Mpc]
    REAL(mcp) :: H_fid,DA_fid,k_val,k_step
    INTEGER   :: inum,k_index,num_k
    REAL(mcp) rombint
    external rombint

    allocate(P0_multipole_theory(this%cmass_mp_overlap%k_num_conv))
    allocate(P2_multipole_theory(this%cmass_mp_overlap%k_num_conv))
    allocate(P4_multipole_theory(this%cmass_mp_overlap%k_num_conv))

    allocate(theory_vec(this%cmass_mp_overlap%size_convolution))
    allocate(P0P2P4_Conv(this%cmass_mp_overlap%size_convolution))
    allocate(k_values_obs(this%cmass_mp_overlap%k_num_obs))
    allocate(k_values_conv(this%cmass_mp_overlap%k_num_conv))

        !----------------------------------------------------------------------------------------------------------
    ! Compute k vectors (observed values and values for convolution matrix)
    !----------------------------------------------------------------------------------------------------------

    do k_index = 1,this%cmass_mp_overlap%k_num_obs
        k_val = this%cmass_mp_overlap%k_min_obs + this%cmass_mp_overlap%k_spacing_obs*(k_index - 1.0)
        k_values_obs(k_index) = k_val
    end do

    do k_index = 1,this%cmass_mp_overlap%k_num_conv
        k_val = this%cmass_mp_overlap%k_min_theory + this%cmass_mp_overlap%dk_theory*(k_index - 1.0)
        k_values_conv(k_index) = k_val
    end do

    !----------------------------------------------------------------------------------------------------------
    ! Setup AP parameters -- can turn on or off
    !----------------------------------------------------------------------------------------------------------

    H_fid = H_fid_no_h/h_fiducial
    DA_fid = DA_fid_no_h*h_fiducial

    !write (*,*) 'CMB',CMB%H0, CMB%h, CMB%h0, this%cmass_mp_overlap%z_eff
    !write (*,*) Hofz(this%cmass_mp_overlap%z_eff), AngularDiameterDistance(this%cmass_mp_overlap%z_eff)
    !write (*,*) H_fid, H_fid_no_h, h_fiducial, HofzhUnit(this%cmass_mp_overlap%z_eff,CMB,c)
    !write (*,*) DA_fid, DA_fid_no_h, D_AhUnit(this%cmass_mp_overlap%z_eff,CMB)
    
    a_perp_cmass_overlap  = D_AhUnit(this%cmass_mp_overlap%z_eff,CMB)/DA_fid
    a_par_cmass_overlap   = H_fid/HofzhUnit(this%cmass_mp_overlap%z_eff,CMB,c)


    !----------------------------------------------------------------------------------------------------------
    ! start theory P0,P2,P4 calculation
    !----------------------------------------------------------------------------------------------------------

    tol_error = 0.000001
    mu_max = 1.0
    mu_min = -1.0

    scale_multipole_L0 = 1.0/2.0
    scale_multipole_L2 = 5.0/2.0
    scale_multipole_L4 = 9.0/2.0

    AP_scaling = 1.0/(a_par_cmass_overlap*(a_perp_cmass_overlap**2))

    do k_index = 1,this%cmass_mp_overlap%k_num_conv


        k_val_overlap_cmass = this%cmass_mp_overlap%k_min_theory + this%cmass_mp_overlap%dk_theory*(k_index - 1.0)
        !write(*,*) k_index, 'P0', k_val_overlap_cmass, this%cmass_mp_overlap%k_min_theory, this%cmass_mp_overlap%dk_theory

        P0_multipole_theory(k_index) = AP_scaling*scale_multipole_L0*rombint(P0_multipole_FN_cmass_overlap,mu_min,mu_max,tol_error)
       !write(*,*) k_index, 'P2'
        P2_multipole_theory(k_index) = AP_scaling*scale_multipole_L2*rombint(P2_multipole_FN_cmass_overlap,mu_min,mu_max,tol_error)
       !write(*,*) k_index, 'P4'
        P4_multipole_theory(k_index) = AP_scaling*scale_multipole_L4*rombint(P4_multipole_FN_cmass_overlap,mu_min,mu_max,tol_error)
        
    end do
      
!    do k_index = 1,this%cmass_mp_overlap%k_num_conv
!                k_val_overlap_cmass = this%cmass_mp_overlap%k_min_theory + this%cmass_mp_overlap%dk_theory*(k_index - 1.0)
!    end do

!    do k_index = 1,this%cmass_mp_overlap%k_num_conv
!                k_val_overlap_cmass = this%cmass_mp_overlap%k_min_theory + this%cmass_mp_overlap%dk_theory*(k_index - 1.0)
!    end do

!    do k_index = 1,this%cmass_mp_overlap%k_num_conv
!                k_val_overlap_cmass = this%cmass_mp_overlap%k_min_theory + this%cmass_mp_overlap%dk_theory*(k_index - 1.0)
!    end do

    !------------------------------------------------------------------------------------------------------------------------
    ! Stack into single vector and Compute convolved multipoles -- done using convolution matrix,
    ! note,     P0P2P4_Conv(:) is the convolved vector!
    !-------------------------------------------------------------------------------------------------------------------------


    theory_vec(1:10)   = P0_multipole_theory
    theory_vec(11:20)  = P2_multipole_theory
    theory_vec(21:30)  = P4_multipole_theory


!    CALL Matrix_MulVec(this%cmass_mp_overlap%Convolution_matrix(:,:),theory_vec(:),P0P2P4_Conv(:))
    !! TT check:
    P0P2P4_Conv = matmul(this%cmass_mp_overlap%Convolution_matrix,theory_vec)

    !write (*,*) 'MatMul worked'
    !----------------------------------------------------------------------------------------------------------
    ! Stack observations into single vector P = [P0,P2,P4]
    !----------------------------------------------------------------------------------------------------------

    num_k = this%cmass_mp_overlap%k_num_obs

        !----------------------------------------------------------------------------------------------------------
    ! Setup interpolation for (convolved) Monopole and quadrupole theory predictions and evaluation at require values.
    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------

    !write (*,*) 'P0P2 splines...'

    call P0_theory_spline_cmass_overlap%Init(k_values_conv,P0P2P4_Conv(1:10),this%cmass_mp_overlap%k_num_conv)
    call P2_theory_spline_cmass_overlap%Init(k_values_conv,P0P2P4_Conv(11:20),this%cmass_mp_overlap%k_num_conv)

    !write (*,*) 'k obs:',this%cmass_mp_overlap%k_num_obs, k_values_obs
    DO inum=1,this%cmass_mp_overlap%k_num_obs
        k_val = k_values_obs(inum)
        k_step = this%cmass_mp_overlap%k_num_obs
        P0P2_final_cmass_overlap(inum) = P0_theory_spline_cmass_overlap%value(k_val)
        P0P2_final_cmass_overlap(inum + k_step) = P2_theory_spline_cmass_overlap%value(k_val)
    ENDDO

    DEALLOCATE(P0_multipole_theory)
    DEALLOCATE(P2_multipole_theory)
    DEALLOCATE(P4_multipole_theory)
    DEALLOCATE(theory_vec)
    DEALLOCATE(P0P2P4_Conv)
    DEALLOCATE(k_values_obs)
    DEALLOCATE(k_values_conv)

    end subroutine

!----------------------------------------------------------------------------------------------------------
! Theory calculation LOWZ - overlap region
!----------------------------------------------------------------------------------------------------------

    subroutine Get_LOWZ_overlap_TheoryMPK(P0P2_final_lowz_overlap)
    REAL(mcp), allocatable, dimension(:) :: P0P2_final_lowz_overlap
    REAL(mcp), allocatable, dimension(:)   :: P0_multipole_theory,P2_multipole_theory,P4_multipole_theory
    REAL(mcp), allocatable, dimension(:)   :: P0P2P4_Conv
    REAL(mcp), allocatable, dimension(:)   :: theory_vec,k_values_obs,k_values_conv
    REAL(mcp) :: mu_max,mu_min,mu_value,d_mu,nu_total,tol_error
    REAL(mcp) :: qnum,wnum
    REAL(mcp) :: scale_multipole_L0,scale_multipole_L2,scale_multipole_L4,AP_scaling
    REAL(mcp), parameter :: H_fid_no_h=82.53, DA_fid_no_h=960.2, h_fiducial=0.7  !! Fiducial parameters at z= 0.32 units [Mpc]
    REAL(mcp) :: H_fid,DA_fid,k_val,k_step
    INTEGER   :: inum,k_index,num_k
    REAL(mcp) rombint
    external rombint

    allocate(P0_multipole_theory(this%lowz_mp_overlap%k_num_conv))
    allocate(P2_multipole_theory(this%lowz_mp_overlap%k_num_conv))
    allocate(P4_multipole_theory(this%lowz_mp_overlap%k_num_conv))

    allocate(theory_vec(this%lowz_mp_overlap%size_convolution))
    allocate(P0P2P4_Conv(this%lowz_mp_overlap%size_convolution))

    allocate(k_values_obs(this%lowz_mp_overlap%k_num_obs))
    allocate(k_values_conv(this%lowz_mp_overlap%k_num_conv))

    !----------------------------------------------------------------------------------------------------------
    ! Compute k vectors (observed values and values for convolution matrix)
    !----------------------------------------------------------------------------------------------------------

    do k_index = 1,this%lowz_mp_overlap%k_num_obs
        k_val = this%lowz_mp_overlap%k_min_obs + this%lowz_mp_overlap%k_spacing_obs*(k_index - 1.0)
        k_values_obs(k_index) = k_val
    end do

    do k_index = 1,this%lowz_mp_overlap%k_num_conv
        k_val = this%lowz_mp_overlap%k_min_theory + this%lowz_mp_overlap%dk_theory*(k_index - 1.0)
        k_values_conv(k_index) = k_val
    end do

    !----------------------------------------------------------------------------------------------------------
    ! Setup AP parameters -- can turn on or off
    !----------------------------------------------------------------------------------------------------------

    H_fid = H_fid_no_h/h_fiducial
    DA_fid = DA_fid_no_h*h_fiducial

    a_perp_lowz_overlap  = D_AhUnit(this%lowz_mp_overlap%z_eff,CMB)/DA_fid
    a_par_lowz_overlap   = H_fid/HofzhUnit(this%lowz_mp_overlap%z_eff,CMB,c)


    !----------------------------------------------------------------------------------------------------------
    ! start theory P0,P2,P4 calculation
    !----------------------------------------------------------------------------------------------------------

    tol_error = 0.000001
    mu_max = 1.0
    mu_min = -1.0

    scale_multipole_L0 = 1.0/2.0
    scale_multipole_L2 = 5.0/2.0
    scale_multipole_L4 = 9.0/2.0
    scale_multipole_L4 = 9.0/2.0

    AP_scaling = 1.0/(a_par_lowz_overlap*(a_perp_lowz_overlap**2))

    do k_index = 1,this%lowz_mp_overlap%k_num_conv

        k_val_overlap_lowz = this%lowz_mp_overlap%k_min_theory + this%lowz_mp_overlap%dk_theory*(k_index - 1.0)

        P0_multipole_theory(k_index) = AP_scaling*scale_multipole_L0*rombint(P0_multipole_FN_lowz_overlap,mu_min,mu_max,tol_error)
        P2_multipole_theory(k_index) = AP_scaling*scale_multipole_L2*rombint(P2_multipole_FN_lowz_overlap,mu_min,mu_max,tol_error)
        P4_multipole_theory(k_index) = AP_scaling*scale_multipole_L4*rombint(P4_multipole_FN_lowz_overlap,mu_min,mu_max,tol_error)

    end do

!    do k_index = 1,this%lowz_mp_overlap%k_num_conv
!            k_val_overlap_lowz = this%lowz_mp_overlap%k_min_theory + this%lowz_mp_overlap%dk_theory*(k_index - 1.0)
!    end do

!    do k_index = 1,this%lowz_mp_overlap%k_num_conv
!            k_val_overlap_lowz = this%lowz_mp_overlap%k_min_theory + this%lowz_mp_overlap%dk_theory*(k_index - 1.0)
!    end do

!    do k_index = 1,this%lowz_mp_overlap%k_num_conv
!            k_val_overlap_lowz = this%lowz_mp_overlap%k_min_theory + this%lowz_mp_overlap%dk_theory*(k_index - 1.0)
!    end do

    !------------------------------------------------------------------------------------------------------------------------
    ! Stack into single vector and Compute convolved multipoles -- done using convolution matrix,
    ! note,     P0P2P4_Conv(:) is the convolved vector!
    !-------------------------------------------------------------------------------------------------------------------------

    theory_vec(1:10)   = P0_multipole_theory
    theory_vec(11:20)  = P2_multipole_theory
    theory_vec(21:30)  = P4_multipole_theory

!    CALL Matrix_MulVec(this%lowz_mp_overlap%Convolution_matrix(:,:),theory_vec(:),P0P2P4_Conv(:))
    !!TT Check:
    P0P2P4_Conv = matmul(this%lowz_mp_overlap%Convolution_matrix,theory_vec)

    !----------------------------------------------------------------------------------------------------------
    ! Stack observations into single vector P = [P0,P2,P4]
    !----------------------------------------------------------------------------------------------------------

    num_k = this%lowz_mp_overlap%k_num_obs

    !----------------------------------------------------------------------------------------------------------
    ! Setup interpolation for (convolved) Monopole and quadrupole theory predictions and evaluation at require values.
    !----------------------------------------------------------------------------------------------------------

    call P0_theory_spline_lowz_overlap%Init(k_values_conv,P0P2P4_Conv(1:10),this%lowz_mp_overlap%k_num_conv)
    call P2_theory_spline_lowz_overlap%Init(k_values_conv,P0P2P4_Conv(11:20),this%lowz_mp_overlap%k_num_conv)

    DO inum=1,this%lowz_mp_overlap%k_num_obs
        k_val = k_values_obs(inum)
        k_step = this%lowz_mp_overlap%k_num_obs
        P0P2_final_lowz_overlap(inum) = P0_theory_spline_lowz_overlap%value(k_val)
        P0P2_final_lowz_overlap(inum + k_step) = P2_theory_spline_lowz_overlap%value(k_val)
    ENDDO


    DEALLOCATE(P0_multipole_theory)
    DEALLOCATE(P2_multipole_theory)
    DEALLOCATE(P4_multipole_theory)
    DEALLOCATE(theory_vec)
    DEALLOCATE(P0P2P4_Conv)
    DEALLOCATE(k_values_obs)
    DEALLOCATE(k_values_conv)

    end subroutine



!----------------------------------------------------------------------------------------------------------
! Theory calculation 2dfloz - overlap region
!----------------------------------------------------------------------------------------------------------

    subroutine Get_2dfloz_overlap_TheoryMPK(P0P2_final_2dfloz_overlap)
    REAL(mcp), allocatable, dimension(:) :: P0P2_final_2dfloz_overlap
    REAL(mcp), allocatable, dimension(:)   :: P0_multipole_theory,P2_multipole_theory,P4_multipole_theory
    REAL(mcp), allocatable, dimension(:)   :: P0P2P4_Conv
    REAL(mcp), allocatable, dimension(:)   :: theory_vec,k_values_obs,k_values_conv
    REAL(mcp) :: mu_max,mu_min,mu_value,d_mu,nu_total,tol_error
    REAL(mcp) :: qnum,wnum
    REAL(mcp) :: scale_multipole_L0,scale_multipole_L2,scale_multipole_L4,AP_scaling
    REAL(mcp), parameter :: H_fid_no_h=82.07, DA_fid_no_h=939.7, h_fiducial=0.7  !! Fiducial parameters at z= 0.31 units [Mpc]
    REAL(mcp) :: H_fid,DA_fid,k_val,k_step
    INTEGER   :: inum,k_index,num_k
    REAL(mcp) rombint
    external rombint

    allocate(P0_multipole_theory(this%twodfloz_mp_overlap%k_num_conv))
    allocate(P2_multipole_theory(this%twodfloz_mp_overlap%k_num_conv))
    allocate(P4_multipole_theory(this%twodfloz_mp_overlap%k_num_conv))

    allocate(theory_vec(this%twodfloz_mp_overlap%size_convolution))
    allocate(P0P2P4_Conv(this%twodfloz_mp_overlap%size_convolution))

    allocate(k_values_obs(this%twodfloz_mp_overlap%k_num_obs))
    allocate(k_values_conv(this%twodfloz_mp_overlap%k_num_conv))

    !----------------------------------------------------------------------------------------------------------
    ! Compute k vectors (observed values and values for convolution matrix)
    !----------------------------------------------------------------------------------------------------------

    do k_index = 1,this%twodfloz_mp_overlap%k_num_obs
        k_val = this%twodfloz_mp_overlap%k_min_obs + this%twodfloz_mp_overlap%k_spacing_obs*(k_index - 1.0)
        k_values_obs(k_index) = k_val
    end do

    do k_index = 1,this%twodfloz_mp_overlap%k_num_conv
        k_val = this%twodfloz_mp_overlap%k_min_theory + this%twodfloz_mp_overlap%dk_theory*(k_index - 1.0)
        k_values_conv(k_index) = k_val
    end do

    !----------------------------------------------------------------------------------------------------------
    ! Setup AP parameters -- can turn on or off
    !----------------------------------------------------------------------------------------------------------

    H_fid = H_fid_no_h/h_fiducial
    DA_fid = DA_fid_no_h*h_fiducial

    a_perp_2dfloz_overlap  = D_AhUnit(this%twodfloz_mp_overlap%z_eff,CMB)/DA_fid
    a_par_2dfloz_overlap   = H_fid/HofzhUnit(this%twodfloz_mp_overlap%z_eff,CMB,c)


    !----------------------------------------------------------------------------------------------------------
    ! start theory P0,P2,P4 calculation
    !----------------------------------------------------------------------------------------------------------

    tol_error = 0.000001
    mu_max = 1.0
    mu_min = -1.0

    scale_multipole_L0 = 1.0/2.0
    scale_multipole_L2 = 5.0/2.0
    scale_multipole_L4 = 9.0/2.0
    scale_multipole_L4 = 9.0/2.0

    AP_scaling = 1.0/(a_par_2dfloz_overlap*(a_perp_2dfloz_overlap**2))

    do k_index = 1,this%twodfloz_mp_overlap%k_num_conv

        k_val_overlap_2dfloz = this%twodfloz_mp_overlap%k_min_theory + this%twodfloz_mp_overlap%dk_theory*(k_index - 1.0)

        P0_multipole_theory(k_index) = AP_scaling*scale_multipole_L0*rombint(P0_multipole_FN_2dfloz_overlap,mu_min,mu_max,tol_error)
        P2_multipole_theory(k_index) = AP_scaling*scale_multipole_L2*rombint(P2_multipole_FN_2dfloz_overlap,mu_min,mu_max,tol_error)
        P4_multipole_theory(k_index) = AP_scaling*scale_multipole_L4*rombint(P4_multipole_FN_2dfloz_overlap,mu_min,mu_max,tol_error)

    end do

!    do k_index = 1,this%twodfloz_mp_overlap%k_num_conv
!            k_val_overlap_2dfloz = this%twodfloz_mp_overlap%k_min_theory + this%twodfloz_mp_overlap%dk_theory*(k_index - 1.0)
!    end do

!    do k_index = 1,this%twodfloz_mp_overlap%k_num_conv
!            k_val_overlap_2dfloz = this%twodfloz_mp_overlap%k_min_theory + this%twodfloz_mp_overlap%dk_theory*(k_index - 1.0)
!    end do

!    do k_index = 1,this%twodfloz_mp_overlap%k_num_conv
!            k_val_overlap_2dfloz = this%twodfloz_mp_overlap%k_min_theory + this%twodfloz_mp_overlap%dk_theory*(k_index - 1.0)
!    end do

    !------------------------------------------------------------------------------------------------------------------------
    ! Stack into single vector and Compute convolved multipoles -- done using convolution matrix,
    ! note,     P0P2P4_Conv(:) is the convolved vector!
    !-------------------------------------------------------------------------------------------------------------------------

    theory_vec(1:10)   = P0_multipole_theory
    theory_vec(11:20)  = P2_multipole_theory
    theory_vec(21:30)  = P4_multipole_theory

!!    CALL Matrix_MulVec(this%twodfloz_mp_overlap%Convolution_matrix(:,:),theory_vec(:),P0P2P4_Conv(:))
    !!TT Check:
    P0P2P4_Conv = matmul(this%twodfloz_mp_overlap%Convolution_matrix,theory_vec)
    !----------------------------------------------------------------------------------------------------------
    ! Stack observations into single vector P = [P0,P2,P4]
    !----------------------------------------------------------------------------------------------------------

    num_k = this%twodfloz_mp_overlap%k_num_obs

    !----------------------------------------------------------------------------------------------------------
    ! Setup interpolation for (convolved) Monopole and quadrupole theory predictions and evaluation at require values.
    !----------------------------------------------------------------------------------------------------------

    call P0_theory_spline_2dfloz_overlap%Init(k_values_conv,P0P2P4_Conv(1:10),this%twodfloz_mp_overlap%k_num_conv)
    call P2_theory_spline_2dfloz_overlap%Init(k_values_conv,P0P2P4_Conv(11:20),this%twodfloz_mp_overlap%k_num_conv)

    DO inum=1,this%twodfloz_mp_overlap%k_num_obs
        k_val = k_values_obs(inum)
        k_step = this%twodfloz_mp_overlap%k_num_obs
        P0P2_final_2dfloz_overlap(inum) = P0_theory_spline_2dfloz_overlap%value(k_val)
        P0P2_final_2dfloz_overlap(inum + k_step) = P2_theory_spline_2dfloz_overlap%value(k_val)
    ENDDO


    DEALLOCATE(P0_multipole_theory)
    DEALLOCATE(P2_multipole_theory)
    DEALLOCATE(P4_multipole_theory)
    DEALLOCATE(theory_vec)
    DEALLOCATE(P0P2P4_Conv)
    DEALLOCATE(k_values_obs)
    DEALLOCATE(k_values_conv)

    end subroutine



!----------------------------------------------------------------------------------------------------------
! Theory calculation 2dfhiz - overlap region
!----------------------------------------------------------------------------------------------------------

    subroutine Get_2dfhiz_overlap_TheoryMPK(P0P2_final_2dfhiz_overlap)
    REAL(mcp), allocatable, dimension(:) :: P0P2_final_2dfhiz_overlap
    REAL(mcp), allocatable, dimension(:)   :: P0_multipole_theory,P2_multipole_theory,P4_multipole_theory
    REAL(mcp), allocatable, dimension(:)   :: P0P2P4_Conv
    REAL(mcp), allocatable, dimension(:)   :: theory_vec,k_values_obs,k_values_conv
    REAL(mcp) :: mu_max,mu_min,mu_value,d_mu,nu_total,tol_error
    REAL(mcp) :: qnum,wnum
    REAL(mcp) :: scale_multipole_L0,scale_multipole_L2,scale_multipole_L4,AP_scaling
    REAL(mcp), parameter :: H_fid_no_h=94.92, DA_fid_no_h=1334.3,h_fiducial=0.7    !! Fiducial parameters at z= 0.56 units [Mpc]
    REAL(mcp) :: H_fid,DA_fid,k_val,k_step
    INTEGER   :: inum,k_index,num_k
    REAL(mcp) rombint
    external rombint

    allocate(P0_multipole_theory(this%twodfhiz_mp_overlap%k_num_conv))
    allocate(P2_multipole_theory(this%twodfhiz_mp_overlap%k_num_conv))
    allocate(P4_multipole_theory(this%twodfhiz_mp_overlap%k_num_conv))

    allocate(theory_vec(this%twodfhiz_mp_overlap%size_convolution))
    allocate(P0P2P4_Conv(this%twodfhiz_mp_overlap%size_convolution))

    allocate(k_values_obs(this%twodfhiz_mp_overlap%k_num_obs))
    allocate(k_values_conv(this%twodfhiz_mp_overlap%k_num_conv))

    !----------------------------------------------------------------------------------------------------------
    ! Compute k vectors (observed values and values for convolution matrix)
    !----------------------------------------------------------------------------------------------------------

    do k_index = 1,this%twodfhiz_mp_overlap%k_num_obs
        k_val = this%twodfhiz_mp_overlap%k_min_obs + this%twodfhiz_mp_overlap%k_spacing_obs*(k_index - 1.0)
        k_values_obs(k_index) = k_val
    end do

    do k_index = 1,this%twodfhiz_mp_overlap%k_num_conv
        k_val = this%twodfhiz_mp_overlap%k_min_theory + this%twodfhiz_mp_overlap%dk_theory*(k_index - 1.0)
        k_values_conv(k_index) = k_val
    end do

    !----------------------------------------------------------------------------------------------------------
    ! Setup AP parameters -- can turn on or off
    !----------------------------------------------------------------------------------------------------------

    H_fid = H_fid_no_h/h_fiducial
    DA_fid = DA_fid_no_h*h_fiducial

    a_perp_2dfhiz_overlap  = D_AhUnit(this%twodfhiz_mp_overlap%z_eff,CMB)/DA_fid
    a_par_2dfhiz_overlap   = H_fid/HofzhUnit(this%twodfhiz_mp_overlap%z_eff,CMB,c)

    !----------------------------------------------------------------------------------------------------------
    ! start theory P0,P2,P4 calculation
    !----------------------------------------------------------------------------------------------------------

    tol_error = 0.000001
    mu_max = 1.0
    mu_min = -1.0

    scale_multipole_L0 = 1.0/2.0
    scale_multipole_L2 = 5.0/2.0
    scale_multipole_L4 = 9.0/2.0
    scale_multipole_L4 = 9.0/2.0

    AP_scaling = 1.0/(a_par_2dfhiz_overlap*(a_perp_2dfhiz_overlap**2))

    do k_index = 1,this%twodfhiz_mp_overlap%k_num_conv

        k_val_overlap_2dfhiz = this%twodfhiz_mp_overlap%k_min_theory + this%twodfhiz_mp_overlap%dk_theory*(k_index - 1.0)

        P0_multipole_theory(k_index) = AP_scaling*scale_multipole_L0*rombint(P0_multipole_FN_2dfhiz_overlap,mu_min,mu_max,tol_error)
        P2_multipole_theory(k_index) = AP_scaling*scale_multipole_L2*rombint(P2_multipole_FN_2dfhiz_overlap,mu_min,mu_max,tol_error)
        P4_multipole_theory(k_index) = AP_scaling*scale_multipole_L4*rombint(P4_multipole_FN_2dfhiz_overlap,mu_min,mu_max,tol_error)

    end do

!    do k_index = 1,this%twodfhiz_mp_overlap%k_num_conv
!            k_val_overlap_2dfhiz = this%twodfhiz_mp_overlap%k_min_theory + this%twodfhiz_mp_overlap%dk_theory*(k_index - 1.0)
!    end do

!    do k_index = 1,this%twodfhiz_mp_overlap%k_num_conv
!            k_val_overlap_2dfhiz = this%twodfhiz_mp_overlap%k_min_theory + this%twodfhiz_mp_overlap%dk_theory*(k_index - 1.0)
!    end do

!    do k_index = 1,this%twodfhiz_mp_overlap%k_num_conv
!            k_val_overlap_2dfhiz = this%twodfhiz_mp_overlap%k_min_theory + this%twodfhiz_mp_overlap%dk_theory*(k_index - 1.0)
!    end do

    !------------------------------------------------------------------------------------------------------------------------
    ! Stack into single vector and Compute convolved multipoles -- done using convolution matrix,
    ! note,     P0P2P4_Conv(:) is the convolved vector!
    !-------------------------------------------------------------------------------------------------------------------------

    theory_vec(1:10)   = P0_multipole_theory
    theory_vec(11:20)  = P2_multipole_theory
    theory_vec(21:30)  = P4_multipole_theory

!!    CALL Matrix_MulVec(this%twodfhiz_mp_overlap%Convolution_matrix(:,:),theory_vec(:),P0P2P4_Conv(:))
    !!TT Check:
    P0P2P4_Conv = matmul(this%twodfhiz_mp_overlap%Convolution_matrix,theory_vec)
    !----------------------------------------------------------------------------------------------------------
    ! Stack observations into single vector P = [P0,P2,P4]
    !----------------------------------------------------------------------------------------------------------

    num_k = this%twodfhiz_mp_overlap%k_num_obs

    !----------------------------------------------------------------------------------------------------------
    ! Setup interpolation for (convolved) Monopole and quadrupole theory predictions and evaluation at require values.
    !----------------------------------------------------------------------------------------------------------

    call P0_theory_spline_2dfhiz_overlap%Init(k_values_conv,P0P2P4_Conv(1:10),this%twodfhiz_mp_overlap%k_num_conv)
    call P2_theory_spline_2dfhiz_overlap%Init(k_values_conv,P0P2P4_Conv(11:20),this%twodfhiz_mp_overlap%k_num_conv)

    DO inum=1,this%twodfhiz_mp_overlap%k_num_obs
        k_val = k_values_obs(inum)
        k_step = this%twodfhiz_mp_overlap%k_num_obs
        P0P2_final_2dfhiz_overlap(inum) = P0_theory_spline_2dfhiz_overlap%value(k_val)
        P0P2_final_2dfhiz_overlap(inum + k_step) = P2_theory_spline_2dfhiz_overlap%value(k_val)
    ENDDO

    DEALLOCATE(P0_multipole_theory)
    DEALLOCATE(P2_multipole_theory)
    DEALLOCATE(P4_multipole_theory)
    DEALLOCATE(theory_vec)
    DEALLOCATE(P0P2P4_Conv)
    DEALLOCATE(k_values_obs)
    DEALLOCATE(k_values_conv)

    end subroutine


        !----------------------------------------------------------------------------------------------------------
        ! functions for calculating theory multipoles for CMASS overlap region-- using 2D power spectrum (local functions)
        !----------------------------------------------------------------------------------------------------------

        function P0_multipole_FN_cmass_overlap(mu)
        REAL(mcp) :: mu,P0_multipole_FN_cmass_overlap

        P0_multipole_FN_cmass_overlap = &
        LinearPower_2D_s_cmass(k_val_overlap_cmass,mu,z_cmass,sigma_v_cmass,a_perp_cmass_overlap,a_par_cmass_overlap,b1_cmass,N_shot_cmass)*Legendre_Pn(mu,0)
        END function P0_multipole_FN_cmass_overlap

        function P2_multipole_FN_cmass_overlap(mu)
        REAL(mcp) :: mu,P2_multipole_FN_cmass_overlap

        P2_multipole_FN_cmass_overlap = &
        LinearPower_2D_s_cmass(k_val_overlap_cmass,mu,z_cmass,sigma_v_cmass,a_perp_cmass_overlap,a_par_cmass_overlap,b1_cmass,N_shot_cmass)*Legendre_Pn(mu,2)
        end function P2_multipole_FN_cmass_overlap

        function P4_multipole_FN_cmass_overlap(mu)
        REAL(mcp) :: mu,P4_multipole_FN_cmass_overlap

        P4_multipole_FN_cmass_overlap = &
        LinearPower_2D_s_cmass(k_val_overlap_cmass,mu,z_cmass,sigma_v_cmass,a_perp_cmass_overlap,a_par_cmass_overlap,b1_cmass,N_shot_cmass)*Legendre_Pn(mu,4)
        end function P4_multipole_FN_cmass_overlap

        !----------------------------------------------------------------------------------------------------------
        ! functions for calculating theory multipoles for LOWZ overlap region-- using 2D power spectrum (local functions)
        !----------------------------------------------------------------------------------------------------------

        function P0_multipole_FN_lowz_overlap(mu)
        REAL(mcp) :: mu,P0_multipole_FN_lowz_overlap

        P0_multipole_FN_lowz_overlap = &
        LinearPower_2D_s_cmass(k_val_overlap_lowz,mu,z_lowz,sigv_lowz,a_perp_lowz_overlap,a_par_lowz_overlap,b1_lowz,N_shot_lowz)*Legendre_Pn(mu,0)
        END function P0_multipole_FN_lowz_overlap

        function P2_multipole_FN_lowz_overlap(mu)
        REAL(mcp) :: mu,P2_multipole_FN_lowz_overlap

        P2_multipole_FN_lowz_overlap = &
        LinearPower_2D_s_cmass(k_val_overlap_lowz,mu,z_lowz,sigv_lowz,a_perp_lowz_overlap,a_par_lowz_overlap,b1_lowz,N_shot_lowz)*Legendre_Pn(mu,2)
        end function P2_multipole_FN_lowz_overlap

        function P4_multipole_FN_lowz_overlap(mu)
        REAL(mcp) :: mu,P4_multipole_FN_lowz_overlap

        P4_multipole_FN_lowz_overlap = &
        LinearPower_2D_s_cmass(k_val_overlap_lowz,mu,z_lowz,sigv_lowz,a_perp_lowz_overlap,a_par_lowz_overlap,b1_lowz,N_shot_lowz)*Legendre_Pn(mu,4)
        end function P4_multipole_FN_lowz_overlap


        !----------------------------------------------------------------------------------------------------------
        ! functions for calculating theory multipoles for 2dfloz overlap region-- using 2D power spectrum (local functions)
        !----------------------------------------------------------------------------------------------------------

        function P0_multipole_FN_2dfloz_overlap(mu)
        REAL(mcp) :: mu,P0_multipole_FN_2dfloz_overlap

        P0_multipole_FN_2dfloz_overlap = &
        LinearPower_2D_s_cmass(k_val_overlap_2dfloz,mu,z_2dfloz,sigv_2dfloz,a_perp_2dfloz_overlap,a_par_2dfloz_overlap,b1_2dfloz,N_shot_2dfloz)*Legendre_Pn(mu,0)
        END function P0_multipole_FN_2dfloz_overlap

        function P2_multipole_FN_2dfloz_overlap(mu)
        REAL(mcp) :: mu,P2_multipole_FN_2dfloz_overlap

        P2_multipole_FN_2dfloz_overlap = &
        LinearPower_2D_s_cmass(k_val_overlap_2dfloz,mu,z_2dfloz,sigv_2dfloz,a_perp_2dfloz_overlap,a_par_2dfloz_overlap,b1_2dfloz,N_shot_2dfloz)*Legendre_Pn(mu,2)
        end function P2_multipole_FN_2dfloz_overlap

        function P4_multipole_FN_2dfloz_overlap(mu)
        REAL(mcp) :: mu,P4_multipole_FN_2dfloz_overlap

        P4_multipole_FN_2dfloz_overlap = &
        LinearPower_2D_s_cmass(k_val_overlap_2dfloz,mu,z_2dfloz,sigv_2dfloz,a_perp_2dfloz_overlap,a_par_2dfloz_overlap,b1_2dfloz,N_shot_2dfloz)*Legendre_Pn(mu,4)
        end function P4_multipole_FN_2dfloz_overlap


        !----------------------------------------------------------------------------------------------------------
        ! functions for calculating theory multipoles for 2dfhiz overlap region-- using 2D power spectrum (local functions)
        !----------------------------------------------------------------------------------------------------------

        function P0_multipole_FN_2dfhiz_overlap(mu)
        REAL(mcp) :: mu,P0_multipole_FN_2dfhiz_overlap

        P0_multipole_FN_2dfhiz_overlap = &
        LinearPower_2D_s_cmass(k_val_overlap_2dfhiz,mu,z_2dfhiz,sigv_2dfhiz,a_perp_2dfhiz_overlap,a_par_2dfhiz_overlap,b1_2dfhiz,N_shot_2dfhiz)*Legendre_Pn(mu,0)
        END function P0_multipole_FN_2dfhiz_overlap

        function P2_multipole_FN_2dfhiz_overlap(mu)
        REAL(mcp) :: mu,P2_multipole_FN_2dfhiz_overlap

        P2_multipole_FN_2dfhiz_overlap = &
        LinearPower_2D_s_cmass(k_val_overlap_2dfhiz,mu,z_2dfhiz,sigv_2dfhiz,a_perp_2dfhiz_overlap,a_par_2dfhiz_overlap,b1_2dfhiz,N_shot_2dfhiz)*Legendre_Pn(mu,2)
        end function P2_multipole_FN_2dfhiz_overlap

        function P4_multipole_FN_2dfhiz_overlap(mu)
        REAL(mcp) :: mu,P4_multipole_FN_2dfhiz_overlap

        P4_multipole_FN_2dfhiz_overlap = &
        LinearPower_2D_s_cmass(k_val_overlap_2dfhiz,mu,z_2dfhiz,sigv_2dfhiz,a_perp_2dfhiz_overlap,a_par_2dfhiz_overlap,b1_2dfhiz,N_shot_2dfhiz)*Legendre_Pn(mu,4)
        end function P4_multipole_FN_2dfhiz_overlap


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function sjgrowtha(reda) !integrate to obtain growth function D(0)/D(z)
        REAL(mcp) :: reda,sjgrowtha

        !sjgrowtha = (1.0d0/reda)*Theory%growth_z%Value(1.0d0/reda-1.0d0)/Theory%sigma8_z%Value(1.0d0/reda-1.0d0)
        !TT Check
        sjgrowtha = (1.0d0/reda)*GrowthRate%Value(reda)
    END function sjgrowtha

    !non-weight component for trapezoid integration
    function sjclsobjnoweight(obj,zlens)
        type(my_type) :: obj
        REAL(mcp) :: zlens,sjclsobjnoweight,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez,ckms
        real(mcp) :: kminsj,kmaxsj

        ckms = 299792.458d0
        distz = ComovingRadialDistance(zlens)
        obj%distlensflatomp = distz
        distz = f_K(distz) !with curvature        
        kval = (obj%ellarromp(obj%ellgenomp)+1.0d0/2.0d0)/((obj%h0omp/100.0d0)*distz)

        kminsj = 1.0d-5
        kmaxsj = 100.0d0
        deltafinal = 0.0d0
        if(kval < kminsj) obj%kline = obj%kline + 1
        if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = NL_MPKPowerAt(kval,zlens)
	hubblez = Hofz(zlens)
        weightpart = 3.0d0/2.0d0*(obj%omdmomp+obj%ombomp)*(obj%h0omp/ckms)**2.0d0*distz*(1.0d0+zlens)
        sjclsobjnoweight = 1.0d0/(distz**2.0d0)/hubblez*deltafinal/(obj%homp)**3.0d0

    END function sjclsobjnoweight

    !weight for trapezoid integration
    function sjclsobjonlyweight(obj,zlens)
        type(my_type) :: obj
        REAL(mcp) :: zlens,sjclsobjonlyweight,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez,ckms
        real(mcp) :: kminsj,kmaxsj

        ckms = 299792.458d0
        distz = ComovingRadialDistance(zlens)
        obj%distlensflatomp = distz !assuming flatness
        distz = f_K(distz) !with curvature        
        weightpart = 3.0d0/2.0d0*(obj%omdmomp+obj%ombomp)*(obj%h0omp/ckms)**2.0d0*distz*(1.0d0+zlens)
        obj%wchooseomp = 1 !doesn't matter here because imposed momp(1) = momp(2)
        sjclsobjonlyweight = weightpart*rombint_obj(obj,weightobjcubic,zlens,obj%mumaxwomp,obj%tol_erroromp) !note changed mumax to mumax1

    END function sjclsobjonlyweight

    !IA-weight for trapezoid integration
    function sjclsiiobjonlyweight(obj,zlens)
        type(my_type) :: obj
        REAL(mcp) :: zlens,sjclsiiobjonlyweight,kval,deltafinal,hubblez,growthsf,growthnorm
        real(mcp) :: kminsj,kmaxsj

        growthnorm = obj%gcubicomp%Value(zlens)	!this is D(0)/D(z)
	hubblez = Hofz(zlens)
        obj%wchooseomp = 1 !doesn't matter here
        sjclsiiobjonlyweight = psourceobjcubic(obj,zlens)*hubblez*(-obj%ampiayesomp*5.0d-14*2.77536627d11*(obj%omdmomp+obj%ombomp)*growthnorm*((1.0d0+zlens)/(1.0d0+0.3d0))**obj%redziayesomp)*(obj%lumarromp(obj%momp(1)))**obj%lumiayesomp !enforced momp(1) = momp(2), so either is fine

    END function sjclsiiobjonlyweight


    !WEAK LENSING INTEGRAND (GG) ---- wrt scale factor instead
    function sjclsobjsf(obj,alens)
        type(my_type) :: obj
        REAL(mcp) :: alens,zlens,sjclsobjsf,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez,ckms
        real(mcp) :: kminsj,kmaxsj

        zlens = 1.0d0/alens - 1.0d0
        ckms = 299792.458d0
        distz = ComovingRadialDistance(zlens)
        obj%distlensflatomp = distz !assuming flatness
        distz = f_K(distz) !with curvature        
        kval = (obj%ellarromp(obj%ellgenomp)+1.0d0/2.0d0)/((obj%h0omp/100.0d0)*distz)
        kminsj = 1.0d-5
        kmaxsj = 100.0d0
        deltafinal = 0.0d0
        if(kval < kminsj) obj%kline = obj%kline + 1
        if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = NL_MPKPowerAt(kval,zlens)
	hubblez = Hofz(zlens)
        weightpart = 3.0d0/2.0d0*(obj%omdmomp+obj%ombomp)*(obj%h0omp/ckms)**2.0d0*distz*(1.0d0+zlens)
        obj%wchooseomp = 1
        weight1 = weightpart*rombint_obj(obj,weightobjcubic,zlens,obj%mumax1omp,obj%tol_erroromp) !note changed mumax to mumax1
        if(obj%bnumomp == 1 .or. obj%bnumomp == 5 .or. obj%bnumomp == 8 .or. obj%bnumomp == 10) then !generalize later
            weight2 = weight1
        else
            obj%wchooseomp = 2
            weight2 = weightpart*rombint_obj(obj,weightobjcubic,zlens,obj%mumax2omp,obj%tol_erroromp)
        end if
        sjclsobjsf = ((1.0d0+zlens)**2.0d0)*weight1*weight2/(distz**2.0d0)/hubblez*deltafinal/(obj%homp)**3.0d0

    END function sjclsobjsf


    !INTRINSIC ALIGNMENT INTEGRAND (II) ---- wrt scale factor instead
    function sjclsiiobjsf(obj,alens)
        type(my_type) :: obj
        REAL(mcp) :: alens,zlens,sjclsiiobjsf,growthnorm,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez,ckms,growthsf
        real(mcp) :: kminsj,kmaxsj

        zlens = 1.0d0/alens - 1.0d0
        ckms = 299792.458d0
        distz = ComovingRadialDistance(zlens)
        obj%distlensflatomp = distz
        distz = f_K(distz) !with curvature        
        kval = (obj%ellarromp(obj%ellgenomp)+1.0d0/2.0d0)/((obj%h0omp/100.0d0)*distz)
        kminsj = 1.0d-5
        kmaxsj = 100.0d0
        deltafinal = 0.0d0
        if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = NL_MPKPowerAt(kval,zlens)
        growthnorm = obj%gcubicomp%Value(zlens)
        deltafinal = deltafinal*((obj%lumarromp(obj%momp(1))*obj%lumarromp(obj%momp(2)))**obj%lumiayesomp)*(-obj%ampiayesomp*5.0d-14*2.77536627d11*(obj%omdmomp+obj%ombomp)*growthnorm*((1.0d0+zlens)/(1.0d0+0.3d0))**obj%redziayesomp)**2.0d0
	hubblez = Hofz(zlens)
        obj%wchooseomp = 1
        weight1 = psourceobjcubic(obj,zlens)*hubblez
        if(obj%bnumomp == 1 .or. obj%bnumomp == 5 .or. obj%bnumomp == 8 .or. obj%bnumomp == 10) then !generalize later
            weight2 = weight1
        else
            obj%wchooseomp = 2
            weight2 = psourceobjcubic(obj,zlens)*hubblez
        end if
        sjclsiiobjsf = ((1.0d0+zlens)**2.0d0)*weight1*weight2/(distz**2.0d0)/hubblez*deltafinal/(obj%homp)**3.0d0

    END function sjclsiiobjsf


    !LENSING - INTRINSIC ALIGNMENT INTEGRAND (GI) --- wrt scale factor instead
    function sjclsgiobjsf(obj,alens)
        type(my_type) :: obj
        REAL(mcp) :: alens,zlens,sjclsgiobjsf,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez,ckms,weight3,weight4,growthsf,growthnorm
        real(mcp) :: kminsj,kmaxsj

        zlens = 1.0d0/alens - 1.0d0
        ckms = 299792.458d0
        distz = ComovingRadialDistance(zlens)
        obj%distlensflatomp = distz !assuming flatness
        distz = f_K(distz) !with curvature        
        kval = (obj%ellarromp(obj%ellgenomp)+1.0d0/2.0d0)/((obj%h0omp/100.0d0)*distz)
        kminsj = 1.0d-5
        kmaxsj = 100.0d0
        deltafinal = 0.0d0
        if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = NL_MPKPowerAt(kval,zlens)
        growthnorm = obj%gcubicomp%Value(zlens)	!this is D(0)/D(z)
        deltafinal = deltafinal*(-obj%ampiayesomp*5.0d-14*2.77536627d11*(obj%omdmomp+obj%ombomp)*growthnorm*((1.0d0+zlens)/(1.0d0+0.3d0))**obj%redziayesomp)
	hubblez = Hofz(zlens)
        weightpart = 3.0d0/2.0d0*(obj%omdmomp+obj%ombomp)*(obj%h0omp/ckms)**2.0d0*distz*(1.0d0+zlens)
        obj%wchooseomp = 1
        weight1 = weightpart*rombint_obj(obj,weightobjcubic,zlens,obj%mumax1omp,obj%tol_erroromp) !note changed mumax to mumax1
        weight3 = psourceobjcubic(obj,zlens)*hubblez
        if(obj%bnumomp == 1 .or. obj%bnumomp == 5 .or. obj%bnumomp == 8 .or. obj%bnumomp == 10) then !generalize later
            weight2 = weight1
            weight4 = weight3
        else
            obj%wchooseomp = 2
            weight2 = weightpart*rombint_obj(obj,weightobjcubic,zlens,obj%mumax2omp,obj%tol_erroromp) !note changed mumax to mumax2
            weight4 = psourceobjcubic(obj,zlens)*hubblez
        end if
        sjclsgiobjsf = ((1.0d0+zlens)**2.0d0)*(weight1*weight4*(obj%lumarromp(obj%momp(2)))**obj%lumiayesomp + weight2*weight3*(obj%lumarromp(obj%momp(1)))**obj%lumiayesomp)/(distz**2.0d0)/hubblez*deltafinal/(obj%homp)**3.0d0

    END function sjclsgiobjsf


    !LENSING-GALAXY INTEGRAND (Gg)
    function sjclscrossobjnoweight(obj,zlens)
        type(my_type) :: obj
        REAL(mcp) :: zlens,sjclscrossobjnoweight,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez
        real(mcp) :: kminsj,kmaxsj

        hubblez = Hofz(zlens)
        distz = ComovingRadialDistance(zlens)
        obj%distlensflatomp = distz
        distz = f_K(distz) !with curvature
        kval = (obj%ellarromp(obj%ellgenomp)+1.0d0/2.0d0)/((obj%h0omp/100.0d0)*distz)

        kminsj = 1.0d-5
        kmaxsj = 100.0d0
        deltafinal = 0.0d0
        if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = NL_MPKPowerAt(kval,zlens)
        sjclscrossobjnoweight = 1.0d0/(distz**2.0d0)/hubblez*deltafinal/(obj%homp**3.0d0)

    END function sjclscrossobjnoweight


    !LENSING-GALAXY INTEGRAND (Gg)
    function sjclscrossobjonlyweight(obj,zlens)
        type(my_type) :: obj
        REAL(mcp) :: zlens,sjclscrossobjonlyweight,weight2,hubblez
        real(mcp) :: kminsj,kmaxsj

        hubblez = Hofz(zlens)

        if(obj%uselensomp == 0) then
            weight2 = plensobjcubiclowz(obj,zlens) !single spec bin
        end if
        if(obj%uselensomp == 1) then
            weight2 = plensobjcubiccmass(obj,zlens) !single spec bin
        end if
        if(obj%uselensomp == 2) then
            weight2 = plensobjcubic2dfloz(obj,zlens) !single spec bin
        end if
        if(obj%uselensomp == 3) then
            weight2 = plensobjcubic2dfhiz(obj,zlens) !single spec bin
        end if
      	sjclscrossobjonlyweight = weight2*obj%bbiazconstomp*hubblez !constant bias within lens-bin

    END function sjclscrossobjonlyweight


    !LENSING-GALAXY INTEGRAND (Gg)
    function sjclscrossobj(obj,zlens)
        type(my_type) :: obj
        REAL(mcp) :: zlens,sjclscrossobj,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez,ckms
        real(mcp) :: kminsj,kmaxsj

        ckms = 299792.458d0
        distz = ComovingRadialDistance(zlens)
        obj%distlensflatomp = distz
        distz = f_K(distz) !with curvature

        kval = (obj%ellarromp(obj%ellgenomp)+1.0d0/2.0d0)/((obj%h0omp/100.0d0)*distz)

        kminsj = 1.0d-5
        kmaxsj = 100.0d0
        deltafinal = 0.0d0
        if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = NL_MPKPowerAt(kval,zlens)

        weightpart = 3.0d0/2.0d0*(obj%omdmomp+obj%ombomp)*(obj%h0omp/ckms)**2.0d0*distz*(1.0d0+zlens)
        obj%wchooseomp = 1 !doesn't matter here because m1 = m2
        weight1 = weightpart*rombint_obj(obj,weightobjcubic,zlens,obj%mumaxwomp,obj%tol_erroromp) !note mumaxw

        if(obj%uselensomp == 0) then
            weight2 = plensobjcubiclowz(obj,zlens) !single spec bin
        end if
        if(obj%uselensomp == 1) then
            weight2 = plensobjcubiccmass(obj,zlens) !single spec bin
        end if
        if(obj%uselensomp == 2) then
            weight2 = plensobjcubic2dfloz(obj,zlens) !single spec bin
        end if
        if(obj%uselensomp == 3) then
            weight2 = plensobjcubic2dfhiz(obj,zlens) !single spec bin
        end if

!the hubble term that should be multiplying plensobjcubic cancels the 1/hubble term that should be below
        sjclscrossobj = weight1*weight2/(distz**2.0d0)*deltafinal/(obj%homp**3.0d0)*obj%bbiazconstomp !constant bias within lens-bin

    END function sjclscrossobj


    !LENSING-GALAXY INTEGRAND (Gg) ---- wrt scale factor instead
    function sjclscrossobjsf(obj,alens)
        type(my_type) :: obj
        REAL(mcp) :: alens,zlens,sjclscrossobjsf,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez,ckms
        real(mcp) :: kminsj,kmaxsj

        zlens = 1.0d0/alens - 1.0d0
        ckms = 299792.458d0
        distz = ComovingRadialDistance(zlens)
        obj%distlensflatomp = distz
        distz = f_K(distz) !with curvature

        kval = (obj%ellarromp(obj%ellgenomp)+1.0d0/2.0d0)/((obj%h0omp/100.0d0)*distz)

        kminsj = 1.0d-5
        kmaxsj = 100.0d0
        deltafinal = 0.0d0
        if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = NL_MPKPowerAt(kval,zlens)

        weightpart = 3.0d0/2.0d0*(obj%omdmomp+obj%ombomp)*(obj%h0omp/ckms)**2.0d0*distz*(1.0d0+zlens)
        obj%wchooseomp = 1 !doesn't matter here because m1 = m2
        weight1 = weightpart*rombint_obj(obj,weightobjcubic,zlens,obj%mumaxwomp,obj%tol_erroromp) !note mumaxw

        if(obj%uselensomp == 0) then
            weight2 = plensobjcubiclowz(obj,zlens) !single spec bin
        end if
        if(obj%uselensomp == 1) then
            weight2 = plensobjcubiccmass(obj,zlens) !single spec bin
        end if
        if(obj%uselensomp == 2) then
            weight2 = plensobjcubic2dfloz(obj,zlens) !single spec bin
        end if
        if(obj%uselensomp == 3) then
            weight2 = plensobjcubic2dfhiz(obj,zlens) !single spec bin
        end if

!the hubble term that should be multiplying plensobjcubic cancels the 1/hubble term that should be below
        sjclscrossobjsf = ((1.0d0+zlens)**2.0d0)*weight1*weight2/(distz**2.0d0)*deltafinal/(obj%homp**3.0d0)*obj%bbiazconstomp !constant bias within lens bin

    END function sjclscrossobjsf


     !INTRINSIC-GALAXY INTEGRAND (Ig) ---- wrt scale factor instead
     function sjclscrossgiobjsf(obj,alens)
         type(my_type) :: obj
         REAL(mcp) :: alens,zlens,sjclscrossgiobjsf,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez,ckms,growthnorm
         real(mcp) :: kminsj,kmaxsj
 
         zlens = 1.0d0/alens - 1.0d0
         ckms = 299792.458d0
         distz = ComovingRadialDistance(zlens)
         obj%distlensflatomp = distz
         distz = f_K(distz) !with curvature
 
         kval = (obj%ellarromp(obj%ellgenomp)+1.0d0/2.0d0)/((obj%h0omp/100.0d0)*distz)
 
         kminsj = 1.0d-5
         kmaxsj = 100.0d0
         deltafinal = 0.0d0
         if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = NL_MPKPowerAt(kval,zlens)
 
         growthnorm = obj%gcubicomp%Value(zlens) !this is D(0)/D(z)
         deltafinal = deltafinal*(-obj%ampiayesomp*5.0d-14*2.77536627d11*(obj%omdmomp+obj%ombomp)*growthnorm*((1.0d0+zlens)/(1.0d0+0.3d0))**obj%redziayesomp)*(obj%lumarromp(obj%momp(1)))**obj%lumiayesomp
         hubblez = Hofz(zlens)
 
         obj%wchooseomp = 1 !doesn't matter here because m1 = m2
         weight1 = psourceobjcubic(obj,zlens)*hubblez
 
         if(obj%uselensomp == 0) then
             weight2 = plensobjcubiclowz(obj,zlens) !single spec bin
         end if
         if(obj%uselensomp == 1) then
             weight2 = plensobjcubiccmass(obj,zlens) !single spec bin
         end if
         if(obj%uselensomp == 2) then
             weight2 = plensobjcubic2dfloz(obj,zlens) !single spec bin
         end if
         if(obj%uselensomp == 3) then
             weight2 = plensobjcubic2dfhiz(obj,zlens) !single spec bin
         end if
 
 !the hubble term that should be multiplying plensobjcubic cancels the 1/hubble term that should be below
         sjclscrossgiobjsf = ((1.0d0+zlens)**2.0d0)*weight1*weight2/(distz**2.0d0)*deltafinal/(obj%homp**3.0d0)*obj%bbiazconstomp !constant bias within lens bin
 
     END function sjclscrossgiobjsf


    !lensing weight, cubic spline of source distribution
    function weightobjcubic(obj,zs)
        type(my_type) :: obj
        REAL(mcp) :: zs,weightobjcubic,chis,diff,ckms

        chis = ComovingRadialDistance(zs)
        diff = f_K(chis-obj%distlensflatomp) !with curvature

        weightobjcubic=0.0d0
        if((zs-obj%aphotarromp(obj%momp(obj%wchooseomp))) >= obj%mumin2vomp)  weightobjcubic = obj%psourcetypeomp(obj%momp(obj%wchooseomp))%Value(zs-obj%aphotarromp(obj%momp(obj%wchooseomp)))
        weightobjcubic = weightobjcubic*diff/f_K(chis)

    END function weightobjcubic


    !cubic spline of source distribution
    function psourceobjcubic(obj,zs)
        type(my_type) :: obj
        REAL(mcp) :: zs,psourceobjcubic

        psourceobjcubic=0.0d0
        if((zs-obj%aphotarromp(obj%momp(obj%wchooseomp))) >= obj%mumin2vomp)  psourceobjcubic = obj%psourcetypeomp(obj%momp(obj%wchooseomp))%Value(zs-obj%aphotarromp(obj%momp(obj%wchooseomp)))

    END function psourceobjcubic


    !cubic spline of lowz lens distribution
    function plensobjcubiclowz(obj,zs)
        type(my_type) :: obj
        REAL(mcp) :: zs,plensobjcubiclowz

        plensobjcubiclowz = obj%plenstypeomplowz%Value(zs)

    END function plensobjcubiclowz

    !cubic spline of cmass lens distribution
    function plensobjcubiccmass(obj,zs)
        type(my_type) :: obj
        REAL(mcp) :: zs,plensobjcubiccmass

        plensobjcubiccmass = obj%plenstypeompcmass%Value(zs)

    END function plensobjcubiccmass

    !cubic spline of 2dfloz lens distribution
    function plensobjcubic2dfloz(obj,zs)
        type(my_type) :: obj
        REAL(mcp) :: zs,plensobjcubic2dfloz

        plensobjcubic2dfloz = obj%plenstypeomp2dfloz%Value(zs)

    END function plensobjcubic2dfloz

    !cubic spline of 2dfhiz lens distribution
    function plensobjcubic2dfhiz(obj,zs)
        type(my_type) :: obj
        REAL(mcp) :: zs,plensobjcubic2dfhiz

        plensobjcubic2dfhiz = obj%plenstypeomp2dfhiz%Value(zs)

      END function plensobjcubic2dfhiz

      function Matrix_QuadForm(Mat,vec)
        !Get vec^T*Mat*vec where Mat is symmetric
        real(mcp) Matrix_QuadForm
        real(mcp), intent(in) :: vec(:)
        real(mcp), intent(in) :: Mat(:,:)
        real(mcp) :: out(size(vec))
        
        !n=size(vec)
        !allocate(out(n))
        !call Matrix_MulVecSymm(Mat,vec,out)
        out = matmul(Mat,vec)
        Matrix_QuadForm = dot_product(out, vec)
        !deallocate(out)

      end function Matrix_QuadForm

      



    END function CosmoLSS_LnLike


    !----------------------------------------------------------------------------------------------------------
    ! functions for calculating  2D (linear) power spectrum in redshift space with AP distortions:
    !
    ! 	 see W. E. Ballinger 1996 (http://arxiv.org/pdf/astro-ph/9605017v1.pdf) for derivation or Beutler 2014 for formula used
    ! We perform the scaling from observed to true k and mu, before calcualting the 2D power spectrum.
    !
    !----------------------------------------------------------------------------------------------------------

    function LinearPower_2D_s_cmass(k_obs,mu_obs,z,sigma_v_cmass,a_perp_cmass,a_par_cmass,b1_cmass,N_shot_cmass)
    implicit none
    !class(TCosmoTheoryPredictions), target :: Theory
    real(mcp),   intent(in) :: z
    REAL(mcp), intent(in) :: k_obs,mu_obs,b1_cmass
    REAL(mcp), intent(in) :: a_par_cmass,N_shot_cmass,sigma_v_cmass,a_perp_cmass
    REAL(mcp) :: LinearPower_2D_s_cmass,F_AP,beta,growth_rate
    REAL(mcp) :: mu_TRUE,k_TRUE,scaling_factor,FoG_term

    F_AP = a_par_cmass/a_perp_cmass
    scaling_factor = 1.0 + (mu_obs**2)*(1.0/F_AP**2 - 1.0)
    mu_TRUE = (mu_obs/F_AP)*(scaling_factor)**(-0.5)
    k_TRUE =  (k_obs/a_perp_cmass)*(scaling_factor)**(0.5)
    !growth_rate = Theory%growth_z%Value(z)/Theory%sigma8_z%Value(z)     !! AL gives f\sigma8, not f
    !TT Check: We want f here, right?
    growth_rate = GrowthRate%Value(1d0/(1d0+z))
    beta =  growth_rate/b1_cmass                                        !! Can be updated to have scale-dependence.
    FoG_term = exp(-(k_TRUE*mu_TRUE*sigma_v_cmass)**2)
    LinearPower_2D_s_cmass = &
    b1_cmass**2*(MPKPowerAt(k_TRUE,z) + N_shot_cmass)*(1.0 + beta*mu_TRUE**2)**2*FoG_term
    !write (*,*) z, k_obs,k_TRUE, MPKPowerAt(k_TRUE,z)
    !write (*,*) F_AP, mu_obs,a_par_cmass,a_perp_cmass
    !write (*,*) '...'
    
    END function LinearPower_2D_s_cmass


    !----------------------------------------------------------------------------------------------------------
    ! Legendre polynomials needed
    !----------------------------------------------------------------------------------------------------------

    function Legendre_Pn(mu,n)
    implicit none
    integer :: n
    real(mcp) :: Legendre_Pn,mu

    if(n .EQ. 0) Legendre_Pn = 1.0
    if(n .EQ. 2) Legendre_Pn = (1.0/2.0)*(3.0*mu**2 - 1.0)
    if(n .EQ. 4) Legendre_Pn = (1.0/8.0)*(35.0*mu**4 - 30.0*mu**2 + 3.0)

    end function Legendre_Pn

    !----------------------------------------------------------------------------------------------------------
    ! Functions needed for AP calculation
    !----------------------------------------------------------------------------------------------------------

    function D_AhUnit(z,CMB)
        !class(CMBParams) :: CMB
        type(myCMBparam) :: CMB
        real(mcp) :: z,D_AhUnit

        D_AhUnit = AngularDiameterDistance(z)*(CMB%h0/100d0)

    end function D_AhUnit

    function HofzhUnit(z,CMB,c)
        !class(CMBParams) :: CMB
      type(myCMBparam), intent(in) :: CMB 
      real(mcp), intent(in) :: z, c
      real(mcp) :: HofzhUnit

      HofzhUnit = (c*Hofz(z)/1.d3)/(CMB%h0/100d0)

    end function HofzhUnit

    function numcat_local(S, num)
        character(LEN=*) S
        character(LEN=1024) numcat_local, numstr
        integer num

        write (numstr, *) num
        numcat_local = trim(S) // trim(adjustl(numstr))

    end function numcat_local


end module LogLikeCosmoLSS_module
