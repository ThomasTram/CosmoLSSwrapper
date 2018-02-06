module LogLikeCosmoLSS_module
  use Interpolation
  use omp_lib
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
    if (k<1d-6) then
       MPKPowerAt = 0d0
    else
       logk = max(log10(k),LinearPowerspectrum%x(1))
       MPKPowerAt =LinearPowerspectrum%Value(logk,z)
    end if
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
