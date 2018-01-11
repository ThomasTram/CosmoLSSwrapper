module LogLikeCosmoLSS1_interface
  use iso_c_binding, only: c_double, c_int, c_bool
  use LogLikeCosmoLSS_module, only: Interpolation_Init_All, Interpolation_Free_All, NonlinearPowerspectrum, this, CMB, set_mp_overlap, CPr, Curvature, CosmoLSS_LnLike
  implicit none
contains
  subroutine c_Interpolation_Init_All(z, H, conf_dist, D, f, size_z, kvec, zvec, pk, pk_nl, size_kvec, size_zvec) bind(c)
    real(c_double), dimension(size_z), intent(in) :: z, H, conf_dist, D, f
    integer(c_int), intent(in) :: size_z, size_kvec, size_zvec
    real(c_double), dimension(size_kvec), intent(in) :: kvec
    real(c_double), dimension(size_zvec), intent(in) :: zvec
    real(c_double), dimension(size_kvec,size_zvec), intent(in) :: pk, pk_nl

    call Interpolation_Init_All(z, H, conf_dist, D, f, size_z, kvec, zvec, pk, pk_nl, size_kvec, size_zvec)
  end subroutine c_Interpolation_Init_All

  subroutine c_Interpolation_Free_All() bind(c)
    call Interpolation_Free_All()
  end subroutine c_Interpolation_Free_All

  subroutine c_set_stuff(CosmoParams) bind(c)
    real(c_double), intent(in), dimension(7) :: CosmoParams
    CPr = CosmoParams(1)
    CMB%H0 = CosmoParams(3)*100!This is H0 in km/s/Mpc   CosmoParams(2)
    CMB%h = CosmoParams(3)
    CMB%omdm = CosmoParams(4)
    CMB%omb = CosmoParams(5)
    CMB%omk = CosmoParams(6)
    if (CosmoParams(7)>0.5d0) then
       Curvature = 1
    elseif (CosmoParams(7)<-0.5d0) then
       Curvature = -1
    else
       Curvature = 0
    end if
  end subroutine c_set_stuff

  subroutine c_set_mp_overlap(ConMat, size_ConMat, intParams, realParams, which_sample) bind(c)
    integer(c_int), intent(in) :: size_ConMat,  which_sample
    real(c_double), dimension(size_ConMat,size_ConMat), intent(in) :: ConMat
    integer(c_int), intent(in), dimension(3) :: intParams
    real(c_double), intent(in), dimension(5) :: realParams
    
    select case (which_sample)
    case (1)
       call set_mp_overlap(this%cmass_mp_overlap, ConMat, size_ConMat, intParams, realParams)
    case (2)
       call set_mp_overlap(this%lowz_mp_overlap, ConMat, size_ConMat, intParams, realParams)
    case (3)
       call set_mp_overlap(this%twodfloz_mp_overlap, ConMat, size_ConMat, intParams, realParams)
    case (4)
       call set_mp_overlap(this%twodfhiz_mp_overlap, ConMat, size_ConMat, intParams, realParams)
    case default
       write (*,*) "Problem!"
    end select
  end subroutine c_set_mp_overlap

  subroutine c_set_this(lenslowz, lens2dfloz, lenscmass, lens2dfhiz, xipm, invcovxipm, sizcov, ellini, maskelements, size_maskelements, bes0arr,bes4arr,bes2arr, intParams, logParams) bind(c)
     real(c_double), dimension(2,29), intent(in) :: lenslowz, lens2dfloz
     real(c_double), dimension(2,28), intent(in) :: lenscmass, lens2dfhiz
     integer(c_int), intent(in) :: sizcov, size_maskelements
     real(c_double), dimension(sizcov), intent(in) :: xipm
     real(c_double), dimension(sizcov,sizcov), intent(in) :: invcovxipm 
     real(c_double), dimension(58999), intent(in) :: ellini
     real(c_double), dimension(size_maskelements), intent(in) :: maskelements
     real(c_double), dimension(9,58999), intent(in) :: bes0arr,bes4arr,bes2arr
     integer(c_int), dimension(5), intent(in) :: intParams
     integer(c_int), dimension(7), intent(in) :: logParams
     integer :: zinit
     
     this%sizcov = sizcov
     this%sizcovpremask = size_maskelements
     this%size_cov = intParams(1)
     this%size_covmask = intParams(2)
     this%klinesum = intParams(3)
     this%set_scenario = intParams(4)
     this%size_covallmask = intParams(5)

     this%use_morell = logParams(1)
     this%use_rombint = logParams(2)
     this%use_conservative = logParams(3)
     this%use_cmass_overlap = logParams(4)
     this%use_lowz_overlap = logParams(5)
     this%use_2dfloz_overlap = logParams(6)
     this%use_2dfhiz_overlap = logParams(7)

     this%arraysjlenslowz = lenslowz
     this%arraysjlens2dfloz = lens2dfloz
     this%arraysjlenscmass = lenscmass
     this%arraysjlens2dfhiz = lens2dfhiz

     this%xipm = xipm
     this%invcovxipm = invcovxipm
     this%ellgentestarrini = ellini
     this%maskelements = maskelements

     this%bes0arr = bes0arr
     this%bes4arr = bes4arr
     this%bes2arr = bes2arr
     
     this%num_z = 37
     allocate(this%exact_z(this%num_z)) !allocate array for matter power spectrum redshifts
     !assign redshifts for matter power spectrum
     this%exact_z(1) = 0.0d0
     this%exact_z(2) = 0.025d0
     do zinit=3,37
        this%exact_z(zinit) = this%exact_z(zinit-1) + 0.154d0*this%exact_z(zinit-1)
     end do
     this%exact_z(37) = 3.474999d0
   end subroutine c_set_this

   subroutine c_set_sources(af1, af2, af3, af4) bind(c)
     real(c_double), dimension(2,70), intent(in) :: af1,af2,af3,af4

     this%arraysjfull(:,:,1) = af1
     this%arraysjfull(:,:,2) = af2
     this%arraysjfull(:,:,3) = af3
     this%arraysjfull(:,:,4) = af4
     
   end subroutine c_set_sources

   subroutine c_CosmoLSS_LnLike(DataParams, loglkl) bind(c)
     real(c_double), dimension(19), intent(in) :: DataParams
     real(c_double), intent(out) :: loglkl
     loglkl = CosmoLSS_LnLike(DataParams)
   end subroutine c_CosmoLSS_LnLike
  
end module LogLikeCosmoLSS1_interface
