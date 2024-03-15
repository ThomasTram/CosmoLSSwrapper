module LogLikeCosmoLSS1_interface
  use iso_c_binding, only: c_double, c_int32_t, c_bool
  use LogLikeCosmoLSS_module, only: Interpolation_Init_All, Interpolation_Free_All, NonlinearPowerspectrum, this, CMB, set_mp_overlap, CPr, Curvature, CosmoLSS_LnLike
  implicit none
contains
  subroutine c_Interpolation_Init_All(z, H, conf_dist, D, f, size_z, kvec, zvec, pk, pk_nl, size_kvec, size_zvec) bind(c)
    real(c_double), dimension(size_z), intent(in) :: z, H, conf_dist, D, f
    integer(c_int32_t), intent(in) :: size_z, size_kvec, size_zvec
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
    integer(c_int32_t), intent(in) :: size_ConMat,  which_sample
    real(c_double), dimension(size_ConMat,size_ConMat), intent(in) :: ConMat
    integer(c_int32_t), intent(in), dimension(3) :: intParams
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

  subroutine c_set_this(lenslowz, lens2dfloz, lenscmass, lens2dfhiz, xipm, invcovxipm, sizcov, ellarr, prefacarrz, size_ellarr, ellini, maskelements, size_maskelements, bes0arr,bes4arr,bes2arr, intParams, logParamsInt) bind(c)
     real(c_double), dimension(2,29), intent(in) :: lenslowz, lens2dfloz
     real(c_double), dimension(2,28), intent(in) :: lenscmass, lens2dfhiz
     integer(c_int32_t), intent(in) :: sizcov, size_maskelements, size_ellarr
     real(c_double), dimension(sizcov), intent(in) :: xipm
     real(c_double), dimension(size_ellarr), intent(in) :: ellarr, prefacarrz
     real(c_double), dimension(sizcov,sizcov), intent(in) :: invcovxipm 
     real(c_double), dimension(58999), intent(in) :: ellini
     real(c_double), dimension(size_maskelements), intent(in) :: maskelements
     real(c_double), dimension(9,58999), intent(in) :: bes0arr,bes4arr,bes2arr
     integer(c_int32_t), dimension(5), intent(in) :: intParams
     integer, parameter :: N_logParams = 15
     integer(c_int32_t), dimension(N_logParams), intent(in) :: logParamsInt
     logical, dimension(N_logParams) :: logParams
     integer :: zinit, zim, setscenario
     
     this%sizcov = sizcov
     this%sizcovpremask = size_maskelements
     this%size_cov = intParams(1)
     this%size_covmask = intParams(2)
     this%klinesum = intParams(3)
     this%set_scenario = intParams(4)
     this%size_covallmask = intParams(5)

     ! Convert the integer array logParamsInt to logical array logParams
     logParams = logParamsInt .ne. 0
     this%use_morell = logParams(1)
     this%use_rombint = logParams(2)
     this%use_conservative = logParams(3)
     this%use_cmass_overlap = logParams(4)
     this%use_lowz_overlap = logParams(5)
     this%use_2dfloz_overlap = logParams(6)
     this%use_2dfhiz_overlap = logParams(7)
     this%use_analyticcov = logParams(8)
     this%use_largescales = logParams(9)
     this%use_bootstrapnz = logParams(10)
     this%write_theoryfiles = logParams(11)
     this%use_accuracyboost = logParams(12)
     this%print_timelike = logParams(13)
     this%use_cholesky = logParams(14)
     this%print_parameters = logParams(15)

     this%arraysjlenslowz = lenslowz
     this%arraysjlens2dfloz = lens2dfloz
     this%arraysjlenscmass = lenscmass
     this%arraysjlens2dfhiz = lens2dfhiz

     this%xipm = xipm
     this%invcovxipm = invcovxipm
     this%ellarr = ellarr
     this%prefacarrz = prefacarrz
     this%nellbins = size_ellarr
     this%ellgentestarrini = ellini
     this%maskelements = maskelements

     this%bes0arr = bes0arr
     this%bes4arr = bes4arr
     this%bes2arr = bes2arr
     
     setscenario = this%set_scenario
     if(setscenario >= 0) then
       do zinit=3,this%wtrapmax
         this%exact_z(zinit) = this%exact_z(zinit-1) + 0.0198d0*this%exact_z(zinit-1)
      end do
      this%exact_z(this%wtrapmax) = 3.474999d0
      else if(setscenario == -2) then
      do zinit=3,this%wtrapmax
          this%exact_z(zinit) = this%exact_z(zinit-1) + 0.021d0*this%exact_z(zinit-1)
      end do
      this%exact_z(this%wtrapmax) = 5.724999d0
      else if(setscenario == -3) then
         do zinit=3,this%wtrapmax
          this%exact_z(zinit) = this%exact_z(zinit-1) + 0.021d0*this%exact_z(zinit-1)
         end do
         this%exact_z(this%wtrapmax) = 5.674999d0
      end if

      !!!trapezoid z-steps for lowz
      this%wtrapmaxlens = 170 !assuming wtrapmaxlens < wtrapmax
      !allocate(this%exact_zlenslowz(this%wtrapmaxlens))
      this%exact_zlenslowz(1) = 0.0d0
      this%exact_zlenslowz(2) = 0.15d0
      do zim=3,this%wtrapmaxlens
         this%exact_zlenslowz(zim) = 0.15d0 + (0.43d0-0.15d0)/168.0d0*(zim-2)
      end do
      this%exact_zlenslowz(this%wtrapmaxlens) = 0.43d0 !0.429999d0

      !!!trapezoid z-steps for cmass
      !allocate(this%exact_zlenscmass(this%wtrapmaxlens))
      this%exact_zlenscmass(1) = 0.0d0
      this%exact_zlenscmass(2) = 0.43d0
      do zim=3,this%wtrapmaxlens
         this%exact_zlenscmass(zim) = 0.43 + (0.7d0-0.43d0)/168.0d0*(zim-2)
      end do
      this%exact_zlenscmass(this%wtrapmaxlens) = 0.7d0 !0.699999d0

   end subroutine c_set_this

   subroutine c_set_sources(sources_for_scenario, dim1, dim2, dim3) bind(c)
      integer(c_int32_t), intent(in) :: dim1, dim2, dim3
      real(c_double), dimension(dim1, dim2, dim3), intent(in) :: sources_for_scenario
      if (this%set_scenario >= 0) then
         this%arraysjfull4bins = sources_for_scenario
      else if (this%set_scenario == -2) then
         this%arraysjfullkv450 = sources_for_scenario
      else if (this%set_scenario == -3) then
         this%arraysjfull = sources_for_scenario
      end if
   end subroutine c_set_sources

   subroutine c_CosmoLSS_LnLike(DataParams, loglkl) bind(c)
     real(c_double), dimension(27), intent(in) :: DataParams
     real(c_double), intent(out) :: loglkl
     loglkl = CosmoLSS_LnLike(DataParams)
   end subroutine c_CosmoLSS_LnLike

   subroutine c_test_array(array, dim1, dim2, dim3) bind(c)
     integer(c_int32_t), intent(in) :: dim1, dim2, dim3
     real(c_double), dimension(dim1, dim2, dim3), intent(in) :: array
     real(c_double), dimension(3, 7, 13) :: fixed_array
     integer :: i, j, k
     fixed_array = array
      do k=1,13
        do j=1,7
           do i=1,3
              write (*,*) "Element ", i, j, k, " of array: ", fixed_array(i, j, k)
           end do
        end do
     end do
   end subroutine c_test_array
  
   subroutine c_test_array2(array) bind(c)
      real(c_double), dimension(3, 7, 13), intent(in) :: array
      real(c_double), dimension(3, 7, 13) :: fixed_array
      integer :: i, j, k
      fixed_array = array
       do k=1,13
         do j=1,7
            do i=1,3
               write (*,*) "Element ", i, j, k, " of array: ", fixed_array(i, j, k)
            end do
         end do
      end do
    end subroutine c_test_array2
end module LogLikeCosmoLSS1_interface
