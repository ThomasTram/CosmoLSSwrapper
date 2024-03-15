"""
.. module:: CosmoLSS
    :synopsis: CosmoLSS likelihood from https://arxiv.org/pdf/1707.06627.pdf, https://github.com/sjoudaki/CosmoLSS

.. moduleauthor:: Thomas Tram <thomas.tram@aias.au.dk>
Last updated December 20, 2017. Based on the CosmoMC module.
"""
import math
import numpy as np
import scipy.linalg as la
import scipy.special
import montepython.io_mp as io_mp
import os
from montepython.likelihood_class import Likelihood_sn
from .pyLogLikeCosmoLSS import interpolation_init_all, interpolation_free_all, set_this, loglkl_from_fortran, set_stuff, set_mp_overlap, set_sources

T_CMB = 2.7255     #CMB temperature
h = 6.62606957e-34     #Planck's constant
kB = 1.3806488e-23     #Boltzmann constant
Ghz_Kelvin = h/kB*1e9  #GHz Kelvin conversion

        
class CosmoLSS(Likelihood_sn):
    """ Class for the MontePython version of the CosmoLSS likelihood"""

    class MultipoleOverlap:
        """ Class for dealing with multipole overlaps."""
        def __init__(self, k_min_theory=0.0, dk_theory=0.05, z_eff=0.57, size_convolution=30, k_num_conv=10,k_spacing_obs=0.05,k_min_obs=0.075, k_fit=125):
            self.k_min_theory = k_min_theory
            self.dk_theory = dk_theory
            self.z_eff = z_eff
            self.size_convolution = size_convolution
            self.k_num_conv = k_num_conv
            self.k_spacing_obs = k_spacing_obs
            self.k_min_obs = k_min_obs
            self.fit_k_0075 = True if k_fit==75 else False
            self.fit_k_0125 = True if k_fit==125 else False
            self.fit_k_0175 = True if k_fit==175 else False

            if self.fit_k_0175:
                self.k_num_obs = 3
                self.size_cov = 9
            elif self.fit_k_0125:
                self.k_num_obs = 2
                self.size_cov = 4
            elif self.fit_k_0075:
                self.k_num_obs = 1
                self.size_cov = 1
            else:
                self.k_num_obs = 0
                self.size_cov = 0
        
        def ReadConvolutionMatrix(self, datafile):
            self.ConvolutionMatrix = np.loadtxt(datafile,usecols=(2,)).reshape(self.size_convolution, self.size_convolution)

    def __init__(self, path, data, command_line):
        # Unusual construction, since the data files are not distributed
        # alongside CosmoLSS (size problems)
        try:
            # Read the .dataset file specifying the data.
            super(CosmoLSS, self).__init__(path, data, command_line)
        except IOError:
            raise io_mp.LikelihoodError('The CosmoLSS data files were not found.')
        
        arguments = {
            'output': 'mPk',
            'non linear':'HALOFIT',
            'z_pk':'0.0, 100.0',
            'P_k_max_h/Mpc':100.0
            }
        self.need_cosmo_arguments(data, arguments)

        # All this info should be inferred from the data files
        if False:
            if (self.set_scenario == 3) and not self.use_conservative:
                self.size_covallmask = 210 #!postmask, fiducial
            else:
                self.size_covallmask = 186 #!postmask, conservative
            
            if (self.set_scenario == 1): #!fiducial KiDS masking
                sizcovishpremask = 180 #!pre-masking
                sizcovish = self.size_covmask #!post-masking
            elif ((self.set_scenario == 3) and not self.use_conservative):
                sizcovishpremask = 340 #!pre-masking
                sizcovish = self.size_covallmask #!post-masking
            elif ((self.set_scenario == 3) and self.use_conservative):
                sizcovishpremask = 332 #!pre-masking
                sizcovish = self.size_covallmask #!post-masking
            #What if scenario is not 1 or 3?

            sizcovishsq = sizcovish**2
            self.sizcov = sizcovish
            self.sizcovpremask = sizcovishpremask

        # !!!!!Lens redshifts!!!!!!!
        self.lens_redshifts = {}
        for samplename in ['cmass','lowz','2dfloz','2dfhiz']:
            self.lens_redshifts[samplename] = np.loadtxt(self.data_directory+'/lensingrsdfiles/nz_'+samplename+'_modelsj.dat')
            self.lens_redshifts[samplename][:,1] /= (np.sum(self.lens_redshifts[samplename][:,1])*0.01)
                                                                    
        # !!!Reading in source distributions
        data_list = []
        if self.use_bootstrapnz and self.set_scenario >= 0:
            for i in range(1, 5):
                filename = self.data_directory + f'/lensingrsdfiles/nz_z{i}_kids_binned_hendrik.dat'
                data = np.loadtxt(filename)
                data_list.append(data)
        elif not self.use_bootstrapnz and self.set_scenario >= 0:
            for i in range(1, 5):
                filename = self.data_directory + f'/lensingrsdfiles/nz_z{i}_kids_binned.dat'
                data = np.loadtxt(filename)
                data_list.append(data)
        elif self.set_scenario == -2:
            for i in range(1, 6):
                filename = self.data_directory + f'/lensingrsdfiles/nz_z{i}_kv450_meannew_blindb.dat'
                data = np.loadtxt(filename)
                data_list.append(data)
        elif self.set_scenario == -3:
            for i in range(1, 6):
                filename = self.data_directory + f'/lensingrsdfiles/kids1000nz{i}.dat'
                data = np.loadtxt(filename)
                data_list.append(data)
        
        # Normalize data_list
        for data in data_list:
            data[:, 1] /= np.sum(data[:, 1]) * 0.05  # Normalize the second column
        
        # Concatenate data from list of arrays
        self.sources_for_scenario = np.array(data_list, order='C')
        # We must call set_sources at the end, after this%set_scenario has been set
        
        #!!!Reading in measurements and masks and covariances
        if self.set_scenario == 1:
            if self.use_analyticcov:
                if self.use_large_scales:
                    xipmtemp = np.loadtxt(self.data_directory+'/lensingrsdfiles/xipmcut_kids_blind1_planckcut.dat',usecols=(1,))
                    covxipmtemp = np.loadtxt(self.data_directory+'/lensingrsdfiles/xipmcutcov_kids_analytic_inc_m_blind1_planckcut.dat',usecols=(2,))
                else:
                    xipmtemp = np.loadtxt(self.data_directory+'/lensingrsdfiles/xipmcut_kids_blind1.dat',usecols=(1,))
                    covxipmtemp = np.loadtxt(self.data_directory+'/lensingrsdfiles/xipmcutcov_kids_analytic_inc_m_blind1.dat',usecols=(2,))
            else:
                xipmtemp = np.loadtxt(self.data_directory+'/lensingrsdfiles/xipmcut_kids_regcomb_blind2_swinburnesj.dat',usecols=(1,))
                covxipmtemp = np.loadtxt(self.data_directory+'/lensingrsdfiles/xipmcutcov_kids_regcomb_blind2_swinburnesj.dat',usecols=(2,))
            if self.use_large_scales:
                masktemp =  np.loadtxt(self.data_directory+'/lensingrsdfiles/xipm_kids4tom_selectsj_planckcut.dat',usecols=(1,))
            else:
                masktemp =  np.loadtxt(self.data_directory+'/lensingrsdfiles/xipm_kids4tom_selectsj.dat',usecols=(1,))
        elif self.set_scenario == 2:
            if self.use_conservative:
                xipmtemp = np.loadtxt(self.data_directory+'/lensingrsdfiles/xipmgtlarge3_kids_regcomb_blind2sj.dat',usecols=(1,))
                masktemp =  np.loadtxt(self.data_directory+'/lensingrsdfiles/xipmgtlarge3_kids4tom_selectsj.dat',usecols=(1,))
                covxipmtemp = np.loadtxt(self.data_directory+'/lensingrsdfiles/xipmgtlarge3cov_kids_regcomb_blind2sj.dat',usecols=(2,))
            else:
                xipmtemp = np.loadtxt(self.data_directory+'/lensingrsdfiles/xipmgtlarge4_kids_regcomb_blind2sj.dat',usecols=(1,))
                masktemp =  np.loadtxt(self.data_directory+'/lensingrsdfiles/xipmgtlarge4_kids4tom_selectsj.dat',usecols=(1,))
                covxipmtemp = np.loadtxt(self.data_directory+'/lensingrsdfiles/xipmgtlarge4cov_kids_regcomb_blind2sj.dat',usecols=(2,))
        elif self.set_scenario == 3:
            if self.use_conservative:
                xipmtemp = np.loadtxt(self.data_directory+'/lensingrsdfiles/xipmgtlarge7_kids_regcomb_blind2sj.dat',usecols=(1,))
                masktemp =  np.loadtxt(self.data_directory+'/lensingrsdfiles/xipmgtlarge7_kids4tom_selectsj.dat',usecols=(1,))
                covxipmtemp = np.loadtxt(self.data_directory+'/lensingrsdfiles/xipmgtlarge7cov_kids_regcomb_blind2sj.dat',usecols=(2,))
            else:
                xipmtemp = np.loadtxt(self.data_directory+'/lensingrsdfiles/xipmgtlarge4_kids_regcomb_blind2sj.dat',usecols=(1,))
                masktemp =  np.loadtxt(self.data_directory+'/lensingrsdfiles/xipmgtlarge4_kids4tom_selectsj.dat',usecols=(1,))
                covxipmtemp = np.loadtxt(self.data_directory+'/lensingrsdfiles/xipmgtlarge4cov_kids_regcomb_blind2sj.dat',usecols=(2,))
        elif self.set_scenario == -2:
            xipmtemp = np.loadtxt(self.data_directory+'/lensingrsdfiles/kv450datanew.dat')
            masktemp =  np.loadtxt(self.data_directory+'/lensingrsdfiles/xipm_kids5tom_selectsj.dat',usecols=(1,))
            covxipmtemp = np.loadtxt(self.data_directory+'/lensingrsdfiles/thps_cov_aug13_blindB_list_nospaces.dat',usecols=(2,))                
        elif (self.set_scenario == -3) :
            xipmtemp = np.loadtxt(self.data_directory+'/lensingrsdfiles/kids1000xipm_desformat.dat')
            if self.use_cholesky:
                self.somdzcholesky = np.loadtxt(self.data_directory+'/lensingrsdfiles/kids1000_SOM_cov_multiplied_cholesky.dat')
            masktemp =  np.loadtxt(self.data_directory+'/lensingrsdfiles/kids1000maskxipm.dat',usecols=(1,))
            covxipmtemp = np.loadtxt(self.data_directory+'/lensingrsdfiles/kids1000cov_xipmsubspace.dat')

        #!Else missing
        if (self.set_scenario == -2 or self.set_scenario == -3):
            if self.set_scenario == -2:
                multarr_values = [-0.0128, -0.0104, -0.0114, 0.0072, 0.0061]
            elif self.set_scenario == -3:
                multarr_values = [-0.009, -0.011, -0.015, 0.002, 0.007] 

            # Construct division factors
            division_factors = [(1.0 + multarr_values[i]) / (1.0 + multarr_values[j]) for i in range(5) for j in range(i+1)]

            # Apply division factors to xipmtemp4ss
            for i in range(15):
                xipmtemp[18*i:18*(i + 1)] /= division_factors[i]


        #Store sizes:
        self.sizcov = math.isqrt(covxipmtemp.shape[0])
        self.sizcovpremask = masktemp.shape[0]
        self.size_covallmask = self.sizcov


        self.xipm = xipmtemp
        #print masktemp.shape
        self.maskelements = masktemp
        #print covxipmtemp.shape, sizcovish
        self.covxipm = covxipmtemp.reshape(self.sizcov, self.sizcov)
        # Invert covariance matrix
        self.invcovxipm = la.inv(self.covxipm)

        #!converted from arcmin to degrees
        if self.set_scenario >= 0:
            thetacfhtini_deg = np.array([0.71336, 1.45210, 2.95582, 6.01675, 12.24745, 24.93039, 50.74726, 103.29898, 210.27107]) / 60.0
        elif self.set_scenario == -2:
            thetacfhtini_deg = np.array([0.7588895, 1.544764, 3.144456, 6.400723, 13.02904, 26.52137, 53.98579, 109.8912, 223.6899]) / 60.0
        elif self.set_scenario == -3:
            thetacfhtini_deg = np.array([0.7588895, 1.544764, 3.144456, 6.400723, 13.02904, 26.52137, 53.98579, 109.8912, 223.6899]) / 60.0

        self.thetaradcfhtini = thetacfhtini_deg * np.pi / 180.0
        nell = 58999
        self.ellgentestarrini = np.arange(2, nell + 2, dtype='float64')

        #!Generate array containing l-values where C(l) is evaluated
        if self.use_morell:
            self.nellbins = 101
            morellfac = 0.1
        else:
            self.nellbins = 31
            morellfac = 0.5
        self.ellarr = np.zeros(self.nellbins)
        self.prefacarrz = np.zeros(self.nellbins)
        for ellgenini in range(1, self.nellbins + 1):
            if ellgenini < 10:
                self.ellarr[ellgenini - 1] = ellgenini + 1
            elif ellgenini > 10 or ellgenini == 10:
                self.ellarr[ellgenini - 1] = self.ellarr[ellgenini - 2] + morellfac*self.ellarr[ellgenini - 2]

        self.ellarr = np.floor(self.ellarr)
        self.prefacarrz = ((self.ellarr + 2.0) * (self.ellarr + 1.0) * self.ellarr * (self.ellarr - 1.0)) ** 0.5 / (self.ellarr + 0.5) ** 2.0

        #!compute Bessel functions
        # I want a column order Fortran contiguous memory block of the form (iii, jjj), so most convenient to form the argument first.
        # We compute the outer product and flatten the matrix
        besarg = np.outer(self.ellgentestarrini,self.thetaradcfhtini).flatten(order='C')
        self.bes0arr = scipy.special.jn(0,besarg)
        self.bes2arr = scipy.special.jn(2,besarg)
        self.bes4arr = scipy.special.jn(4,besarg)

        # Initialise multipole overlaps
        # Default settings for MultipoleOverlap Class, can be overwritten at initialisation:
        # k_min_theory=0.0, dk_theory=0.05, z_eff=0.57, size_convolution=30, k_num_conv=10,k_spacing_obs=0.05,k_min_obs=0.075,k_0175 = False, fit_k_0125 = False,fit_k_0075 = False
        CosmoLSS.mp_overlaps = {'cmass':self.MultipoleOverlap(z_eff=0.57, k_fit=self.kfit_cmass),
                                    'lowz': self.MultipoleOverlap(z_eff=0.32,k_fit=self.kfit_lowz),
                                    '2dfloz': self.MultipoleOverlap(z_eff=0.31,k_fit=self.kfit_2dfloz),
                                    '2dfhiz':self.MultipoleOverlap(z_eff=0.56,k_fit=self.kfit_2dfhiz)}

        
        # Read convolution matrix for the used overlaps from files and push to Fortran
        name_to_number = {'cmass':1,'lowz':2,'2dfloz':3,'2dfhiz':4}
        for key, value in self.mp_overlaps.items():
            if key=='lowz':
                fname_from_dataset = self.LOWZ_overlap_conv
            elif key=='cmass':
                fname_from_dataset = self.CMASS_overlap_conv
            elif key=='2dfloz':
                fname_from_dataset = self.twodfloz_overlap_conv
            elif key=='2dfhiz':
                fname_from_dataset = self.twodfhiz_overlap_conv
            fname_from_dataset = fname_from_dataset.strip(r'%DATASETDIR%')
            fname = os.path.join(self.data_directory,fname_from_dataset)
            value.ReadConvolutionMatrix(fname)
            # Push data pointers to Fortran
            intParams = np.array([value.k_num_obs,value.size_convolution,value.k_num_conv],dtype='int32')
            realParams = np.array([value.z_eff, value.k_min_theory, value.dk_theory, value.k_spacing_obs,value.k_min_obs])
            set_mp_overlap(value.ConvolutionMatrix, value.size_convolution, intParams, realParams, name_to_number[key])
            
        # Set the logk and z arrays which will de used in the power spectrum interpolation.
        self.Nk = 203
        self.Nz = 38
        self.logkh_for_pk = np.linspace(-4,2.2,self.Nk)
        self.kh_for_pk = 10**self.logkh_for_pk
        self.z_for_pk = np.linspace(0,6.1,self.Nz)

        # Push some fields in the python class to the corresponding derived type in Fortran, called this:
        intParams = np.array([self.size_cov, self.sizcov, self.klinesum, self.set_scenario, self.size_covallmask],dtype='int32')
        logParams = np.array([self.use_morell, self.use_rombint, self.use_conservative,
                                  self.use_cmass_overlap, self.use_lowz_overlap, self.use_2dfloz_overlap, self.use_2dfhiz_overlap,
                                  self.use_analyticcov, self.use_largescales, self.use_bootstrapnz, self.write_theoryfiles, self.use_accuracyboost, 
                                  self.print_timelike, self.use_cholesky, self.print_parameters],dtype='int32')

        set_this(self.lens_redshifts['lowz'], self.lens_redshifts['2dfloz'], self.lens_redshifts['cmass'],self.lens_redshifts['2dfhiz'],
                     self.xipm, self.invcovxipm, self.sizcov, self.ellarr, self.prefacarrz, self.nellbins, 
                     self.ellgentestarrini, self.maskelements, len(self.maskelements),
                     self.bes0arr, self.bes4arr, self.bes2arr, intParams, logParams)
        
        # Set sources
        set_sources(self.sources_for_scenario)
        
    def loglkl(self, cosmo, data):
        """
        Construct and pass relevant pointers to Fortran and call the fortran likelihood
        """
        #from pyLogLikeCosmoLSS import interpolation_init_all
        # Recover background evolution
        bg = cosmo.get_background()

        # Recover cosmological parameters and pass them
        H0 = cosmo.Hubble(0.)
        h = cosmo.h()
        Omega0_k = cosmo.Omega0_k()
        if np.abs(Omega0_k)==0.:
            CPr = 1.0
        else:
            CPr = 1./np.sqrt(np.abs(Omega0_k)*H0**2)
        omcdm = cosmo.Omega0_cdm()
        omb = cosmo.Omega_b()
        set_stuff(CPr, H0, h, omcdm, omb, Omega0_k*h*h, np.sign(-Omega0_k))

         # Recover linear and nonlinear powerspectra from CLASS using the fast methods:
        pk_l = h**3*cosmo.get_pk_array(self.kh_for_pk*h,self.z_for_pk,self.Nk,self.Nz,nonlinear=False).reshape(self.Nz,self.Nk)
        pk_nl = h**3*cosmo.get_pk_array(self.kh_for_pk*h,self.z_for_pk,self.Nk,self.Nz,nonlinear=True).reshape(self.Nz,self.Nk)

       
         # Initialise interpolations. Remember that interpolation x-arrays must be increasing, and passing a reverse view using [::-1]
        # will not work since we are relying on the .data pointer inside cython. So use a for the background instead.
        interpolation_init_all(1./(1.+bg['z']), bg['H [1/Mpc]'], bg['comov. dist.'], bg['gr.fac. D'],  bg['gr.fac. f'], len(bg['z']),
                                   self.logkh_for_pk, self.z_for_pk, pk_l, pk_nl, self.Nk, self.Nz)

        # Update bootstrap sources with integer corresponding to bootstrap realisation 0-999
        # Update: I think we are no longer using bootstraps, so we will set the sources in the _init_ method
        #self.update_sources(np.random.randint(0,999))
       
        #Finally, recover nuisance parameters and call the Fortran likelihood:
        nuisance = {}
        for name in self.use_nuisance:
            nuisance[name] = data.mcmc_parameters[name]['current']*data.mcmc_parameters[name]['scale']
        
        loglkl = loglkl_from_fortran(**nuisance)
        # Free interpolation structures
        interpolation_free_all()
        # Add priors (Please check!)
        loglkl += (nuisance['dcamp']/0.00023)**2
        # Remember correct sign
        return -loglkl