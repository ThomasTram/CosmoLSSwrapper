from numpy import linspace, empty
from numpy cimport ndarray as ar
import numpy as np

cdef extern from "pyLogLikeCosmoLSS.h":
    void c_interpolation_init_all(double *z, double *H, double * conf_dist, double *D,  double *f, int *size_z, double *kvec, double *zvec,  double *pk, double *pk_nl, int *size_kvec, int*size_zvec)
    void c_interpolation_free_all()
    void c_set_stuff(double *CosmoParams)
    void c_set_mp_overlap(double *ConMat, int *size_ConMat, int *intParams, double *realParams, int *which_sample);
    void c_set_this(double *lenslowz, double *lens2dfloz, double *lenscmass, double *lens2dfhiz, double *xipm, double *invcovxipm, int *sizcov, 
        double* ellarr, double* prefacarrz, int* size_ellarr, double *ellini, double *maskelements, int *size_maskelements, double *bes0arr,double *bes4arr,double *bes2arr, int *intParams, int *logParams)
    void c_cosmolss_lnlike(double * DataParams, double *loglkl)
    void c_set_sources(double *source_data, int* dim1, int* dim2, int* dim3)
    void c_test_array(double *source_data, int* dim1, int* dim2, int* dim3)
    void c_test_array2(double *source_data)

#subroutine c_set_this(double[:] lenslowz, double[:] lens2dfloz, double[:] lenscmass, double[:] lens2dfhiz, double[:] xipm, double[:, :] invcovxipm, int sizcov, double[:] ellarr, double[:] prefacarrz, int size_ellarr, double[:] ellini, double[:] maskelements, int size_maskelements, double[:] bes0arr, double[:] bes4arr, double[:] bes2arr, int[:] intParams, int[:] logParams)
    #void c_set_sources(double *af1, double *af2, double *af3, double *af4)
    
def interpolation_init_all(ar[double,ndim=1] z, ar[double,ndim=1] H, ar[double,ndim=1] conf_dist, ar[double,ndim=1] D, ar[double,ndim=1] f, int size_z,
                               ar[double,ndim=1] kvec, ar[double,ndim=1] zvec, ar[double,ndim=2, mode="c"] pk, ar[double,ndim=2, mode="c"] pk_nl, int size_kvec, int size_zvec):
    c_interpolation_init_all(<double*> z.data, <double*> H.data, <double*> conf_dist.data,  <double*> D.data,  <double*> f.data, &size_z,
                                 <double*> kvec.data, <double*> zvec.data,  <double*> pk.data, <double*> pk_nl.data, &size_kvec, &size_zvec)
    return

def interpolation_free_all():
    c_interpolation_free_all()
    return

def set_stuff(double CPr, double H0, double h, double omdm, double omb, double omk, double K):
    cdef:
        ar[double] pars = np.array([CPr,H0,h,omdm,omb,omk,K])
    c_set_stuff(<double*> pars.data)
    return

def set_mp_overlap(ar[double,ndim=2] ConMat, int size_ConMat, ar[int,ndim=1] intParams, ar[double,ndim=1] realParams, int which_sample):
    c_set_mp_overlap(<double*> ConMat.data, &size_ConMat, <int*> intParams.data, <double*> realParams.data, &which_sample)
    return

def set_this(ar[double,ndim=2] lenslowz, ar[double,ndim=2] lens2dfloz, ar[double,ndim=2] lenscmass, ar[double,ndim=2] lens2dfhiz,
                 ar[double,ndim=1] xipm, ar[double,ndim=2] invcovxipm, int sizcov,
                 ar[double,ndim=1] ellarr, ar[double,ndim=1] prefacarrz, int size_ellarr,
                 ar[double,ndim=1] ellini, ar[double,ndim=1] maskelements, int size_maskelements,
                 ar[double,ndim=1] bes0arr,  ar[double,ndim=1] bes4arr,  ar[double,ndim=1] bes2arr,
                 ar[int,ndim=1] intParams, ar[int,ndim=1] logParams):
    
    c_set_this(<double *> lenslowz.data, <double *> lens2dfloz.data, <double *> lenscmass.data, <double *> lens2dfhiz.data,
                   <double *> xipm.data, <double *> invcovxipm.data, &sizcov,
                   <double *> ellarr.data, <double *> prefacarrz.data, &size_ellarr,
                   <double *> ellini.data, <double *> maskelements.data, &size_maskelements,
                   <double *> bes0arr.data, <double *> bes4arr.data, <double *> bes2arr.data,
                   <int *> intParams.data, <int*> logParams.data)
    return

def set_sources(double[:, :, :] source_data):
    cdef:
        int dim1 = source_data.shape[0]
        int dim2 = source_data.shape[1]
        int dim3 = source_data.shape[2]
    c_set_sources(&source_data[0,0,0], &dim3, &dim2, &dim1)
    #ar[double,ndim=2] af1, ar[double,ndim=2] af2, ar[double,ndim=2] af3, ar[double,ndim=2] af4):
    #c_set_sources(<double *> af1.data, <double *> af2.data, <double *> af3.data, <double *> af4.data)
    return

def test_array(double[:,:,:] source_data):
    cdef:
        int dim1 = source_data.shape[0]
        int dim2 = source_data.shape[1]
        int dim3 = source_data.shape[2]
    print(dim1, dim2, dim3)
    cdef double *data_ptr = &source_data[0,0,0]
    c_test_array(&source_data[0,0,0], &dim3, &dim2, &dim1)

def test_array2(ar[double,ndim=3] source_data):
    c_test_array2(&source_data[0,0,0])
    return

def loglkl_from_fortran(sigma_v_cmass=-1,b1_cmass=-1,N_shot_cmass=-1,sigv_lowz=-1, b1_lowz=-1,N_shot_lowz=-1,b1_2dfloz=-1,
                            b1_2dfhiz=-1,sigv_2dfloz=-1,sigv_2dfhiz=-1,N_shot_2dfloz=-1,N_shot_2dfhiz=-1,ampia=-1,
                            redzia=-1,lumia=-1,a1phot=-1,a2phot=-1,a3phot=-1,a4phot=-1):
    cdef:
        double loglkl
        ar[double,ndim=1] DataParams = np.array([sigma_v_cmass,b1_cmass,N_shot_cmass,sigv_lowz, b1_lowz,N_shot_lowz,b1_2dfloz,
                            b1_2dfhiz,sigv_2dfloz,sigv_2dfhiz,N_shot_2dfloz,N_shot_2dfhiz,ampia,
                            redzia,lumia,a1phot,a2phot,a3phot,a4phot])
    c_cosmolss_lnlike(<double *> DataParams.data, &loglkl)
    return loglkl
