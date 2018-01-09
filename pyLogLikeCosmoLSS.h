extern void c_interpolation_init_all(double *z, double *H, double *conf_dist, double *D,  double *f, int *size_z, double *kvec, double *zvec, double *pk, double *pk_nl, int *size_kvec, int*size_zvec);
extern void c_set_stuff(double *CosmoParams);
extern void c_set_mp_overlap(double *ConMat, int *size_ConMat, int *intParams, double *realParams, int *which_sample);
extern void c_set_this(double *af1, double *af2, double *af3, double *af4, double *lenslowz, double *lens2dfloz, double *lenscmass, double *lens2dfhiz, double *xipm, double *invcovxipm, int *sizcov, double *ellini, double *maskelements, int *size_maskelements, double *bes0arr,double *bes4arr,double *bes2arr, int *intParams, int* logParams);
extern void c_cosmolss_lnlike(double *dataparams, double *loglkl);

