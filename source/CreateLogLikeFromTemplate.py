
# coding: utf-8

# In[1]:

import numpy as np
import shutil


# In[2]:

outpath = '/Users/au192734/Projects/montepython_CosmoLSS/montepython/likelihoods/CosmoLSS/'
inpath = '/Users/au192734/Downloads/cosmolss_onlyupdatedfiles/'


# In[ ]:


    


# In[15]:

with open(inpath+'source/CosmoLSS.f90') as fid:
    all_lines = fid.readlines()
    
comments = [r'TCosmoTheoryPredictions',r'obj%theoryomp = Theory',r'obj%exact_z_index = this%exact_z_index',
           r'allocate(obj%theoryomp)',r'allocate(obj%exact_z_index',r'class(CosmoLSSLikelihood)']
substitutions = {r'Theory%MPK%ny':r'NonlinearPowerspectrum%ny',
                 r'Theory%MPK%y':r'NonlinearPowerspectrum%y',
                 r'Theory%MPK%PowerAt':r'MPKPowerAt',
                 r'FIRSTPRIVATE(ellgen,obj)':r'PRIVATE(ellgen) FIRSTPRIVATE(obj)',
                 r'FIRSTPRIVATE(ttt,obj)':r'PRIVATE(ttt) FIRSTPRIVATE(obj)',
                 r'FIRSTPRIVATE(wttt,obj)':r'PRIVATE(wttt) FIRSTPRIVATE(obj)',
                 r'FIRSTPRIVATE(ellgen,obj)':r'PRIVATE(ellgen) FIRSTPRIVATE(obj)',
                 r'FIRSTPRIVATE(obj,ellgen)':r'PRIVATE(ellgen) FIRSTPRIVATE(obj)',
                 r'CALL Matrix_MulVec(this%':r'P0P2P4_Conv = matmul(theory_vec, this%',
                 r'(:,:),theory_vec(:),P0P2P4_Conv(:))':')',
                 r'LinearPower_2D_s_cmass(Theory,':r'LinearPower_2D_s_cmass(',
                r'Theory%growth_z%Value(1.0d0/reda-1.0d0)/Theory%sigma8_z%Value(1.0d0/reda-1.0d0)':r'GrowthRate%Value(reda)',
                r'obj%theoryomp%NL_MPK%':r'NL_MPK',
                r'Theory%growth_z%Value(z)/Theory%sigma8_z%Value(z)':r'GrowthRate%Value(1d0/(1d0+z))',
                r'class(CMBParams)':r'type(myCMBparam)',
                 r'end module CosmoLSS':r'end module LogLikeCosmoLSS_module',
                 r'function CosmoLSS_LnLike(this,CMB,Theory,DataParams)':r'function CosmoLSS_LnLike(DataParams)',
                 r'== .':r'.eqv. .',
                 r'CMB%h0/100':r'CMB%h0/100d0',
                 r'real(mcp) :: DataParams(:)':'real(mcp), intent(in) :: DataParams(19) \n real(mcp), parameter :: c = 2.99792458d8',
                }

firstCMBparam = True
for i, line in enumerate(all_lines):
    if firstCMBparam and r'class(CMBParams)' in line:
        all_lines[i] = r'!'+line
        firstCMBparam = False
        continue
    
    iscomment = False
    for substring in comments:
        if substring in line:
            iscomment = True
            break
    if iscomment:
        all_lines[i] = r'!'+line
    else:
        for key, val in substitutions.iteritems():
            line = line.replace(key,val)
        all_lines[i] = line


# In[16]:

shutil.copy(outpath+'source/LogLikeTemplate.f90',outpath+'source/LogLikeCosmoLSS.f90')
with open(outpath+'source/LogLikeCosmoLSS.f90','a') as fid:
    fid.writelines(all_lines[701:])


# In[ ]:



