from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system

ext_modules = [Extension(# module name:
                         'pyLogLikeCosmoLSS',
                         # source file:
                         ['pyLogLikeCosmoLSS.pyx'],
                         # other compile args for gcc
                         include_dirs=[get_include()],
                         extra_compile_args=['-fPIC', '-O3', '-lgfortran','-lquadmath'],
                        #library_dirs = [".","/usr/local/Cellar/gcc/6.3.0_1/lib/gcc/6"],
                        library_dirs = ["."],
                         # other files to link to
                         extra_link_args=['libCosmoLSS.a',
                                              #'-L/usr/local/Cellar/gcc/6.3.0_1/lib/gcc/6',
                                              '-lgomp','-lgfortran','-lquadmath',
                                              ])]

setup(name = 'pyLogLikeCosmoLSS',
      cmdclass = {'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = ext_modules)
