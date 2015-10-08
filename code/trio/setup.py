#Rick Dewey 8.14.12
#Last modified 8.14.12 Rick Dewey
#setup utility for cython extension module hmmUtils
#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

setup(cmdclass = {'build_ext':build_ext}, name = "hmmUtils", ext_modules = [Extension("hmmUtils", ["hmmUtils.pyx"])])
