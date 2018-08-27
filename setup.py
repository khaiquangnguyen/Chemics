# from distutils.core import setup
# from Cython.Build import cythonize
#
# setup(
#     ext_modules = cythonize("FastDpCalculator.pyx")
# )

#
from setuptools import setup
from setuptools import Extension
from Cython.Distutils import build_ext
import numpy


ext_modules=[ Extension("FastDpCalculator",
              ["FastDpCalculator.pyx"],
              extra_compile_args = ["-ffast-math"],
              include_dirs = [numpy.get_include()])]

setup(
    name = "FastDpCalculator",
    cmdclass = {"build_ext": build_ext},
    ext_modules = ext_modules,
    include_dirs = [numpy.get_include()],
    install_requires=['matplotlib', 'PySide', 'scipy', 'peakutils', 'numpy',
                                                          'pandas']
)

