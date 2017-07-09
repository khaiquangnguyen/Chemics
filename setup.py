# from distutils.core import setup
# from Cython.Build import cythonize
#
# setup(
#     ext_modules = cythonize("FastDpCalculator.pyx")
# )


from setuptools import setup
from setuptools import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext


ext_modules=[ Extension("FastDpCalculator",
              ["FastDpCalculator.pyx"],
              extra_compile_args = ["-ffast-math"])]

setup(
  name = "FastDpCalculator",
  cmdclass = {"build_ext": build_ext},
  ext_modules = ext_modules)
