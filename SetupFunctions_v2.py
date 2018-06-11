from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

ext_modules=[ Extension("CalculateLLReadDataFunctions_v2",
              ["CalculateLLReadDataFunctions_v2.pyx"],
              libraries=["m"],
              extra_compile_args = ["-ffast-math"])]

setup(
  name = 'Hello world app',
  ext_modules = cythonize("CalculateLLReadDataFunctions_v2.pyx"),
)
