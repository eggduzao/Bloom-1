
import numpy
from setuptools import setup
from Cython.Build import cythonize


setup(
    name="c",
    packages=["c"],
    ext_modules=cythonize("c/cython_algorithm.pyx", compiler_directives={"language_level": "3"}),
    include_dirs=[numpy.get_include()]
)

# Compilation
# python cython_algorithm.py build_ext --inplace

