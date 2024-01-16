"""
Setup file for package 'GASTLI'.
"""
from setuptools import find_packages
from numpy.distutils.core import Extension, setup
import os



extra_compile_args = []


mod_dimensions = Extension(\
    name='GASTLI.dimensions',\
    sources=['GASTLI/mod_dimensions.f95'],\
    f2py_options=['only:', 'subroutine_name', ':'])


extensions = [mod_dimensions]


def setup_function():

    setup(
        name='GASTLI',
        version='1.0',
        description='GAS gianT modeL for Interiors',
        long_description=open(os.path.join(
          os.path.dirname(__file__), 'README.rst')).read(),
        author='Lorena Acuña',
        author_email='acuna@mpia.de',
        packages=find_packages(),
        include_package_data=True,
        install_requires=[
            'numpy',
            'pandas',
            'scipy',
            'h5py',
            'matplotlib'
        ],
        zip_safe=False,
        ext_modules=extensions,
    )


setup_function()
