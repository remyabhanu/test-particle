import os
import glob
# Need to use the enhanced version of distutils packaged with
# numpy so that we can compile fortran extensions
from setuptools.command import install as _install
from setuptools import find_packages
from numpy.distutils.core import Extension, setup
from numpy.distutils import exec_command

# Output debugging information while installing
os.environ['DISTUTILS_DEBUG'] = "1"

#############################################################################
# First, check to make sure we are executing
# 'python setup.py install' from the same directory
# as setup.py (root davitpy directory)
#############################################################################
path = os.getcwd()
assert('setup.py' in os.listdir(path)), \
       "You must execute 'python setup.py install' from within the \
davitpy root directory."

#############################################################################
# define a read function for using README.md for long_description
#############################################################################


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


#############################################################################
# Now we must define all of our C and Fortran extensions
#############################################################################
# Fortran extensions
#############################################################################

tsyg = Extension('tsygFort',
                 sources=['tsyganenko/T96.f',
                          'tsyganenko/TS04c.f',
                          'tsyganenko/geopack08.for',
                          'tsyganenko/geopack08.pyf'])

weim = Extension('weimFort',
                 sources=['weimer/w2k.f',
                          'weimer/w2k.pyf'])

#############################################################################
# Now execute the setup
#############################################################################
setup(name='testpartpy',
      version="0.8",
      description="Space Science Toolkit",
      author="VT SuperDARN Lab and friends",
      author_email="german.farinas@gmail.com",
      url="",
      download_url="",
      packages=find_packages(),
      long_description=read('README.md'),
      zip_safe=False,
      ext_modules=[tsyg,weim],
      install_requires=[],
      classifiers=[
            "Development Status :: 7 - Beta",
            "Topic :: Scientific/Engineering",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU General Public License (GPL)",
            "Natural Language :: English",
            "Programming Language :: Python"
            ],
      )

if os.environ['DISTUTILS_DEBUG'] == "1":
    print('Sources', find_packages())
