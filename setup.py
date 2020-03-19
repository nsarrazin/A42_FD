from setuptools import setup, find_packages
from distutils.util import convert_path
import sys


# ver_path = convert_path('sim_common/version.py')
# with open(ver_path) as ver_file:
#     ns = {}
#     exec(ver_file.read(), ns)
#     version = ns['version']

setup(
    name='b42fd',
    version='0.1.0',
    description='Package for B42 FD',
    author='B42',
    author_email='sarrazin.nathan@gmail.com',
    # url='',

    install_requires=['numpy',
                      'numba',
                      'scipy',
                      'pandas',
                      'pytest',
                      'matplotlib'],

    packages=find_packages('.', exclude=["test"]),
    classifiers=[
        'Programming Language :: Python :: 3 :: Only',
        'Development Status :: 2 - Pre-Alpha']
)