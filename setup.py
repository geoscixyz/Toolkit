#!/usr/bin/env python
from __future__ import print_function
"""Geosci Toolkit

Geophysical Toolkit for Geologists
http://toolkit.geosci.xyz.
"""

from distutils.core import setup
from setuptools import find_packages


CLASSIFIERS = [
    'Development Status :: 1 - Beta',
    'Intended Audience :: Education',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Mathematics',
    'Topic :: Scientific/Engineering :: Physics',
    'Operating System :: Microsoft :: Windows',
    'Operating System :: POSIX',
    'Operating System :: Unix',
    'Operating System :: MacOS',
    'Natural Language :: English',
]

with open('README.md') as f:
    LONG_DESCRIPTION = ''.join(f.readlines())

setup(
    name='GeoToolkit',
    version='0.0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'tqdm',
        'matplotlib',
        'geopandas',
        'scikit-image',
    ],
    author='GeoSci Developers',
    author_email='fourndo@gmail.com',
    description='Geophysical toolkit for geologists',
    long_description=LONG_DESCRIPTION,
    keywords='geophysics, geologists',
    url='http://toolkit.geosci.xyz',
    download_url='https://github.com/geoscixyz/Toolkit',
    classifiers=CLASSIFIERS,
    platforms=['Windows', 'Linux', 'Solaris', 'Mac OS-X', 'Unix'],
    license='MIT License',
    use_2to3=False,
)
