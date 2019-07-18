#!/usr/bin/env python
from setuptools import setup, find_packages

import re
with open('biopython_convert/__version.py') as version_file:
    __versionstr__ = '.'.join(re.search(r"__version__ = \[(\d+), (\d+), (\d+)\]", version_file.read(), re.M).group()[1:])

with open('README.rst') as readme:
    setup(
        name='biopython.convert',
        version=__versionstr__,
        python_requires='>=3.7',
        packages=find_packages(),
        long_description=readme.read(),
        install_requires=['biopython>=1.73', 'gffutils>=0.9', 'jmespath>=0.9.4'],
        scripts=['bin/biopython.convert'],
        url='https://github.com/brinkmanlab/biopython-convert',
        license='MIT Amended',
        author='Nolan Woods',
        author_email='nolan_w@sfu.ca',
        description='Interconvert various file formats supported by biopython. Supports querying records with JMESPath.',
        include_package_data=True,
        test_suite="tests",
        classifiers=[
            "Environment :: Console",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Operating System :: POSIX",
            "Programming Language :: Python :: 3.7",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
        ]
    )
