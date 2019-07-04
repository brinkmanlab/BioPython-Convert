from setuptools import setup, find_packages
from biopython_convert.__version import __versionstr__

with open('README.rst') as readme:
    setup(
        name='biopython.convert',
        version=__versionstr__,
        packages=find_packages(),
        long_description=readme.read(),
        install_requires=['biopython', 'gffutils', 'jmespath'],
        scripts=['bin/biopython.convert'],
        url='https://github.com/brinkmanlab/biopython-convert',
        license='MIT Amended',
        author='Nolan Woods',
        author_email='nolan_w@sfu.ca',
        description='Interconvert various file formats supported by biopython. Supports querying records with JMESPath.',
        include_package_data=True
    )
