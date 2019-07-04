from setuptools import setup, find_packages
from biopython.convert.__version import __versionstr__

with open('README.rst') as readme:
    setup(
        name='BioPython-Convert',
        version=__versionstr__,
        packages=find_packages(),
        long_description=readme.read(),
        install_requires=['biopython', 'gffutils', 'jmespath'],
        scripts=['feature_merge/feature_merge.py'],
        url='https://github.com/brinkmanlab/biopython-convert',
        license='MIT',
        author='Nolan Woods',
        author_email='nolan_w@sfu.ca',
        description='Interconvert various file formats supported by biopython. Supports querying records with JMESPath.',
        include_package_data=True
    )
