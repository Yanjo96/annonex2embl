import os
import glob
import unittest
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def my_test_suite():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('tests', pattern='*_test.py')
    return test_suite

setup(
    name='annonex2embl',
    version='0.4.6',
    description='Converts annotated DNA sequence alignments in NEXUS format to ENA submission files',
    long_description=read('README.md'),
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'License :: OSI Approved :: GNU General Public License (GPLv3)',
        'Programming Language :: Python :: 2.7',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ],
    keywords='DNA sequence submission to ENA',
    url='https://github.com/michaelgruenstaeudl/annonex2embl',
    author='Michael Gruenstaeudl',
    author_email='m.gruenstaeudl@fu-berlin.de',
    license='GPLv3',
    packages=['annonex2embl'], # So that the subfolder 'annonex2embl' is read immediately.
    #packages = find_packages(),
    install_requires=['biopython', 'unidecode', 'termcolor'],
    scripts=glob.glob('scripts/*'),
    test_suite='setup.my_test_suite',
    include_package_data=True,
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'annonex2embl = annonex2embl.CLIOps:start_annonex2embl'
        ],
    },
)
