from setuptools import setup, find_packages

setup(
    name='genomics',
    version='0.1.0',
    description='Collection of data structures and '
                'utility functions for genomics',
#   url='NA',
    author='Bulak Arpat',
    author_email='bulak.arpat@gmail.com',
    license='GPLv3',
    packages=find_packages(),
    entry_points={
        'console_scripts': ['region = genomics.region:main']
    },
    install_requires=['numpy'],
    dependency_links=['http://github.com/gatfieldlab/pypackages#egg=extended&subdirectory=extended'],
    zip_safe=False
)
