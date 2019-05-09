import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='qc3C',
    description='Hi-C quality control analysis',
    long_description=long_description,
    version='0.2',
    author='Matthew Z DeMaere',
    author_email='matt.demaere@gmail.com',
    platforms='Linux-86_x64',
    packages=setuptools.find_packages(),
    url='https://github.com/cerebis/qc3C',
    license='GNU Affero General Public License v3',

    install_requires=[
        'numpy',
        'pandas',
        'pysam',
        'tqdm',
        'biopython',
        'recordclass',
    ],

    classifiers=[
        'Programming Language :: Python :: 3.6',
        'License :: OSI Approved :: GNU Affero General Public License v3',
        'Operating System :: POSIX :: Linux',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Development Status :: 5 - Alpha'
    ],

    entry_points={
        'console_scripts': ['qc3C=qc3C.command_line:main'],
    }
)
