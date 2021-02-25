import re
import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

version_str = None
VERSION_FILE = "qc3C/_version.py"
with open(VERSION_FILE, "rt") as vh:
    for _line in vh:
        mo = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", _line, re.M)
        if mo:
            version_str = mo.group(1)
            break

if version_str is None:
    raise RuntimeError("Unable to find version string in {}".format(VERSION_FILE))

setuptools.setup(
    name='qc3C',
    description='Hi-C quality control analysis',
    long_description=long_description,
    version=version_str,
    author='Matthew Z DeMaere',
    author_email='matt.demaere@gmail.com',
    platforms='Linux-86_x64',
    packages=setuptools.find_packages(),
    url='https://github.com/cerebis/qc3C',
    license='GNU Affero General Public License v3',

    install_requires=[
        'astropy<4.1',
        'biopython==1.76',
        'leven',
        'json2html',
        'numba',
        'numpy',
        'pandas',
        'psutil',
        'pysam',
        'pytest',
        'recordclass',
        'scipy',
        'simplejson',
        'tqdm'
    ],

    classifiers=[
        'Programming Language :: Python :: 3.7',
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
