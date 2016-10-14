#!/usr/bin/env python3
"""
Run 'python3 setup.py install' to install Unicycler.
"""

# Make sure this is being run with Python 3.4 or later.
import sys
if sys.version_info.major != 3 or sys.version_info.minor < 4:
    print('Error: you must execute setup.py using Python 3.4 or later')
    sys.exit(1)

import os
import shutil
from distutils.command.build import build
from distutils.core import Command
import subprocess
import multiprocessing
import fnmatch
import importlib.util

# Install setuptools if not already present.
if not importlib.util.find_spec("setuptools"):
    import ez_setup
    ez_setup.use_setuptools()

from setuptools import setup
from setuptools.command.install import install

# Get the program version from another file.
exec(open('unicycler/version.py').read())

with open('README.md', 'rb') as readme:
    LONG_DESCRIPTION = readme.read().decode()


class UnicycleBuild(build):
    """
    The build process runs the Makefile to build the C++ functions into a shared library.
    """

    def run(self):
        build.run(self)  # Run original build code

        clean_cmd = ['make', 'clean']
        try:
            make_cmd = ['make', '-j', str(min(8, multiprocessing.cpu_count()))]
        except NotImplementedError:
            make_cmd = ['make']

        def clean_cpp():
            subprocess.call(clean_cmd)

        def compile_cpp():
            subprocess.call(make_cmd)

        self.execute(clean_cpp, [], 'Cleaning previous compilation: ' + ' '.join(clean_cmd))
        self.execute(compile_cpp, [], 'Compiling Unicycler: ' + ' '.join(make_cmd))


class UnicycleInstall(install):
    """
    The install process copies the C++ shared library to the install location.
    """

    def run(self):
        install.run(self)  # Run original install code
        shutil.copyfile(os.path.join('unicycler', 'cpp_functions.so'),
                        os.path.join(self.install_lib, 'unicycler', 'cpp_functions.so'))
        gene_data_dir = os.path.join(self.install_lib, 'unicycler', 'gene_data')
        if not os.path.exists(gene_data_dir):
            os.makedirs(gene_data_dir)
        shutil.copyfile(os.path.join('unicycler', 'gene_data', 'start_genes.fasta'),
                        os.path.join(gene_data_dir, 'start_genes.fasta'))
        shutil.copyfile(os.path.join('unicycler', 'gene_data', 'lambda_phage.fasta'),
                        os.path.join(gene_data_dir, 'lambda_phage.fasta'))


class UnicycleClean(Command):
    """
    Custom clean command that really cleans up everything, except for:
      - the compiled *.so file needed when running the programs
      - setuptools-*.egg file needed when running this script
    """
    user_options = []

    def initialize_options(self):
        self.cwd = None

    def finalize_options(self):
        self.cwd = os.getcwd()

    def run(self):
        assert os.getcwd() == self.cwd, 'Must be in package root: %s' % self.cwd

        delete_directories = []
        for root, dir_names, filenames in os.walk(self.cwd):
            for dir_name in fnmatch.filter(dir_names, '*.egg-info'):
                delete_directories.append(os.path.join(root, dir_name))
            for dir_name in fnmatch.filter(dir_names, 'build'):
                delete_directories.append(os.path.join(root, dir_name))
            for dir_name in fnmatch.filter(dir_names, '__pycache__'):
                delete_directories.append(os.path.join(root, dir_name))
        for delete_directory in delete_directories:
            print('Deleting directory:', delete_directory)
            shutil.rmtree(delete_directory)

        delete_files = []
        for root, dir_names, filenames in os.walk(self.cwd):
            for filename in fnmatch.filter(filenames, 'setuptools*.zip'):
                delete_files.append(os.path.join(root, filename))
            for filename in fnmatch.filter(filenames, '*.o'):
                delete_files.append(os.path.join(root, filename))
            for filename in fnmatch.filter(filenames, '*.pyc'):
                delete_files.append(os.path.join(root, filename))
        for delete_file in delete_files:
            print('Deleting file:', delete_file)
            os.remove(delete_file)


setup(name='unicycler',
      version=__version__,
      description='bacterial genome assembler for hybrid read sets',
      long_description=LONG_DESCRIPTION,
      url='http://github.com/rrwick/unicycler',
      author='Ryan Wick',
      author_email='rrwick@gmail.com',
      license='GPL',
      packages=['unicycler'],
      entry_points={"console_scripts": ['unicycler = unicycler.unicycler:main',
                                        'unicycler_align = unicycler.unicycler_align:main',
                                        'unicycler_check = unicycler.unicycler_check:main',
                                        'pacbio_hybrid_polish = unicycler.pacbio_hybrid_polish:main']},
      zip_safe=False,
      cmdclass={'build': UnicycleBuild,
                'install': UnicycleInstall,
                'clean': UnicycleClean}
      )
