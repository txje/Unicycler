'''
Run 'python3 setup.py install' to install Unicycler.
'''

import sys
import os
import shutil
from distutils.command.build import build
from distutils.core import Command
import subprocess
import multiprocessing
import glob
import fnmatch
import imp

# Install setuptools if not already present
try:
    imp.find_module('setuptools')
except ImportError:
    import ez_setup
    ez_setup.use_setuptools()
from setuptools import setup
from setuptools.command.install import install

# Get the program version from another file.
exec(open('unicycler/version.py').read())

with open('README.md', 'rb') as readme:
    LONG_DESCR = readme.read().decode('utf-8')


class UnicycleBuild(build):
    '''
    The build process runs the Makefile to build the C++ functions into a shared library.
    '''
    def run(self):
        build.run(self) # Run original build code
        try:
            cmd = ['make', '-j', str(max(8, multiprocessing.cpu_count()))]
        except NotImplementedError:
            cmd = ['make']
        def compile():
            subprocess.call(cmd)
        self.execute(compile, [], 'Compiling Unicycler: ' + ' '.join(cmd))


class UnicycleInstall(install):
    '''
    The install process copies the C++ shared library to the install location.
    '''
    def run(self):
        install.run(self) # Run original install code
        shutil.copyfile('unicycler/cpp_functions.so',
                        os.path.join(self.install_lib, 'unicycler', 'cpp_functions.so'))


class UnicycleClean(Command):
    '''
    Custom clean command that really cleans up everything, except for:
      - the compiled *.so file needed when running the programs
      - setuptools-*.egg file needed when running this script
    '''
    user_options = []
    def initialize_options(self):
        self.cwd = None
    def finalize_options(self):
        self.cwd = os.getcwd()
    def run(self):
        assert os.getcwd() == self.cwd, 'Must be in package root: %s' % self.cwd

        delete_directories = []
        for root, dirnames, filenames in os.walk(self.cwd):
            for dirname in fnmatch.filter(dirnames, '*.egg-info'):
                delete_directories.append(os.path.join(root, dirname))
            for dirname in fnmatch.filter(dirnames, 'build'):
                delete_directories.append(os.path.join(root, dirname))
            for dirname in fnmatch.filter(dirnames, '__pycache__'):
                delete_directories.append(os.path.join(root, dirname))
        for delete_directory in delete_directories:
            print('Deleting directory:', delete_directory)
            shutil.rmtree(delete_directory)

        delete_files = []
        for root, dirnames, filenames in os.walk(self.cwd):
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
      long_description=LONG_DESCR,
      url='http://github.com/rrwick/unicycler',
      author='Ryan Wick',
      author_email='rrwick@gmail.com',
      license='GPL',
      packages=['unicycler'],
      entry_points={"console_scripts": ['unicycler = unicycler.unicycler:main',
                                        'antpath = unicycler.antpath:main',
                                        'scrutinate = unicycler.scrutinate:main']},
      zip_safe=False,
      cmdclass={'build': UnicycleBuild,
                'install': UnicycleInstall,
                'clean': UnicycleClean}
)
