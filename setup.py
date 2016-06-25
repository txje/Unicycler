'''
Run 'python3 setup.py install' to install Unicycler.
'''

import os
import shutil
from distutils.command.build import build
import subprocess
from multiprocessing import cpu_count
from setuptools import setup
from setuptools.command.install import install

exec(open('unicycler/version.py').read())

with open('README.md', 'rb') as readme:
    LONG_DESCR = readme.read().decode('utf-8')

# BASEPATH = os.path.dirname(os.path.abspath(__file__))
# UNICYCLER_PATH = os.path.join(BASEPATH, 'unicycler')

class UnicycleBuild(build):
    '''
    The build process runs the Makefile to build the C++ functions into a shared library.
    '''
    def run(self):
        build.run(self) # Run original build code
        try:
            cmd = ['make', '-j', str(max(8, cpu_count()))]
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
                'install': UnicycleInstall}
)
