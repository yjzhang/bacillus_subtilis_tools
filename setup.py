from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install

from setuptools import Command
# https://stackoverflow.com/questions/20288711/post-install-script-with-python-setuptools

import os
import subprocess

PATH = os.path.dirname(__file__)
DB_PATH = os.path.join(PATH, 'subtiwiki', 'data')
DB_FILES = ['subtiwiki.db']
GZ_DB_FILES = ['subtiwiki.db.gz']


class DevelopWithGunzip(develop):

    def run(self):
        print('unzipping db files')
        for f in GZ_DB_FILES:
            subprocess.run('gunzip {0}'.format(os.path.join(DB_PATH, f)), shell=True)
        develop.run(self)


class InstallWithGunzip(install):

    def run(self):
        print('unzipping db files')
        for f in GZ_DB_FILES:
            subprocess.run('gunzip {0}'.format(os.path.join(DB_PATH, f)), shell=True)
        install.run(self)

class GunzipDBFiles(Command):
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass

    def run(self):
        for f in GZ_DB_FILES:
            subprocess.run('gunzip {0}'.format(os.path.join(DB_PATH, f)), shell=True)

class GzipDBFiles(Command):
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass

    def run(self):
        for f in DB_FILES:
            subprocess.run('gzip {0}'.format(os.path.join(DB_PATH, f)), shell=True)



setup(
    name='bacillus_subtilis_tools',
    version='0.0.1',
    author='Yue Zhang',
    author_email='yjzhang@cs.washington.edu',
    url='https://github.com/yjzhang/bacillus_subtilis_tools',
    license='MIT',
    install_requires=['backports.functools_lru_cache'],
    packages=find_packages("."),
    cmdclass={
        'develop': DevelopWithGunzip,
        'install': InstallWithGunzip,
        'unzip': GunzipDBFiles,
        'zip': GzipDBFiles,
    },
    # this is for including the data dir in the package.
    zip_safe=False,
    package_data={'cellmesh': ['data/subtiwiki.db']},
    include_package_data=True,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
    ],

)
