# setup.py

from setuptools import setup, find_packages

setup(
  name = "athena",
  version = "1.1",
  packages = filter(lambda x: x.find('athena') == 0, find_packages()),
  py_modules = ['main'],
  entry_points = {
    'console_scripts': [ 'athena-meta = main:main' ]
  },
  install_requires=[
    'bx-python>=0.7.3',
    'ipython-cluster-helper>=0.5.2',
    'numpy>=1.11.0',
    'pysam>=0.9.0',
  ],
  author='Alex Bishara',
  description='athena assembler',
)
