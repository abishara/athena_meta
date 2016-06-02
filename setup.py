# setup.py

from setuptools import setup, find_packages

setup(
  name = "athena",
  version = "0.1",
  packages = filter(lambda x: x.find('athena') == 0, find_packages()),
  py_modules = ['main'],
  entry_points = {
    'console_scripts': [
      'main = main:main',
    ]
  },

  install_requires=[
  ],
  author='Alex Bishara',
  description='athena assembler',
)
