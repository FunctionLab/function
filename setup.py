#!/usr/bin/env python

from setuptools import setup, find_packages

reqs = []
with open('requirements.txt') as f:
    for l in f.read().splitlines():
        if not l.strip().startswith('#'):
            reqs.append(l.strip())

# Keyword order: https://packaging.python.org/distributing
setup(
    # These 9 fields are inserted into PKG-INFO. Unspecified keys are set to
    # UNKNOWN values.
    name='flib',
    version='0.1.0',
    description="Function Lab Python Library",
  # long_description="",
    url='https://github.com/FunctionLab/function',
    author='Function Lab',
  # author_email='',
  # license='',
  # platform=''

    packages=find_packages(),
    install_requires=reqs,
)
