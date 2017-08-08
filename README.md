# Function Library

[ ![Codeship Status for FunctionLab/function](https://app.codeship.com/projects/793cc3c0-47e5-0135-b732-66f05eb232be/status?branch=master)](https://app.codeship.com/projects/231720)

flib is short for Functionlab Library. These are a set of assorted scripts and python classes that have been useful. These are provided without support in case someone else finds them useful.

See the [wiki](https://github.com/FunctionLab/function/wiki) for code recipes!

## Installation with pip

        $ pip install git+https://github.com/FunctionLab/function.git@master

## Development setup

Create a virtualenv and activate it:

        $ virtualenv fenv
        $ source fenv/bin/activate

### Clone from GitHub

        $ git clone git@github.com:FunctionLab/function.git

### Install dependencies

        $ cd function
        $ pip install -r requirements.txt


## Unit tests

        $ python -m unittest discover flib/tests/

# Contributing to flib

flib operates on a pull request model. If you have code that you would like to contribute, please file a pull request. We aim for pep8 compliance in future code. If you touch old code, we'd appreciate it if you could make it pep8 compliant as well.
