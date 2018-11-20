# Function Library

[ ![Codeship Status for FunctionLab/function](https://app.codeship.com/projects/793cc3c0-47e5-0135-b732-66f05eb232be/status?branch=master)](https://app.codeship.com/projects/231720)

flib is short for Functionlab Library. These are a set of assorted scripts and python classes that have been useful. These are provided without support in case someone else finds them useful.
These scripts should run in python2 and python3, just use different requirements files to set up environment. 

See the [wiki](https://github.com/FunctionLab/function/wiki) for code recipes!

## Installation with pip

        $ pip install git+https://github.com/FunctionLab/function.git@master

## Development setup

Create a virtualenv and activate it:

        $ virtualenv fenv
        $ source fenv/bin/activate

### Clone from GitHub

        $ git clone git@github.com:FunctionLab/function.git

### Install dependencies python2

        $ cd function
        $ pip install -r requirements.txt


### Install dependencies python3

        $ cd function
        $ pip install -r requirements3.txt


## Unit tests
### Unit test for python2

        $ python -m unittest discover flib/tests/
### Unit test for python3        
To run all the test units:

        $ python -m unittest discover flib/tests/
        
To run all in one file:

        $ python -m unittest flib.tests.test_goa.TestOBO
        
To run a specific test e.g. `test_get_value`

        $ python -m unittest flib.tests.test_dab.TestDab.test_get_value
        
## Running a single script

You can run a single script e.g. `obo.py`
Remember to export directory where flib directory resides if there is problem with finding the module

    $ export PYTHONPATH=$HOME/dev/workspace/function
    $ cd $HOME/dev/workspace/function/flib/core
    $ python <script_name>
    

# Contributing to flib

flib operates on a pull request model. If you have code that you would like to contribute, please file a pull request. We aim for pep8 compliance in future code. If you touch old code, we'd appreciate it if you could make it pep8 compliant as well.
