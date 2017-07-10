# Function Library

[ ![Codeship Status for FunctionLab/function](https://app.codeship.com/projects/793cc3c0-47e5-0135-b732-66f05eb232be/status?branch=master)](https://app.codeship.com/projects/231720)

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
