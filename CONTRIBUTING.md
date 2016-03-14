# Contributing

LingPy development follows the model and workflow described as the 
["fork-branch-and-pull-request style"](https://gun.io/blog/how-to-github-fork-branch-and-pull-request/).
In order to keep the code transparent even for multiple contributors, we have set up a list 
of [coding conventions](https://github.com/lingpy/lingpy/blob/master/CONVENTIONS.md).


## Getting started

To install LingPy to hack on it, fork the repository, open a terminal and type:
```bash
$ git clone https://github.com/<your-github-user>/lingpy/
$ cd lingpy
$ python setup.py develop
```
This will install LingPy in ["develpment mode"](http://pythonhosted.org//setuptools/setuptools.html#development-mode),
i.e. you will be able edit the sources in the cloned repository and import the altered code just as the regular python package.

LingPy development makes use of several additional python libraries to help with testing and
quality assurance:

- `mock`
- `nose` as test aggregator and runner
- `coverage` for test coverage measurement and reporting
- `flake8` for coding style checks
- `tox` to run tests across different supported python versions.


## Running tests

To run our tests we use the [nose library](https://nose.readthedocs.org/en/latest/),
just enter the `lingpy` directory and call `nosetests` or `nosetests.exe` on the command line. 

To run tests on all supported python versions, you need to run `tox`, which will read the
configuration in `tox.ini`. At the moment this requires commands `python2.7` and `python3.4`
being available in your system, which will be the case on Ubuntu 14.04, for example.

*Please make sure all tests pass without failure or error when run with tox before 
submitting a pull request.*


## Creating the Documentation

Currently, documentation is created using the following steps:

* whenever code is added to LingPy, the contributors add documentation inline in their code, following the style also used in other projects such as SciPy, NumPy, and Networkx.
* the general website structure is added around the code, you can find its content by browsing the lingpy/doc/sources/ directory.
* before compiling the code, we create a full reference containing links to all code in LingPy, using the sphinx-apidoc command:

  ```
  $ cd lingpy/ # got to relevenant main folder of lingpy
  $ sphinx-apidoc -o doc/sources/reference/ lingpy
  ```
  
  This creates specific documentation for all the LingPy package, structured with respect to the modules in LingPy.
* run the sphinx build as described on the sphinx website at http://sphinx-doc.org, e.g. on Linux or Mac OSX running

  ```
  $ cd doc
  $ make html
  ```