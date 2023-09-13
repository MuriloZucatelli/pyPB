# Python support for CFD

<[![Join the chat at https://gitter.im/robertsawko/pyfd](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/robertsawko/pyfd?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)>

<[![Build Status](https://travis-ci.org/robertsawko/pyfd.svg?branch=master)](https://travis-ci.org/robertsawko/pyfd)>


.. image:: https://app.travis-ci.com/anuga-community/anuga_core.svg?branch=main
    :target: https://app.travis-ci.com/anuga-community/anuga_core
    :alt: travis ci status
   
.. image:: https://ci.appveyor.com/api/projects/status/x5airjv7eq2u805w/branch/main?svg=true
    :target: https://ci.appveyor.com/project/stoiver/anuga-core-nwgr0
    :alt: appveyor status

.. image:: https://img.shields.io/pypi/v/anuga.svg
    :target: https://pypi.python.org/pypi/anuga/
    :alt: Latest PyPi Version

.. image:: https://img.shields.io/pypi/dm/anuga.svg
    :target: https://pypistats.org/packages/anuga
    :alt: PyPi download statistics

.. image:: https://img.shields.io/conda/vn/conda-forge/anuga.svg
    :target: https://anaconda.org/conda-forge/anuga
    :alt: Latest Conda Version
 
.. image:: https://img.shields.io/conda/dn/conda-forge/anuga.svg
    :target: https://anaconda.org/conda-forge/anuga
    :alt: Conda Forge download statistics

.. image:: https://readthedocs.org/projects/anuga/badge/?version=latest
    :target: https://anuga.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
    
This is a collection of Python modules for routine CFD calculation 

 * Colebrook-White friction factor check with single phase pipe
 * Calculation of turbulence intensities for inlet BC specification
 * Wheeler and PD moment inversion algorithms
 * Method of classes method

## General

### Pressure drop

Simple pressure drop calculation based on Colebrook-White. Tested with efunda calculator.

## Turbulence

### Turbulent BC

## Population balance modelling
### Methods of moments
#### Wheeler's moment inversion

Also known as Modified Chebyshev algorithm (Gautschi 2004). Obtains a set of weights and abcissas from the set of moments.

#### Realizability check
Checks whether Hadamard matrices have non-negative determinants


