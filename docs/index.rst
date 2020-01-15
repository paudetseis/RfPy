
.. figure:: ../rfpy/examples/picture/RfPy_logo.png
   :align: center

Documentation
=============

``RfPy`` is a package containing Python tools for calculating teleseismic
receiver functions. Methods are available to plot and post-process 
receiver function data. Post processing includes: H-k method,
Harmonic decomposition and CCP stacking. The code uses 
the ``StDb`` package for querying and building a station database 
used in command-line scripts.

.. image:: https://travis-ci.com/paudetseis/RfPy.svg?branch=master
    :target: https://travis-ci.com/paudetseis/RfPy
.. image:: https://codecov.io/gh/paudetseis/RfPy/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/paudetseis/RfPy

.. note::

   ``RfPy`` was written independently of the Python software `rf <https://rf.readthedocs.io>`_,
   although it is possible that some classes and methods defined here might be applied
   to ``rf`` objects, since both are heavily based on `obspy <http://www.obspy.org>`_. The main
   differences between ``RfPy`` and ``rf`` are as follows:

   * ``RfPy`` only calculates `P` receiver functions, whereas ``rf`` can also calculate
     `S` receiver functions. 
   * ``RfPy`` uses a Wiener spectral deconvolution technique that minimizes noise based on
     observed seismic noise on all components, whereas ``rf`` uses either a water-level 
     spectral deconvolution or a time-domain deconvolution.
   * ``RfData`` objects are used to calculate single-station and single-event receiver 
     functions, whereas ``rf`` can handle multiple stations at once. 

.. toctree::
   :maxdepth: 1
   :caption: Quick Links

   links

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   rfpy

.. toctree::
   :maxdepth: 2
   :caption: API Documentation

   api

.. toctree::
   :maxdepth: 2
   :caption: Scripts & Tutorials

   scripts