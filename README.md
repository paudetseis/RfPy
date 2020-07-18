
![](./rfpy/examples/picture/RfPy_logo.png)

## Teleseismic receiver function calculation and post-processing 

RfPy is a software to calculate single event-station receiver functions from the spectral deconvolution technique. Methods are available to post-process the receiver function data to calculate H-k stacks, back-azimuth harmonics and common-conversion-point (CCP) imaging. The code uses the ``StDb`` package for querying and building a station database and can be used through command-line scripts.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3905414.svg)](https://doi.org/10.5281/zenodo.3905414)
[![Build Status](https://travis-ci.com/paudetseis/RfPy.svg?branch=master)](https://travis-ci.com/paudetseis/RfPy)
[![codecov](https://codecov.io/gh/paudetseis/RfPy/branch/master/graph/badge.svg)](https://codecov.io/gh/paudetseis/RfPy)

Installation, Usage, API documentation and scripts are described at 
https://paudetseis.github.io/RfPy/.

Authors: [`Pascal Audet`](https://www.uogeophysics.com/authors/admin/) (Developer and Maintainer) and [`Jeremy Gosselin`](https://www.uogeophysics.com/authors/gosselin/) (Contributor)
<!-- #### Citing

If you use `SplitPy` in your work, please cite the 
[`Zenodo DOI`](https://zenodo.org/badge/latestdoi/211722700).
 -->
#### Contributing

All constructive contributions are welcome, e.g. bug reports, discussions or suggestions for new features. You can either [open an issue on GitHub](https://github.com/paudetseis/RfPy/issues) or make a pull request with your proposed changes. Before making a pull request, check if there is a corresponding issue opened and reference it in the pull request. If there isn't one, it is recommended to open one with your rationale for the change. New functionality or significant changes to the code that alter its behavior should come with corresponding tests and documentation. If you are new to contributing, you can open a work-in-progress pull request and have it iteratively reviewed.

Examples of straightforward contributions include editing the documentation or adding notebooks that describe published examples of teleseismic receiver functions. Suggestions for improvements (speed, accuracy, flexibility, etc.) are also welcome.

#### References

- Audet, P. (2010) Temporal Variations in Crustal Scattering Structure near Parkfield, California, Using Receiver Functions, Bulletin of the Seismological Society of America (2010) 100 (3): 1356-1362. https://doi.org/10.1785/0120090299

- Audet, P. (2015) Layered crustal anisotropy around the San Andreas Fault near Parkfield, California, J. Geophys. Res. Solid Earth, 120, 3527-3543, https://doi.org/10.1002/2014JB011821

- Cossette, E., Audet, P., and Schneider, D.A. (2016) Structure and anisotropy of the crust in the Cyclades, Greece, using receiver functions constrained by in situ rock textural data, J. Geophys. Res. Solid Earth, 121, 2661-2678, https://doi.org/10.1002/2015JB012460

- Tarayoun, A., P. Audet, S. Mazzotti, and A. Ashoori (2017) Architecture of the crust and uppermost mantle in the northern Canadian Cordillera from receiver functions, J. Geophys. Res. Solid Earth, 122, 5268–5287, https://doi.org/10.1002/2017JB014284.

#### Use cases

- Audet, P., Schutt, D., Schaeffer, A.J., Estève, C., Aster, R., and Cubley, J. (2020). Moho variations across the northern Canadian Cordillera, Seism. Res. Lett., accepted.