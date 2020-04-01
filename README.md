# MarineEcosystemNotebooks

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/gaelforget/Cbiomes2019Notebooks/master)
[![DOI](https://zenodo.org/badge/185446209.svg)](https://zenodo.org/badge/latestdoi/185446209)

[Jupyter](https://jupyter.org) / [Julia](https://julialang.org) notebooks related to marine ecosystem models and observations.

1. how ocean colour data and model ouptut can be used jointly in [CBIOMES](https://https://github.com/CBIOMES)
2. how model output and data available online are easily accessed in `julia`

<p align="center">
  <img width="300" src="https://raw.githubusercontent.com/gaelforget/Cbiomes2019Notebooks/master/figs/cbiomes-01.png">
</p>

_Notes:_

- _Each `.ipynb` notebook is paired with a `.jl` file via `jupytext`_
- _An interactive version can readily be spun up via the `launch binder` badge_
- _Please use the [repository issue tracker](https://guides.github.com/features/issues/) for queries, bug reports, new contributions, etc._
- _The `src/` folder contains helper functions & scripts._

## Ocean Color

- `01. OceanColourAlgorithms.ipynb` provides simple recipes to compare [CBIOMES model output](https://github.com/gaelforget/CBIOMES) and [ocean color data](https://www.oceancolour.org).
- `02. ModelReflectanceMap.ipynb` uses `Plots.jl` to map out [CBIOMES model output](https://github.com/gaelforget/CBIOMES) and [ocean color data](https://www.oceancolour.org).
- `03. Classifications.ipynb` 2D array demo of [OC-CCI](https://www.oceancolour.org) / [Jackson et al 2017](http://doi.org/10.1016/j.rse.2017.03.036) classifier.

## Accessing Data

- `01. DarwinModelOutput.ipynb` from either (1) the [MIT-CBIOMES opendap](http://engaging-opendap.mit.edu:8080/las/) server or (2) the [Simons CMAP](https://cmap.readthedocs.io/en/latest/) data base to access model output from the [CBIOMES](https://cbiomes.org) project.
- `02. GradientsCruiseData.ipynb` uses [Simons' CMAP](https://cmap.readthedocs.io/en/latest/) to download [SCOPE-Gradients](http://scope.soest.hawaii.edu/data/gradients/data/) which is then read and plotted in `julia` using the `Plots.jl` package.


