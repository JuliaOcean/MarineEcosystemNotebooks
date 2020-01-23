# MarineEcosystemNotebooks

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/gaelforget/Cbiomes2019Notebooks/master)
[![DOI](https://zenodo.org/badge/185446209.svg)](https://zenodo.org/badge/latestdoi/185446209)

[Jupyter](https://jupyter.org) notebooks related to marine ecosystem models and observations in [Julia](https://julialang.org).

1. how ocean colour data and model ouptut can be used jointly in [CBIOMES](https://https://github.com/CBIOMES)
2. how the [Simons' CMAP](https://cmap.readthedocs.io/en/latest/) database can be queried from within `julia`.

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
- `05. OceanColourClassifiers.ipynb` 2D array demo of [OC-CCI](https://www.oceancolour.org) / [Jackson et al 2017](http://doi.org/10.1016/j.rse.2017.03.036) classifier.

## Data Bases

- `03. DarwinCmapJulia.ipynb` uses [Simons' CMAP](https://cmap.readthedocs.io/en/latest/) to download output of the [CBIOMES model output](https://github.com/gaelforget/CBIOMES) as a file or dataframe using `julia`'s `pycall` package to call the `pycmap` python package.
- `04. DarwinCmapPython.ipynb` uses [Simons' CMAP](https://cmap.readthedocs.io/en/latest/) to map output of the Darwin model using the `pycmap` python package (this is a **python rather than julia** notebook).


