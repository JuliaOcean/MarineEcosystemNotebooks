# Cbiomes2019Notebooks

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/gaelforget/Cbiomes2019Notebooks/master)
[![DOI](https://zenodo.org/badge/185446209.svg)](https://zenodo.org/badge/latestdoi/185446209)

The `Notebooks/` folder illustrate:

1. how ocean colour data and model ouptut can be used jointly in [CBIOMES](https://cbiomes.org) 
2. how [Simons' CMAP](https://cmap.readthedocs.io/en/latest/) can be queried from within `julia` or `python` via [pycmap](https://github.com/simonscmap/pycmap)

<img src="https://raw.githubusercontent.com/gaelforget/Cbiomes2019Notebooks/master/figs/cbiomes-01.png" alt="Drawing" height="50"/>

_Notes:_

- To start an interactive version of the notebooks, please click on the `launch binder` badge.
- Included notebooks are written in the [Julia](https://julialang.org) language.
- Notebooks are paired with corresponding `.jl file` that run directly in `julia`.
- The `src/` folder contains `julia` utility scripts.

## Ocean color

- `01. OceanColourAlgorithms.ipynb` provides simple recipes to compare [CBIOMES model output](https://github.com/gaelforget/CBIOMES) and [ocean color data](https://www.oceancolour.org).
- `02. ModelReflectanceMap.ipynb` uses `Plots.jl` to map out [CBIOMES model output](https://github.com/gaelforget/CBIOMES) and [ocean color data](https://www.oceancolour.org).
- `05. OceanColourClassifiers.ipynb` 2D array demo of [OC-CCI](https://www.oceancolour.org) / [Jackson et al 2017](http://doi.org/10.1016/j.rse.2017.03.036) classifier.

## Query CMAP

- `03. DarwinCmapJulia.ipynb` uses [Simons' CMAP](https://cmap.readthedocs.io/en/latest/) to download output of the [CBIOMES model output](https://github.com/gaelforget/CBIOMES) as a file or dataframe using `julia`'s `pycall` package to call the `pycmap` python package.
- `04. DarwinCmapPython.ipynb` uses [Simons' CMAP](https://cmap.readthedocs.io/en/latest/) to map output of the Darwin model using the `pycmap` python package (this is a **python rather than julia** notebook).


