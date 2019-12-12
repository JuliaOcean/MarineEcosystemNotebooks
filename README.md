# Cbiomes2019Notebooks

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/gaelforget/Cbiomes2019Notebooks/master)
[![DOI](https://zenodo.org/badge/185446209.svg)](https://zenodo.org/badge/latestdoi/185446209)

The included notebooks illustrate:

- 1) how ocean colour data and model ouptut can be used jointly in [CBIOMES](https://cbiomes.org) using the [Julia](https://julialang.org) language
- 2) how [Simons' CMAP](https://cmap.readthedocs.io/en/latest/) can be queried from within `julia` or `python` via [pycmap](https://github.com/simonscmap/pycmap)

You can start an interactive version of these notebooks using the above `launch binder` badge. Each Julia notebook is paired with a *light script* `.jl file` that can be run directly in julia.

<img src="https://raw.githubusercontent.com/gaelforget/Cbiomes2019Notebooks/master/figs/cbiomes-01.png" alt="Drawing" style="height: 100px;"/>

## `Notebooks/`

- `01. OceanColourAlgorithms.ipynb` provides simple recipes to compare [CBIOMES model output](https://github.com/gaelforget/CBIOMES) and [ocean color data](https://www.oceancolour.org).
- `02. ModelReflectanceMap.ipynb` uses `Plots.jl` to map out [CBIOMES model output](https://github.com/gaelforget/CBIOMES) and [ocean color data](https://www.oceancolour.org).
- `03. DarwinCmapJulia.ipynb` uses [Simons' CMAP](https://cmap.readthedocs.io/en/latest/) to download output of the [CBIOMES model output](https://github.com/gaelforget/CBIOMES) as a file or dataframe using `julia`'s `pycall` package to call the `pycmap` python package.
- `04. DarwinCmapPython.ipynb` uses [Simons' CMAP](https://cmap.readthedocs.io/en/latest/) to map output of the Darwin model using the `pycmap` python package (this is a **python rather than julia** notebook).
- `05. OceanColourClassifiers.ipynb` 2D array demo of [OC-CCI](https://www.oceancolour.org) / [Jackson et al 2017](http://doi.org/10.1016/j.rse.2017.03.036) classifier.

## `src/`

Julia utility scripts.

