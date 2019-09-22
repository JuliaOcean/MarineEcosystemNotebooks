# Cbiomes2019Notebooks

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/gaelforget/Cbiomes2019Notebooks/master)
[![DOI](https://zenodo.org/badge/185446209.svg)](https://zenodo.org/badge/latestdoi/185446209)

The included notebooks illustrate:

- 1) how ocean colour data and model ouptut can be used jointly in [CBIOMES](https://cbiomes.org) using the [Julia](https://julialang.org) language
- 2) how [Simons' CMAP](https://cmap.readthedocs.io/en/latest/) can be queried from within `julia` or `python` via [pycmap](https://github.com/simonscmap/pycmap)

<img src="https://raw.githubusercontent.com/gaelforget/Cbiomes2019Notebooks/master/figs/cbiomes-01.png" alt="Drawing" style="height: 100px;"/>

## `Notebooks/`

- `01. OceanColourAlgorithms.ipynb` provides simple recipes to compare model output and ocean color data.
- `02. ModelReflectanceMap.ipynb` uses `Plots.jl` to map output of the Darwin model and OC-CCI data.
- `03. DarwinCmapJulia.ipynb` uses [Simons' CMAP](https://cmap.readthedocs.io/en/latest/) to download output of the `Darwin` model as a file or dataframe using `julia`'s `pycall` package to call the `pycmap` python package.
- `04. DarwinCmapPython.ipynb` uses [Simons' CMAP](https://cmap.readthedocs.io/en/latest/) to map output of the Darwin model using the `pycmap` python package (this is a **python rather than julia** notebook).

## `src/`

Julia utility scripts.

