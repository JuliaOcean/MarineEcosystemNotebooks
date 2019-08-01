# Cbiomes2019Notebooks

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/gaelforget/Cbiomes2019Notebooks/master)

These notebooks illustrate how ocean colour data and model ouptut can be used jointly in context of the [CBIOMES](https://cbiomes.org) project using the [Julia](https://julialang.org) language.

## `Notebooks/`

- `OceanColourAlgorithms.ipynb` provides simple recipes to compare model output and ocean color data.
- `ModelReflectanceMap.ipynb` uses `Plots.jl` to map output of the Darwin model and OC-CCI data.
- `DarwinCmapJulia.ipynb` uses [Simons' CMAP](https://cmap.readthedocs.io/en/latest/) to download output of the Darwin model as a file or dataframe using julia's pycall package to call the opedia python package.
- `DarwinCmapRegional.ipynb` uses [Simons' CMAP](https://cmap.readthedocs.io/en/latest/) to map output of the Darwin model using the opedia python package (this is a **python rather than julia** notebook).

## `src/`

Julia utility scripts.

## Source

![CBIOMES logo](https://raw.githubusercontent.com/gaelforget/Cbiomes2019Notebooks/master/figs/cbiomes-01.png)
