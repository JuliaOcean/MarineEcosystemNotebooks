# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Query Simons' `CMAP` from `python`
#
# **This python notebook** uses [Simons' CMAP](https://cmap.readthedocs.io/en/latest/) to extract Darwin model output from the [CBIOMES](https://cbiomes.org) project.
#
# <img src="../figs/cbiomes-01.png" alt="Drawing" style="height: 50px;"/>
#
# **Pre-requisites**
#
# - install the [PyCmap](https://github.com/simonscmap/pycmap) python package and its dependencies using `pip`.

import pycmap

# **Generate regional map of Chlorophyll in the Darwin model**

# +
from pycmap import viz

tables = ['tblDarwin_Nutrient']    # see catalog.csv  for the complete list of tables and variable names
variables = ['FeT']                            # see catalog.csv  for the complete list of tables and variable names   
startDate = '1994-01-03'
endDate = '1994-01-03'
lat1, lat2 = 10, 70
lon1, lon2 = -180, -80
depth1, depth2 = 0, 10
exportDataFlag = False
showMapFlag = True

viz.plot_map(tables, variables, startDate, endDate, 
             lat1, lat2, lon1, lon2, depth1, depth2, exportDataFlag, showMapFlag)
# -


