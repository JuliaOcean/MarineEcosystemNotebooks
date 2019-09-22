# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.1
#   kernelspec:
#     display_name: Julia 1.1.0
#     language: julia
#     name: julia-1.1
# ---

# This julia notebook uses [Simons' CMAP](https://cmap.readthedocs.io/en/latest/) to depict output of the Darwin model in context of the [CBIOMES](https://cbiomes.org) project.
#
# <img src="../figs/cbiomes-01.png" alt="Drawing" style="height: 50px;"/>
#
# **Pre-requisites:**
#
# - 1) install the [PyCmap python package and its dependencies](https://github.com/simonscmap/pycmap) using `pip`.
# - 2) compile [PyCall.jl](https://github.com/simonscmap/pycmap) using external python distribution that installed `PyCmap`.
#
# _Note:_ the commented out `pip install pycmap` call does not need to be repeated every time. However, if the `pyimport` call later returns an error then please try uncommenting the `pip install pycmap` call and start over.

# +
#run(`pip install pycmap`) #pycmap is used via PyCall later

run(pipeline(`which python`,"whichpython.txt")) #external python path
ENV["PYTHON"]=readline("whichpython.txt")
import Pkg; Pkg.add("PyCall"); Pkg.build("PyCall")
# -

# **Now we can import `pycmap` inside `julia` using `PyCall.jl`.**

using PyCall
PyCmap = pyimport("pycmap")

# #### Get the cmap data catalog
#
# This will create a `data/` folder should contain `catalog.csv`. Another folder named `embed/` gets created below.

getCatalog = pyimport("opedia.getCatalog")

# #### Import opedia functionalities

opediareg = pyimport("opedia.plotRegional")
opediasub = pyimport("opedia.subset")

# #### Download a regional Chlorophyll data set from the Darwin model
#
# This will download a file called `data/RM_tblDarwin_Nutrient_3day_FeT_darwin_3day.csv`.

# +
tables = ["tblDarwin_Nutrient_3day"] # see catalog.csv  for the complete list of tables and variable names
variables = ["FeT_darwin_3day"] # see catalog.csv  for the complete list of tables and variable names   
startDate = "1994-01-03"
endDate = "1994-01-03"
lat1, lat2 = 10, 70
lon1, lon2 = -180, -80
depth1, depth2 = 0, 10
fname = "regional"
exportDataFlag = true # false if you you do not want to download data

opediareg.regionalMap(tables, variables, startDate, endDate, lat1, lat2, lon1, lon2, depth1, depth2, fname, exportDataFlag)
# -

# #### Load the same data directly as a dataframe

df = opediasub.spaceTime(tables[1], variables[1], startDate, endDate, lat1, lat2, lon1, lon2, depth1, depth2)


