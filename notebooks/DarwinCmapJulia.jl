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

# This julia notebook uses [Simons' CMAP](https://cmap.readthedocs.io/en/latest/) to extract Darwin model output from the [CBIOMES](https://cbiomes.org) project.
#
# <img src="../figs/cbiomes-01.png" alt="Drawing" style="height: 100px;"/>
#
# **Pre-requisites:**
#
# - 1) install the [PyCmap](https://github.com/simonscmap/pycmap) python package and its dependencies using `pip`.
# - 2) compile [PyCall.jl](https://github.com/simonscmap/pycmap) using external python distribution that installed `PyCmap`.
#
# **Note:** 
#
# The commented out `pip install pycmap` call does not need to be repeated every time. However, if the `pyimport` call later returns an error then please try uncommenting the `pip install pycmap` call and start over.

# +
#run(`pip install pycmap`) #pycmap is used via PyCall later

run(pipeline(`which python`,"whichpython.txt")) #external python path
ENV["PYTHON"]=readline("whichpython.txt")
import Pkg; Pkg.add("PyCall"); Pkg.build("PyCall")
# -

# **Now we can import `pycmap` inside `julia` using `PyCall.jl`.**

using PyCall
PyCmap = pyimport("pycmap")
cmap = PyCmap.API(token="929df5e0-d501-11e9-9d93-a1d831a9e2b2")

# #### Get the cmap data catalog
#
# This `df.to_csv` call creates the `catalog.csv` with the content of `df`.

df = cmap.get_catalog()
df.to_csv("catalog.csv")

# #### Download a regional Chlorophyll data set from the Darwin model (1/2)
#
# The `plot_map` function can in principle export data to file when `exportDataFlag=true`.

# +
tables = ["tblDarwin_Nutrient"] # see catalog.csv  for the complete list of tables and variable names
variables = ["FeT"] # see catalog.csv  for the complete list of tables and variable names   
startDate = "1994-01-03"
endDate = "1994-01-03"
lat1, lat2 = 10, 70
lon1, lon2 = -180, -80
depth1, depth2 = 0, 10
fname = "regional"
exportDataFlag = false # set to true if you want to download data

viz = pyimport("pycmap.viz")
viz.plot_map(tables, variables, startDate, endDate, 
    lat1, lat2, lon1, lon2, depth1, depth2, fname, exportDataFlag)
# -

# #### Download a regional Chlorophyll data set from the Darwin model (2/2)
#
# The `space_time` function returns the same data to memory as a dataframe.

df = cmap.space_time(tables[1], variables[1], startDate, endDate, 
    lat1, lat2, lon1, lon2, depth1, depth2)
