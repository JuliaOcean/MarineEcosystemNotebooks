# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.4
#   kernelspec:
#     display_name: Julia 1.3.1
#     language: julia
#     name: julia-1.3
# ---

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# # CBIOMES Global Model Output
#
# Here we retrieve model output from the [CBIOMES](https://cbiomes.org) project via two methods:
#
# - 1. in array format from the [MIT-CBIOMES opendap](http://engaging-opendap.mit.edu:8080/las/) (e.g. visit [this page](http://engaging-opendap.mit.edu:8080/las/UI.vm#panelHeaderHidden=false;differences=false;autoContour=false;xCATID=3C6AA795DF3E9F4E1208CEFE8341F298;xDSID=id-ab2a4e0c65;varid=MXLDEPTH-id-cdfa319965;imageSize=auto;over=xy;compute=Nonetoken;tlo=15-Jan-1992%2000:00;thi=15-Jan-1992%2000:00;catid=3C6AA795DF3E9F4E1208CEFE8341F298;dsid=id-ab2a4e0c65;varid=MXLDEPTH-id-cdfa319965;avarcount=0;xlo=-180;xhi=180;ylo=-90;yhi=90;operation_id=Plot_2D_XY_zoom;view=xy))
# - 2. in tabular format from the [Simons CMAP data base](https://cmap.readthedocs.io/en/latest/) (go to [this page](https://cmap.readthedocs.io/en/latest/catalog/datasets/Darwin_clim.html#darwin-clim))
#
#
# <img src="../figs/cbiomes-01.png" alt="Drawing" style="height: 200px;"/>

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## Using Opendap
#  
# The `NCDatasets.jl` package readily supports lazy access to `opendap` data sets. The function below for example access Iron concentration for a chosen depth level and month. Only this two-dimensional slice is transferred over the network.
#
# _warning: this method has failed on mybinder in the past_
#
# ### ➥ `Array` / gridded map
#
#

# + {"slideshow": {"slide_type": "subslide"}}
using NCDatasets, Plots

function test_opendapp(k,t)
    srv="http://engaging-opendap.mit.edu:8080/thredds/dodsC/las/"
    fil="id-2e0ea5ca2c/data_usr_local_tomcat_content_cbiomes_20200206_17_Nutrients_FeT.0001.nc.jnl"
    test=Dataset(srv*fil)
    lon=test["LON_C"]
    lat=test["LAT_C"]
    tmp=test["FeT"][:,:,k,t]
    tmp[findall(tmp.<=-0.99e34)].=NaN
    println( test["FeT"].attrib["long_name"]*" in "*test["FeT"].attrib["units"] )    
    return vec(lon)[1:2:end],vec(lat)[1:2:end],transpose(tmp[1:2:end,1:2:end])
end

# + {"slideshow": {"slide_type": "subslide"}}
heatmap(test_opendapp(1,1))

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## Using pycmap
#
# The `space_time` API call below returns the same data as before but as a python table object which is easily saved to a `csv` file.
#
# ### ➥ `Table` / data cloud

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# _Pre-requisites:_
#
# - 1. _install the [PyCmap](https://github.com/simonscmap/pycmap) python package and its dependencies using `pip`_
# - 2. _compile [PyCall.jl](https://github.com/simonscmap/pycmap) using external python distribution that installed `PyCmap`_
# - 3. _obtain your own API key from the [SimonsCMAP website](https://simonscmap.com) (free; takes `<30s`)_
# - 4. _import `pycmap` via `pycall`_

# + {"slideshow": {"slide_type": "skip"}}
if false
    run(`pip install pycmap`) #pycmap is used via PyCall later
    run(pipeline(`which python`,"whichpython.txt")) #external python path
    ENV["PYTHON"]=readline("whichpython.txt")
    import Pkg; Pkg.build("PyCall")
end

# + {"slideshow": {"slide_type": "skip"}, "cell_type": "markdown"}
# _Import `pycmap` Into julia:_
#
# [PyCmap](https://github.com/simonscmap/pycmap) is the `Python` API that we will use in `julia` (to query the `CMAP` data base) via the [PyCall.jl](https://github.com/JuliaPy/PyCall.jl) package. 
#
# _You may need to replace `your-own-API-key` (as outline below) with your own API key from the [SimonsCMAP website](https://simonscmap.com) and uncomment the command below._

# + {"slideshow": {"slide_type": "skip"}}
using PyCall
PyCmap = pyimport("pycmap")
#cmap = PyCmap.API(token="your-own-API-key")
cmap = PyCmap.API()

# + {"slideshow": {"slide_type": "skip"}, "cell_type": "markdown"}
# ### Get Data Catalog
#
# _See [SimonsCMAP website](https://simonscmap.com) for more information. The commented `df.to_csv` call below creates `catalog.csv` and writes the data from `df` to this file._

# + {"slideshow": {"slide_type": "skip"}}
df = cmap.get_catalog();
#df.to_csv("catalog.csv")

# + {"slideshow": {"slide_type": "subslide"}}
tables = ["tblDarwin_Nutrient"] # see catalog.csv  for the complete list of tables and variable names
variables = ["FeT"] # see catalog.csv  for the complete list of tables and variable names   
startDate = "1994-01-03"
endDate = "1994-01-03"
lat1, lat2 = 10, 70
lon1, lon2 = -180, -80
depth1, depth2 = 0, 10

df = cmap.space_time(tables[1], variables[1], startDate, endDate, 
    lat1, lat2, lon1, lon2, depth1, depth2)

pth="../samples/gradients/"
!isdir("$pth") ? mkdir("$pth") : nothing
df.to_csv("$pth"*"FeT.csv")

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# The data can then easily be read, manipulated, and plotted using Julia's `DataFrames.jl` package.
# -

using CSV, DataFrames
df = CSV.File("$pth"*"FeT.csv") |> DataFrame!

# + {"slideshow": {"slide_type": "subslide"}}
jan=findall( (df.time .== df.time[1]).&( (!ismissing).(df.FeT)) )
jan=rand(jan,1000)
scatter(df.lon[jan],df.lat[jan],zcolor=df.FeT[jan],title="random selection")
# -


