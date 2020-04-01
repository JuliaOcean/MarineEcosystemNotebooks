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

# # Download & Plot _Gradients_ Cruise Data
#
# This julia notebook downloads and plots [SCOPE-Gradients](http://scope.soest.hawaii.edu/data/gradients/data/) cruise data from the [Simons' CMAP](https://cmap.readthedocs.io/en/latest/) data base.
#
# <img src="../figs/cbiomes-01.png" alt="Drawing" style="height: 200px;"/>

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# _Pre-requisites:_
#
# - 1. _install the [PyCmap](https://github.com/simonscmap/pycmap) python package and its dependencies using `pip`_
# - 2. _compile [PyCall.jl](https://github.com/simonscmap/pycmap) using external python distribution that installed `PyCmap`_
# - 3. _obtain your own API key from the [SimonsCMAP website](https://simonscmap.com) (free, just need email, takes `<30s`)_
# -

if false
    run(`pip install pycmap`) #pycmap is used via PyCall later
    run(pipeline(`which python`,"whichpython.txt")) #external python path
    ENV["PYTHON"]=readline("whichpython.txt")
    import Pkg; Pkg.build("PyCall")
end

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ### Import _pycmap_ Into _julia_
#
# [PyCmap](https://github.com/simonscmap/pycmap) is the `Python` API that we will use in `julia` (to query the `CMAP` data base) via the [PyCall.jl](https://github.com/JuliaPy/PyCall.jl) package. 
#
# _You may need to replace `your-own-API-key` (as outline below) with your own API key from the [SimonsCMAP website](https://simonscmap.com) and uncomment the command below._
# -

using PyCall
PyCmap = pyimport("pycmap")
#cmap = PyCmap.API(token="your-own-API-key")
cmap = PyCmap.API()

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ### Get Data Catalog
#
# _See [SimonsCMAP website](https://simonscmap.com) for more information. The commented `df.to_csv` call below creates `catalog.csv` and writes the data from `df` to this file._

# + {"slideshow": {"slide_type": "-"}}
df = cmap.get_catalog()
#df.to_csv("catalog.csv")

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ### Download Data
#
# The lists provided below contain `CMAP table` names associated with the [SCOPE-Gradients](http://scope.soest.hawaii.edu/data/gradients/data/) cruise data. 
#
# _Uncomment one code bloc at a time to download more files._

# + {"slideshow": {"slide_type": "subslide"}}
pth="../samples/gradients/"
!isdir("$pth") ? mkdir("$pth") : nothing

list0=[]
if !isfile("$pth"*"tblKM1906_Gradients3_uway_optics.csv")
    list0=["tblKM1906_Gradients3","tblKM1906_Gradients3_uway_optics","tblKM1906_Gradients3_uwayCTD","tblKM1906_Gradients3_uw_tsg",
        "tblMGL1704_Gradients2_CTD","tblMGL1704_Gradients2_uway_optics","tblKOK1606_Gradients1_CTD","tblKOK1606_Gradients1_uway_optics"]
end

#list0=["tblMGL1704_Gradients2_Nutrients","tblMGL1704_Gradients2_Diazotroph",
#    "tblMGL1704_Gradients2_TargetedMetabolites","tblMGL1704_Gradients2_Trace_Metals"]

#list0=["tblKOK1606_Gradients1_Nutrients","tblKOK1606_Gradients1_Dissolved_Gasses",
#    "tblKOK1606_Gradients1_TargetedMetabolites","tblKOK1606_Gradients1_Diazotroph"]

for i in list0
    df=cmap.get_dataset(i)
    df.to_csv("$pth$i.csv")
end

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ### Read And Plot Data From File
#
# As a use case example below we read the `LISST` data collected during the `Gradients 2` cruise and plot a subset of the data.

# + {"slideshow": {"slide_type": "fragment"}}
using CSV, DataFrames
df = CSV.File("$pth"*"tblKM1906_Gradients3_uway_optics.csv") |> DataFrame! ;
# -

cmap.get_dataset_metadata("tblKM1906_Gradients3_uway_optics")
cmap.get_metadata("tblKM1906_Gradients3_uway_optics","LISST_small")
cmap.get_unit("tblKM1906_Gradients3_uway_optics", "LISST_small")


# + {"slideshow": {"slide_type": "subslide"}}
using Plots

tmp=df[1:10:end,:]
scatter(tmp[!,:lat],tmp[!,:LISST_small],marker = 2,label="1.25-2.0 micron",
    xlabel="°N",ylabel="umol C/L",title="Particles (LISST C) in Gradients 3")
scatter!(tmp[!,:lat],tmp[!,:LISST_medium],marker = 1,label="2.0-20 micron")
scatter!(tmp[!,:lat],tmp[!,:LISST_large],marker = 1,label="20-100.0 micron")
# -

