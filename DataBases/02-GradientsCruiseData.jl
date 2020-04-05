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

# # _Gradients_ Cruise Data
#
# Here we retrieve [SCOPE-Gradients](http://scope.soest.hawaii.edu/data/gradients/data/) cruise data from the [Simons' CMAP](https://cmap.readthedocs.io/en/latest/) data base.
#
# <img src="../figs/cbiomes-01.png" alt="Drawing" style="height: 100px;"/>

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# _Pre-requisites:_
#
# - 1. _install the [PyCmap](https://github.com/simonscmap/pycmap) python package and its dependencies using `pip`_
# - 2. _compile [PyCall.jl](https://github.com/simonscmap/pycmap) using external python distribution that installed `PyCmap`_
# - 3. _obtain your own API key from [Simons' CMAP](https://simonscmap.com) (free; takes `<30s`)_
# - 4. _import `pycmap` (via `pycall`), `Plots`, and `helper functions`_
# -

if false
    run(`pip install pycmap`) #pycmap is used via PyCall later
    run(pipeline(`which python`,"whichpython.txt")) #external python path
    ENV["PYTHON"]=readline("whichpython.txt")
    import Pkg; Pkg.build("PyCall")
end

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ### Import _pycmap_ , _Plots_ , and  _helper functions_
#
# [PyCmap](https://github.com/simonscmap/pycmap) is the `Python` API that we will use in `julia`, via the [PyCall.jl](https://github.com/JuliaPy/PyCall.jl) package, to query the `CMAP` data base. `Plots.jl` is one of `julia`'s plotting packages.
#
# _You may need to replace `your-own-API-key` (as outline below) with your own API key from [Simons' CMAP](https://simonscmap.com) and uncomment the command below._

# +
using PyCall
PyCmap = pyimport("pycmap")
#cmap = PyCmap.API(token="your-own-API-key")
cmap = PyCmap.API()

using Plots
include("helper_functions.jl")

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ### Get Data Catalog
#
# _The commented `df.to_csv` command writes the content of `df` to a new `catalog.csv` file. Alternatively, `Pandas.jl` can be used as also shown._

# + {"slideshow": {"slide_type": "-"}}
df = cmap.get_catalog()
#df.to_csv("catalog.csv")

#df=Pandas.DataFrame(cmap.get_catalog())
#to_csv(df,"catalog.csv")

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ### Download & Reload Data Set
#
# The lists provided by `gradients_list()` contain `CMAP table` names associated with the [SCOPE-Gradients](http://scope.soest.hawaii.edu/data/gradients/data/) cruise data. 

# + {"slideshow": {"slide_type": "-"}}
if false
    pth="../samples/gradients/"
    !isdir("$pth") ? mkdir("$pth") : nothing
    
    list0=gradients_list("main")
    for i in list0
        df=cmap.get_dataset(i)
        df.to_csv("$pth$i.csv")
    end
    
    using CSV, DataFrames
    df = CSV.File("$pth"*"tblKM1906_Gradients3_uway_optics.csv") |> DataFrame!
end

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ### Read Data & Meta-Data
#
# As an example below we read the `LISST` data collected during the `Gradients 3` cruise and then plot a subset of the data.
# -

s=cmap_get("tblKM1906_Gradients3_uway_optics","LISST_small")
m=cmap_get("tblKM1906_Gradients3_uway_optics","LISST_medium")
l=cmap_get("tblKM1906_Gradients3_uway_optics","LISST_large")

# + {"slideshow": {"slide_type": "subslide"}}
t=1:10:length(s["lat"])
Plots.plot(s["lat"][t],s["val"][t],marker = 2,label=s["Long_Name"],
    xlabel="°N",ylabel=s["Unit"], title="Particles (LISST C) in Gradients 3")
Plots.plot!(m["lat"][t],m["val"][t],marker = 2,label=m["Long_Name"])
Plots.plot!(l["lat"][t],l["val"][t],marker = 2,label=l["Long_Name"])
# -

t=1:5:length(s["lat"])
Plots.plot(s["val"][t],marker = 2,label=s["Long_Name"],
    ylabel=s["Unit"], title="Particles (LISST C) in Gradients 3")
Plots.plot!(m["val"][t],marker = 2,label=m["Long_Name"])
Plots.plot!(l["val"][t],marker = 2,label=l["Long_Name"])
Plots.plot!(s["lat"][t]/10,linecolor=:black,linestyle = :dash,label="°N / 10")

# +
#cmap.get_dataset_metadata("tblMGL1704_Gradients2_uway_optics")

# +
s=cmap_get("tblMGL1704_Gradients2_uway_optics","LISST_small")
m=cmap_get("tblMGL1704_Gradients2_uway_optics","LISST_medium")
l=cmap_get("tblMGL1704_Gradients2_uway_optics","LISST_large")

t=1:5:length(s["lat"])
Plots.plot(s["val"][t],marker = 2,label=s["Long_Name"],
    ylabel=s["Unit"], title="Particles (LISST C) in Gradients 2")
Plots.plot!(m["val"][t],marker = 2,label=m["Long_Name"])
Plots.plot!(l["val"][t],marker = 2,label=l["Long_Name"])
Plots.plot!(s["lat"][t]/10,linecolor=:black,linestyle = :dash,label="°N / 10")

# +
s=cmap_get("tblKOK1606_Gradients1_uway_optics","LISST_small")
m=cmap_get("tblKOK1606_Gradients1_uway_optics","LISST_medium")
l=cmap_get("tblKOK1606_Gradients1_uway_optics","LISST_large")

t=1:5:length(s["lat"])
Plots.plot(s["val"][t],marker = 2,label=s["Long_Name"],
    ylabel=s["Unit"], title="Particles (LISST C) in Gradients 1")
Plots.plot!(m["val"][t],marker = 2,label=m["Long_Name"])
Plots.plot!(l["val"][t],marker = 2,label=l["Long_Name"])
Plots.plot!(s["lat"][t]/10,linecolor=:black,linestyle = :dash,label="°N / 10")
# -


