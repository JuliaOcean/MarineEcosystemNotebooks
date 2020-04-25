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

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ### Set up tools
#
# [PyCmap](https://github.com/simonscmap/pycmap) is the Python API that we will use in Julia, via [PyCall.jl](https://github.com/JuliaPy/PyCall.jl), to query the CMAP data base. [Plots.jl](http://docs.juliaplots.org/latest/) is a common `Julia` plotting package and `helper_functions.jl` adds a few convenience functions.

# + {"slideshow": {"slide_type": "skip"}, "cell_type": "markdown"}
# _Pre-requisites:_
#
# - 1. _install the [PyCmap](https://github.com/simonscmap/pycmap) python package and its dependencies using `pip`_
# - 2. _compile [PyCall.jl](https://github.com/simonscmap/pycmap) using external python distribution that installed `PyCmap`_
# - 3. _obtain your own API key from [Simons' CMAP](https://simonscmap.com) (free; takes `<30s`)_
# - 4. _import `pycmap` (via `pycall`), `Plots`, and `helper functions`_
#
# _You may need to replace `your-own-API-key` (as outline below) with your own API key from [Simons' CMAP](https://simonscmap.com) and uncomment the command below._

# + {"cell_style": "split", "slideshow": {"slide_type": "skip"}}
if false
    run(`pip install pycmap`) #pycmap is used via PyCall later
    run(pipeline(`which python`,"whichpython.txt")) #external python path
    ENV["PYTHON"]=readline("whichpython.txt")
    import Pkg; Pkg.build("PyCall")
end

# + {"cell_style": "split", "slideshow": {"slide_type": "skip"}}
using PyCall
PyCmap = pyimport("pycmap")
#cmap = PyCmap.API(token="your-own-API-key")
cmap = PyCmap.API()

using Plots
include("helper_functions.jl")

# + {"slideshow": {"slide_type": "skip"}, "cell_type": "markdown"}
# ### Get Data Catalog
#
# _Downloading the catalog is a simple way to verify that cmap is all set-up. The commented `df.to_csv` command writes the content of `df` to a new `catalog.csv` file. Alternatively, `Pandas.jl` can be used as also shown._

# + {"slideshow": {"slide_type": "skip"}}
df = cmap.get_catalog();
#df.to_csv("catalog.csv")

#df=Pandas.DataFrame(cmap.get_catalog());
#to_csv(df,"catalog.csv")

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ### Download & File Data
#
# A simple method is to download data from `CMAP` and store it to a `CSV` file which any software can then reload. Below, `cmap_helpers.tables` provide CMAP `table` lists for [SCOPE-Gradients](http://scope.soest.hawaii.edu/data/gradients/data/) cruise data, which `cmap.get_dataset` downloads one at a time.
#
# Alternatively one can use the computer's memory (next slides).

# + {"slideshow": {"slide_type": "-"}}
pth="../samples/gradients/"
!isdir("$pth") ? mkdir("$pth") : nothing

Γ=1:3
γ=filter(x -> occursin("tblKM1906",x), readdir(pth))
!isempty(γ) ? Γ = [] : nothing

for g in Γ
    list0=cmap_helpers.tables("G$g")
    for i in list0
        df=cmap.get_dataset(i)
        df.to_csv("$pth$i.csv")
    end
end

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ### Read Data + Meta-Data
#
# As an example below we read, and then plot, the `LISST` data collected during the `Gradients 3` cruise.
# -

s=cmap_helpers.get("tblKM1906_Gradients3_uway_optics","LISST_small")
m=cmap_helpers.get("tblKM1906_Gradients3_uway_optics","LISST_medium")
l=cmap_helpers.get("tblKM1906_Gradients3_uway_optics","LISST_large")

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ### Get Ancillary Data
#
# Here we interpolate a `sea surface height` climatology estimate ([Forget et al 2015](http://doi.org/10.5194/gmd-8-3071-2015)) along the ship track. _Any number of commonly available methods can readily be used to interpolate gridded estimates to observed locations. For `MITgcm` output in our example, we use `MeshArrays.jl` method._
# -

ssh=cbiomes_helpers.myinterp(pth,"SSH",s["lon"],s["lat"])
ssh=merge(ssh,Dict("Unit" => "m", "Variable" => "SSH", "Long_Name" => "Sea Surface Height (Mean Dynamic Topography)"))

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ### Plot Data vs Latitude
#
# Here we plot the LISST data collected in Gradients 3 and the sea surface height climatology.
# -

t=1:5:length(s["lat"])
scatter(s["lat"][t],s["val"][t],marker = 1.5,label=s["Long_Name"],
    xlabel="°N",ylabel=s["Unit"], title="Gradients 3")
scatter!(m["lat"][t],m["val"][t],marker = 1.5,label=m["Long_Name"])
scatter!(l["lat"][t],l["val"][t],marker = 1.5,label=l["Long_Name"])
plot!(ssh["lat"][t],(1.0.-ssh["val"][t]).*5.0,linecolor=:black,label="(1.- ssh in m) * 5")

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ### Plot Data vs Station ID
#
# Here we plot `Gradients 1`, `Gradients 2`, and then `Gradients 3` one after the other.
#
# _Note how strongly everything covaries with the ship moving back and forth across the gyre_

# + {"cell_style": "center", "slideshow": {"slide_type": "subslide"}}
s1=cmap_helpers.get("tblKOK1606_Gradients1_uway_optics","LISST_small")
m1=cmap_helpers.get("tblKOK1606_Gradients1_uway_optics","LISST_medium")
l1=cmap_helpers.get("tblKOK1606_Gradients1_uway_optics","LISST_large")
ssh1=cbiomes_helpers.myinterp(pth,"SSH",s1["lon"],s1["lat"])

t=1:1:length(s1["lat"])
scatter(s1["val"][t],marker = 2,label=s["Long_Name"],
    ylabel=s["Unit"], title="Gradients 1")
scatter!(m1["val"][t],marker = 2,label=m["Long_Name"])
scatter!(l1["val"][t],marker = 2,label=l["Long_Name"])
plot!((1.0.-ssh1["val"][t]).*20.0,linecolor=:black,label="(1.- ssh in m) * 20")
plot!(s1["lat"][t]./10.,linecolor=:black,linestyle = :dash,label="°N / 10")

# + {"cell_style": "center", "slideshow": {"slide_type": "subslide"}}
s2=cmap_helpers.get("tblMGL1704_Gradients2_uway_optics","LISST_small")
m2=cmap_helpers.get("tblMGL1704_Gradients2_uway_optics","LISST_medium")
l2=cmap_helpers.get("tblMGL1704_Gradients2_uway_optics","LISST_large")
ssh2=cbiomes_helpers.myinterp(pth,"SSH",s2["lon"],s2["lat"])

t=1:1:length(s2["lat"])
scatter(s2["val"][t],marker = 2,label=s["Long_Name"],
    ylabel=s["Unit"], title="Gradients 2")
scatter!(m2["val"][t],marker = 2,label=m["Long_Name"])
scatter!(l2["val"][t],marker = 2,label=l["Long_Name"])
plot!((1.0.-ssh2["val"][t]).*5.0,linecolor=:black,label="(1.- ssh in m) * 5")
plot!(s2["lat"][t]./10.0,linecolor=:black,linestyle = :dash,label="°N / 10")

# + {"slideshow": {"slide_type": "subslide"}}
t=1:2:length(s["lat"])
plot(s["val"][t],marker = 2,label=s["Long_Name"],
    ylabel=s["Unit"], title="Gradients 3")
plot!(m["val"][t],marker = 2,label=m["Long_Name"])
plot!(l["val"][t],marker = 2,label=l["Long_Name"])
plot!((1.0.-ssh["val"][t]).*5.0,linecolor=:black,label="(1.- ssh in m) * 5")
plot!(s["lat"][t]./10.0,linecolor=:black,linestyle = :dash,label="(lat in °N) / 5")
# -


