# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.3
#   kernelspec:
#     display_name: Julia 1.9.1
#     language: julia
#     name: julia-1.9
# ---

# # Model And Observed Maps
#
# Here we (1) read `CBIOMES model output` and `ocean color data` from file, and (2) illustrate these arrays as global maps in `julia`.
#
# <img src="../figs/cbiomes-01.png" alt="Drawing" style="height: 100px;"/>

# + [markdown] {"slideshow": {"slide_type": "slide"}}
# ### Read variables from file
#
# Two-dimensional arrays for longitude (`lon`), latitude (`lat`), irradiance reflectance from a model (`drwn3_Rirr` at `wv_drwn3`), and remotely sensed reflectance from satellite data (`cci_Rrs_490` at 490nm) are read from files in the `samples/` folder.

# +
using Downloads

OceanColor_samples=joinpath(tempdir(),"OceanColor_samples")
url="https://github.com/JuliaOcean/MarineEcosystemNotebooks/raw/master/samples/"
files=("lon.bin","lat.bin","cci_Rrs_490.bin","drwn3_Rirr.bin")
if !isdir(OceanColor_samples)
    mkdir(OceanColor_samples)
    [Downloads.download(url*f,joinpath(OceanColor_samples,f)) for f in files]
end
        
"Done downloading files"
# -

function read_sample(p,f;n=1)
    n==1 ? a = zeros(Float32,(720,360)) : a = zeros(Float32,(720,360,n))
    hton.(read!(joinpath(p,f),a))
end

# + {"slideshow": {"slide_type": "-"}}
lon=read_sample(OceanColor_samples,files[1])
lat=read_sample(OceanColor_samples,files[2])
cci_Rrs_490=read_sample(OceanColor_samples,files[3])
drwn3_Rirr=read_sample(OceanColor_samples,files[4],n=13)

cci_Rrs_490[findall(cci_Rrs_490.==0)].=NaN
drwn3_Rirr[findall(drwn3_Rirr.==0)].=NaN

"Done reading files"

# + [markdown] {"slideshow": {"slide_type": "subslide"}}
# ### Model and data wavebands
#
# Currently, the `OC-CCI` [satellite data set](https://esa-oceancolour-cci.org) provides remotely sensed reflectance at 6 wavelengths (`wv_cci` in `nm`) while the `CBIOMES-global` [ocean model](https://cbiomes.readthedocs.io/) outputs irradiance reflectance at 13 wavelengths (`wv_drwn3` in `nm`). 
# -

wv_cci=[412, 443, 490, 510, 555, 670]
wv_drwn3=[400,425,450,475,500,525,550,575,600,625,650,675,700];

# + [markdown] {"slideshow": {"slide_type": "slide"}}
# ### Display model and data maps
#
# The `Plots.jl` package provides convenient plotting recipes like `heatmap`. the following command assumes that `Plots.jl` has aleary been installed using `julia`'s package manager (documentation available [here](https://docs.julialang.org/en/)). 
# -

using Plots

# + {"slideshow": {"slide_type": "subslide"}}
ii=3
heatmap(vec(lon[:,1]),vec(lat[1,:]), transpose(drwn3_Rirr[:,:,ii]), clims=(0,0.08))
title!("irradiance reflectance from a model (at $(wv_drwn3[ii])nm)")

# + {"slideshow": {"slide_type": "subslide"}}
heatmap(vec(lon[:,1]),vec(lat[1,:]), transpose(cci_Rrs_490), clims=(0,0.01))
title!("remotely sensed reflectance from satellite (at 490nm)")

# + [markdown] {"slideshow": {"slide_type": "slide"}}
# ### Your turn!
#
# Here is one idea in case you want to take this a bit further:
#
# - Apply these to all points in `drwn3_Rirr` using for loops or broadcast.
# - Plot the resulting map of `Rrs` to compare with the data map (`cci_Rrs_490`).
