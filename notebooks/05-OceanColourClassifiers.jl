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
#     display_name: Julia 1.1.0
#     language: julia
#     name: julia-1.1
# ---

# # Simple ocean color data classifiers
#
# This notebook applies simple classifications to `ocean color data` in context of the [CBIOMES](https://cbiomes.org) project. It is written in [Julia](https://julialang.org) and can be used interactively via [binder](https://mybinder.org/v2/gh/gaelforget/Cbiomes2019Notebooks/master). 
#
# <img src="../figs/cbiomes-01.png" alt="Drawing" style="height: 50px;"/>

# ### Activate packages for later use
#
# It is assumed that listed packages have aleary been installed using `julia`'s package manager (documentation available [here](https://docs.julialang.org/en/)). 

using Plots, Distributions, NetCDF, NCDatasets

# ### Optical classification using reflectances
#
# `Fuzzy logic` classifiers defined in [Moore et al 2009](https://doi.org/10.1016/j.rse.2009.07.016) and [Jackson et al 2017](http://dx.doi.org/10.1016/j.rse.2017.03.036) can be used to assign optical class memberships from an `Rrs` vector. While Moore et al define `n=8` classes using an in-situ database, Jackson et al instead define `n=14` classes using a satellite database. The latter benefits from better data coverage across all of the ecological provinces of the global ocean and is used in `OC-CCI`. 
#
# In both cases the classifier is encoded in a mean reflectance spectra (`M[i][1:6]`) and a covariance matrix (`S[i][1:6,1:6]`) provided for each optical class (`i` in `1:n`). Class memberships are then derived by computing the squared Mahalanobis distance to each `M[i]` and passing the result to cumulative chi-squared distribution function (Equations 11 and 12 in [Moore et al 2011](https://doi.org/10.1109/36.942555)).

# +
include("../samples/M09.jl")

M09=Dict("M" => M, "S" => S, "Sinv" => inv.(S))
plot(wv_cci,M,w=3); xlabel!("nm"); ylabel!("Rrs")

# +
#Jackson et al 2017:
tmpM = ncread("../samples/J17.nc", "cluster_means")
tmpSinv = ncread("../samples/J17.nc", "inverse_covariance")

M=Array{Any,1}(undef,14)
Sinv=Array{Any,1}(undef,14)
for ii=1:length(M)
    M[ii]=vec(tmpM[ii,:])
    Sinv[ii]=tmpSinv[1:6,1:6,ii]
end

J17=Dict("M" => M, "Sinv" => Sinv, "S" => inv.(Sinv))
plot(wv_cci,M,w=3); xlabel!("nm"); ylabel!("Rrs")
# -

# ### Class membership function

function fcm(M,Sinv,Rrs)
    f=Array{Any,1}(undef,length(M))
    for ii=1:length(M)
        X=vec(Rrs)-M[ii]
        Z=transpose(X)*Sinv[ii]*X
        f[ii]=ccdf(Chisq(6),Z)
    end
    f
end

# ### Apply J17 classifier to 2D region

# Read file and display one waveband

# +
dir0="../samples/"
fil=dir0*"ESACCI-OC-RRS-sample-fv4.0.nc"
ds = Dataset(fil)

Rrs_412=ds["Rrs_412"]
Rrs_443=ds["Rrs_443"]
Rrs_490=ds["Rrs_490"]
Rrs_510=ds["Rrs_510"]
Rrs_555=ds["Rrs_555"]
Rrs_670=ds["Rrs_670"]

heatmap(ds["Rrs_490"][:,:,1])
# -

# Find points that have a full set of input

tmp=fill(false,size(Rrs_412))
for ii in eachindex(Rrs_412)
    !ismissing(Rrs_412[ii]) ? tmp[ii]=!ismissing(Rrs_443[ii].*Rrs_490[ii].*Rrs_510[ii].*Rrs_555[ii].*Rrs_670[ii]) : nothing
end
ii=findall(tmp)

# Compute memberships

# +
mbrshp=Array{Float64,3}(undef,(1440,960,14))
for jj=1:length(ii); 
    kk=ii[jj]
    Rrs_tmp=[Rrs_412[kk] Rrs_443[kk] Rrs_490[kk] Rrs_510[kk] Rrs_555[kk] Rrs_670[kk]]
    mbrshp[kk[1],kk[2],:]=fcm(J17["M"],J17["Sinv"],Rrs_tmp)
end

heatmap(mbrshp[:,:,10])
# -


