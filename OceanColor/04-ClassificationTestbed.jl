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

# # Prepare data set for clustering project
#
# Code below does the bulk of a workflow that may look something like:
#
# - read variable + metadata
# - interpolate to chosen grid
# - add column to dataframe
# - save dataframe to csv
# - save metadata to md
#
# _Side note: use PARF_

# + {"slideshow": {"slide_type": "slide"}}
using MeshArrays
γ=GridSpec("LatLonCap","GRID_LLC90/")
Γ=GridLoad(γ)

lon=[i for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]
lat=[j for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]
(f,i,j,w,_,_,_)=InterpolationFactors(Γ,vec(lon),vec(lat))
λ=Dict("f" => f,"i" => i,"j" => j,"w" => w);

# + {"slideshow": {"slide_type": "slide"}}
using MITgcmTools

fileName="nctiles_climatology/THETA/THETA"
THETA=read_nctiles(fileName,"THETA",γ)
show(THETA)

# + {"slideshow": {"slide_type": "slide"}}
n=length(lon[:])
tmp=vec(fill(NaN,(n*12,1)))
for m=1:12
    k=collect(1:n) .+(m-1)*n
    tmp[k] .= Interpolate(THETA[:,1,m],λ["f"],λ["i"],λ["j"],λ["w"])
end

all=Dict("SST" => tmp)

# + {"slideshow": {"slide_type": "slide"}}
using Plots
SST1=reshape(all["SST"][1:n,1],size(lon))
SST7=reshape(all["SST"][(1:n) .+ n*7,1],size(lon))
heatmap([transpose(SST1) transpose(SST7)])
# -


