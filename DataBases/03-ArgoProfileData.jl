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

# # Argo Profile Data : _a North Pacific Example_
#
# Here we plot Argo temperature anomalies over the years in a North Pacific box, where we see possible hints of `NP blob 2.0` & its origin.
#
# _Algorithm:_
#
# - 1. get Argo data set from `GDAC` (ftp) 
# - 2. apply QC & interpolate to standard levels (`MITprof`)
# - 3. collocate with climatology & uncertainty field (`MITprof`)
# - 4. subtract climatology & bin average over region (this notebook)
#
# _Data Source & Format:_
#
# - Roemmich, et al, 2019: On the Future of Argo: A Global, Full-Depth, Multi-Disciplinary Array. Frontiers in Marine Science, 6. https://doi.org/10.3389/fmars.2019.00439
# - Forget, G., J.-M. Campin, P. Heimbach, C. N. Hill, R. M. Ponte, and C. Wunsch, 2015: ECCO version 4: an integrated framework for non-linear inverse modeling and global ocean state estimation. Geoscientific Model Development, 8, 3071-3104, <https://doi.org/10.5194/gmd-8-3071-2015>
#

# + {"slideshow": {"slide_type": "skip"}, "cell_type": "markdown"}
# ### Set up tools
#

# + {"cell_style": "center", "slideshow": {"slide_type": "skip"}}
using ArgoData, Plots, MAT, Dates, Statistics
argo_T = matread("argo_T.mat")
argo_S = matread("argo_S.mat")

# + {"slideshow": {"slide_type": "skip"}, "cell_type": "markdown"}
# ### Depth-Time Anomaly Function

# + {"slideshow": {"slide_type": "skip"}}
years=collect(2006:2019); ny=length(years); nd=Int(argo_T["nr"])
tim=[y+m/12 for m in 1:12, y in years][:]
dep=-argo_T["prof_depth"]


δt=Millisecond(86400. *1000. *30.) # 30 for +/-30 days in miliseconds
δl=5 # 5 for +/-5 degree ranges
argo_time = julian2datetime.( datetime2julian(DateTime(0)) .+ argo_T["prof_date"] )

function TimeDepthMedian(x,λ=40.)
    anom=Array{Float64,2}(undef,ny*12,nd)
    for y=1:ny, m=1:12, d=1:nd
        tt=findall( (abs.(argo_time.-DateTime(years[y],m,15)).<δt) .&
                    (abs.(argo_T["prof_lat"].-λ).<δl) )
        tt=[tt[i][1] for i in 1:length(tt)]
        xx=(x)[tt,d]
        xx=xx[findall((!isnan).(xx))]
        ~isempty(xx) ? xm=median(xx) : xm=NaN
        anom[m+(y-1)*12,d]=xm
    end

    return anom,λ
end

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ### Depth-Time Anomaly Plots
#

# + {"cell_style": "center", "slideshow": {"slide_type": "slide"}}
anom,λ=TimeDepthMedian(argo_T["prof"]-argo_T["monclim"],50.); v="Tobs-Tclim"
contourf(tim,vec(dep)[end:-1:1],transpose(anom[:,end:-1:1]),title="median( $v ) in $λ +/- $δl °N",xlim=(2014.,2020.),ylim=(-300.,0.),clim=(-1.5,1.5))

# + {"cell_style": "center", "slideshow": {"slide_type": "slide"}}
anom,λ=TimeDepthMedian(argo_S["prof"]-argo_S["monclim"],50.); v="Sobs-Sclim"
contourf(tim,vec(dep)[end:-1:1],transpose(anom[:,end:-1:1]),title="median( $v ) in $λ +/- $δl °N",xlim=(2014.,2020.),ylim=(-300.,0.),clim=(-0.3,0.3))
# -


