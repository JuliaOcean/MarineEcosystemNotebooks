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

# # Single Particle Simulation Example
#
# <table><tr>
# <td> <img src="SolidBodyRotation.gif" alt="Drawing" style="width: 400px;"/> </td>
# <td> <img src="ConvergingSpiral.gif" alt="Drawing" style="width: 400px;"/> </td>
# </tr></table>
#
# Here we simulate the trajectory of a particle drifting in an idealized flow which is a crude representation for an ocean eddy. 
#
# To start, we use solid body rotation around a central point. Exercises are provided at the end (e.g. to add a convergence / divergence term). 
#
# This demonstrated the use of differential equation solvers in Julia which are broadly applicable to many other problems.
#
# - 1. setup the software and initialize example
# - 2. simulate trajectories & plot results
# - 3. experiment with parameters (user)

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ### For More Documentation 
#
# On Julia :
# - <https://julialang.org>
# - <https://docs.julialang.org>
#
# On this notebook :
# - <https://docs.juliadiffeq.org/latest> 
# - <https://en.wikipedia.org/wiki/Displacement_(vector)>
# - <https://juliaclimate.github.io/IndividualDisplacements.jl/dev>
# - <https://juliaclimate.github.io/MeshArrays.jl/dev>

# + {"slideshow": {"slide_type": "slide"}, "cell_style": "center", "cell_type": "markdown"}
# ## 1.1 Import Software
#
# The following `Julia` Packages are used for solving differential equations and plotting results.

# + {"cell_style": "center"}
using OrdinaryDiffEq, Plots
using IndividualDisplacements, MeshArrays
include("helper_functions.jl")

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 1.2  Gridded Domain Example
#
# The `SetPeriodicDomain` helper function sets up a basic grid of size `np x np`.

# + {"slideshow": {"slide_type": "-"}}
np=16
Γ=SetupPeriodicDomain(np);

# +
#?SetupPeriodicDomain
#show(Γ["XC"])
#scatter(Γ["XC"][1],Γ["YC"][1],leg=:none,title="grid points")

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 1.3 Define Time Period & Velocity Fields
#
# For convenience, at the end we store all parameters in `𝑃` (a dictionary). 
#
# _Skipped in initial presentation mode_

# + {"slideshow": {"slide_type": "skip"}}
#time range
t0=0.0
t1=0.98*2*pi
#t1=3.0*2*pi
𝑇 = (t0,t1)

#solid-body rotation around central location
i=Int(np/2+1)
u=-(Γ["YG"].-Γ["YG"][1][i,i])
v=(Γ["XG"].-Γ["XG"][1][i,i])

#add some convergence to / divergence from central location
d=0.0 
#d=-0.1
u=u+d*(Γ["XG"].-Γ["XG"][1][i,i])
v=v+d*(Γ["YG"].-Γ["YG"][1][i,i])

#store everything in a dictionnary
𝑃=Dict("u0" => u, "u1" => u, "v0" => v, "v1" => v, "t0" => t0, "t1" => t1)
𝑃=merge(𝑃,Γ)

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1.4 Define initial position and time period
# -

u0=np*[1/3,1/3]
du=fill(0.0,2);

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2.1 solve for particle trajectory
#
# - `ODEProblem` formulates the differential equation along with the time period `𝑇`, parameters `𝑃`
# - `solve` then performs the integration over `𝑇`, starting from `u0`, using the `Tsit5` solver
#
# _Tsit5 - Tsitouras 5/4 Runge-Kutta method. (free 4th order interpolant)._ ([from the docs](https://docs.sciml.ai/dev/solvers/ode_solve/#))
#
# _Interested in additional documentation? Try `?ODEProblem` or `?solve`_
# -

prob = ODEProblem(⬡,u0,𝑇,𝑃)
sol = solve(prob,Tsit5(),reltol=1e-8);

# + {"slideshow": {"slide_type": "-"}, "cell_type": "markdown"}
# ## 2.2 Post-Process Output
#
# _(Not much in the case of this example)_
# -

x,y=sol[1,:],sol[2,:]
nt=length(x)

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2.3 plot result (particle trajectory)
#
# - define `myplot` convenience function
# - generate animation using `myplot`

# + {"slideshow": {"slide_type": "-"}}
myplot(i)=plot(x[1:i],y[1:i],linewidth=2,arrow = 2,
    title="Solid body rotation / Spiral example",leg=false,
    xaxis="x",yaxis="y",xlims=(0,np),ylims=(0,np))


# + {"slideshow": {"slide_type": "subslide"}}
p=Int(ceil(nt/100))
anim = @animate for i ∈ 1:p:nt
    myplot(i)
end
pth=tempdir()*"/"
gif(anim, pth*"SolidBodyRotation.gif", fps = 15)
# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ### Want a single plot?
#
# Try uncommenting one line at a time in the next cell below.


# + {"slideshow": {"slide_type": "-"}}
#plt=myplot(nt)
#scatter!(plt,[u0[1]],[u0[2]])
#savefig(plt,pth*"SolidBodyRotation.png")

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ### Want to change other parameters?
#
# Try uncommenting, in the code above, these lines:
#
# - `#?SetupPeriodicDomain`
# - `#show(Γ["XC"])`
# - `#t1=3.0*2*pi`
# - `#d=-0.1`
# -



