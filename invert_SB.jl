import Pkg; Pkg.activate(".")
using Revise
using TMI

fname = "/Users/zhihua/Documents/Work/Research/Projects/SBtracer/data/SB_ctrl.nc"
TMIfname = replace(fname, ".nc" => "_TMI.nc")
γ =  Grid(fname, "maskC", "XC", "YC", "Z")

y = (θ = readfield(fname,"THETA",γ),
     S = readfield(fname,"SALT",γ))

w = (θ = 0.001,
     S = 0.001)

@time m̃ = massfractions(y, w)
Ã = watermassmatrix(m̃, γ)

# a first guess: observed surface boundary conditions are perfect.
# set surface boundary condition to the observations.
bθ = getsurfaceboundary(y.θ)

## reconstruct temperature
using LinearAlgebra: lu
Ãlu = lu(Ã)
θ̃ = steadyinversion(Ãlu,bθ,γ)

# compare to c.θ
println("difference max: ", Base.maximum(y.θ - θ̃))
println("difference min: ", Base.minimum(y.θ - θ̃))

# ocean volume that has originated from each surface box
volume = volumefilled("",Ãlu,γ)
writefield(TMIfname, volume)
