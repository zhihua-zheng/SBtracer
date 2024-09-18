import Pkg; Pkg.activate(".")
using Revise
using TMI
using Test
# using GeoPythonPlot
using Interpolations
using Statistics
using LinearAlgebra
import Pkg; 
Pkg.activate("./scripts")
Pkg.add("Plots")
Pkg.add("PyPlot")
Pkg.instantiate()

using Plots

#======================================================================================================#
# What is the mean MOT and SST during modern time and LGM from random samples of N=20?

TMIversion = "LGM_90x45x33_OG18"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);
N=20
N_surf =4050 #20
N_ocean = 4050*33
# take synthetic, noisy observations
y, W⁻, ctrue, ytrue, locs, wis = TMI.surface_synthetic_observations(TMIversion,"θ",γ,N)

y_o, W_o⁻, ctrue_o, ytrue_o, locs_o, wis_o =  TMI.ocean_synthetic_observations(TMIversion,"θ",γ,N)
mean_y = round(mean(y),digits=2)
mean_y_o = round(mean(y_o),digits=0)

pyplot() 
plot(y, xlabel="Samples", ylabel="Temp", label="Surface (Modern)")
# Add a line plot for y_o to the same figure
plot!(y_o, xlabel="Samples", ylabel="Temp", label="Ocean (Modern)")

# Set the title with the mean values
title!("Mean SST: $mean_y,  MOT: $mean_y_o")



#####===============================================================##########
#####================Bootstrapping==================================##########
Bootst = 1000  # Number of bootstrap samples

# Function to perform bootstrap resampling and calculate means and standard deviations
function bootstrap_stats(TMIversion, γ, N, B)
    mean_ys = Float64[]
    std_ys = Float64[]
    mean_y_os = Float64[]
    std_y_os = Float64[]

    for i in 1:Bootst
        y, _, _, _, _, _ = TMI.surface_synthetic_observations(TMIversion, "θ", γ, N,σ=0.0)
        y_o, _, _, _, _, _ = TMI.ocean_synthetic_observations(TMIversion, "θ", γ, N,σ=0.0)
        
        push!(mean_ys, mean(y))
        push!(std_ys, std(y))
        push!(mean_y_os, mean(y_o))
        push!(std_y_os, std(y_o))
    end

    return mean_ys, std_ys, mean_y_os, std_y_os
end

# Perform bootstrap resampling
mean_ys, std_ys, mean_y_os, std_y_os = bootstrap_stats(TMIversion, γ, N, B)

# Compute overall means
mean_y = mean(mean_ys)
mean_y_o = mean(mean_y_os)

# Compute standard deviations of the means
std_y = std(mean_ys)
std_y_o = std(mean_y_os)




# Generate plot with standard deviation shading
x = 1:Bootst

plot(x, mean_ys, ribbon=std_ys, label="Surface Mean ± Std Dev", xlabel="Samples", ylabel="Temp")
plot!(x, mean_y_os, ribbon=std_y_os, label="Ocean Mean ± Std Dev", title="Bootstrap Mean SST and MOT with Standard Deviation Shading")

scatter(mean_ys,mean_y_os,xlims=(0,20), ylims=(0,20), xlabel="Mean Surface Temperature \$\\bar{\\theta}_s\$",
ylabel="Mean Ocean Temperature \$\\bar{\\theta}\$",title="Scatter Plot of Mean Temperatures",label="mean from random sampling of 20 points")

######===========================================================##################
# Hottest and coldest on the scatter plot
#y_all, _, _, _, _, _ = surface_synthetic_observations(TMIversion, "θ", γ, N_surf,σ=0.0)
#y_o_all, _, _, _, _, _ = ocean_synthetic_observations(TMIversion, "θ", γ, N_ocean,0.0)
y_all, _, _ = surface_synthetic_observations(TMIversion, "θ", γ,σ=0.0)
y_o_all, _, _ = synthetic_observations(TMIversion, "θ", γ,σ=0.0)
n = 100;
y_hot20 = sort(y_all)[1:n]
y_o_hot20 = sort(y_o_all)[1:n]

y_cool20 = sort(y_all,rev=false)[1:n]
y_o_cool20 = sort(y_o_all,rev=false)[1:n]

scatter(mean_ys, mean_y_os, xlims=(-3,30), ylims=(-3,30), 
        xlabel="Mean Surface Temperature \$\\bar{\\theta}_s\$", 
        ylabel="Mean Ocean Temperature \$\\bar{\\theta}\$", 
        title="LGM Scatter Plot of Mean Temperatures", 
        label="mean from random sampling of 20 points",
        legend=:topleft)

# Add more scatter plots to the same figure
scatter!(y_hot20, y_o_hot20, label="Hot 20 points")
scatter!(y_cool20, y_o_cool20, label="Cool 20 points")
scatter!([mean(y_all,γ)], [mean(y_o_all)], label="True mean LGM")

