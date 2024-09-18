import Pkg; Pkg.activate(".")
using Revise
using TMI
using Test
using Interpolations
using Statistics
using LinearAlgebra
Pkg.activate("./scripts")
Pkg.add("Plots")
Pkg.add("PyPlot")
Pkg.instantiate()

using Plots

Bootst = 2000 
function bootstrap_stats(TMIversion, γ, N, B)
    mean_ys = Float64[]
    std_ys = Float64[]
    mean_y_os = Float64[]
    std_y_os = Float64[]

    for i in 1:Bootst
        y, _, _, _, _, _ = surface_synthetic_observations(TMIversion, "θ", γ, N,σ=0.0)
        y_o, _, _, _, _, _ = ocean_synthetic_observations(TMIversion, "θ", γ, N,σ=0.0)
        
        push!(mean_ys, mean(y))
        push!(std_ys, std(y))
        push!(mean_y_os, mean(y_o))
        push!(std_y_os, std(y_o))
    end

    return mean_ys, std_ys, mean_y_os, std_y_os
end


#======================================================================================================#
# What is the mean MOT and SST during LGM from random samples of N=20?
TMIversion_lgm = "LGM_90x45x33_G14"
A_lgm, Alu_lgm, γ_lgm, TMIfile_lgm, L_lgm, B_lgm = config_from_nc(TMIversion_lgm);
N=20

# Perform bootstrap resampling on LGM data
mean_ys_lgm, std_ys_lgm, mean_y_os_lgm, std_y_os_lgm = bootstrap_stats(TMIversion_lgm, γ_lgm, N, B_lgm)

# Hottest, coldest and true mean  temp. during LGM on the scatter plot
y_all_lgm, _, _ = surface_synthetic_observations(TMIversion_lgm, "θ", γ_lgm,σ=0.0)
y_o_all_lgm, _, _ = synthetic_observations(TMIversion_lgm, "θ", γ_lgm,σ=0.0)

mean_y_all_lgm = round(mean(y_all_lgm, γ_lgm), digits=2)
mean_y_o_all_lgm = round(mean(y_o_all_lgm), digits=2)

n = 100;
y_hot20_lgm = sort(y_all_lgm)[1:n]
y_o_hot20_lgm = sort(y_o_all_lgm)[1:n]

y_cool20_lgm = sort(y_all_lgm,rev=false)[1:n]
y_o_cool20_lgm = sort(y_o_all_lgm,rev=false)[1:n]



###==============================================================================
# What is the mean MOT and SST during modern time from random samples of N=20?

TMIversion_mod = "modern_90x45x33_G14_v2"
A_mod, Alu_mod, γ_mod, TMIfile_mod, L_mod, B_mod = config_from_nc(TMIversion_mod);
mean_ys_mod, std_ys_mod, mean_y_os_mod, std_y_os_mod = bootstrap_stats(TMIversion_mod, γ_mod, N, B_mod)
# Hottest, coldest and true mean  temp. during LGM on the scatter plot
y_all_mod, _, _ = surface_synthetic_observations(TMIversion_mod, "θ", γ_mod,σ=0.0)
y_o_all_mod, _, _ = synthetic_observations(TMIversion_mod, "θ", γ_mod,σ=0.0)

mean_y_all_mod = round(mean(y_all_mod, γ_mod), digits=2)
mean_y_o_all_mod = round(mean(y_o_all_mod), digits=2)
n = 100;
y_hot20_mod = sort(y_all_mod)[1:n]
y_o_hot20_mod = sort(y_o_all_mod)[1:n]

y_cool20_mod = sort(y_all_mod,rev=false)[1:n]
y_o_cool20_mod = sort(y_o_all_mod,rev=false)[1:n]




scatter(mean_ys_lgm, mean_y_os_lgm, xlims=(-3,30), ylims=(-3,30), 
        xlabel="Mean Surface Temperature \$\\bar{\\theta}_s(^\\circ C)\$", 
        ylabel="Mean Ocean Temperature \$\\bar{\\theta}(^\\circ C)\$", 
        title="LGM vs. Mod. Ocean Temp.", 
        label="mean from random sampling of 20 LGM points",
        legend=:topleft,markerstrokewidth=0.3,markersize=3,
        legendfontsize=8, dpi=1000)

# Add more scatter plots to the same figure
scatter!(mean_ys_mod, mean_y_os_mod,
        label="mean from random sampling of 20 mod. points",marker=(:diamond, 3),
        markerstrokewidth=0.3,opacity=0.5,
        legendfontsize=8)


scatter!(y_hot20_lgm, y_o_hot20_lgm, label="Hot 100 points",markerstrokewidth=0,
        legendfontsize=8)
scatter!(y_hot20_mod, y_o_hot20_mod, label="Hot 100 points",marker=(:diamond, 3),markerstrokewidth=0,
        legendfontsize=8)

scatter!(y_cool20_lgm, y_o_cool20_lgm, label="Cool 100 points",markerstrokewidth=0,
        legendfontsize=8)
scatter!(y_cool20_mod, y_o_cool20_mod, label="Cool 100 points",marker=(:diamond, 3),markerstrokewidth=0,
        legendfontsize=8)

scatter!([mean_y_all_lgm], [mean_y_o_all_lgm], 
        label="True mean LGM ($mean_y_all_lgm Surf. and $mean_y_o_all_lgm Global)",
        markerstrokewidth=2, markersize=8,
        legendfontsize=8)

scatter!([mean_y_all_mod], [mean_y_o_all_mod], 
         label="True mean Mod. ($mean_y_all_mod Surf. and $mean_y_o_all_mod Global)",
         marker=(:diamond, 5), markerstrokewidth=2, opacity=0.9,
        legendfontsize=8)

savefig("Sampling.png",)