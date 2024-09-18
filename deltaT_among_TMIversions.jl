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

# Example version list
vlist = TMI.versionlist()[1:end-1]

# Filter modern and LGM versions
modern_versions = collect(filter(v -> startswith(v, "modern"), vlist))
lgm_versions = collect(filter(v -> startswith(v, "LGM"), vlist)[1:end-1])


modern_versions = filter!(x -> x != "modern_90x45x33_G14", modern_versions)
lgm_versions = filter!(x -> x != "LGM_90x45x33_OG18", lgm_versions)

# Initialize matrices to store differences
surface_temp_diffs = zeros(length(modern_versions), length(lgm_versions))
ocean_temp_diffs = zeros(length(modern_versions), length(lgm_versions))


surface_temp_mean = zeros(length(modern_versions), length(lgm_versions))
ocean_temp_mean = zeros(length(modern_versions), length(lgm_versions))


surface_temp_mean_mod = zeros(length(modern_versions), length(lgm_versions))
surface_temp_mean_lgm = zeros(length(modern_versions), length(lgm_versions))
ocean_temp_mean_mod = zeros(length(modern_versions), length(lgm_versions))
ocean_temp_mean_lgm = zeros(length(modern_versions), length(lgm_versions))


# Iterate over pairs of TMI versions and compute differences
for i in 1:length(modern_versions)
    for j in 1:length(lgm_versions)
        TMIversion_mod = modern_versions[i]
        TMIversion_lgm = lgm_versions[j]

        # Load configurations and compute means for the modern version
        try
            A_mod, Alu_mod, γ_mod, TMIfile_mod, L_mod, B_mod = config_from_nc(TMIversion_mod)
            y_all_mod, _, _ = surface_synthetic_observations(TMIversion_mod, "θ", γ_mod, σ=0.0)
            y_o_all_mod, _, _ = synthetic_observations(TMIversion_mod, "θ", γ_mod, σ=0.0)
            mean_y_all_mod = round(mean(y_all_mod, γ_mod), digits=2)
            mean_y_o_all_mod = round(mean(y_o_all_mod), digits=2)

            # Load configurations and compute means for the LGM version
            A_lgm, Alu_lgm, γ_lgm, TMIfile_lgm, L_lgm, B_lgm = config_from_nc(TMIversion_lgm)
            y_all_lgm, _, _ = surface_synthetic_observations(TMIversion_lgm, "θ", γ_lgm, σ=0.0)
            y_o_all_lgm, _, _ = synthetic_observations(TMIversion_lgm, "θ", γ_lgm, σ=0.0)
            mean_y_all_lgm = round(mean(y_all_lgm, γ_lgm), digits=2)
            mean_y_o_all_lgm = round(mean(y_o_all_lgm), digits=2)

            # Compute differences and store them
            surface_temp_diffs[i, j] = mean_y_all_mod - mean_y_all_lgm
            ocean_temp_diffs[i, j] = mean_y_o_all_mod - mean_y_o_all_lgm

            surface_temp_mean[i, j] = mean_y_all_mod - mean_y_all_lgm
            ocean_temp_mean[i, j] = mean_y_o_all_mod - mean_y_o_all_lgm


            surface_temp_mean_mod[i, j] = mean_y_all_mod
            surface_temp_mean_lgm[i, j] = mean_y_all_lgm
            ocean_temp_mean_mod[i, j] = mean_y_o_all_mod
            ocean_temp_mean_lgm[i, j] = mean_y_o_all_lgm


        catch e
            println("Can't grab TMIVersion: ", TMIversion_mod, " or ", TMIversion_lgm)
            println("Error: ", e)
        end
    end
end





# Function to shorten the version names
function shorten_version_name(version)
    parts = split(version, '_')
    prefix = parts[1] == "modern" ? "M" : "L"
    grid_size = join(split(parts[2], "x")[1:2], "x")
    suffix = parts[end]
    return "$prefix\\_$grid_size\\_$suffix"
end

# Shorten the version names
short_modern_versions = [shorten_version_name(v) for v in modern_versions]
short_lgm_versions = [shorten_version_name(v) for v in lgm_versions]


heatmap1 = heatmap(
    surface_temp_diffs,
    xlabel="LGM Versions",
    ylabel="Modern Versions",
    xticks=(1:length(short_lgm_versions), short_lgm_versions),
    yticks=(1:length(short_modern_versions), short_modern_versions),
    title="Surface Temperature Differences",
    color=:thermal,
    xrotation=45,
    yrotation=55,
    tickfontsize=5,
    titlefontsize=10,
    guidefontsize=10
)

# Create heatmap for ocean temperature differences
heatmap2 = heatmap(
    ocean_temp_diffs,
    xlabel="LGM Versions",
    xticks=(1:length(short_lgm_versions), short_lgm_versions),
    yticks=nothing,
    title="Ocean Temperature Differences",
    color=:thermal,
    xrotation=45,
    yrotation=55,
    tickfontsize=5,
    titlefontsize=10,
    guidefontsize=10
)

plot(heatmap1, heatmap2, layout=2,left_margin=2Plots.mm,bottom_margin=3Plots.mm,dpi=1000)
savefig("heatmap.png",)