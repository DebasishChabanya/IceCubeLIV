# analysis/Extension.jl

using Plots
using Plots.PlotMeasures
using QuadGK
using Statistics
using Random


# Settings used for the report
const ENERGY_FRACTION = 0.30
const Z_MIN = 1e-4
const MC_SEED = 2026
const RESULTS_DIR = joinpath(@__DIR__, "..", "results")
mkpath(RESULTS_DIR)

# A bounded range, not a Gaussian 1-sigma error
sample_energy(E) = E * (1 - ENERGY_FRACTION + 2 * ENERGY_FRACTION * rand())

function mc_summary(values)
    return (mean(values), std(values), median(values),
            quantile(values, 0.025), quantile(values, 0.975))
end

function save_summary(filename, rows)
    open(joinpath(RESULTS_DIR, filename), "w") do io
        println(io, "case,mean,sd,median,p025,p975")
        for (name, values) in rows
            println(io, name, ",", join(mc_summary(values), ","))
        end
    end
end


include("../src/GRB_RedshiftSampler.jl")
using .GRB_RedshiftSampler

# 13 selected candidates; for neutrino 107, retaining the smallest angular separation
E_TeV  = [79.76, 78.45, 64.05, 168.32, 70.79, 136.87, 62.98, 83.37, 261.23, 138.47, 191.39, 133.85, 95.12]
dt_days = [2.14, 1.86, 0.85, 2.65, 1.25, 1.57, 0.27, 2.50, 1.35, 1.00, 0.72, 2.88, 1.40]
dt_sec = dt_days .* 86400.0 # Convert days to seconds

# 2 measured redshifts, 11 missing (-999.0)
measured_z_array = [-999.0, -999.0, 1.38, -999.0, -999.0, -999.0, -999.0, -999.0, 0.645, -999.0, -999.0, -999.0, -999.0]

# T90 parameters (2 short GRBs <= 2.0s, 11 long GRBs > 2.0s or missing)
t90_array = [38.66, 33.28, 25.40, 90.50, 4.61, 0.13, 2.43, 0.10, 36.10, 127.75, 25.86, 70.40, 24.58]

# QG constants (Planck 2018)
const H0_km_s_Mpc = 67.4  
const MPC_TO_KM = 3.08567758e19
const H0_s = H0_km_s_Mpc / MPC_TO_KM 
const M_p_TeV = 1.2209e16 

function D_z(z::Float64)
    omega_m = 0.315
    omega_l = 0.685
    integrand(zp) = (1.0 + zp) / sqrt(omega_m * (1.0 + zp)^3 + omega_l)
    res, _ = quadgk(integrand, 0.0, z)
    return res / H0_s
end

const D1 = D_z(1.0)
const QG_CONVERSION = D1 / M_p_TeV 

function make_fit_line(m_sec_per_tev, max_x)
    x_line = range(0, max_x * 1.1, length=100)
    y_line = m_sec_per_tev .* x_line
    return x_line, y_line
end

# 10,000-RUN MONTE CARLO (Sampled redshift + error)
println("--- Extended data: Sampling (10,000 MC Runs) ---")
N_sims = 10000
Random.seed!(MC_SEED)
clamp_values = sort(unique([1e-6, 1e-4, 1e-3, 1e-2, Z_MIN]))
eta_by_clamp = zeros(N_sims, length(clamp_values))
mc_etas = zeros(N_sims)

for i in 1:N_sims
    sim_z = zeros(13)
    for j in 1:13
        if measured_z_array[j] != -999.0
            sim_z[j] = measured_z_array[j]
        else
            sim_z[j] = GRB_RedshiftSampler.get_final_physical_redshift(-999.0, t90_array[j])
        end
    end
    
    # CLAMP: Prevent division by zero if MC rolls exactly z=0
    raw_sim_z = copy(sim_z)
    sim_z = max.(sim_z, Z_MIN)
    
    sim_E = sample_energy.(E_TeV)
    sim_K = [sim_E[j] * (D_z(sim_z[j]) / D1) for j in eachindex(E_TeV)]
    
    # Origin-constrained WLS with relative weights proportional to 1/K^2
    sim_slopes = dt_sec ./ sim_K
    mc_etas[i] = mean(sim_slopes) / QG_CONVERSION

    for (k, cutoff) in enumerate(clamp_values)
        cutoff_z = max.(raw_sim_z, cutoff)
        cutoff_K = [sim_E[j] * (D_z(cutoff_z[j]) / D1) for j in eachindex(E_TeV)]
        eta_by_clamp[i, k] = mean(dt_sec ./ cutoff_K) / QG_CONVERSION
    end
end

mean_eta_mc = mean(mc_etas)
std_eta_mc = std(mc_etas)

println("Extended eta = $(round(mean_eta_mc, digits=1)) ± $(round(std_eta_mc, digits=1)) (Dimensionless)")

# One redshift example, shown at recorded energies with the assumed energy range
Random.seed!(2027)
phys_z_plot = zeros(13)
for i in 1:13
    if measured_z_array[i] != -999.0
        phys_z_plot[i] = measured_z_array[i]
    else
        phys_z_plot[i] = GRB_RedshiftSampler.get_final_physical_redshift(-999.0, t90_array[i])
    end
end
phys_z_plot = max.(phys_z_plot, Z_MIN) 

K_phys_plot = [E_TeV[i] * (D_z(phys_z_plot[i]) / D1) for i in 1:13]
K_phys_err_plot = ENERGY_FRACTION .* K_phys_plot # 30% instrumental error

# Plot Generation
p1 = scatter(K_phys_plot, dt_sec, xerror=K_phys_err_plot, color=:darkgreen, markersize=6, 
             label="13 selected candidates", legend=:bottomright)
x_line, y_line = make_fit_line(mean_eta_mc * QG_CONVERSION, maximum(K_phys_plot))

plot!(p1, x_line, y_line, color=:black, linewidth=2.5, 
      label="MC mean (η = $(round(mean_eta_mc, digits=1)))")

x_med, y_med = make_fit_line(median(mc_etas) * QG_CONVERSION, maximum(K_phys_plot))
plot!(p1, x_med, y_med, color=:black, linestyle=:dash, linewidth=2,
      label="MC median (η = $(round(median(mc_etas), digits=1)))")

plot_opts = (size=(800, 650), xlabel="Kinematic Factor K (TeV)", ylabel="Time Delay Δt (seconds)", 
             margin=10mm, bottom_margin=15mm, left_margin=15mm,
             guidefontsize=12, tickfontsize=10, legendfontsize=9,
             ylims=(0, 260000)) # Clamped to keep the 3-day data points visible

plot!(p1; plot_opts...)
savefig(p1, joinpath(RESULTS_DIR, "4_Extended_13_Events.png"))

println("\nSuccessfully rendered! Saved Extended_13_Events.png")

save_summary("extension_summary.csv", [("redshift_and_energy", mc_etas)])

# The same energy and redshift draws are used for each cutoff
open(joinpath(RESULTS_DIR, "extension_cutoff_sensitivity.csv"), "w") do io
    println(io, "z_min,mean,sd,median,p025,p975")
    for (k, cutoff) in enumerate(clamp_values)
        println(io, cutoff, ",", join(mc_summary(eta_by_clamp[:, k]), ","))
    end
end
open(joinpath(RESULTS_DIR, "extension_mc_draws.csv"), "w") do io
    println(io, "draw,eta")
    for i in eachindex(mc_etas)
        println(io, i, ",", mc_etas[i])
    end
end

println("MC median = ", round(median(mc_etas), digits=2))
println("Central 95% of MC draws = ", round.(quantile(mc_etas, [0.025, 0.975]), digits=2))
println("The MC standard deviation is not a fitted confidence interval.")

# Show the tail rather than hiding it behind a mean +/- standard deviation
eta_bins = 10 .^ range(log10(minimum(mc_etas)), log10(maximum(mc_etas)) + 1e-8, length=61)
eta_counts = [count(x -> eta_bins[k] <= x < eta_bins[k+1], mc_etas) for k in 1:length(eta_bins)-1]
eta_centres = sqrt.(eta_bins[1:end-1] .* eta_bins[2:end])
occupied = eta_counts .> 0
p_dist = bar(eta_centres[occupied], eta_counts[occupied],
             bar_width=diff(eta_bins)[occupied], xscale=:log10, yscale=:log10,
             fillrange=0.5, ylims=(0.5, maximum(eta_counts)*2),
             xlims=(minimum(mc_etas)*0.8, maximum(mc_etas)*1.2),
             color=:darkgreen, alpha=0.75, label="10,000 draws",
             xlabel="Eta (log scale)", ylabel="Number of draws (log scale)",
             size=(800, 550), margin=10mm, legend=:topright)
vline!(p_dist, [mean_eta_mc], color=:black, linewidth=2, label="Mean")
vline!(p_dist, [median(mc_etas)], color=:black, linestyle=:dash, linewidth=2, label="Median")
savefig(p_dist, joinpath(RESULTS_DIR, "6_Eta_distribution_13_events.png"))
