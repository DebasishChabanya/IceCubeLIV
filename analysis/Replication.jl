# analysis/Replication.jl

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

E_TeV  = [61.7, 134.2, 98.5, 86.5, 86.1, 186.6, 66.7]
dt_sec = [73690.0, 135731.0, 15446.0, 160909.0, 187050.0, 229039.0, 23286.0]

# Proxies used in AC2023
ac2023_z = [1.38, 0.60, 1.38, 1.38, 1.38, 1.38, 1.38] 

# Physical inputs for the z sampler
measured_z_array = [1.38, -999.0, -999.0, -999.0, -999.0, -999.0, -999.0]
t90_array        = [999.0, 0.5, 5.0, 5.0, 5.0, 5.0, 5.0]

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

# Dividing the slope by D(1)/M_p to obtain dimensionless eta
const QG_CONVERSION = D1 / M_p_TeV 

function make_fit_line(m_sec_per_tev, max_x)
    x_line = range(0, max_x * 1.1, length=100)
    y_line = m_sec_per_tev .* x_line
    return x_line, y_line
end

# 1: AC2023 base (Ordinary Least Squares)
K_ac2023 = [E_TeV[i] * (D_z(ac2023_z[i]) / D1) for i in 1:7]
m_ols = sum(K_ac2023 .* dt_sec) / sum(K_ac2023.^2)
residuals_ols = dt_sec .- (m_ols .* K_ac2023)
se_m_ols = sqrt(sum(residuals_ols.^2) / (6 * sum(K_ac2023.^2))) 

eta_ols = m_ols / QG_CONVERSION
eta_ols_err = se_m_ols / QG_CONVERSION

println("--- Panel 1: AC2023 Baseline ---")
println("Eta = $(round(eta_ols, digits=1)) ± $(round(eta_ols_err, digits=1)) (Dimensionless)")

p1 = scatter(K_ac2023, dt_sec, color=:firebrick, markersize=7, 
             label="AC2023 Data", title="", legend=:topleft)
x_line, y_line = make_fit_line(m_ols, maximum(K_ac2023))
plot!(p1, x_line, y_line, color=:black, linestyle=:dash, linewidth=2.5, 
      label="OLS Fit (η = $(round(eta_ols, digits=1)))")

# 2: Energy sampling with fixed redshift proxies
Random.seed!(2025)
energy_only_etas = [
    begin
        sim_E = sample_energy.(E_TeV)
        sim_K = sim_E .* (K_ac2023 ./ E_TeV)
        mean(dt_sec ./ sim_K) / QG_CONVERSION
    end
    for _ in 1:10000
]
eta_wls = mean(energy_only_etas)
eta_wls_err = std(energy_only_etas)
m_wls = eta_wls * QG_CONVERSION

println("--- Panel 2: Instrumental error ---")
println("Eta = $(round(eta_wls, digits=1)) ± $(round(eta_wls_err, digits=1)) (Dimensionless)")

K_err = ENERGY_FRACTION .* K_ac2023
p2 = scatter(K_ac2023, dt_sec, xerror=K_err, color=:mediumblue, markersize=7, 
             label="±30% assumed energy range", title="", legend=:topleft)
x_line, y_line = make_fit_line(m_wls, maximum(K_ac2023))
plot!(p2, x_line, y_line, color=:black, linewidth=2.5, 
      label="MC mean (η = $(round(eta_wls, digits=1)))")

# 3: 10,000-RUN MONTE CARLO (Sampled redshift + error)
println("--- Panel 3: Sampling (10,000 MC Runs) ---")
N_sims = 10000
Random.seed!(MC_SEED)
clamp_values = sort(unique([1e-6, 1e-4, 1e-3, 1e-2, Z_MIN]))
eta_by_clamp = zeros(N_sims, length(clamp_values))
mc_etas = zeros(N_sims)

for i in 1:N_sims
    sim_z = [measured_z_array[j] != -999.0 ? measured_z_array[j] : GRB_RedshiftSampler.get_final_physical_redshift(-999.0, t90_array[j]) for j in 1:7]
    
    # CLAMP: Prevent division by zero if MC rolls exactly z=0
    raw_sim_z = copy(sim_z)
    sim_z = max.(sim_z, Z_MIN)
    
    sim_E = sample_energy.(E_TeV)
    sim_K = [sim_E[j] * (D_z(sim_z[j]) / D1) for j in eachindex(E_TeV)]
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

println("Eta = $(round(mean_eta_mc, digits=1)) ± $(round(std_eta_mc, digits=1)) (Dimensionless)")

# One redshift example, shown at recorded energies with the assumed energy range
Random.seed!(2027)
phys_z_plot = [measured_z_array[i] != -999.0 ? measured_z_array[i] : GRB_RedshiftSampler.get_final_physical_redshift(-999.0, t90_array[i]) for i in 1:7]
phys_z_plot = max.(phys_z_plot, Z_MIN) # Apply the same clamp for plotting

K_phys_plot = [E_TeV[i] * (D_z(phys_z_plot[i]) / D1) for i in 1:7]
K_phys_err_plot = ENERGY_FRACTION .* K_phys_plot

p3 = scatter(K_phys_plot, dt_sec, xerror=K_phys_err_plot, color=:purple, markersize=7, 
             label="Example redshift draw", title="", legend=:topleft)
x_line3, y_line3 = make_fit_line(mean_eta_mc * QG_CONVERSION, maximum(K_phys_plot))
plot!(p3, x_line3, y_line3, color=:black, linewidth=2.5, 
      label="MC mean (η = $(round(mean_eta_mc, digits=1)))")

x_med, y_med = make_fit_line(median(mc_etas) * QG_CONVERSION, maximum(K_phys_plot))
plot!(p3, x_med, y_med, color=:black, linestyle=:dash, linewidth=2,
      label="MC median (η = $(round(median(mc_etas), digits=1)))")

# Save individual plots
plot_opts = (size=(800, 650), xlabel="Kinematic Factor K (TeV)", ylabel="Time Delay Δt (seconds)", 
             margin=10mm, bottom_margin=15mm, left_margin=15mm,
             guidefontsize=12, tickfontsize=10, titlefontsize=14, legendfontsize=9)

plot!(p1; plot_opts...)
savefig(p1, joinpath(RESULTS_DIR, "1_Replication.png"))

plot!(p2; plot_opts...)
savefig(p2, joinpath(RESULTS_DIR, "2_Energy_uncertainty.png"))

plot!(p3; ylims=(0, 250000), plot_opts...)
savefig(p3, joinpath(RESULTS_DIR, "3_Redshift_and_energy.png"))

println("\nSuccessfully rendered! Saved 3 individual plots.")

# Baseline standard error and MC spread are different quantities
open(joinpath(RESULTS_DIR, "replication_baseline.csv"), "w") do io
    println(io, "case,eta,ols_standard_error")
    println(io, "baseline,", eta_ols, ",", eta_ols_err)
end
save_summary("replication_summary.csv", [("energy_only", energy_only_etas), ("redshift_and_energy", mc_etas)])

# The same energy and redshift draws are used for each cutoff
open(joinpath(RESULTS_DIR, "replication_cutoff_sensitivity.csv"), "w") do io
    println(io, "z_min,mean,sd,median,p025,p975")
    for (k, cutoff) in enumerate(clamp_values)
        println(io, cutoff, ",", join(mc_summary(eta_by_clamp[:, k]), ","))
    end
end
open(joinpath(RESULTS_DIR, "replication_mc_draws.csv"), "w") do io
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
             color=:purple, alpha=0.75, label="10,000 draws",
             xlabel="Eta (log scale)", ylabel="Number of draws (log scale)",
             size=(800, 550), margin=10mm, legend=:topright)
vline!(p_dist, [mean_eta_mc], color=:black, linewidth=2, label="Mean")
vline!(p_dist, [median(mc_etas)], color=:black, linestyle=:dash, linewidth=2, label="Median")
savefig(p_dist, joinpath(RESULTS_DIR, "5_Eta_distribution_7_events.png"))
