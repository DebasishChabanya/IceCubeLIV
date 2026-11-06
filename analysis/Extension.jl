# analysis/Extension.jl

using Plots
using Plots.PlotMeasures
using QuadGK
using Statistics
using Random

const ENERGY_FRACTION = 0.30
const Z_MIN = 1e-4
const MC_SEED = 2026
const RESULTS_DIR = joinpath(@__DIR__, "..", "results")
mkpath(RESULTS_DIR)

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

# Correlation.jl loads the catalogue readers and redshift sampler, then creates matches.
include("Correlation.jl")
using .GRB_RedshiftSampler

# Keeping one GRB per neutrino. If a neutrino has several candidate GRBs, uses the one with the smallest angular separation.
selected_by_neutrino = Dict{Int32, Any}()
for m in matches
    id = m.nu.id
    if !haskey(selected_by_neutrino, id) ||
       m.d_psi < selected_by_neutrino[id].d_psi ||
       (m.d_psi == selected_by_neutrino[id].d_psi && m.grb.id < selected_by_neutrino[id].grb.id)
        selected_by_neutrino[id] = m
    end
end

selected_matches = sort(collect(values(selected_by_neutrino)), by=m -> m.nu.id)
isempty(selected_matches) && error("No GRB-neutrino associations passed the matching cuts")

N_events = length(selected_matches)
E_TeV = Float64[m.E_TeV for m in selected_matches]
dt_days = Float64[m.delta_t_days for m in selected_matches]
dt_sec = dt_days .* 86400.0
measured_z_array = Float64[m.grb.redshift for m in selected_matches]
t90_array = Float64[m.grb.t90 for m in selected_matches]

println("Selected $(N_events) neutrinos from $(length(matches)) candidate pairs")
for m in selected_matches
    println("Neutrino $(m.nu.id): $(m.grb.id), E=$(round(m.E_TeV, digits=3)) TeV, delay=$(round(m.delta_t_days, digits=4)) days, separation=$(round(m.d_psi, digits=3)) deg")
end

open(joinpath(RESULTS_DIR, "extension_selected_events.csv"), "w") do io
    println(io, "neutrino_id,grb_id,energy_TeV,delay_days,redshift,t90,separation_deg")
    for m in selected_matches
        println(io, join((m.nu.id, m.grb.id, m.E_TeV, m.delta_t_days,
                          m.grb.redshift, m.grb.t90, m.d_psi), ","))
    end
end

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

println("--- Selected events: sampling (10,000 MC runs) ---")
N_sims = 10000
Random.seed!(MC_SEED)
clamp_values = sort(unique([1e-6, 1e-4, 1e-3, 1e-2, Z_MIN])) #putting different clamps to see robustness as it changes quite substancially
eta_by_clamp = zeros(N_sims, length(clamp_values))
mc_etas = zeros(N_sims)

for i in 1:N_sims
    sim_z = [GRB_RedshiftSampler.get_final_physical_redshift(measured_z_array[j], t90_array[j])
             for j in 1:N_events]

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
println("Extended eta = $(round(mean_eta_mc, digits=1)) ± $(round(std_eta_mc, digits=1))")

Random.seed!(2027)
phys_z_plot = [GRB_RedshiftSampler.get_final_physical_redshift(measured_z_array[i], t90_array[i])
               for i in 1:N_events]
phys_z_plot = max.(phys_z_plot, Z_MIN)

K_phys_plot = [E_TeV[i] * (D_z(phys_z_plot[i]) / D1) for i in 1:N_events]
K_phys_err_plot = ENERGY_FRACTION .* K_phys_plot

p1 = scatter(K_phys_plot, dt_sec, xerror=K_phys_err_plot, color=:darkgreen, markersize=6,
             label="$(N_events) selected candidates", legend=:bottomright)
x_line, y_line = make_fit_line(mean_eta_mc * QG_CONVERSION, maximum(K_phys_plot))
plot!(p1, x_line, y_line, color=:black, linewidth=2.5,
      label="MC mean (η = $(round(mean_eta_mc, digits=1)))")

x_med, y_med = make_fit_line(median(mc_etas) * QG_CONVERSION, maximum(K_phys_plot))
plot!(p1, x_med, y_med, color=:black, linestyle=:dash, linewidth=2,
      label="MC median (η = $(round(median(mc_etas), digits=1)))")

plot_opts = (size=(800, 650), xlabel="Kinematic Factor K (TeV)", ylabel="Time Delay Δt (seconds)",
             margin=10mm, bottom_margin=15mm, left_margin=15mm,
             guidefontsize=12, tickfontsize=10, legendfontsize=9,
             ylims=(0, 260000))

plot!(p1; plot_opts...)
savefig(p1, joinpath(RESULTS_DIR, "4_Extended_$(N_events)_Events.png"))
println("Saved 4_Extended_$(N_events)_Events.png")

save_summary("extension_summary.csv", [("redshift_and_energy", mc_etas)])

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
savefig(p_dist, joinpath(RESULTS_DIR, "6_Eta_distribution_$(N_events)_events.png"))
