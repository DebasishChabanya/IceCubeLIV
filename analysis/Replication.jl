# analysis/Replication.jl

using Plots
using Plots.PlotMeasures
using QuadGK
using Statistics

include("../src/GRB_RedshiftSampler.jl")
using .GRB_RedshiftSampler

E_TeV  = [61.7, 134.2, 98.5, 86.5, 86.1, 186.6, 66.7]
dt_sec = [73690.0, 135731.0, 15446.0, 160909.0, 187050.0, 229039.0, 23286.0]

# Proxies used in AC2023
ac2023_z = [1.38, 0.60, 1.38, 1.38, 1.38, 1.38, 1.38] 

# Physical inputs for the z sampler
measured_z_array = [1.38, -999.0, -999.0, -999.0, -999.0, -999.0, -999.0]
t90_array        = [999.0, 0.5, 5.0, 5.0, 5.0, 5.0, 5.0]

# Quantum Gravity Physics Constants (Planck 2018)
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

# Conversion factor: seconds/TeV -> Dimensionless η
const QG_CONVERSION = D1 / M_p_TeV 

function make_fit_line(m_sec_per_tev, max_x)
    x_line = range(0, max_x * 1.1, length=100)
    y_line = m_sec_per_tev .* x_line
    return x_line, y_line
end

# 1: AC2023 BASELINE (Ordinary Least Squares)
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

# 2: Instrumental error (Weighted Least Squares)
slopes_wls = dt_sec ./ K_ac2023
m_wls = mean(slopes_wls)
se_m_wls = std(slopes_wls) / sqrt(7)

eta_wls = m_wls / QG_CONVERSION
eta_wls_err = se_m_wls / QG_CONVERSION

println("--- Panel 2: Instrumental Error ---")
println("Eta = $(round(eta_wls, digits=1)) ± $(round(eta_wls_err, digits=1)) (Dimensionless)")

K_err = 0.30 .* K_ac2023
p2 = scatter(K_ac2023, dt_sec, xerror=K_err, color=:mediumblue, markersize=7, 
             label="±30% Instrumental Error", title="", legend=:topleft)
x_line, y_line = make_fit_line(m_wls, maximum(K_ac2023))
plot!(p2, x_line, y_line, color=:black, linewidth=2.5, 
      label="WLS Fit (η = $(round(eta_wls, digits=1)))")

# 3: 10,000-RUN MONTE CARLO (Sampled redshift + error)
println("--- Panel 3: Dynamic Sampling (10,000 MC Runs) ---")
N_sims = 10000
mc_etas = zeros(N_sims)

for i in 1:N_sims
    sim_z = [measured_z_array[j] != -999.0 ? measured_z_array[j] : GRB_RedshiftSampler.get_final_physical_redshift(-999.0, t90_array[j]) for j in 1:7]
    
    # CLAMP: Prevent division by zero if MC rolls exactly z=0
    sim_z = max.(sim_z, 1e-4)
    
    sim_K = [E_TeV[j] * (D_z(sim_z[j]) / D1) for j in 1:7]
    sim_slopes = dt_sec ./ sim_K
    mc_etas[i] = mean(sim_slopes) / QG_CONVERSION
end

mean_eta_mc = mean(mc_etas)
std_eta_mc = std(mc_etas)

println("Eta = $(round(mean_eta_mc, digits=1)) ± $(round(std_eta_mc, digits=1)) (Dimensionless)")

# one visual example for the scatter plot
phys_z_plot = [measured_z_array[i] != -999.0 ? measured_z_array[i] : GRB_RedshiftSampler.get_final_physical_redshift(-999.0, t90_array[i]) for i in 1:7]
phys_z_plot = max.(phys_z_plot, 1e-4) # Apply the same clamp for plotting

K_phys_plot = [E_TeV[i] * (D_z(phys_z_plot[i]) / D1) for i in 1:7]
K_phys_err_plot = 0.30 .* K_phys_plot

p3 = scatter(K_phys_plot, dt_sec, xerror=K_phys_err_plot, color=:purple, markersize=7, 
             label="1 Example MC Universe", title="", legend=:topleft)
x_line3, y_line3 = make_fit_line(mean_eta_mc * QG_CONVERSION, maximum(K_phys_plot))
plot!(p3, x_line3, y_line3, color=:black, linewidth=2.5, 
      label="10k MC Avg Fit (η = $(round(mean_eta_mc, digits=1)))")

# Save individual plots
plot_opts = (size=(600, 500), xlabel="Kinematic Factor K (TeV)", ylabel="Time Delay Δt (seconds)", 
             margin=10mm, bottom_margin=15mm, left_margin=15mm,
             guidefontsize=12, tickfontsize=10, titlefontsize=14, legendfontsize=11)

plot!(p1; plot_opts...)
savefig(p1, "1. Replication.png")

plot!(p2; plot_opts...)
savefig(p2, "2. Instrumental error included.png")

plot!(p3; ylims=(0, 250000), plot_opts...)
savefig(p3, "3. Error + Dynamic redshift (avg).png")

println("\nSuccessfully rendered! Saved 3 individual plots.")