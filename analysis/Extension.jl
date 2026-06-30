# analysis/Extension.jl

using Plots
using Plots.PlotMeasures
using QuadGK
using Statistics

include("../src/GRB_RedshiftSampler.jl")
using .GRB_RedshiftSampler

# Extracted data for the 13 unique valid matches (Neutrino 107 double-counting resolved)
E_TeV  = [79.76, 78.45, 64.05, 168.32, 70.79, 136.87, 62.98, 83.37, 261.23, 138.47, 191.39, 133.85, 95.12]
dt_days = [2.14, 1.86, 0.85, 2.65, 1.25, 1.57, 0.27, 2.50, 1.35, 1.00, 0.72, 2.88, 1.40]
dt_sec = dt_days .* 86400.0 # Convert days to seconds

# 2 measured redshifts, 11 missing (-999.0)
measured_z_array = [-999.0, -999.0, 1.38, -999.0, -999.0, -999.0, -999.0, -999.0, 0.645, -999.0, -999.0, -999.0, -999.0]

# T90 parameters (2 short GRBs <= 2.0s, 11 long GRBs > 2.0s or missing)
t90_array = [38.66, 33.28, 25.40, 90.50, 4.61, 0.13, 2.43, 0.10, 36.10, 127.75, 25.86, 70.40, 24.58]

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
const QG_CONVERSION = D1 / M_p_TeV 

function make_fit_line(m_sec_per_tev, max_x)
    x_line = range(0, max_x * 1.1, length=100)
    y_line = m_sec_per_tev .* x_line
    return x_line, y_line
end

# 10,000-RUN MONTE CARLO (Sampled redshift + error)
println("--- Extended dataset: Dynamic sampling (10,000 MC Runs) ---")
N_sims = 10000
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
    sim_z = max.(sim_z, 1e-4)
    
    sim_K = [E_TeV[j] * (D_z(sim_z[j]) / D1) for j in 1:13]
    
    # Weighted Least Squares formulation
    sim_slopes = dt_sec ./ sim_K
    mc_etas[i] = mean(sim_slopes) / QG_CONVERSION
end

mean_eta_mc = mean(mc_etas)
std_eta_mc = std(mc_etas)

println("Extended eta = $(round(mean_eta_mc, digits=1)) ± $(round(std_eta_mc, digits=1)) (Dimensionless)")

# Generate visual example for the scatter plot
phys_z_plot = zeros(13)
for i in 1:13
    if measured_z_array[i] != -999.0
        phys_z_plot[i] = measured_z_array[i]
    else
        phys_z_plot[i] = GRB_RedshiftSampler.get_final_physical_redshift(-999.0, t90_array[i])
    end
end
phys_z_plot = max.(phys_z_plot, 1e-4) 

K_phys_plot = [E_TeV[i] * (D_z(phys_z_plot[i]) / D1) for i in 1:13]
K_phys_err_plot = 0.30 .* K_phys_plot # 30% instrumental error

# Plot Generation
p1 = scatter(K_phys_plot, dt_sec, xerror=K_phys_err_plot, color=:darkgreen, markersize=6, 
             label="13-Event MC", legend=:bottomright)
x_line, y_line = make_fit_line(mean_eta_mc * QG_CONVERSION, maximum(K_phys_plot))

plot!(p1, x_line, y_line, color=:black, linewidth=2.5, 
      label="10k MC Avg Fit (η = $(round(mean_eta_mc, digits=1)))")

plot_opts = (size=(650, 500), xlabel="Kinematic Factor K (TeV)", ylabel="Time Delay Δt (seconds)", 
             margin=10mm, bottom_margin=15mm, left_margin=15mm,
             guidefontsize=12, tickfontsize=10, legendfontsize=11,
             ylims=(0, 260000)) # Clamped to keep the 3-day data points visible

plot!(p1; plot_opts...)
savefig(p1, "4. Extended_13_Events.png")

println("\nSuccessfully rendered! Saved Extended_13_Events.png")