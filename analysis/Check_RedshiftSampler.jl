using Plots
using Plots.PlotMeasures 
using Statistics


include("../src/GRB_RedshiftSampler.jl")
using .GRB_RedshiftSampler


function lgrb_theory(z)
    z < 3.0 ? (1.0 + z)^2.1 : (z <= 6.0 ? 128.0 * (1.0 + z)^-1.4 : 0.0)
end

function sgrb_theory(z)
    z <= 0.9 ? exp((z - 0.9) / 0.39) : (z <= 6.0 ? exp(-(z - 0.9) / 0.26) : 0.0)
end

z_grid = range(0.0, 6.0, length=1000)
dz = z_grid[2] - z_grid[1]


lgrb_norm_pdf = lgrb_theory.(z_grid) ./ sum(lgrb_theory.(z_grid) .* dz)
sgrb_norm_pdf = sgrb_theory.(z_grid) ./ sum(sgrb_theory.(z_grid) .* dz)



println("Running Inverse CDF Sampler for 100,000 Long GRBs and 100,000 Short GRBs...")
N_samples = 100000


lgrb_samples = [GRB_RedshiftSampler.get_final_physical_redshift(-999.0, 5.0) for _ in 1:N_samples]


sgrb_samples = [GRB_RedshiftSampler.get_final_physical_redshift(-999.0, 0.5) for _ in 1:N_samples]




p1 = histogram(lgrb_samples, bins=100, normalize=:pdf, alpha=0.6, color=:steelblue, 
               label="LGRB MC Samples", title="Long GRBs (T90 > 2s)", legend=:topright)
plot!(p1, z_grid, lgrb_norm_pdf, linewidth=3, color=:red, label="W&P Theory")


p2 = histogram(sgrb_samples, bins=100, normalize=:pdf, alpha=0.6, color=:darkorange, 
               label="SGRB MC Samples", title="Short GRBs (T90 <= 2s)", legend=:topright)
plot!(p2, z_grid, sgrb_norm_pdf, linewidth=3, color=:red, label="W&P Theory")


plot(p1, p2, layout=(1, 2), size=(1000, 500), 
     xlabel="Redshift (z)", ylabel="Density (dN/dz)",
     margin=8mm, bottom_margin=12mm, left_margin=12mm) # This pushes the plot inward so text fits

savefig("InverseCDF_Visual.png")
println("Plot saved as InverseCDF_Visual.png")