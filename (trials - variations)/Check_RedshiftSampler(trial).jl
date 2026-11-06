using Plots
using Plots.PlotMeasures 
using Statistics
using Random

Random.seed!(2027)


include("../src/GRB_RedshiftSampler.jl")
using .GRB_RedshiftSampler


z_grid = GRB_RedshiftSampler.Z_GRID
lgrb_norm_pdf = GRB_RedshiftSampler.LGRB_CDF.pdf
sgrb_norm_pdf = GRB_RedshiftSampler.SGRB_CDF.pdf

println("Running Inverse CDF Sampler for 100,000 Long GRBs and 100,000 Short GRBs...")
N_samples = 100000


lgrb_samples = [GRB_RedshiftSampler.get_final_physical_redshift(-999.0, 5.0) for _ in 1:N_samples]


sgrb_samples = [GRB_RedshiftSampler.get_final_physical_redshift(-999.0, 0.5) for _ in 1:N_samples]




p1 = histogram(lgrb_samples, bins=100, normalize=:pdf, alpha=0.6, color=:steelblue, 
               label="LGRB MC Samples", title="Long GRBs (T90 > 2s)", legend=:topright)
plot!(p1, z_grid, lgrb_norm_pdf, linewidth=3, color=:red, label="W&P + volume / (1+z)")


p2 = histogram(sgrb_samples, bins=100, normalize=:pdf, alpha=0.6, color=:darkorange, 
               label="SGRB MC Samples", title="Short GRBs (T90 <= 2s)", legend=:topright)
plot!(p2, z_grid, sgrb_norm_pdf, linewidth=3, color=:red, label="W&P + volume / (1+z)")


plot(p1, p2, layout=(1, 2), size=(1000, 500), 
     xlabel="Redshift (z)", ylabel="Density (dN/dz)",
     margin=8mm, bottom_margin=12mm, left_margin=12mm) # This pushes the plot inward so text fits

savefig("InverseCDF_Visual.png")
println("Plot saved as InverseCDF_Visual.png")
