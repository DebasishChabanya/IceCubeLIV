# analysis/Correlation.jl

using Statistics

# Including the files
include("../src/Neutrino.jl")
include("../src/GRB.jl")
include("../src/GRB_RedshiftSampler.jl")

function haversine_deg(ra1, dec1, ra2, dec2)
    ra1_rad, dec1_rad = deg2rad(ra1), deg2rad(dec1)
    ra2_rad, dec2_rad = deg2rad(ra2), deg2rad(dec2)
    
    d_ra = ra2_rad - ra1_rad
    d_dec = dec2_rad - dec1_rad
    
    a = sin(d_dec/2)^2 + cos(dec1_rad) * cos(dec2_rad) * sin(d_ra/2)^2
    c = 2 * asin(sqrt(min(a, 1.0)))
    
    return rad2deg(c)
end

function run_correlation()
    println("Loading the catalogs...")
    
    neutrino_path = joinpath(@__DIR__, "../data/Neutrino.npy")
    grb_path = joinpath(@__DIR__, "../data/GRB.txt")
    
    # Load from the modules
    neutrinos = Main.Neutrino.load_neutrino_catalog(neutrino_path) 
    grbs = Main.GRB.load_grb_catalog(grb_path)               
    
    # Handle missing GRB positional errors via median
    valid_errors = [g.pos_error for g in grbs if g.pos_error > 0.0 && g.pos_error != 999.0 && g.pos_error != -999.0]
    median_grb_err = median(valid_errors)
    println("Median GRB positional error calculated: $(round(median_grb_err, digits=2)) degrees")

    matches = []
    
    println("Sweeping $(length(neutrinos)) neutrinos against $(length(grbs)) GRBs...")
    
    for nu in neutrinos
        # 1. ENERGY: Convert GeV to TeV
        E_TeV = nu.energy / 1000.0
        
        if E_TeV < 60.0 || E_TeV > 500.0
            continue
        end
        
        # Convert IceCube Radians to Degrees
        nu_ra_deg = rad2deg(nu.ra)
        nu_dec_deg = rad2deg(nu.decl)
        
       # 2. NEUTRINO RESOLUTION ERROR
        recon_type = lowercase(nu.reconstruction)
        sigma_nu = (recon_type == "shower" || recon_type == "cascade") ? 15.0 : 1.0
        
        for grb in grbs
            # 3. TEMPORAL CUT (Strictly Late-Arriving)
            delta_t_days = nu.mjd - grb.mjd
            
            if delta_t_days < 0.0 || delta_t_days > 3.0
                continue
            end
            
            # 4. SPATIAL (Spherical Geometry)
            sigma_grb = (grb.pos_error > 0.0 && grb.pos_error != 999.0 && grb.pos_error != -999.0) ? grb.pos_error : median_grb_err
            
            # Use the unit-corrected neutrino degrees
            d_psi = haversine_deg(nu_ra_deg, nu_dec_deg, grb.ra, grb.decl)
            
            sigma_total = sqrt(sigma_nu^2 + sigma_grb^2)
            boundary = 3.0 * sigma_total
            
            # 5. MATCH CONFIRMATION
            if d_psi <= boundary
                push!(matches, (nu=nu, E_TeV=E_TeV, grb=grb, delta_t_days=delta_t_days, d_psi=d_psi, boundary=boundary))
            end
        end
    end
    
    println("\n=============================================")
    println("CROSS-MATCHING COMPLETE (AC2023 METHODOLOGY)")
    println("Total Neutrinos Scanned: $(length(neutrinos))")
    println("Total GRBs Scanned: $(length(grbs))")
    println("VALID MATCHES FOUND: $(length(matches))")
    println("=============================================")
    
    return matches
end

matches = run_correlation()