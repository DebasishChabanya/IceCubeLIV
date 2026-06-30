module GRB_RedshiftSampler

export get_final_physical_redshift, get_final_ac2023_redshift


# LGRB model (Wanderman & Piran 2015) - Collapsars

function lgrb_pdf(z::Float64)
    if z < 0.0
        return 0.0
    elseif z < 3.0
        return (1.0 + z)^2.1
    elseif z <= 6.0
        # Added the 128.0 continuity constant (4^3.5)
        return 128.0 * (1.0 + z)^-1.4
    else
        return 0.0
    end
end

# SGRB exponential model (Wanderman & Piran 2015) - Neutron Star Mergers

function sgrb_pdf(z::Float64)
    if z < 0.0
        return 0.0
    elseif z <= 0.9
        return exp((z - 0.9) / 0.39)
    elseif z <= 6.0
        return exp(-(z - 0.9) / 0.26)
    else
        return 0.0
    end
end

# 2. INVERSE CDF

# Building the math grid from 0 to 6

const Z_GRID = collect(range(0.0, 6.0, length=10000))

# Integrating the PDF to find the normalized CDF

function build_cdf(pdf_func)
    pdf_vals = pdf_func.(Z_GRID)
    cdf_vals = cumsum(pdf_vals) 
    cdf_normalized = cdf_vals ./ cdf_vals[end] # Forcing the max to exactly 1.0
    return cdf_normalized
end

const LGRB_CDF = build_cdf(lgrb_pdf)
const SGRB_CDF = build_cdf(sgrb_pdf)

# Maping the uniform random draw 

function invert_cdf(U::Float64, cdf_array)
    idx = searchsortedfirst(cdf_array, U)
    
    # Grid boundaries

    if idx == 1
        return Z_GRID[1]
    elseif idx > length(cdf_array)
        return Z_GRID[end]
    end
    
    
    U_low = cdf_array[idx-1]
    U_high = cdf_array[idx]
    z_low = Z_GRID[idx-1]
    z_high = Z_GRID[idx]
    
    fraction = (U - U_low) / (U_high - U_low)
    return z_low + fraction * (z_high - z_low)
end

# 3. T90 LOGIC

# Check duration. Missing data defaults to Long GRB since they are the majority.

function is_long_grb(t90::Float64)
    if isnan(t90) || t90 < 0.0
        return true
    elseif t90 > 2.0
        return true
    else
        return false
    end
end

# 4. continuous

function sample_physical_redshift(t90::Float64)
    U = rand()
    if is_long_grb(t90)
        return invert_cdf(U, LGRB_CDF)
    else
        return invert_cdf(U, SGRB_CDF)
    end
end

# A flat proxy 

function sample_ac2023_proxy_redshift(t90::Float64)
    if is_long_grb(t90)
        return 1.38 # AC2023 arbitrary Long GRB distance
    else
        return 0.60 # AC2023 arbitrary Short GRB distance
    end
end


function get_final_physical_redshift(measured_z::Float64, t90::Float64)

    # Catch missing catalog flags like -999.0 or NaN

    if isnan(measured_z) || measured_z < 0.0
        return sample_physical_redshift(t90)
    else

        # We have real telescope data, keep it

        return measured_z
    end
end


function get_final_ac2023_redshift(measured_z::Float64, t90::Float64)
    if isnan(measured_z) || measured_z < 0.0
        return sample_ac2023_proxy_redshift(t90)
    else
        return measured_z
    end
end

end