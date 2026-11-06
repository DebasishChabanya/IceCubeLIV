module GRB_RedshiftSampler

export get_final_physical_redshift, get_final_ac2023_redshift

function lgrb_pdf(z::Float64)
    if z < 0.0
        return 0.0
    elseif z < 3.0
        return (1.0 + z)^2.1
    elseif z <= 6.0

        return 128.0 * (1.0 + z)^-1.4
    else
        return 0.0
    end
end

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

const Z_GRID = collect(range(0.0, 6.0, length=10000))

function build_cdf(pdf_func)
    p = pdf_func.(Z_GRID)
    c = vcat(0.0, cumsum(diff(Z_GRID) .* (p[1:end-1] .+ p[2:end]) ./ 2))
    return (pdf=p ./ c[end], cdf=c ./ c[end])
end

using QuadGK: quadgk

function volume_time_weight(z)
    h(x) = sqrt(0.315 * (1 + x)^3 + 0.685)
    distance, _ = quadgk(x -> 1 / h(x), 0.0, z; rtol=1e-10)
    return distance^2 / (h(z) * (1 + z))
end

const LGRB_CDF = build_cdf(z -> lgrb_pdf(z) * volume_time_weight(z))
const SGRB_CDF = build_cdf(z -> sgrb_pdf(z) * volume_time_weight(z))

function invert_cdf(U::Float64, distribution)
    p, c = distribution.pdf, distribution.cdf
    U <= 0.0 && return Z_GRID[1]
    U >= 1.0 && return Z_GRID[end]
    i = clamp(searchsortedlast(c, U), 1, length(c) - 1)
    h = Z_GRID[i+1] - Z_GRID[i]
    slope = (p[i+1] - p[i]) / h
    area = U - c[i]
    area == 0.0 && return Z_GRID[i]
    offset = 2 * area / (p[i] + sqrt(max(0.0, p[i]^2 + 2 * slope * area)))
    return Z_GRID[i] + clamp(offset, 0.0, h)
end

function is_long_grb(t90::Float64)
    if isnan(t90) || t90 < 0.0
        return true
    elseif t90 > 2.0
        return true
    else
        return false
    end
end

function sample_physical_redshift(t90::Float64)
    U = rand()
    if is_long_grb(t90)
        return invert_cdf(U, LGRB_CDF)
    else
        return invert_cdf(U, SGRB_CDF)
    end
end

function sample_ac2023_proxy_redshift(t90::Float64)
    if is_long_grb(t90)
        return 1.38
    else
        return 0.60
    end
end

function get_final_physical_redshift(measured_z::Float64, t90::Float64)

    if isnan(measured_z) || measured_z < 0.0
        return sample_physical_redshift(t90)
    else

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

