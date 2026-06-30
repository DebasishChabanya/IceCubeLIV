module GRB

export GRBEvent, load_grb_catalog


struct GRBEvent
    id::String
    mjd::Float64
    ra::Float64
    decl::Float64
    pos_error::Float64
    t90::Float64
    redshift::Float64
end

function load_grb_catalog(filepath::String)
    lines = readlines(filepath)
    grb_catalog = GRBEvent[]
    
    for line in lines
        # Skips the header block and completely empty lines
        if startswith(strip(line), "#") || isempty(strip(line))
            continue
        end
        
        # split() handles the uneven spaces between the 15 columns
        cols = split(line)
        
        # Extracting the 7 parameters exactly as they are
        id = cols[1]                           # Column 1: GRB_name
        ra = parse(Float64, cols[4])           # Column 4: ra
        decl = parse(Float64, cols[5])         # Column 5: decl
        pos_err = parse(Float64, cols[6])      # Column 6: pos_error
        t90 = parse(Float64, cols[7])          # Column 7: T90
        redshift = parse(Float64, cols[12])    # Column 12: redshift
        mjd = parse(Float64, cols[15])         # Column 15: mjd (T0)
        
        push!(grb_catalog, GRBEvent(id, mjd, ra, decl, pos_err, t90, redshift))
    end
    
    return grb_catalog
end

end