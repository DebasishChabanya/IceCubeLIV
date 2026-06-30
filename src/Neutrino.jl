module Neutrino

export NeutrinoEvent, load_neutrino_catalog


struct NeutrinoEvent
    id::Int32
    mjd::Float64
    ra::Float64
    decl::Float64
    energy::Float64
    reconstruction::String
end

function load_neutrino_catalog(filepath::String)
    neutrino_catalog = NeutrinoEvent[]
    
    open(filepath, "r") do f
        
        seek(f, 384)
        
        
        while !eof(f)
            mjd = read(f, Float64)         # <f8 (8 bytes)
            ra = read(f, Float32)          # <f4 (4 bytes)
            decl = read(f, Float32)        # <f4 (4 bytes)
            
           
            skip(f, 64)
            
            
            recon_raw = read(f, 10)        # |S10 (10 bytes)
            recon_str = strip(String(recon_raw), ['\0', ' '])
            
            energy = read(f, Float32)      # <f4 (4 bytes)
            id = read(f, Int32)            # <i4 (4 bytes)
            
           
            skip(f, 4)
            
            
            push!(neutrino_catalog, NeutrinoEvent(id, mjd, Float64(ra), Float64(decl), Float64(energy), recon_str))
        end
    end
    
    return neutrino_catalog
end

end