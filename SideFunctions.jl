

function U_LJ(r, U0)
    return 4*U0*(r^(-12) -(r)^-6)
end


function ∇rU_LJ(r, U0)
    return -48 * U0 *((r)^-12 - 0.5(r)^-6) / r
end

 
function U_soft(r, U0, σi, σj)
    if r>(σi+σj)
        return 0
    else
        return U0 *(1-r/(σi+σj))^2
    end
end

function ∇rU_soft(r, U0, σi, σj)
    if r>(σi+σj)
        return 0
    else
        return - 2 * U0 /(σi+σj) * (1-r/(σi+σj))
    end
end

function U(r, U0,σi, σj, name )
    if name=="LJ"
        return(U_LJ(r, U0))
    elseif name=="soft"
        return(U_soft(r, U0, σi, σj))
    end
end

function ∇rU(r, U0,σi, σj, name )
    if name=="LJ"
        return(∇rU_LJ(r, U0))
    elseif name=="soft"
        return(∇rU_soft(r, U0, σi, σj))
    end
end




function distance_periodic(ri, rj; Lx=1.0, Ly=1.0)
    #dx, dy  = abs.(ri- rj)

    dx = rj[1]-ri[1] #> Lx/2 ? rj[1]-ri[1] - Lx : rj[1]-ri[1] < -Lx/2 ? rj[1]-ri[1] + Lx : rj[1]-ri[1]
    dy = rj[2]-ri[2]

     # Apply periodic boundary conditions for x-axis
     if dx > Lx / 2
        dx -= Lx
    elseif dx < -Lx / 2
        dx += Lx
    end

    # Apply periodic boundary conditions for y-axis
    if dy > Ly / 2
        dy -= Ly
    elseif dy < -Ly / 2
        dy += Ly
    end
    
   #=  if dx>(Lx-dx)
        dx = -(Lx-dx)
    end
    if dy>(Ly-dy)
        dy = -(Ly-dy)
    end =#
    #dx = min(dx, Lx - dx)
    #dy = min(dy, Ly - dy)
    return sqrt(dx^2 + dy^2), dx, dy
end

function rescale_periodic!(positions ; limits= [1.0,1.0], dimension=2)
    for d in 1:dimension
       @. positions[d]= positions[d] - floor(positions[d] /limits[d]) * limits[d]
    end
end    


function evaluate_kinetic_energy(velocities, m)
    return(m/2*sum(velocities[1].^2+ velocities[1].^2))
end




function wave!(name,phases, t, x, y, parameters )
    if name=="planar"

    #function planar!(phases, t, x, y, parameters)
        k = parameters["k"]
        ω = parameters["ω"]
        phases.=  (1 .+ cos.(ω*t .- k .* x)) ./ 2
    elseif name=="circular"


    #function circular!(phases, t, x, y, parameters)
    k = parameters["k"]
    ω = parameters["ω"]
    x0 = parameters["x0"]
    y0 = parameters["y0"]
    phases.=  (1 .+ cos.(ω*t .- k .* sqrt.((x .- x0).^2 .+ (y .- y0).^2))) ./ 2
    end
end

function compute_packing_fraction(sizes ; Lx=1.0, Ly=1.0)
    return sum(π.* sizes.^2)/(Lx*Ly)
end


function time_string_to_minutes(time_string::AbstractString)::Float64
    days_split = split(time_string, "-")
    if length(days_split)>1
        days = parse(Int, days_split[1])
        rest = days_split[2]
    else
        days=0
        rest = days_split[1]
    end
    parts = split(rest, ":")
    
    if length(parts) == 1
        # If only ss format, assume 00 hours and 00 minutes
        hours = 0
        minutes = 0
        seconds = parse(Int, parts[1])
    elseif length(parts) == 2
        # If mm:ss format, assume 00 hours
        hours = 0
        minutes = parse(Int, parts[1])
        seconds = parse(Int, parts[2])
    elseif length(parts) == 3
        # If hh:mm:ss format
        hours = parse(Int, parts[1])
        minutes = parse(Int, parts[2])
        seconds = parse(Int, parts[3])
    else
        throw(ArgumentError("Invalid time string format. Expected 'hh:mm:ss', 'mm:ss', or 'ss'"))
    end
    
    total_minutes = days*24*60 +  hours * 60 + minutes + seconds / 60.0
    
    return total_minutes
end

# Function to calculate the area of a polygon
function polygon_area(polygon)
    n = length(polygon)
    area = 0.0
    for i in 1:n
        x1, y1 = polygon[i]
        x2, y2 = polygon[mod(i, n) + 1]
        area += x1 * y2 - x2 * y1
    end
    return abs(area) / 2.0
end
