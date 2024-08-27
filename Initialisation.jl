using Distributions

macro printvar(var)
    ex = esc(var)
    varname = string(var)
    quote
        println(string($varname), " = ", $ex)
        $ex
    end
end



function random_positions(n ;  Lx=1.0, Ly=1.0)
    return [rand(n)*Lx, rand(n)*Ly]
end


function initialize_polydisperse_sizes(n, σ0, δσ)
    return rand(Normal(σ0, δσ), n)
end


function initialize_sizes_from_phases(phases, σ0, δσ, amplitude)
    n = length(phases)
    sizes =  rand(Normal(σ0, δσ), n) .+ amplitude/2 .* (1 .- cos.(phases * 2π ))
    return sizes
end

function initialize_random_phases(n)
    return(rand(n))
end

function initialize_random_phases_gaussian_positive(n, average_phase, δphase )
    phases = (rand(Normal(average_phase, δphase), n))
    for i in 1:length(phases)
        if phases[i]<0
            phases[i]=-phases[i]
        elseif phases[i]>1
            phases[i]=2-phases[i]
        end
    end

    return(phases)



end  


function initialize_uniform_phases(n)
    return(zeros(n))
end

function initialize_candidate_division(n)
    return(trues(n))
end

function position_lattice_noisy(n ; Lx=1.0, Ly=1.0, η=0.02)
    nx = Int(floor(sqrt(n))) 
    ny = n÷nx +Int(n%nx!=0)
    empty = nx*ny-n

    positions = [[],[]]
    for i in 1:n
        x = 1/(2*nx) + ((i-1) % nx)*1/nx + (rand()-0.5)*η
        y = 1/(2*ny) + ((i-1) ÷ nx)*1/ny + (rand()-0.5)*η
        push!(positions[1], x)
        push!(positions[2], y)
        
    end
rescale_periodic!(positions)
return(positions)
end




function initialize_forces(n; dimension=2)
    return [zeros(n) for _ in 1:dimension]
end

function initialize_norm_force(n)  
    return zeros(n)
end


function initialize_zero_velocities(n ; dimension=2)
    return [zeros(n) for _ in 1:dimension]
end


function initialize_velocities(n, kinetic_energy, m; dimension=2)
    velocities = [rand(Normal(0,1), n) for d in 1:dimension]
    for d in 1:dimension
        velocities[d].-= mean(velocities[d])
    end
    current_kinetic_energy = evaluate_kinetic_energy(velocities, m)

    if current_kinetic_energy!=0
        rescale = sqrt( kinetic_energy / current_kinetic_energy )
        
        for d in 1:dimension
            velocities[d] .*= rescale  
        end
    end
    return velocities
end


function initialize_verlet_list(n)
    verlet_list = [Set{Int64}() for k in 1:n]
end



function initialize_pairwise_distance(positions ; Lx=1.0, Ly=1.0)
    n = length(positions[1])
    distances = zeros(n,n)
    for i in 1:n
        for j in 1:i-1
            d, dx, dy = distance_periodic([positions[1][i],positions[2][i]],
            [positions[1][j],positions[2][j]], Lx=Lx, Ly=Ly)
            distances[i,j] = d
        end
    end
    return(distances)
end


function initialize_cell_list(positions, rx, ry ; Lx=1.0, Ly=1.0 )
    nx, ny = Int(Lx÷rx)+Int(Lx%rx>1e-8), Int(Ly÷ry)+Int(Ly%ry>1e-8)
    n = length(positions[1])
    
    cell_to_index = Dict{Tuple{Int64,Int64}, Set{Int64}}()
    index_to_cell = Vector{Tuple{Int64, Int64}}()

    # cell_to index maps the coordinates of the cell to the index of the particle inside
    for i in 1:nx
        for j in 1:ny
            cell_to_index[(i,j)] = Set()
        end
    end

    for k in 1:n
        x,y = positions[1][k], positions[2][k]
        i = Int(x÷rx)+1
        j = Int(y÷ry)+1
        push!(cell_to_index[(i,j)], k)
        push!(index_to_cell, (i,j))
    end
    return nx, ny, cell_to_index, index_to_cell


end



