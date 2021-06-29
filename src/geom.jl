export generate_config, generate_config2D, writeGranu, plot2d

function generate_config(N, L, radii)
    X = L * rand(JVecF, N)
    Z = rand(1:length(radii), N)
    at = Atoms(:X, X)
    at.Z = Z 
    set_cell!(at, [L 0 0; 0 L 0; 0 0 L])
    set_pbc!(at, true)
    return at
end

function generate_config2D(N, L, radii)
    X = L * rand(3, N)
    X[3, :] .= 0.0
    Z = rand(1:length(radii), N)
    at = Atoms(:X, X)
    at.Z = Z 
    set_cell!(at, [L 0 0; 0 L 0; 0 0 1])
    set_pbc!(at, (true, true, false))
    inplane!(at)
    return at
end

function writeGranu(filename, config, L)
    dl = " "
    nu = length(config)
    L = convert(Float64, L)
    open(filename, "w") do f
        println(f, "ITEM: TIMESTEP")
        println(f, 1)
        println(f, "ITEM: NUMBER OF ATOMS")
        println(f, nu)
        println(f, "ITEM: BOX BOUNDS pp pp pp")
        println(f, 0.00, dl, L)
        println(f, 0.00, dl, L)
        println(f, 0.00, dl, L)
        println(f, "ITEM: ATOMS id type x y z fx fy fx radius")
        for i = 1:nu
            u = config.X[i]
            r = config.M[i]
            z = config.Z[i]
            println(f, i, dl, z, dl, u[1], dl, u[2], dl, u[3], dl, 0.00, dl, 0.00, dl, 0.00, dl, r)
        end
    end
end

using JuLIP, PyPlot

function plot2d(at; ttl = nothing)
    x, y, _ = xyz(at)
    N = length(at.Z)
    for i = 1:N
        z = at.Z[i];
        if z == 1
            color = "b.";
        else
            color = "r.";
        end
        plot(mod(x[i],10), mod(y[i],10), color, markersize=20) #set to be mod(x, Lx), mod(y, Ly)
    end
    if !isnothing(ttl)
        title(ttl)
    end 
    PyPlot.draw()
    PyPlot.pause(0.0001)    
end 
