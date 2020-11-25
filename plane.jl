using JuLIP
using LinearAlgebra
using Printf
using ProgressMeter
using SparseArrays
using PyPlot

include("./utils.jl")

r0 = 1
R = 50;
L = 2*R #+ 5*r0
X = longhex_2d(; K0 = 0, K1 = L)
r = [ norm(X[:,i],Inf) for i=1:size(X,2)]
Idom = findall(r .< R)
X = X[:, Idom]
x = X[1,:]
y = X[2,:]
X = [X; zeros(Float64, 1, size(X,2))]

@printf("------------- constructing configuration ------------------\n")
X0 = deepcopy(X)
at = Atoms(:X, X)
r = [ norm(x, Inf) for x in positions(at) ]
Ifree = findall(r .< R-3r0)
## check attributes of `at`
lj = lennardjones(; r0=r0)
set_calculator!(at, lj)
set_pbc!(at, [false, false, false])
# set_pbc!(at, false)
mask = fill(false, 3, length(at))
mask[[1 2], Ifree] .= true
set_mask!(at, mask)

r = [ norm(x, 2) for x in positions(at) ]
@printf("------------- parameter set ------------------\n")
T = 0.01
Δt = 0.001
n = convert(Int64, T/Δt) + 1
N = length(Ifree)
stol = 1e-13

σ = 15
H = σ/4
A = 0.015
b = 0.1
rc = σ
uc = A*exp(-5)^2

ux = zeros(Float64, N)
uy = zeros(Float64, N)

vx = zeros(Float64, N)
vy = zeros(Float64, N)

ax = zeros(Float64, N)
ay = zeros(Float64, N)

uxp = zeros(Float64, N)
uyp = zeros(Float64, N)

vxp = zeros(Float64, N)
vyp = zeros(Float64, N)

axp = zeros(Float64, N)
ayp = zeros(Float64, N)

@printf("------------- initial conditions ------------------\n")
IcIdx = findall(r .≤ rc)
# initial value

uxp[IcIdx] = A/(A-uc) * (1 .+ b .*cos.(2π.*r[IcIdx]./H)) .*( A*exp.( -(r[IcIdx]./σ).^2 ) .-uc ) .* x[IcIdx]./r[IcIdx]
uyp[IcIdx] = A/(A-uc) * (1 .+ b .*cos.(2π.*r[IcIdx]./H)) .*( A*exp.( -(r[IcIdx]./σ).^2 ) .-uc ) .* y[IcIdx]./r[IcIdx]
# ux[1,1] = uy[1,1] = NaN
uxp[1] = 0.0
uyp[1] = 0.0

X[1,Ifree] = X0[1,Ifree] + uxp[:]
X[2,Ifree] = X0[2,Ifree] + uyp[:]
set_positions!(at, X)

f = forces(at) |> mat
f = f[:, Ifree]

axp[:] = f[1,:]
ayp[:] = f[2,:]

@printf("--------------- solving -------------------\n")
Ix = findall(x -> x≠0, uxp)
@assert length(Ix) < N/3
Jx = ones(Int64,size(Ix))
Vx = uxp[Ix]
Txidx = Array{Int64, 1}()
Iy = findall(x -> x≠0, uyp)
@assert length(Ix) < N/3
Jy = ones(Int64,size(Iy))
Vy = uyp[Iy]
Tyidx = Array{Int64, 1}()

nUd = 500
Udx = zeros(Float64, N, nUd)
Udy = zeros(Float64, N, nUd)
# global uxcnt = 1
# const global uycnt = 1
# sparse not help much in this case

p = Progress(n-1)
for t = 2:n
    ux[:] = uxp[:] .+ vxp[:] .* Δt .+ axp[:] .* (Δt)^2/2
    uy[:] = uyp[:] .+ vyp[:] .* Δt .+ ayp[:] .* (Δt)^2/2
    X[1,Ifree] = X0[1,Ifree] + ux[:]
    X[2,Ifree] = X0[2,Ifree] + uy[:]
    set_positions!(at, X)

    f = forces(at) |> mat
    f = f[:,Ifree]
    ax[:] = f[1,:] 
    ay[:] = f[2,:]
        
    vx[:] = vxp[:] .+ axp[:] .* (Δt)/2 .+ ax[:] .* (Δt)/2
    vy[:] = vyp[:] .+ ayp[:] .* (Δt)/2 .+ ay[:] .* (Δt)/2

    # store displacement
    Ixtmp = findall(x -> abs(x) > 1e-13, ux)
    ltx = length(Txidx) + 1
    if length(Ixtmp) < N/3
        # @printf("--------------- saving sparse x-------------------\n")
        append!(Ix, Ixtmp)
        append!(Jx, t*ones(Int64,size(Ixtmp)))
        append!(Vx, ux[Ixtmp])
    else
        # @printf("--------------- saving dense x-------------------\n")
        push!(Txidx, t)
        Udx[:, ltx] = ux
    end

    Iytmp = findall(x -> abs(x) > 1e-13, uy)
    lty = length(Tyidx) + 1
    if length(Iytmp) < N/3
        # @printf("--------------- saving sparse y-------------------\n")
        append!(Iy, Iytmp)
        append!(Jy, t*ones(Int64,size(Iytmp)))
        append!(Vy, uy[Iytmp])
    else
        # @printf("--------------- saving dense y-------------------\n")
        Udy[:, lty] = uy
        push!(Tyidx, t)
    end

    copyto!(uxp, ux)
    copyto!(uyp, uy)
    copyto!(vxp, vx)
    copyto!(vyp, vy)
    copyto!(axp, ax)
    copyto!(ayp, ay)
    
    ProgressMeter.update!(p, t)
end

Usx = sparse(Ix, Jx, Vx, N, n)
Usy = sparse(Iy, Jy, Vy, N, n)

m = 998
if m in Txidx
    ux = Udx[:, findall(x -> x==m, Txidx)]
else
    ux = Usx[:, m]
end
if m in Tyidx
    uy = Udy[:, findall(x -> x==m, Tyidx)]
else
    uy = Usy[:, m]
end
u = [norm([ux[i]; uy[i]],2) for i = 1:N]

fig = plt.figure()
scatter(x[Ifree], y[Ifree], c=u, cmap="coolwarm", s = 5)
plt.colorbar()
plt.xticks([])
plt.yticks([])
# plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
axis("equal")
# axis("off")
# plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)