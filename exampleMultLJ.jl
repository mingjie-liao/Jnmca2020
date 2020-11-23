using JuLIP
using PyPlot

include("multLJ.jl")
include("geom.jl")

using JuLIP, PyPlot

N, L = 80, 10
r0 = 0.5
ratio = 1.4
r1 = r0*ratio
σ0 = [2*r0 r0+r1; r0+r1 2*r1]./(2^(1/6))
ϵ0 = [1.0 1.0; 1.0 1.0]
calc = MultiLJ.MLJ([1 2], ϵ0, σ0)
radii = [1, 1.4]

cfg = generate_config2D(N, L, radii)
set_calculator!(cfg, calc)

cfg0 = deepcopy(cfg)
x0, y0, z0 = xyz(cfg)
x0 = mod.(x0, L)
y0 = mod.(y0, L)


SR = 1200
scatter(x0, y0, s = SR*radii[cfg.Z], marker="o", edgecolors="k")
idx = findall(x -> x == 2, cfg.Z)
scatter(x0[idx], y0[idx], s = SR*radii[cfg.Z][idx], marker="o", edgecolors="k")
axis("square")
plt.xticks([])
plt.yticks([])
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.subplots_adjust(top=1, bottom=0, left=0, right=1, hspace=0, wspace=0)
plt.margins(0, 0)
savefig("grad2dinit.pdf")
plt.clf()

minimise!(cfg; verbose=2)

x, y, z = xyz(cfg)
x = mod.(x, L)
y = mod.(y, L)
scatter(x, y, s = SR*radii[cfg.Z], marker="o", edgecolors="k")
idx = findall(x -> x == 2, cfg.Z)
scatter(x[idx], y[idx], s = SR*radii[cfg.Z][idx], marker="o", edgecolors="k")
axis("square")
plt.xticks([])
plt.yticks([])
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.subplots_adjust(top=1, bottom=0, left=0, right=1, hspace=0, wspace=0)
plt.margins(0, 0)
savefig("grad2dstab.pdf")
plt.clf()