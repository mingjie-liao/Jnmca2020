using JuLIP
using PyPlot

include("multLJ.jl")
include("geom.jl")

using JuLIP, PyPlot

N, L = 50, 10
r0 = 1
ratio = 0.3
r1 = r0*ratio
σ0 = [2*r0 r0+r1; r0+r1 2*r1]./(2^(1/6))
ϵ0 = [1.0 1.0; 1.0 1.0]
calc = MultiLJ.MLJ([1 2], ϵ0, σ0)
radii = [1, 0.3]

cfg = generate_config2D(N, L, radii)
set_calculator!(cfg, calc)

cfg0 = deepcopy(cfg)
x0, y0, z0 = xyz(cfg)
x0 = mod.(x0, L)
y0 = mod.(y0, L)

fig, ax = plt.subplots()
FC = ["blue", "yellow"]
for i = 1:length(x)
    sp = cfg.Z[i]
    fc = FC[sp]
    r = radii[sp]
    c = matplotlib.patches.Circle([x0[i],y0[i]], r, facecolor=fc, edgecolor="black")
    ax.add_patch(c)
end
# axis([0, 10, 0, 10])
axis("equal")
plt.xticks([])
plt.yticks([])
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.subplots_adjust(top=1, bottom=0, left=0, right=1, hspace=0, wspace=0)
plt.margins(0, 0)
savefig("grad2dinit.pdf")

idx2 = findall(x -> x == 2, cfg.Z)
println(length(idx2))

minimise!(cfg; verbose=2)

x, y, z = xyz(cfg)
x = mod.(x, L)
y = mod.(y, L)

fig, ax = plt.subplots()
for i = 1:length(x)
    sp = cfg.Z[i]
    fc = FC[sp]
    r = radii[sp]
    c = matplotlib.patches.Circle([x[i],y[i]], r, facecolor=fc, edgecolor="black")
    ax.add_patch(c)
end
# axis([0, 10, 0, 10])
axis("equal")
plt.xticks([])
plt.yticks([])
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.subplots_adjust(top=1, bottom=0, left=0, right=1, hspace=0, wspace=0)
plt.margins(0, 0)
savefig("grad2dstab.pdf")