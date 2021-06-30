using Jnmca2020
using JuLIP
using JuLIPMaterials
radii = 1
X = longhex_2d(; K0 = 0, K1 = 10)
S = size(X,2)
z = zeros(1,S)
X = [X;z]
x = X[1,:]
y = X[2,:]
z = X[3,:]
Mx = maximum(x)
My = maximum(y)
ϵ = 1/(2*Mx)

for i = 1:length(x)
    if y[i] > -1/2*My && y[i] <1/2*My
        x[i] -= ϵ*x[i]
    else x[i] += ϵ*x[i]
    end
end
X[1,:] = x
at = Atoms(:X,X)
set_cell!(at, [2*Mx 0 0; 0 2*My 0; 0 0 1])
set_pbc!(at, (true, true, false))
inplane!(at)
calc = lennardjones()
set_calculator!(at, calc)
minimise!(at)

#plot
fig, ax = plt.subplots()
FC = ["blue", "yellow"]
for i = 1:length(x)
    j = 1:length(y)
    k = 1:length(z)
    sp = at.Z[i]
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
savefig("grad2dinit.pdf")
