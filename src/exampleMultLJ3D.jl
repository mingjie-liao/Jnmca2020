using Jnmca2020
using JuLIP
using PyPlot
import Jnmca2020.MultiLJ: MLJ

N, L =650, 10
r0 = 0.5
ratio = 1.4
r1 = r0*ratio
σ0 = [2*r0 r0+r1; r0+r1 2*r1]./(2^(1/6))
ϵ0 = [1.0 1.0; 1.0 1.0]
calc = MLJ([1 2], ϵ0, σ0)
radii = [1, 1.4]

cfg = generate_config(N, L, radii)
set_calculator!(cfg, calc)

cfg0 = deepcopy(cfg)
x0, y0, z0 = xyz(cfg)
x0 = mod.(x0, L)
y0 = mod.(y0, L)
z0 = mod.(z0, L)
idx = findall(x -> x == 2, cfg.Z)
cfg.M[:] .= radii[1]
cfg.M[idx] .= radii[2]
writeGranu("x0.dump", cfg, L)

minimise!(cfg; verbose=2, g_calls_limit = 2000)
# minimise!(cfg; verbose=2, precond = FF(cfg; innerstab = 0.1, stab = 0.001))
x, y, z = xyz(cfg)
x = mod.(x, L)
y = mod.(y, L)
z = mod.(z, L)
set_positions!(cfg, [x y z]')
writeGranu("x1.dump", cfg, L)

