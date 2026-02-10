using Test

include("../src/HeliosProfiles.jl")
using .HeliosProfiles

fname = "BDT_rho800Ti25a1H2.15M50"
run = HeliosRun(fname)

@test length(run.t) >= 2
@test size(run.r, 1) >= 2

# Exact grid point check.
tp0 = run.t[1]
rp0 = run.r[1, 1]
vals0 = interpolate(run, tp0, rp0)

@test vals0.ni == run.ni[1, 1]
@test vals0.ne == run.ne[1, 1]
@test vals0.Ti == run.Ti[1, 1]
@test vals0.Te == run.Te[1, 1]
@test vals0.u == run.u[1, 1]
@test vals0.rho == run.rho[1, 1]

# Midpoint check in both dimensions.

tp1 = 0.5 * (run.t[1] + run.t[2])
rp1 = 0.5 * (run.r[1, 1] + run.r[2, 1])

v00 = 0.5 * (run.rho[1, 1] + run.rho[2, 1])
r0 = run.r[:, 1]
r1 = run.r[:, 2]

i0 = searchsortedlast(r1, rp1)
i1 = i0 + 1
w = (rp1 - r1[i0]) / (r1[i1] - r1[i0])
v11 = (1 - w) * run.rho[i0, 2] + w * run.rho[i1, 2]

expected_rho = 0.5 * (v00 + v11)
vals1 = interpolate(run, tp1, rp1)

@test isapprox(vals1.rho, expected_rho; rtol=1e-12, atol=0.0)
@test all(isfinite, (vals1.ni, vals1.ne, vals1.Ti, vals1.Te, vals1.u, vals1.rho))
