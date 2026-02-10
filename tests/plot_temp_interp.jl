using Plots

include("../src/HeliosProfiles.jl")
using .HeliosProfiles

fname = length(ARGS) >= 1 ? ARGS[1] : "BDT_rho800Ti25a1H2.15M50"
tp = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.5

run = HeliosRun(fname)

if tp < run.t[1] || tp > run.t[end]
    error("tp is out of bounds: $tp not in [$(run.t[1]), $(run.t[end])]")
end

it0 = searchsortedlast(run.t, tp)
it1 = min(it0 + 1, length(run.t))

# Use the radius grid from the earlier time slice for the interpolated curve.
rp = run.r[:, it0]
Te_interp = [interpolate(run, tp, r).Te for r in rp]

p = plot(rp, Te_interp; label="Te interp @ tp", lw=2, size=(1200, 800))
plot!(run.r[:, it0], run.Te[:, it0]; label="Te @ t=$(run.t[it0])", ls=:dash)
if it1 != it0
    plot!(run.r[:, it1], run.Te[:, it1]; label="Te @ t=$(run.t[it1])", ls=:dashdot)
end

xlabel!("r")
ylabel!("Te")
title!("Temperature vs radius")

display(p)
readline()