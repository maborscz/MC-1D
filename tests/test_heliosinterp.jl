using Test

include("../src/HeliosProfiles.jl")
using .HeliosProfiles

fname = "BDT_rho800Ti25a1H2.15M50"
run = HeliosRun(fname)

@test ndims(run.t) == 1
@test length(run.t) > 0

@test ndims(run.r) == 2
@test size(run.r, 2) == length(run.t)

@test size(run.r) == size(run.ni)
@test size(run.r) == size(run.ne)
@test size(run.r) == size(run.Ti)
@test size(run.r) == size(run.Te)
@test size(run.r) == size(run.u)
@test size(run.r) == size(run.rho)
