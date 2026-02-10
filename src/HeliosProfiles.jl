"""
HeliosProfiles

Utilities for loading Helios HDF5 runs and interpolating radial profiles
in time and radius.
"""
module HeliosProfiles

using HDF5

export HeliosRun, interpolate

# Directory to look for Helios data files, relative to *this* Julia source file.
# Put your HDF5 files in:  data/helios-runs/
const DATA_DIR = normpath(joinpath(@__DIR__, "..", "data", "helios-runs"))


"""
        HeliosRun(fname::AbstractString)

Load a Helios HDF5 run from `data/helios-runs/`. The `fname` argument is the
file stem (without the `.h5` extension).

### Fields
- `t`       -- `N` element vector of times [s].
- `r`:      -- `M x N` matrix of radial positions for each time index [cm]
- `ni`:     -- `M x N` matrix of ion number density [cm^-3]
- `ne`:     -- `M x N` matrix of electron number density [cm^-3]
- `Ti`:     -- `M x N` matrix of ion temperature [eV]
- `Te`:     -- `M x N` matrix of electron temperature [eV]
- `u`:      -- `M x N` matrix of radial velocity [cm/s]
- `rho`:    -- `M x N` matrix of mass density [g/cm^3]
"""
struct HeliosRun
    t::Vector{Float64} # (N,)
    r::Matrix{Float64} # (M,N)
    ni::Matrix{Float64}
    ne::Matrix{Float64}
    Ti::Matrix{Float64}
    Te::Matrix{Float64}
    u::Matrix{Float64}
    rho::Matrix{Float64}

    function HeliosRun(fname::AbstractString)
        required = ("t", "r", "ni", "ne", "Ti", "Te", "u", "rho")
        
        fpath = joinpath(DATA_DIR, fname * ".h5")
        isfile(fpath) || error("HDF5 file not found: $fpath. Please ensure the file exists in the data/helios-runs/ directory.")

        h5open(fpath, "r") do f
            for key in required
                haskey(f, key) || error("Missing dataset '$key' in HDF5 file: $fpath")
            end

            t   = vec(Float64.(read(f["t"])))
            r   = Matrix{Float64}(Float64.(read(f["r"])))
            ni  = Matrix{Float64}(Float64.(read(f["ni"])))
            ne  = Matrix{Float64}(Float64.(read(f["ne"])))
            Ti  = Matrix{Float64}(Float64.(read(f["Ti"])))
            Te  = Matrix{Float64}(Float64.(read(f["Te"])))
            u   = Matrix{Float64}(Float64.(read(f["u"])))
            rho = Matrix{Float64}(Float64.(read(f["rho"])))

            run = new(t, r, ni, ne, Ti, Te, u, rho)
            return run
        end
    end
end

"""
    interpolate(run::HeliosRun, tp::Float64, rp::Float64)

Linearly interpolate all fields in time and radius and return a NamedTuple
with values at the requested point. The radius grid is time-dependent, so
the radius bracket is computed separately at each neighboring time slice.

### Input
- `run`     -- a `HeliosRun` object containing the data to interpolate.
- `tp`      -- the time point to interpolate at [s].
- `rp`      -- the radius point to interpolate at [cm].

### Output
A `NamedTuple` with the following fields:
- `ni`      -- interpolated ion number density [cm^-3]
- `ne`      -- interpolated electron number density [cm^-3]
- `Ti`      -- interpolated ion temperature [eV]
- `Te`      -- interpolated electron temperature [eV]
- `u`       -- interpolated radial velocity [cm/s]
- `rho`     -- interpolated mass density [g/cm^3]
"""
function interpolate(run::HeliosRun, tp::Float64, rp::Float64)
    it0, it1, wt = _bracket(run.t, tp, "t")
    ir00, ir01, wr0 = _bracket(run.r[:, it0], rp, "r")
    ir10, ir11, wr1 = _bracket(run.r[:, it1], rp, "r")

    @inline function interp_field(field::Matrix{Float64})
        y00 = field[ir00, it0]
        y01 = field[ir01, it0]
        y10 = field[ir10, it1]
        y11 = field[ir11, it1]

        v0 = (1 - wr0) * y00 + wr0 * y01
        v1 = (1 - wr1) * y10 + wr1 * y11
        return (1 - wt) * v0 + wt * v1
    end

    return (
        ni = interp_field(run.ni),
        ne = interp_field(run.ne),
        Ti = interp_field(run.Ti),
        Te = interp_field(run.Te),
        u = interp_field(run.u),
        rho = interp_field(run.rho),
    )
end

"""
    _bracket(x::Vector{Float64}, xp::Float64, label::AbstractString)

Return the index pair `(i0, i1)` that brackets `xp` in the sorted vector `x`
and the linear weight `w` in [0, 1].
"""
function _bracket(x::Vector{Float64}, xp::Float64, label::AbstractString)
    if xp < x[1] || xp > x[end]
        error("$label is out of bounds: $xp not in [$(x[1]), $(x[end])]")
    end
    if xp == x[end]
        return length(x), length(x), 0.0
    end
    i0 = searchsortedlast(x, xp)
    i1 = i0 + 1
    dx = x[i1] - x[i0]
    dx == 0 && error("$label grid has duplicate values at indices $i0 and $i1")
    w = (xp - x[i0]) / dx
    return i0, i1, w
end

end