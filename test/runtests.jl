using SafeTestsets

const GROUPS = isempty(ARGS) ? ["core", "slow"] : ARGS

if "core" in GROUPS
    @safetestset "Aqua" include("quality/aqua.jl")
    @safetestset "Data Series" include("dataseries.jl")
    @safetestset "Time Series" include("timeseries.jl")
    @safetestset "Solution" include("solutions.jl")
    @safetestset "Diagnostics" include("diagnostics.jl")
    @safetestset "HDF5" include("hdf5.jl")
end
