using SafeTestsets

@safetestset "Data Series                                                                     " begin
    include("dataseries_tests.jl")
end
@safetestset "Time Series                                                                     " begin
    include("timeseries_tests.jl")
end
@safetestset "Solution                                                                        " begin
    include("solution_tests.jl")
end
@safetestset "Diagnostics                                                                     " begin
    include("diagnostics_tests.jl")
end
@safetestset "HDF5                                                                            " begin
    include("hdf5_tests.jl")
end
