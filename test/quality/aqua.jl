using Aqua
using GeometricSolutions
using Test

# Two ambiguities between `==` on `TimeSeries` and GeometricBase's `==` on `AbstractVariable`.
Aqua.test_all(GeometricSolutions; ambiguities = (; broken = true)) # https://github.com/JuliaGNI/GeometricSolutions.jl/issues/33
