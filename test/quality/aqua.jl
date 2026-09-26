using Aqua
using GeometricSolutions
using Test

# Two ambiguities between `==` on `TimeSeries` and GeometricBase's `==` on `AbstractVariable`:
# https://github.com/JuliaGNI/GeometricSolutions.jl/issues/33
Aqua.test_all(GeometricSolutions; ambiguities = (; broken = true))
