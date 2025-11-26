
# _periodicity(::AT, periodicity::AT) where {AT <: AbstractArray} = periodicity
# _periodicity(q::AbstractArray, ::NullPeriodicity) = zero(q)

function _periodicity(s::NamedTuple, periodicity::NamedTuple)
    NamedTuple{keys(s)}(Tuple(haskey(periodicity, k) ? periodicity[k] : NullPeriodicity()
    for k in keys(s)))
end
