
# _periodicity(::AT, periodicity::AT) where {AT <: AbstractArray} = periodicity
# _periodicity(q::AbstractArray, ::NullPeriodicity) = zero(q)

function _periodicity(s::NamedTuple, periodicity::NamedTuple)
    NamedTuple{keys(s)}(Tuple(haskey(periodicity, k) ? periodicity[k] : NullPeriodicity() for k in keys(s)))
end


function relative_norm_error(sol, ref, p=2)
    norm(sol .- ref, p) / norm(ref, p)
end

function maximum_error(sol, ref)
    maximum(abs.(sol .- ref))
end

function relative_maximum_error(sol, ref)
    maximum_error(sol, ref) / maximum(abs.(ref))
end
