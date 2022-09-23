
abstract type AbstractSolution{dType, tType} end

GeometricBase.ntime(sol::AbstractSolution) = error("ntime() not implemented for ", typeof(sol))

GeometricBase.timesteps(sol::AbstractSolution) = error("timesteps() not implemented for ", typeof(sol))
GeometricBase.eachtimestep(sol::AbstractSolution) = error("eachtimestep() not implemented for ", typeof(sol))

Base.step(sol::AbstractSolution) = error("step() not implemented for ", typeof(sol))
nstore(sol::AbstractSolution) = error("nstore() not implemented for ", typeof(sol))

Base.write(::IOStream, sol::AbstractSolution) = error("Base.write(::IOStream, ::Solution) not implemented for ", typeof(sol))


function test_interface(sol::AbstractSolution)

    @test_nowarn step(sol)
    @test_nowarn ntime(sol)
    @test_nowarn nstore(sol)

    # @test_nowarn lastentry(sol)
    @test_nowarn timesteps(sol)
    @test_nowarn eachtimestep(sol)
    
end
