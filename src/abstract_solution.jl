
abstract type AbstractSolution{dType, tType, N} end

GeometricBase.nsave(sol::AbstractSolution) = error("nsave() not implemented for ", typeof(sol))
GeometricBase.ntime(sol::AbstractSolution) = error("ntime() not implemented for ", typeof(sol))
GeometricBase.nsamples(sol::AbstractSolution) = error("nsamples() not implemented for ", typeof(sol))

GeometricBase.timesteps(sol::AbstractSolution) = error("timesteps() not implemented for ", typeof(sol))
GeometricBase.eachtimestep(sol::AbstractSolution) = error("eachtimestep() not implemented for ", typeof(sol))
GeometricBase.eachsample(sol::AbstractSolution) = 1:nsamples(sol)

counter(sol::AbstractSolution) = error("counter() not implemented for ", typeof(sol))
offset(sol::AbstractSolution) = error("offset() not implemented for ", typeof(sol))
lastentry(sol::AbstractSolution) = error("lastentry() not implemented for ", typeof(sol))

Base.write(::IOStream, sol::AbstractSolution) = error("Base.write(::IOStream, ::Solution) not implemented for ", typeof(sol))


function test_interface(sol::AbstractSolution)

    @test_nowarn nsave(sol)
    @test_nowarn ntime(sol)
    @test_nowarn nsamples(sol)

    @test_nowarn counter(sol)
    @test_nowarn offset(sol)
    @test_nowarn lastentry(sol)
    @test_nowarn timesteps(sol)

    @test_nowarn eachtimestep(sol)
    @test_nowarn eachsample(sol)
    
end
