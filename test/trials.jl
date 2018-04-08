abstract type DSTest <: DataSource end

@testset "Testing trials.jl" begin
    testpath = joinpath(tempname(), "Subject 23", "trialname")
    testdiffsubbasepath = joinpath(tempname(), "Participant 23", "trialname")

    @testset "Trial constructor tests" begin
        @test_nowarn Trial{DSTest}(testpath)
        @test_nowarn Trial{DSTest}(testdiffsubbasepath, "Participant")
        @test_throws DomainError Trial{DSTest}(testdiffsubbasepath)

        @test_nowarn Trial{DSTest}(23, "trialname", testpath)
        @test_nowarn Trial{DSTest}(23, "trialname", testpath, Dict(:nut => :pecan))
        @test_throws ArgumentError Trial{DSTest}(23, "trialname", testpath[4:end])
    end

    @test_throws MethodError readtrial(Trial{DSTest}(testpath))
end

