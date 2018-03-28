abstract type TestTD <: TrialDescriptor end

@testset "Testing trials.jl" begin
    testpath = joinpath(tempname(), "Subject 23", "trialname")
    testdiffsubbasepath = joinpath(tempname(), "Participant 23", "trialname")

    @testset "Trial constructor tests" begin
        @test_nowarn Trial{TestTD}(testpath)
        @test_nowarn Trial{TestTD}(testdiffsubbasepath, "Participant")
        @test_throws DomainError Trial{TestTD}(testdiffsubbasepath)

        @test_nowarn Trial{TestTD}(23, "trialname", testpath)
        @test_nowarn Trial{TestTD}(23, "trialname", testpath, Dict(:nut => :pecan))
        @test_throws ArgumentError Trial{TestTD}(23, "trialname", testpath[4:end])
    end

    @test_throws MethodError readtrial(Trial{TestTD}(testpath))
end

