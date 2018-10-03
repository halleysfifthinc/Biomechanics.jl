abstract type DSTest <: AbstractDataSource end
struct DSTestData <: DSTest
    events::Dict{Symbol,Any}
    data::AbstractArray
end

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

    t = Trial{DSTest}(testpath)

    @testset "Segment tests" begin
        @test_nowarn d = DSTestData(Dict{Symbol,Any}(),rand(10,10))
        @test_nowarn s = Segment(t, Dict{Symbol,Any}(:pecans => 3), d)
    end

    d = DSTestData(Dict{Symbol,Any}(),rand(10,10))
    s = Segment(t, Dict{Symbol,Any}(:pecans => 3), d)

    @testset "AnalyzedSegment tests" begin
        @test_nowarn as = AnalyzedSegment(s, Dict{Symbol,Any}())
    end

    @test_throws MethodError readtrial(Trial{DSTest}(testpath))
end

