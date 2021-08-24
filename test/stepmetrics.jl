using Statistics, Biomechanics.SpatiotemporalMetrics, Random
using Biomechanics: interleave

Random.seed!(0xBEEF)

stride = 1.3
lstep = 0.6*stride
rstep = stride - lstep

rswing_σ = randn(100) .* .001;
lswing_σ = randn(100) .* .001;
rds_σ = randn(100) .* .001;
lds_σ = randn(100) .* .001;

@testset "Walking" begin
    lstance = 0.70*stride
    rstance = 0.75*stride
    lswing = stride - lstance

    rswing = stride - rstance
    lswings = lswing .+ lswing_σ;
    rswings = rswing .+ rswing_σ;
    ldss = lstep .- lswings .+ lds_σ;
    rdss = rstep .- rswings .+ rds_σ;

    lsteps = lswings + ldss;
    rsteps = rswings + rdss;

    strides = lsteps + rsteps;

    ge = cumsum([0; interleave(ldss, lswings, rdss, rswings)]);

    rfs = ge[1:4:end];
    lfo = ge[2:4:end];
    lfs = ge[3:4:end];
    rfo = ge[4:4:end];

    @test stridetimes(rfs) ≈ strides

    @test collect(steptimes(;lfs, rfs)) ≈ [lsteps, rsteps]

    @test collect(swingstance(rfs, rfo; normalize=false)) ≈ [rswings, strides - rswings]

    @test swing(rfs, rfo; normalize=false) ≈ rswings
    @test stance(rfs, rfo; normalize=false) ≈ strides - rswings

    @test hcat(swingstance(rfs, rfo)...) ≈ [rswings strides - rswings] ./ strides

    @test swing(rfs, rfo) ≈ rswings ./ strides
    @test stance(rfs, rfo) ≈ (strides - rswings) ./ strides

    @test singlesupport(;lfs, lfo, rfs, rfo, normalize=false) ≈ (lswings + rswings)[1:99]
    @test doublesupport(;lfs, lfo, rfs, rfo, normalize=false) ≈ (strides - lswings - rswings)[1:100]

    @test singlesupport(;lfs, lfo, rfs, rfo, normalize=true) ≈
        (lswings + rswings)[1:99] ./ strides[1:99]
    @test doublesupport(;lfs, lfo, rfs, rfo, normalize=true) ≈
        (strides - lswings - rswings)[1:100] ./ strides[1:100]
end

@testset "Running" begin
    lstance = 0.30*stride
    rstance = 0.35*stride
    lswing = stride - lstance
    rswing = stride - rstance
    lstances = lstance .+ lswing_σ;
    rstances = rstance .+ rswing_σ;
    lfloats = lstep .- lstances .+ lds_σ;
    rfloats = rstep .- rstances .+ rds_σ;

    lsteps = rstances + lfloats;
    rsteps = lstances + rfloats;

    strides = lsteps + rsteps;

    ge = cumsum([0; interleave(rstances, lfloats, lstances, rfloats)]);

    rfs = ge[1:4:end];
    rfo = ge[2:4:end];
    lfs = ge[3:4:end];
    lfo = ge[4:4:end];

    @test stridetimes(rfs) ≈ strides

    @test collect(steptimes(;lfs, rfs)) ≈ [lsteps, rsteps]

    @test collect(swingstance(rfs, rfo; normalize=false)) ≈ [strides - rstances, rstances]

    @test swing(rfs, rfo; normalize=false) ≈ strides - rstances
    @test stance(rfs, rfo; normalize=false) ≈ rstances

    @test hcat(swingstance(rfs, rfo)...) ≈ [strides - rstances rstances] ./ strides

    @test swing(rfs, rfo) ≈ (strides - rstances) ./ strides
    @test stance(rfs, rfo) ≈ rstances ./ strides

    @test singlesupport(;lfs, lfo, rfs, rfo, normalize=false) ≈ (lstances + rstances)[1:99]
    @test (doublesupport(;lfs, lfo, rfs, rfo, normalize=false) .== 0) |> prod

    @test singlesupport(;lfs, lfo, rfs, rfo, normalize=true) ≈ (lstances + rstances)[1:99] ./ strides[1:99]
    @test (doublesupport(;lfs, lfo, rfs, rfo, normalize=true) .== 0) |> prod
end

