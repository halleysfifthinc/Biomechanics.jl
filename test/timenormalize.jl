@testset "Testing normalizing functions" begin
    @testset "normtime" begin
        nt1 = rand()

        normed1 = normtime(nt1,nt1+1)
        normed2 = normtime(nt1, nt1+1, 65)

        @test length(normed1) === 100
        @test length(normed2) === 65
        @test first(normed1) === nt1
        @test first(normed2) === nt1
        @test last(normed1) < (nt1 + 1)
        @test last(normed1) < (nt1 + 1)
    end

    @testset "normalizeevents" begin
        @test prod(normalizeevents(1, 2, 1.5) .=== [ 51 ])
        @test prod(normalizeevents(1, 2, [1.25, 1.75]) .=== [26, 76])
        @test prod(normalizeevents(1, 2, [1.2, 1.8], 10) .=== [3, 9])
    end

    @testset "timenormalize" begin
        t200 = [ normtime(0, 2*pi, 200); normtime(2*pi, 2*2*pi, 200) ]
        t100 = [ normtime(0, 2*pi, 100); normtime(2*pi, 2*2*pi, 100) ]
        
        x200 = sin.(t200)
        x100 = sin.(t100)

        normx200 = timenormalize(x200, [1, 201, 401], 100)

        @test prod(isapprox.(x100, normx200; atol=1e-10))
    end
end

