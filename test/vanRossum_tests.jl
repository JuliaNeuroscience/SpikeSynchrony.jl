using Test, SpikeSynchrony

# These test are based on the analytic results of van Rossum (2001)
f = collect(1.0:10.0)

# Distance of f with f must be always 0
for τ in 10 .* rand(10)
    @test vanRossum(f, f, τ) == 0.0
    # @test vanRossum_fast(f, f, τ) == 0.0
end

# Distance of two identical vectors but the second has one extra spike,
# it is analytically known to be 1.0 (with our normalization)
# irrespectively of the length of the vectors, for any τ.: eq (7)
global g = copy(f);
push!(g, 11.0)
for τ in 10 .* rand(10)
    @test vanRossum(f, g, τ) ≈ 1.0
    # @test vanRossum_fast(f, g, 1.0) ≈ 1.0
end

# Shifting a spike from f to g produces a known result: eq. (8)
for δt in (0.25, 0.5)
    for τ in 10 .* rand(10)
        local g = copy(f);
        g[end] = f[end] + δt
        # divide by 2 due to our normalization
        @test vanRossum(f, g, τ)^2 / 2 ≈ 1.0 - exp(-δt/τ)
        # @test vanRossum_fast(f, g, τ)^2 / 2 ≈ 1.0 - exp(-δt/τ)
    end
end
