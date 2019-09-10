using Test, SpikeSynchrony

####### ANALYTICAL TEST ########
y1 = [1]
y2 = [2]
ta = [0, 0, 1, 1, 2, 2, 3, 3]
# I calculated S by hand for these trains:
Sa = [0, 0, 5//9, 3//8, 3//8, 5//9, 0, 0]

t, S = SPIKE_distance_profile(y1, y2; tf = 3)
@test t == ta
@test S ≈ Sa    atol = 1e-15

D_S_a = sum( 1//2 * (Sa[2:end] .+ Sa[1:end-1]) .* (ta[2:end] .- ta[1:end-1]))//(ta[end] - ta[1])
D_S = SpikeSynchrony.trapezoid_integral(t, S)/(t[end] - t[1])
@test D_S ≈ D_S_a    atol = 1e-15

####### Second Analytic test ########
y1 = [2.0, 5.0, 8.0]
y2 = [1, 5, 9]
t, S  = SPIKE_distance_profile(y1, y2; t0 = 0, tf = 10)

analytic_t = [0, 1, 2, 5, 8, 9, 10]
analytic_S = [0, 0, 0.555555555555555580,  0.222222222222222210,
11/36, 25/98, 0, 0, 25/98, 11/36, 0.222222222222222210,  0.555555555555555580, 0, 0]

@test S ≈ analytic_S  atol = 1e-15

####### TEST WITH ARBITRARY REAL DATA ########
t0 = 0; tf = 500
y1 = cumsum(rand(10:10:tf, 100))
y2 = cumsum(rand(10:10:tf, 100))
t, S  = SPIKE_distance_profile(y1, y2)
t, S2 = SPIKE_distance_profile(y2, y1)
t, S0 = SPIKE_distance_profile(y1, y1)
# tr, Sr = SPIKE_distance_profile(abs.(reverse(y1) .- y1[end]), abs.(reverse(y1) .- y1[end]),

# Fundamental tests from the definition of S and D_S
@test issorted(t)
@test (0 .≤ S) == trues(length(S))
@test S == S2
@test S0 == zeros(length(S0))

D_S = SpikeSynchrony.trapezoid_integral(t, S)/(t[end]-t[1]) # == SPIKE_distance(y1, y2)
@test 0 ≤ D_S ≤ 1
