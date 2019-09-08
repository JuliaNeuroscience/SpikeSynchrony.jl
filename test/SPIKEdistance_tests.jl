using Test, SpikeSynchrony

# test trains
y1 = cumsum(rand(10:10:500, 100))
y2 = cumsum(rand(10:10:500, 100))

t, S = SPIKE_distance_profile(y1, y2)

# some basic tests:
@test (S .≥ 0) == trues(length(S))
@test issorted(t)
z = findall(isequal(0), S)

t, S = SPIKE_distance_profile(y1, y1)
@test S == zeros(length(S))


y1 = [2, 5, 8]
y2 = [1, 5, 9]
expected_times = [0.0, 1.0, 2.0, 5.0, 8.0, 9.0, 10.0]

t, S = SPIKE_distance_profile(y1, y2)
t, S2 = SPIKE_distance_profile(y2, y1)

@test S == S2
@test t == expected_times
