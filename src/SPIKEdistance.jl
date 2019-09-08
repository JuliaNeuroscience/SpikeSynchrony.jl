#=
Julia implementation of the SPIKE distance, from the description of the Scholarpedia
http://www.scholarpedia.org/article/SPIKE-distance and the paper
Kreuz et al. (2013), [Monitoring spike train synchrony](https://doi.org/10.1152/jn.00873.2012)

Implementation by Geoge Datseris. It assumes sorted input spike trains and
at the moment is only for a single pair.

I am assuming that S(t) is a continuous function... In the Python code this
is not the case. But why?!
=#
export SPIKE_distance_profile, SPIKE_distance

"""
    SPIKE_distance_profile(y1, y2) -> t, S(t)
Calculate the SPIKE distance profile function ``S(t)``, Eq. (19) of [1], given input
spike trains `y1, y2`. Return `t, S(t)`.

The keywords `t0 = min(y1[1], y2[1]) - oneunit(T)` and
`t_end = max(y1[end], y2[end]) + oneunit(T)` add auxilary spikes to both trains
at these times, as demanded by the algorithm.

[1] : Kreuz et al. (2013), [Monitoring spike train synchrony](https://doi.org/10.1152/jn.00873.2012)
"""
function SPIKE_distance_profile(y1::AbstractVector{T}, y2::AbstractVector{T};
    t0 = min(y1[1], y2[1]) - oneunit(T),
    t_end = max(y1[end], y2[end]) + oneunit(T),
    ) where {T}
    @assert issorted(y1) && issorted(y2)
    # Add auxilary spikes:
    y1 = vcat(t0, y1, t_end)
    y2 = vcat(t0, y2, t_end)
    # initialize appropriate data structures:
    tvec, indices = _initialize_spike_data(y1, y2)
    # S are the values of S(t), which is a piecewise linear function
    S = zeros(Float64, length(tvec))
    i1, i2 = 1, 1 # indices of the spikes of two spike trains
    # The following loop calculates the value of S at the k-th overall spike
    # when all spikes are sorted (values at k=1, length(t) are = 0)
    for k in 2:length(tvec)-1
        t = tvec[k]
        # find out which indices to increment (where is current spike from)
        if indices[k] == 3 # both trains share current spike
            i1 += 1; i2 += 1
        elseif indices[k] == 1 # current spike from train 1
            i1 += 1
        else # current spike from train 2
            i2 += 1
        end
        # Because of the way `tvec` and `indices` are defined,
        # this is guaranteed to be correct from the data preperation:
        tP1 = y1[i1]; tP2 = y2[i2]
        tF1 = y1[i1+1]; tF2 = y2[i2+1]
        xISI1, xISI2 = tF1 - tP1, tF2 - tP2
        xF1 = tF1 - t; xF2 = tF2 - t
        # Find Δt:
        ΔtP1 = minimum_Δt(tP1, y2, i2); ΔtP2 = minimum_Δt(tP2, y1, i1)
        ΔtF1 = minimum_Δt(tF1, y2, i2); ΔtF2 = minimum_Δt(tF2, y1, i1)
        # Compute S1, S2
        S1 = (ΔtP1*xF1 + ΔtF1*xP1)/xISI1; S2 = (ΔtP2*xF2 + ΔtF2*xP2)/xISI2
        # Compute S
        S[k] = 0.5*(S1*xISI2 + S2*xISI1)/(xISI1 + xISI2)^2
    end
    return tvec, S
end

"""
Create smart containers for the time vector (that has ALL spikes sorted) as
well as the indices of where each spike comes from. This allows for exceptionally
simplified source code for the SPIKE distance profile ``S(t)``.
"""
function _initialize_spike_data(y1::AbstractVector{T}, y2) where {T}
    tvec = Vector{T}()
    indices = Vector{Int8}()
    i1 = i2 = 1 # indices of the spike trains
    while i1 ≤ length(y1) && i2 ≤ length(y2)
        v1, v2 = y1[i1], y2[i2]
        if v1 == v2
            push!(tvec, v1); push!(indices, 3)
            i1 += 1; i2 += 1;
        elseif v1 < v2
            push!(tvec, v1); push!(indices, 1)
            i1 += 1
        else
            push!(tvec, v2); push!(indices, 2)
            i2 += 1
        end
    end
    # @assert tvec == sort!(unique(vcat(y1, y2)))
    return tvec, indices
end

"""
    minimum_Δt(t, y, i = 1) -> Δt
Return the minimum distance `Δt` between time point `t` and the spikes in `y`.
Start searching from index `i` of `y` and efficiently scan indices starting from
`i` and incrementally increasing and decreasing.
"""
function minimum_Δt(t::T, y, i = 1) where {T}
    dleft, dright = typemax(t), typemax(t)
    for j in i:-1:1 # scan from i until 1
        z = abs(t - y[j])
        z > dleft ? break : (dleft = z)
    end
    for j in i+1:1:length(y) # scan from i+1 until the end
        z = abs(t - y[j])
        z > dright ? break : (dright = z)
    end
    return min(dright, dleft)
end

"""
    SPIKE_distance(y1, y2) -> D_S
Calculate the SPIKE distance ``D_S``, which is the average of
[`SPIKE_distance_profile`](@ref).
"""
function SPIKE_distance(y1, y2; kwargs...)
    t, S = SPIKE_distance(y1, y2; kwargs...)
    return trapezoid_integral(t, S)/(t[end] - t[1])
end

function trapezoid_integral(t, S::AbstractVector{T}) where {T}
    I = zero(T)
    for i in 2:length(t)
        @inbounds I += (S[i] + S[i-1])*(t[i] - t[i-1])
    end
    return T(0.5) * I
end
