pauli0 = [1. 0.; 0. 1.]
pauli = [
    [0. 1.; 1. 0.],
    [0. -1.0im; 1.0im 0.],
    [1. 0.; 0. -1.]
]

antisymmetricproduct(u, v) = u[1] * v[2] - u[2] * v[1]

function sumqtt(qtt)
    return prod(sum(T, dims=2)[:, 1, :] for T in qtt)[1]
end
function maxlinkdim(n::Integer, localdim::Integer=2)
    return 0:n-2, [min(localdim^i, localdim^(n - i)) for i in 1:(n-1)]
end

function sum_quantics_mps(mps)
    m = mps[1] * ITensor(1, siteind(mps, 1))
    for i in 2:length(mps)
        m *= mps[i] * ITensor(1, siteind(mps, i))
    end
    return scalar(m)
end

function mps_to_array(mps)
    result = Vector{Array{Float64, 3}}()
    T1 = Array(mps[1], siteind(mps, 1), linkind(mps, 1))
    push!(result, reshape(T1, 1, size(T1)...))
    for i in 2:length(mps)-1
        push!(result, Array(mps[i], linkind(mps, i-1), siteind(mps, i), linkind(mps, i)))
    end
    Tlast = Array(mps[end], linkind(mps, length(mps)-1), siteind(mps, length(mps)))
    push!(result, reshape(Tlast, size(Tlast)..., 1))
    return result
end

function maxrelerror(f, qtt::Vector{Array{Float64, 3}}, indices::Vector{Vector{Int}})
    return maximum(abs(f(i) - evaluate_qtt(qtt, i)) / abs(f(i)) for i in indices)
 end

 function maxabserror(f, qtt::Vector{Array{Float64, 3}}, indices::Vector{Vector{Int}})
     return maximum(abs(f(i) - evaluate_qtt(qtt, i)) for i in indices)
 end

function scalar(a::Matrix)
    if size(a) == (1, 1)
        return first(a)
    else
        throw(ArgumentError("$a is not a scalar."))
    end
end

function evaluate_qtt(qtt, q::Vector{<:Integer})
    return scalar(prod(T[:, i, :] for (T, i) in zip(qtt, q)))
end
