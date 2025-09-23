
struct cachedfunc{ValueType}
    f::Function
    d::Dict{Vector{Int}, ValueType}

    function cachedfunc(::Type{ValueType}, f::Function) where ValueType
        new{ValueType}(f, Dict())
    end
end

function (cf::cachedfunc{ValueType})(x::Vector{Int})::ValueType where {ValueType}
    if haskey(cf.d, x)
        return cf.d[x]
    else
        val = cf.f(x)
        cf.d[deepcopy(x)] = val
        return val
    end
end

Base.broadcastable(x::cachedfunc) = Ref(x)
