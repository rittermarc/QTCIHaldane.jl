function haldane(k, t2, ϕ, m)
    a::Vector{Vector{Float64}} =
    [
        [1, 0],
        [-0.5, 0.5sqrt(3)],
        [-0.5, -0.5sqrt(3)]
    ]
    b::Vector{Vector{Float64}} = [a[2] - a[3], a[3] - a[1], a[1] - a[2]]

    return 2 * t2 * cos(ϕ) * sum(cos(k' * bi) for bi in b) * pauli0 +    # NNN hopping
        sum(cos(k' * ai) * pauli[1] + sin(k' * ai) * pauli[2] for ai in a) + # NN hopping
        (m - 2 * t2 * sin(ϕ) * sum(sin(k' * bi) for bi in b)) * pauli[3]    # staggered offset
end

function crossinterpolate_chern(
    ::Type{ValueType},
    f,
    localdims::Vector{Int},
    firstpivot::TCI.MultiIndex=ones(Int, length(localdims));
    tolerance::Float64=1e-8,
    maxiter::Int=200,
    sweepstrategy::Symbol=:backandforth,
    pivottolerance::Float64=1e-12,
    normalizeerror=true,
    verbosity::Int=0,
    additionalpivots::Vector{TCI.MultiIndex}=TCI.MultiIndex[],
    evalooserror::Bool=false,
    oosindices::Vector{TCI.MultiIndex}=[rand([1, 2], length(localdims)) for _ in 1:2000],
) where {ValueType}
    tci = TCI.TensorCI1{ValueType}(f, localdims, firstpivot)
    n = length(tci)
    errors = Float64[]
    cherns = Float64[]
    ranks = Int[]
    inserrors = Float64[]
    ooserrors = Float64[]

    for pivot in additionalpivots
        println("Adding pivot $pivot")
        TCI.addglobalpivot!(tci, f, pivot, tolerance)
        println("Rank $(TCI.rank(tci))")
    end

    for iter in TCI.rank(tci)+1:maxiter
        foward_sweep = (
            sweepstrategy == :forward ||
            (sweepstrategy != :backward && isodd(iter))
        )

        if foward_sweep
            TCI.addpivot!.(tci, 1:n-1, f, pivottolerance)
        else
            TCI.addpivot!.(tci, (n-1):-1:1, f, pivottolerance)
        end

        push!(errors, TCI.lastsweeppivoterror(tci))
        push!(ranks, maximum(TCI.rank(tci)))

        if evalooserror
            tt = TCI.tensortrain(tci)
            insindices = collect(setdiff(keys(f.d), oosindices))
            push!(inserrors, maxabserror(f, tt, insindices))
            push!(ooserrors, maxabserror(f, tt, oosindices))
            push!(cherns, sumqtt(tt) / 4 / 2pi)
        end

        if verbosity > 0 && (mod(iter, 10) == 0 || last(errors) < tolerance)
            if evalooserror
                println("rank = $(last(ranks)), error = $(last(errors)), chern = $(last(cherns))")
            else
                println("rank = $(last(ranks)), error = $(last(errors))")
            end
        end

        errornormalization = normalizeerror ? tci.maxsamplevalue : 1.0
        if last(errors) < tolerance * errornormalization
            break
        end
    end

    errornormalization = normalizeerror ? tci.maxsamplevalue : 1.0
    return tci, ranks, errors ./ errornormalization, cherns, inserrors, ooserrors
end

function evaluatechern_haldane(
    deltam::Float64,
    nquantics::Int;
    t2::Float64=1e-1,
    tolerance::Float64=1e-4,
    evalooserror::Bool=false,
    datadirectory=".",
)
    phi = pi/2
    m = 3sqrt(3) * t2 + deltam
    domainboundx = [-4pi/3, 4pi/3]
    domainboundy = [-6pi/(3sqrt(3)), 6pi/(3sqrt(3))]

    ndiscretization = 2^nquantics
    kxvals = range(domainboundx..., length=ndiscretization+1)#[2:end]
    kxvals = 0.5 .* (kxvals[1:ndiscretization] .+ kxvals[2:end])
    kyvals = range(domainboundy..., length=ndiscretization+1)#[2:end]
    kyvals = 0.5 .* (kyvals[1:ndiscretization] .+ kyvals[2:end])

    f(q) = berrycurvature_quantics_dets(
        kindex -> haldane([kxvals[kindex[1]], kyvals[kindex[2]]], t2, phi, m),
        1, q, nquantics)

    localdims = fill(4, nquantics)
    cf = cachedfunc(Float64, f)
    firstpivot = TCI.optfirstpivot(cf, localdims)
    println("$firstpivot, $(cf(firstpivot))")

    cutxstep = div(ndiscretization, 4)
    cutxvals = 1:cutxstep:ndiscretization
    cutystep = div(ndiscretization, min(2^nquantics, 8192))
    oosindices = [
        index_to_quantics_fused((kxi, kyi); numdigits=nquantics)
        for kxi in cutxvals, kyi in 684:cutystep:ndiscretization
    ]

    walltime = @elapsed tci, ranks, errors, cherns, inserrors, ooserrors = crossinterpolate_chern(
        Float64,
        cf,
        localdims,
        firstpivot,
        tolerance=tolerance,
        maxiter=200,
        verbosity=1,
        pivottolerance=1e-16,
        evalooserror=evalooserror,
        oosindices=oosindices[:]
    )

    chernnumber = NaN
    if !evalooserror
        walltimeint = @elapsed chernnumber = sumqtt(TCI.tensortrain(tci)) / 4 / 2pi
        walltime += walltimeint
    else
        chernnumber = last(cherns)
    end

    println("δm = $deltam, R = $nquantics : C = $chernnumber")

    savepath::String = (
        evalooserror
        ? "$datadirectory/nq$(nquantics)_deltam$(deltam)_oos.jld2"
        : "$datadirectory/nq$(nquantics)_deltam$(deltam).jld2"
    )

    jldsave(
        savepath;
        deltam=deltam,
        nquantics=nquantics,
        cherns=cherns,
        chernnumber=chernnumber,
        tci=tci,
        ranks=ranks,
        errors=errors,
        nevals=length(cf.d),
        walltime=walltime,
        inserrors=inserrors,
        ooserrors=ooserrors
    )
end
