module Dominance

using Combinatorics, Stella, DataFrames, FreqTables, GLM, StatsModels, StatsBase, StatsAPI, Printf, PrettyTables

export Domin, dominance

struct Domin
    nobs::Int64
    nreg::Int64
    dep::Symbol
    indeps::Vector{Any}
    covars::Vector{Symbol}
    fit_overall::Float64
    fitstat::AbstractDataFrame
    domstat::Array # cominance statistics
    comstat::Array # Complete Dominance
    constat::Array # conditional dominance
end

canonical = Dict(
    LogitLink => Bernoulli,
    InverseLink => Gamma,
    LogLink => Poisson,
    InverseSquareLink => InverseGaussian,
    NegativeBinomialLink => NegativeBinomial,
    IdentityLink => Normal
)

function dominance(_data::AbstractDataFrame,
    dep::Symbol, # dependent variable name
    indeps::Vector; # independent variables in a vector, sets are allowed in tuples
    covars=[],
    fitstat=:McFadden,
    link=nothing,
    family=nothing)

    # prepare the data set
    df = dropmissing(select(_data, vcat(dep, untuple(indeps), covars)))

    # get all combination of the indeps vector
    nvars = length(indeps)
    vvec = collect(combinations(collect(1:nvars)))
    nreg = length(vvec)
    println("A total of ", nreg, " regressions will be estimated.\n")

    # df for saving fitstats
    fs = DataFrame()
    fs.terms = Vector{Any}(undef, nreg)
    fs.terms_sorted = Vector{Any}(undef, nreg)
    fs.nterms = zeros(Int16, nreg)
    fs.r2m = Vector{Union{Missing,Float64}}(missing, nreg)
    for i in 1:length(indeps)
        fs[!, Symbol(i)] = Vector{Union{Missing,Float64}}(missing, nreg)
    end

    for (i, vindex) in enumerate(vvec)
        vars = indeps[vindex]

        fs[i, :terms] = vars
        fs[i, :terms_sorted] = sort(untuple(vars))
        fs[i, :nterms] = length(vars)

        fm = get_formula(dep, vcat(vars, covars))

        if link == nothing
            tmpreg = lm(fm, df)
            fs[i, :r2m] = r2(tmpreg)
        else
            if family == nothing
                family = canonical[link]
            end
            if !in(fitstat, [:McFadden, :Nagelkerke])
                throw(ArgumentError(fitstat, " is not allowed"))
            end
            tmpreg = glm(fm, df, family(), link())
            fs[i, :r2m] = r2(tmpreg, fitstat)
        end
    end

    # additional contribution of each indep
    for i = 1:nreg
        for j = 1:length(indeps)
            addterms = sort(unique(untuple(vcat(fs[i, :terms], indeps[j]))))
            ind = findfirst(x -> x == addterms, fs.terms_sorted)
            if i != ind
                fs[i, Symbol(j)] = fs[ind, :r2m] - fs[i, :r2m]
            end
        end
    end

    # complete dominance
    complete = zeros(Int8, nvars, nvars)
    for i = 1:nvars
        for j = (i+1):nvars
            tmpfs = dropmissing(fs[:, [Symbol(i), Symbol(j)]])
            fs1 = vcat(fs[i, :r2m], tmpfs[:, 1])
            fs2 = vcat(fs[j, :r2m], tmpfs[:, 2])
            compared = (fs1 .- fs2)
            if Base.all(compared .> 0.0)
                complete[i, j] = 1
                complete[j, i] = -1
            elseif Base.all(compared .< 0.0)
                complete[i, j] = -1
                complete[j, i] = 1
            end
        end
    end

    # conditional dominance
    conditional_dom = zeros(Float64, nvars, nvars)
    conditional_dom[:, 1] = fs.r2m[1:nvars]
    for subdf in groupby(fs, :nterms)
        for i = 1:nvars
            j = subdf[1, :nterms]
            if j < nvars
                conditional_dom[i, j+1] = mean(skipmissing(subdf[:, Symbol(i)]))
            end
        end
    end

    # general dominance
    gen_dom = mean(conditional_dom, dims=2)[:, 1]
    sta_dom = gen_dom ./ sum(gen_dom)
    ranking = ordinalrank(sta_dom, rev=true)
    domstat = hcat(gen_dom, sta_dom, ranking)

    return Dominance.Domin(
        nrow(df),
        nreg,
        dep,
        indeps,
        covars,
        fs[nreg, :r2m],
        select(fs, Not(:terms_sorted)),
        domstat,
        complete,
        conditional_dom
    )
end

function Base.show(io::IO, dom::Domin)

    println(io, "Number of observations      = ", @sprintf("%15s", dom.nobs))
    println(io, "Number of regression models = ", @sprintf("%15.0d", dom.nreg))
    println(io, "Overall Fit Statistic       = ", @sprintf("%15.4f", dom.fit_overall))

    print(io, "\n")

    nvars = length(dom.indeps)

    # define sets
    sets = []
    indepnames = []
    n = 1
    for i = 1:nvars
        if isa(dom.indeps[i], Tuple)
            push!(sets, dom.indeps[i])
            push!(indepnames, string("Set ", n))
            n += 1
        else
            push!(indepnames, string(dom.indeps[i]))
        end
    end

    if length(sets) > 0
        for i = 1:length(sets)
            println(io, "Set ", @sprintf("%4d", i), " = ", sets[i])
        end
    end
    print(io, "\n")

    println(io, "\nDominance Statistics and Ranking:")
    pretty_table(io,
        dom.domstat,
        header=["General Dominance", "Standardized Dominance", "         Ranking"],
        row_labels=indepnames,
        row_label_column_title=string(dom.dep),
        formatters=(ft_printf("%6.4f", 1:2), ft_printf("%4d", 3)),
        hlines=[0, 1, nvars + 1],
        vlines=[1]
    )

    # output complete dominance
    println(io, "\nComplete dominance:")
    pretty_table(io,
        dom.comstat,
        header=indepnames,
        row_labels=indepnames,
        row_label_column_title="dominates?",
        hlines=[0, 1, nvars + 1],
        vlines=[1]
    )

    # conditional dominance
    println(io, "\nConditional dominance:")
    pretty_table(io,
        dom.constat,
        header=collect(1:nvars),
        row_labels=indepnames,
        row_label_column_title="Variables",
        formatters=ft_printf("%6.4f", 1:nvars),
        hlines=[0, 1, nvars + 1],
        vlines=[1]
    )
    print(io, "\n")


end

function untuple(vec)
    vv = []
    for v in vec
        if isa(v, Tuple)
            vv = vcat(vv, v...)
        else
            push!(vv, v)
        end
    end
    return vv
end

function get_r2add(df, terms, indepvars)
    for i = 1:size(df, 1)
        if df[i, :terms_sorted] == indepvars && terms != indepvars
            return df[i, :r2m]
        end
    end
    return missing
end

function get_formula(dep, indeps)
    # construct another vector that untuple all variables
    if indeps == nothing
        return Term(dep) ~ term(1)
    end

    vv = untuple(indeps)

    return Term(dep) ~ sum(term.(vv))
end

end

#-----------------------------------------------------------------------------
"""
definitions from Azen and Traxel pp. 323-324
Dominance analysis considers one predictor (Xt) to completely dominate
another (Xj) if its additional contribution to every possible subset model
(that does not include Xf and Xj) is greater than that of the other
predictor. In cases where complete dominance cannot be established,
conditional or general dominance can also be used. A predictor
is said to conditionally dominate another if its average additional
contribution within each model size is greater than that of the
other predictor. Finally, one predictor generally dominates another
if its average conditional contribution over all model sizes is
greater than that of the other predictor.

https://cran.r-project.org/web/packages/domir/vignettes/domir_basics.html

linear regression - R²
logistic regression - Pseudo R² (McFadden, Nagelkerke, Estrella, Cox-Snell all produce the same dominance relationship)

Joseph Luchman's Introduction to Dominance Analysis on CRAN at https://cran.r-project.org/web/packages/domir/vignettes/domir_basics.html
"""

