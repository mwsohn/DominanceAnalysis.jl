module Dominance

using Combinatorics, Stella, DataFrames, StatsModels, StatsBase, StatsAPI, Printf, PrettyTables, Statistics, GLM

export Domin, dominance, dominance_designations

struct Domin
    nobs::Int64
    nreg::Int64
    dep::Symbol
    indeps::Vector{Any}
    covars::Vector{Symbol}
    fit_overall::Float64
    fitstat::AbstractDataFrame
    domstat::Array # General dominance statistics
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

"""
    dominance(df::AbstractDataFrame, dep::Symbol, indeps::Vector; 
        covars = [], fitstat = :McFadden, link = nothing, family = nothing, verbose = true,
        kwargs...)

Performs dominance analysis. 

## Options:

    - df - input DataFrame
    - dep - dependent variable (Symbol only)
    - indeps - a vector of Symbols or tuple of Symbols. Tuples are used to create sets of variables
    - covars - covariates to be included in all models. These variables are not used in dominance analysis
    - fitstat - R2 variant to be used for GLM models. Currently, :McFadden (default) and :Nagelkerke are available
    - link - link function for GLM models
    - family - distribution family to go with the `link` function (Default: family associated with the canonical link)
    - verbose - a dot for every 10 regressions. Only used for 100 regressions or more
    - multithreads - set it to `true` to turn on multithreading
    
"""
function dominance(_data::AbstractDataFrame,
    dep::Symbol, # dependent variable name
    indeps::Vector; # independent variables in a vector, sets are allowed in tuples
    covars=[],
    fitstat=:McFadden,
    link=nothing,
    family=nothing,
    verbose=true,
    multithreads=true,
    wts=nothing)

    # prepare the data set

    # MF, MM, and y
    allvars = untuple(vcat(dep, indeps, covars))
    df = dropmissing(select(_data, allvars))
    MF = ModelFrame(term(dep) ~ sum(term.(allvars)), df)
    MM = ModelMatrix(MF)
    y = convert(Vector{Float64},response(MF))
    assign = MM.assign
    mm = MM.m
     
    # link and family
    if link != nothing
        if family == nothing
            family = canonical[link]
        end
        if !in(fitstat, [:McFadden, :Nagelkerke])
            throw(ArgumentError(fitstat, " is not allowed"))
        end
    end

    # get all combination of the indeps vector
    nvars = length(indeps)
    vvec = collect(combinations(collect(1:nvars)))
    nreg = length(vvec)
    println("\nA total of ", nreg, " regressions will be estimated.\n")

    # df for saving fitstats
    fs = DataFrame()
    fs.terms = Vector{Any}(undef, nreg)
    fs.terms_sorted = Vector{Any}(undef, nreg)
    fs.nterms = zeros(Int16, nreg)
    fs.r2m = Vector{Union{Missing,Float64}}(missing, nreg)
    for i in 1:length(indeps)
        fs[!, Symbol(i)] = Vector{Union{Missing,Float64}}(missing, nreg)
    end

    # fit stats for null models
    fm1 = get_formula(dep, covars)
    if link == nothing
        fitnull = r2(lm(fm1, df))
    else
        if wts == nothing
            # fitnull = r2(glm(fm1, df, family(), link()), fitstat)
            fitnull = r2(glm(mm[:, vidx(covars,allvars,assign)], y, family(), link()), fitstat)
        else
            fitnull = r2(glm(mm[:, vidx(covars, allvars, assign)], y, family(), link(), wts=wts), fitstat)
        end
    end

    # formulae
    fm = []
    for (i, vindex) in enumerate(vvec)
        vars = indeps[vindex]
        fs[i, :terms] = vars
        fs[i, :terms_sorted] = sort(untuple(vars))
        fs[i, :nterms] = length(vars)
        push!(fm, get_formula(dep, vcat(vars, covars)))
    end

    if multithreads == false
        for i = 1:nreg
            verbose && show_progress(i,nreg)
            # fs[i, :r2m] = get_fitstat(df, fm[i], family=family, link=link, fitstat=fitstat, wts=wts)
            fs[i, :r2m] = get_fitstat(mm[:, vidx(untuple(vcat(fs[i,:terms],covars)), allvars, assign)], y, family=family, link=link, fitstat=fitstat, wts=wts)
        end
    else
        Threads.@threads for i = 1:nreg
            verbose && show_progress(i, nreg)
            # fs[i, :r2m] = get_fitstat(df, fm[i], family=family, link=link, fitstat=fitstat, wts=wts)
            fs[i, :r2m] = get_fitstat(mm[:, vidx(untuple(vcat(fs[i, :terms], covars)), allvars, assign)], y, family=family, link=link, fitstat=fitstat, wts=wts)
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
    for (i,j) in permutations(1:nvars,2)
        tmpfs = dropmissing(fs[:, [Symbol(i), Symbol(j)]])
        fs1 = vcat(fs[i, :r2m], tmpfs[:, 1])
        fs2 = vcat(fs[j, :r2m], tmpfs[:, 2])
        compared = (fs1 .- fs2)
        if Base.all(compared .> 0.0)
            complete[i, j] = 1
        elseif Base.all(compared .< 0.0)
            complete[i, j] = -1
        end
    end

    # conditional dominance
    conditional_dom = zeros(Float64, nvars, nvars)
    conditional_dom[:, 1] = fs.r2m[1:nvars] .- fitnull
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

function show_progress(i, nreg)
    if nreg >= 100
        if mod(i, 20) == 0
            print(".")
        end
        if mod(i, 1600) == 0
            println(" ", @sprintf("%5d", i))
        end
    end
end

function vidx(indeps, allvars, assign)
    n = [1]
    for v in indeps
         i = findfirst(x -> x == v, allvars)
         n = vcat(n, findall(x -> x == i, assign))
    end
    return n
end

function get_fitstat(X, y; family = nothing, link = nothing, fitstat = nothing, wts = nothing)
    if link == nothing
        return r2(lm(X,y))
    end
    if wts == nothing
        return r2(glm(X,y, family(), link()), fitstat)
    end
    return r2(glm(X,y, family(), link(), wts = wts), fitstat)
end
# function get_fitstat(df, fm; family=nothing, link=nothing, fitstat=nothing, wts=nothing)
#     if link == nothing
#         return r2(lm(fm, df))
#     end
#     if wts == nothing
#         return r2(glm(fm, df, family(), link()), fitstat)
#     end
#     return r2(glm(fm, df, family(), link(), wts=wts), fitstat)
# end

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
    ntables = ceil(Int8, nvars / 7)
    if ntables > 1
        for i in 1:ntables
            fr = 1 + 7*(i-1)
            to = min(nvars,7*i)
            pretty_table(io,
                dom.comstat[:, fr:to],
                header=indepnames[fr:to],
                row_labels=indepnames,
                row_label_column_title="dominates?",
                hlines=[0, 1, nvars + 1],
                vlines=[1]
            )
        end
    else
        pretty_table(io,
            dom.comstat,
            header=indepnames,
            row_labels=indepnames,
            row_label_column_title="dominates?",
            hlines=[0, 1, nvars + 1],
            vlines=[1]
        )
    end

    # conditional dominance
    println(io, "\nConditional dominance:")
    if ntables > 1
        for i in 1:ntables
            fr = 1 + 7 * (i - 1)
            to = min(nvars, 7 * i)
            pretty_table(io,
                dom.constat[:,fr:to],
                header=collect(fr:to),
                row_labels=indepnames,
                row_label_column_title="Variables",
                formatters=ft_printf("%6.4f", 1:nvars),
                hlines=[0, 1, nvars + 1],
                vlines=[1]
            )
        end
    else
        pretty_table(io,
            dom.constat,
            header=collect(1:nvars),
            row_labels=indepnames,
            row_label_column_title="Variables",
            formatters=ft_printf("%6.4f", 1:nvars),
            hlines=[0, 1, nvars + 1],
            vlines=[1]
        )
    end

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
    if indeps == []
        return Term(dep) ~ term(1)
    end
    return Term(dep) ~ sum(term.(untuple(indeps)))
end

"""
    dominance_designations(d::Domin)

Prints strongest dominance designations based on the dominance analysis provided as input.
"""
function dominance_designations(dom)
    vars = Dominance.untuple(dom.indeps)
    nvars = length(vars)

    # compute dominance
    davg = mean(dom.constat, dims=2)
    ddesig = zeros(Int8, nvars, nvars)
    for (i, j) in permutations(1:nvars, 2)
        if dom.comstat[i, j] == 1
            ddesig[i, j] = 1 # complete dominance
        elseif sum(dom.constat[i, :] .> dom.constat[j, :]) == nvars
            ddesig[i, j] = 2 # conditional dominance
        elseif davg[i] > davg[j]
            ddesig[i, j] = 3 # general dominance
        end
    end

    # sort by value in reverse order
    # dd = sort(dd,byvalue=true,rev=true)
    println("Strongest dominance desginations:\n")
    # complete dominance
    for (i, j) in permutations(1:nvars, 2)
        if ddesig[i, j] == 1
            println(vars[i], " completely dominates ", vars[j])
        end
    end
    # conditional dominance
    for (i, j) in permutations(1:nvars, 2)
        if ddesig[i, j] == 2
            println(vars[i], " conditionally dominates ", vars[j])
        end
    end
    # general dominance
    for (i, j) in permutations(1:nvars, 2)
        if ddesig[i, j] == 3
            println(vars[i], " generally dominates ", vars[j])
        end
    end
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

Luchman's intro and formulas at 
https://cran.r-project.org/web/packages/domir/vignettes/domir_basics.html
"""

