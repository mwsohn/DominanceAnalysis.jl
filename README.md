# Dominance
 A package for conducting dominance analysis.

 ## Installation

 `] add https://github.com/mwsohn/Dominance.jl`

 ## Syntax

```
    dominance(df,dep,indeps,covars,link,family)
```

### Options:
- df - input data in a DataFrame
- dep - dependent variable (Symbol)
- indeps - a vector of independent variables. A set of variables can be
    specified as a tuple within the vector
- covars - a vector of covariates that will be included in all models.
    They will not be used as independent variables in the dominance analysis
- link - a link function for GLM models. Valid link functions are LogitLink,
    LogLink, IdentifyLink, etc. If `family` option is not specified, `link` will
    be used to determine the type of GLM model
- family - a distribution family. If not specified, a distribution will be chosen

## Example

```
julia> dominance(auto,:price, [:mpg, :trunk, :foreign, :weight])
A total of 15 regressions will be estimated.

Number of observations      =              74
Number of regression models =              15
Overall Fit Statistic       =          0.5082


Dominance Statistics and Ranking:
─────────┬────────────────────────────────────────────────────────────────
   price │ Dominance Statistic  Standardized Dominance            Ranking 
─────────┼────────────────────────────────────────────────────────────────
     mpg │              0.1031                  0.2028                  3
   trunk │              0.0399                  0.0785                  4
 foreign │              0.1192                  0.2346                  2
  weight │              0.2460                  0.4840                  1
─────────┴────────────────────────────────────────────────────────────────

Complete dominance:
────────────┬─────────────────────────────
 dominates? │ mpg  trunk  foreign  weight 
────────────┼─────────────────────────────
        mpg │   0      0        0      -1
      trunk │   0      0        0      -1
    foreign │   0      0        0      -1
     weight │   1      1        1       0
────────────┴─────────────────────────────

Conditional dominance:
───────────┬────────────────────────────────
 Variables │      1       2       3       4 
───────────┼────────────────────────────────
       mpg │ 0.2196  0.1361  0.0564  0.0003
     trunk │ 0.0988  0.0444  0.0078  0.0086
   foreign │ 0.0024  0.1010  0.1636  0.2099
    weight │ 0.2901  0.2553  0.2236  0.2149
───────────┴────────────────────────────────
```



