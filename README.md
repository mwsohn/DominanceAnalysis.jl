# Dominance Analysis
 A Julia package for conducting dominance analysis. An excellent introduction and formulas on `dominance analysis` by Joseph Luchman can be found 
 at https://cran.r-project.org/web/packages/domir/vignettes/domir_basics.html. See below for references as well.[^1][^2][^3][^4][^5][^6]

  ## Installation

 `] add https://github.com/mwsohn/Dominance.jl`

 ## Syntax

```
    dominance(df,dep,indeps; covars = [],link = nothing,family = nothing, verbose = false, multithreads = false, wts=nothing)
```

### - Options
- df - input data in a DataFrame
- dep - dependent variable
- indeps - a vector of independent variables. A `set` can be
    specified as a tuple within the vector
- covars - a vector of covariates that will be included in all models.
    They will not be used as independent variables in the dominance analysis
- link - a link function for GLM models. Valid link functions are LogitLink,
    LogLink, IdentifyLink, etc. If `family` option is not specified, `link` will
    be used to determine the type of GLM model
- family - a distribution family. If not specified, a distribution will be
    automatically chosen based on the link function.
- fitstat - choose a pseudo R² method (:McFadden or :Nagelkerke). The default is :McFadden.
- multithreads - set it to `false` to turn off the use of multithreading. 
    Multithreading with 8 threads increased the execution speed by 2 times.
- verbose - set it to `true` to turn on the verbose model that will display dots to be printed
    as a progress indicator
- wts - specify a weight vector for complex survey data. 

### - Notes

All variables must be in Symbols. A group of variables to be treated as a `set` should be entered as a tuple
in the `indeps` vector. All multi-valued CategoricalArrays will be treated as a set by default. Linear regression
and logistic regression models with or without weights have been tested against the output from Stata domin
program (https://journals.sagepub.com/doi/pdf/10.1177/1536867X211025837).

If you turn on multithreading (`multithreads = true`), this program on Julia v1.11.4 is about 10 times faster than the Stata `domin.ado` on Stata MP v.18.0 for a domianance analysis involving 255 models with a data set of over 200,000 observations. 

For multithreading to be effective, please turn on multithreading when starting up Julia (`--threads=auto` or `"julia.NumThreads" : "auto"` in Visual Studio Editor settings).

When performing dominance analysis with a complex survey data, you only need to specify the probability weights as a
vector for the `wts` option.[^7] 

## Example

I am using `auto.dta` downloaded from http://www.stata-press.com/data/r13/auto.dta, converted to Julia DataFrame using
`read_stata` in Stella.jl.

### 1. Linear regression

```
julia> dominance(auto,:price, [:mpg, :trunk, :foreign, :weight])

A total of 15 regressions will be estimated.

Number of observations      =              74
Number of regression models =              15
Overall Fit Statistic       =          0.5082



Dominance Statistics and Ranking:
─────────┬─────────────────────────────────────────────────────────────
   price │ General Dominance  Standardized Dominance           Ranking 
─────────┼─────────────────────────────────────────────────────────────
     mpg │            0.1031                  0.2028                 3
   trunk │            0.0399                  0.0785                 4
 foreign │            0.1192                  0.2346                 2
  weight │            0.2460                  0.4840                 1
─────────┴─────────────────────────────────────────────────────────────

Conditional dominance:
───────────┬────────────────────────────────
 Variables │      1       2       3       4 
───────────┼────────────────────────────────
       mpg │ 0.2196  0.1361  0.0564  0.0003
     trunk │ 0.0988  0.0444  0.0078  0.0086
   foreign │ 0.0024  0.1010  0.1636  0.2099
    weight │ 0.2901  0.2553  0.2236  0.2149
───────────┴────────────────────────────────

Complete dominance proportions:
────────────┬─────────────────────────────────
 dominates? │    mpg   trunk  foreign  weight 
────────────┼─────────────────────────────────
        mpg │ 0.0000  0.5000   0.5000  0.0000
      trunk │ 0.5000  0.0000   0.2500  0.0000
    foreign │ 0.5000  0.7500   0.0000  0.0000
     weight │ 1.0000  1.0000   1.0000  0.0000
────────────┴─────────────────────────────────
```

### 2. Linear regression with a set of variables

```
julia> dominance(auto, :price, [:mpg, :trunk, :foreign, (:weight, :turn)])

A total of 15 regressions will be estimated.

Number of observations      =              74
Number of regression models =              15
Overall Fit Statistic       =          0.5336

Set    1 = (:weight, :turn)


Dominance Statistics and Ranking:
─────────┬─────────────────────────────────────────────────────────────
   price │ General Dominance  Standardized Dominance           Ranking 
─────────┼─────────────────────────────────────────────────────────────
     mpg │            0.1036                  0.1941                 2
   trunk │            0.0388                  0.0727                 4
 foreign │            0.0881                  0.1651                 3
   Set 1 │            0.3031                  0.5681                 1
─────────┴─────────────────────────────────────────────────────────────

Conditional dominance:
───────────┬────────────────────────────────
 Variables │      1       2       3       4 
───────────┼────────────────────────────────
       mpg │ 0.2196  0.1373  0.0574  0.0000
     trunk │ 0.0988  0.0437  0.0062  0.0064
   foreign │ 0.0024  0.0813  0.1225  0.1462
     Set 1 │ 0.3776  0.3236  0.2711  0.2403
───────────┴────────────────────────────────

Complete dominance proportions:
────────────┬─────────────────────────────────
 dominates? │    mpg   trunk  foreign   Set 1 
────────────┼─────────────────────────────────
        mpg │ 0.0000  0.7500   0.5000  0.0000
      trunk │ 0.2500  0.0000   0.2500  0.0000
    foreign │ 0.5000  0.7500   0.0000  0.0000
      Set 1 │ 1.0000  1.0000   1.0000  0.0000
────────────┴─────────────────────────────────
```

### 3. Logistic regression with a set of variables

```
julia> using Statistics
julia> auto.pricelow = auto.price .< mean(auto.price);
julia> dominance(auto, :pricelow, [:mpg, :trunk, :foreign, (:weight, :turn)], family=Bernoulli, link=LogitLink)

A total of 15 regressions will be estimated.

Number of observations      =              74
Number of regression models =              15
Overall Fit Statistic       =          0.4089

Set    1 = (:weight, :turn)


Dominance Statistics and Ranking:
──────────┬─────────────────────────────────────────────────────────────
 pricelow │ General Dominance  Standardized Dominance           Ranking 
──────────┼─────────────────────────────────────────────────────────────
      mpg │            0.0750                  0.1835                 3
    trunk │            0.0158                  0.0386                 4
  foreign │            0.1419                  0.3470                 2
    Set 1 │            0.1762                  0.4309                 1
──────────┴─────────────────────────────────────────────────────────────

Conditional dominance:
───────────┬────────────────────────────────
 Variables │      1       2       3       4 
───────────┼────────────────────────────────
       mpg │ 0.1253  0.1133  0.0603  0.0012
     trunk │ 0.0325  0.0251  0.0034  0.0022
   foreign │ 0.0202  0.1384  0.1923  0.2167
     Set 1 │ 0.1571  0.1910  0.1845  0.1723
───────────┴────────────────────────────────

Complete dominance proportions:
────────────┬─────────────────────────────────
 dominates? │    mpg   trunk  foreign   Set 1 
────────────┼─────────────────────────────────
        mpg │ 0.0000  0.7500   0.5000  0.0000
      trunk │ 0.2500  0.0000   0.2500  0.0000
    foreign │ 0.5000  0.7500   0.0000  0.5000
      Set 1 │ 1.0000  1.0000   0.5000  0.0000
────────────┴─────────────────────────────────


```

## References

[^1]: Joseph N. Luchman, 2013. "DOMIN: Stata module to conduct dominance analysis," Statistical Software Components S457629, Boston College Department of Economics, revised 07 Jan 2025. 

[^2]: Azen R, Budescu DV. The dominance analysis approach for comparing predictors in multiple regression. Psychol Methods. 2003 Jun;8(2):129-48. https://doi.org/10.1037/1082-989x.8.2.129. PMID: 12924811.

[^3]: Azen, R., & Traxel, N. (2009). Using Dominance Analysis to Determine Predictor Importance in Logistic Regression. Journal of Educational and Behavioral Statistics, 34(3), 319-347. https://doi.org/10.3102/1076998609332754 (Original work published 2009)

[^4]: Azen, Razia, David V Budescu, and Benjamin Reiser. 2001. “Criticality of Predictors in Multiple Regression.” British Journal of Mathematical and Statistical Psychology 54 (2): 201–25. https://doi.org/10.1348/000711001159483.

[^5]: Budescu, David V. 1993. “Dominance Analysis: A New Approach to the Problem of Relative Importance of Predictors in Multiple Regression.” Psychological Bulletin 114 (3): 542–51. https://doi.org/10.1037/0033-2909.114.3.542.

[^6]: Grömping, Ulrike. 2007. “Estimators of Relative Importance in Linear Regression Based on Variance Decomposition.” The American Statistician 61 (2): 139–47. https://doi.org/10.1198/000313007X188252.

[^7]: Luchman, Joseph Nicholas. 2015. “Determining Subgroup Difference Importance with Complex Survey Designs: An Application of Weighted Dominance Analysis.” Survey Practice 8 (5). https://​doi.org/​10.29115/​SP-2015-0022.





