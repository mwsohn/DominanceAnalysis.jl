using CategoricalArrays, DataFrames, GLM, RDatasets, DominanceAnalysis, Stella

mtcars = dataset("datasets","mtcars")
renvars!(mtcars)


@testset "Linear Regression without sets" begin

    dom1 = dominance(mtcars, :mpg, [:am, :cyl, :carb])

    # general dominance compared to the values obtained from domir.R
    @test isapprox(dom1.fit_overall, 0.8113023, rtol = 0.0001)
    @test isapprox(dom1.domstat[1, 1], 0.2156848, rtol=0.0001) # am
    @test isapprox(dom1.domstat[2, 1], 0.4173094, rtol=0.0001) # cyl
    @test isapprox(dom1.domstat[3, 1], 0.1783081, rtol=0.0001) # carb

    # conditional dominance
    @test isapprox(dom1.constat[1, 3], 0.07076149, rtol=0.0001) # am at 3
    @test isapprox(dom1.constat[2, 3], 0.10762967, rtol=0.0001) # cyl at 3
    @test isapprox(dom1.constat[3, 3], 0.05228872, rtol = 0.0001) # carb at 3

end

@testset "Linear Regression with a set of variables" begin

    dom2 = dominance(mtcars, :mpg, [:am, :cyl, (:carb, :wt)])

    @test isapprox(dom2.fit_overall, 0.8502, rtol=0.0001) # overall fit statistic

    # general dominance compared to the values from Stata domin
    @test isapprox(dom2.domstat[1, 1], 0.1307, rtol=0.01) # am
    @test isapprox(dom2.domstat[2, 1], 0.3308, rtol=0.01) # cyl
    @test isapprox(dom2.domstat[3, 1], 0.3887, rtol=0.01) # Set 1 (carb wt)

    # conditional dominance compared to the values from Stata domin
    @test isapprox(dom2.constat[1, 3], 0.0077, rtol=0.01) # am
    @test isapprox(dom2.constat[2, 3], 0.0417, rtol=0.01) # cyl
    @test isapprox(dom2.constat[3, 3], 0.0912, rtol=0.01) # Set 1 (carb wt)
end







