using CategoricalArrays, DataFrames, GLM, RDatasets

mtcars = dataset("datasets","mtcars")


@testset "Linear Regression without sets" begin

    dom1 = dominance(mtcars, :mpg, [:am, :cyl, :carb])

    # general dominance compared to the values obtained from domir.R
    @test isapprox(dom1.fit_overall, 0.8113023)
    @test isapprox(dom1.domstat[1,1], 0.2156848) # am
    @test isapprox(dom1.domstat[2,1], 0.4173094) # cyl
    @test isapprox(dom1.domstat[3,1], 0.1783081) # carb

    # conditional dominance
    @test isapprox(dom1.constat[1, 3], 0.07076149) # am at 3
    @test isapprox(dom1.constat[2, 3], 0.10762967) # cyl at 3
    @test isapprox(dom1.constat[3, 3], 0.05228872) # carb at 3

end

@testset "Linear Regression with a set of variables" begin

    dom2 = dominance(mtcars, :mpg, [:am, :cyl, (:carb, :wt)])

    @test isapprox(dom2.fit_overall, 0.8502) # overall fit statistic

    # general dominance compared to the values from Stata domin
    @test isapprox(dom2.domstat[1, 1], 0.1307) # am
    @test isapprox(dom2.domstat[2, 1], 0.3308) # cyl
    @test isapprox(dom2.domstat[3, 1], 0.3887) # Set 1 (carb wt)

    # conditional dominance compared to the values from Stata domin
    @test isapprox(dom2.constat[1, 3], 0.0077) # am
    @test isapprox(dom2.constat[2, 3], 0.0417) # cyl
    @test isapprox(dom2.constat[3, 3], 0.0912) # Set 1 (carb wt)
end







