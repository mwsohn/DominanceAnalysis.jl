using CategoricalArrays, DataFrames, GLM, RDatasets, DominanceAnalysis, JLD2

mtcars = dataset("datasets","mtcars")
rename!(mtcars, Pair.(names(mtcars), lowercase.(names(mtcars))))

@testset "Linear Regression without sets" begin

    dom1 = dominance(mtcars, :mpg, [:am, :cyl, :carb])

    # general dominance compared to the values obtained from domir.R
    @test isapprox(dom1.fit_overall, 0.8113023, atol = 0.0001)
    @test isapprox(dom1.domstat[1, 1], 0.2156848, atol=0.0001) # am
    @test isapprox(dom1.domstat[2, 1], 0.4173094, atol=0.0001) # cyl
    @test isapprox(dom1.domstat[3, 1], 0.1783081, atol=0.0001) # carb

    # conditional dominance
    @test isapprox(dom1.constat[1, 3], 0.07076149, atol=0.0001) # am at 3
    @test isapprox(dom1.constat[2, 3], 0.10762967, atol=0.0001) # cyl at 3
    @test isapprox(dom1.constat[3, 3], 0.05228872, atol = 0.0001) # carb at 3

end

@testset "Linear Regression with a set of variables" begin

    dom2 = dominance(mtcars, :mpg, [:am, :cyl, (:carb, :wt)])

    @test isapprox(dom2.fit_overall, 0.8502, atol=0.0001) # overall fit statistic

    # general dominance compared to the values from Stata domin
    @test isapprox(dom2.domstat[1, 1], 0.1307, atol=1e-3) # am
    @test isapprox(dom2.domstat[2, 1], 0.3308, atol=1e-3) # cyl
    @test isapprox(dom2.domstat[3, 1], 0.3887, atol=1e-3) # Set 1 (carb wt)

    # conditional dominance compared to the values from Stata domin
    @test isapprox(dom2.constat[1, 3], 0.0077, atol=1e-3) # am
    @test isapprox(dom2.constat[2, 3], 0.0417, atol=1e-3) # cyl
    @test isapprox(dom2.constat[3, 3], 0.0912, atol=1e-3) # Set 1 (carb wt)
end

llcp22 = JLD2.load_object("llcp2022_shortsleep_subset.jld2")
indeps = [:agecat, :male, :race2, :bmicat,
    (:satisfaction, :support, :isolation, :lostemp, :foodstamp,
     :foodinsec, :billpay, :utilpay, :transportation, :stress)
]

@testset "Logistic Regression with a set of variables and without weights" begin

    dom3 = dominance(llcp22, :sleep_short, indeps, family=Bernoulli, link=LogitLink)

    # general dominance compared to the values from Stata domin
    @test isapprox(dom3.domstat[1, 1], 0.0044, atol=1e-3) # agecat
    @test isapprox(dom3.domstat[2, 1], 0.0001, atol=1e-3) # male
    @test isapprox(dom3.domstat[3, 1], 0.0059, atol=1e-3) # race2
    @test isapprox(dom3.domstat[4, 1], 0.0027, atol=1e-3) # bmicat
    @test isapprox(dom3.domstat[5, 1], 0.0776, atol=1e-3) # Set 1 (10 variables)


    # conditional dominance compared to the values from Stata domin
    @test isapprox(dom3.constat[1, 5], 0.0017, atol=1e-3) # agecat
    @test isapprox(dom3.constat[2, 5], 0.0002, atol=1e-3) # male
    @test isapprox(dom3.constat[3, 5], 0.0034, atol=1e-3) # race2
    @test isapprox(dom3.constat[4, 5], 0.0015, atol=1e-3) # bmicat
    @test isapprox(dom3.constat[5, 5], 0.0727, atol=1e-3) # Set 1

end


@testset "Logistic Regression with a set of variables and without weights" begin

    dom3 = dominance(llcp22, :sleep_short, indeps, family=Bernoulli, link=LogitLink)

    # general dominance compared to the values from Stata domin
    @test isapprox(dom3.domstat[1, 1], 0.0044, atol=1e-3) # agecat
    @test isapprox(dom3.domstat[2, 1], 0.0001, atol=1e-3) # male
    @test isapprox(dom3.domstat[3, 1], 0.0059, atol=1e-3) # race2
    @test isapprox(dom3.domstat[4, 1], 0.0027, atol=1e-3) # bmicat
    @test isapprox(dom3.domstat[5, 1], 0.0776, atol=1e-3) # Set 1 (10 variables)


    # conditional dominance compared to the values from Stata domin
    @test isapprox(dom3.constat[1, 5], 0.0017, atol=1e-3) # agecat
    @test isapprox(dom3.constat[2, 5], 0.0002, atol=1e-3) # male
    @test isapprox(dom3.constat[3, 5], 0.0034, atol=1e-3) # race2
    @test isapprox(dom3.constat[4, 5], 0.0015, atol=1e-3) # bmicat
    @test isapprox(dom3.constat[5, 5], 0.0727, atol=1e-3) # Set 1

end

llcp22 = dropmissing(llcp22) 

@testset "Logistic Regression with a set of variables and weights" begin

    dom4 = dominance(llcp22, :sleep_short, indeps, family=Bernoulli, link=LogitLink, wts=llcp22._llcpwt)

    # general dominance compared to the values from Stata domin
    @test isapprox(dom4.domstat[1, 1], 0.0029, atol=1e-3) # agecat
    @test isapprox(dom4.domstat[2, 1], 0.0001, atol=1e-3) # male
    @test isapprox(dom4.domstat[3, 1], 0.0050, atol=1e-3) # race2
    @test isapprox(dom4.domstat[4, 1], 0.0025, atol=1e-3) # bmicat
    @test isapprox(dom4.domstat[5, 1], 0.0775, atol=1e-3) # Set 1 (10 variables)


    # conditional dominance compared to the values from Stata domin
    @test isapprox(dom4.constat[1, 5], 0.0020, atol=1e-3) # agecat
    @test isapprox(dom4.constat[2, 5], 0.0003, atol=1e-3) # male
    @test isapprox(dom4.constat[3, 5], 0.0036, atol=1e-3) # race2
    @test isapprox(dom4.constat[4, 5], 0.0015, atol=1e-3) # bmicat
    @test isapprox(dom4.constat[5, 5], 0.0752, atol=1e-3) # Set 1

end




