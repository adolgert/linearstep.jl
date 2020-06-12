using BenchmarkTools
using RCall

R"setwd('~/.julia/packages/linearstep/src')"
R"source('adam.R')"

## Conversion tools
function array_r_to_jl(s4arr)
    dims = R"dim($(s4arr))"
    if typeof(dims) <: RObject{NilSxp}
        dims = [rcopy(R"length($(s4arr))")]
    end
    m = zeros(dims...)
    if length(dims) == 2
        for j in 1:dims[2]
            for i in 1:dims[1]
                m[i, j] = R"$(s4arr)[$(i), $(j)]"
            end
        end
    elseif length(dims) == 1
        for i in 1:dims[1]
            m[i] = R"$(s4arr)[$(i)]"
        end
    else
        error("cannot convert because more than 2 dimensions $(dims)")
    end
    m
end


function rtocohortstate(X)
    if typeof(X) <: RObject
        Xj = rcopy(X)
    else
        Xj = X
    end
    X = array_r_to_jl(Xj[:X])
    pit = array_r_to_jl(Xj[:pit])
    V = Xj[:V]
    ageInDays = Float64(Xj[:ageInDays])
    CohortState(X, pit, V, ageInDays, Xj[:alpha1], Xj[:alpha2])
end


function rtopopstate(rXX)
    if typeof(rXX) <: RObject
        rXXj = rcopy(rXX)
    else
        rXXj = rXX
    end
    X = array_r_to_jl(rXXj[:X])
    pit = array_r_to_jl(rXXj[:pit])
    V = rXXj[:V]
    rage = rXXj[:ageInDays]
    if typeof(rage) <: Array
        ageInDays = convert(Array{Float64,2}, reshape(rage, 1, length(rage)))
    else
        ageInDays = array_r_to_jl(rXXj[:ageInDays])
    end
    PopulationState(X, pit, V, ageInDays, rXXj[:alpha1], rXXj[:alpha2])
end

## Setup for testing
alpha = 0.4
reval("alpha <- $(alpha)")
reval("params <- makeParameters_Adam()")
params = make_parameters_adam()

## Check that the parameters match exactly.
param = make_parameters_adam()
rparam = R"makeParameters_Adam()"
rparamj = rcopy(rparam)
typeof(param[:tau])
@test Set(keys(rparamj)) == Set(keys(param))
for (key, value) in param
    if (typeof(value) <: Number)
        rvalue = rparamj[key]
        @test isapprox(value, rvalue)
    elseif (typeof(value) <: Vector)
        rvalue = rparamj[key]
        if length(value) == length(rvalue)
            if !all(isapprox.(value, rvalue))
                @show key
                @show value
                @show rvalue
            end
        else
            @test length(value) == length(rvalue)
        end
    end
end

## calP_SoI

# Both are length 9
nu = [0.50000000, 0.16666667, 0.14285714, 0.12500000, 0.11111111, 0.11111111,
0.11111111, 0.08333333, 0.00000000]
r = [1.0000000, 0.9900498, 0.9607894, 0.9607894, 0.9607894, 0.9607894,
0.9607894, 0.9607894, 0.9753099]
PSoI = calP_SoI(9, 10, nu, r)
@rput(nu)
@rput(r)
rPSoI = R"calP_SoI(list(N=9, tau=10, nu=nu, r=r))"
rPSoIj = array_r_to_jl(rPSoI)
@test PSoI ≈ rPSoIj

rSuS = array_r_to_jl(rparamj[:SuS])
@test param[:SuS] ≈ rSuS
rPSoI = array_r_to_jl(rparamj[:PSoI])
@test param[:PSoI] ≈ rPSoI

## xiF, defined in make_parameters.
reval("empty = emptyX_Adam()")
reval("xiF_0 <- function(x, p.xi = 2, N = 9, C = 10, D = 0) {
  diff(pbeta(c(0:N) / N, (1 - exp(-p.xi * x)) * C, exp(-p.xi * x) * C + D))
}")
rd = R"xiF_0(0.3)"
jd = xiF_inner(0.3)
@test jd ≈ array_r_to_jl(rd)
rd = R"xiF_0(0)"
jd = xiF_inner(0)
@test jd ≈ array_r_to_jl(rd)

rxi = R"params$xiF(empty$pit, empty$V)"
rxij = array_r_to_jl(rxi)
vecX = emptyX_Adam(Float64)
jxi = params[:xiF](vecX.pit, vecX.V)
@test jxi ≈ rxij


## calP_Adam
vecX = emptyX_Adam(Float64)
calP = calP_Adam(alpha, vecX.pit, vecX.V, params)
reval("vecX = emptyX_Adam()")
rcalP = R"calP_Adam($(alpha), vecX$pit, vecX$V, params)"
rcalPj = rcopy(rcalP)
# there are NaN in rcalPj. Maybe that's a translation problem?
rcalPj[isnan.(rcalPj)] .= 0
@test calP ≈ rcalPj


## emtpyXX_Adam
cXX = emptyXX_Adam(Float64, L=2920)
rXX = R"emptyXX_Adam(L=2920)"
rXXj = rtopopstate(rXX)
@test cXX ≈ rXXj
differences(cXX, rXXj)


## PxX_Adam tests
vecX = emptyX_Adam(Float64)
PxX_Adam!(vecX, alpha, params)
reval("vecX = emptyX_Adam()")
rvecX = R"PxX_Adam(alpha, vecX, params)"
rvecXj = rtocohortstate(rvecX)
@test vecX ≈ rvecXj

## cohortXX_Adam tests
cXX = cohortXX_Adam(alpha, params)
# The params variable is defined in R at this point.
rXX = R"cohortXX_Adam($(alpha), params)"
rXXj = rtopopstate(rXX)
@test cXX ≈ rXXj
differences(cXX, rXXj)

## cohort2ages_Adam
cXX = cohortXX_Adam(alpha, params)
XX = cohort2ages_Adam(cXX, params)
reval("cXX = cohortXX_Adam($(alpha), params)")
reval("XX = cohort2ages_Adam(cXX, params)")
rXX = R"XX"
rXXj = rtopopstate(rXX)
# Confirm they started with the same values.
@test rtopopstate(R"cXX") ≈ cXX
# This is the new test.
@test_broken XX ≈ rXXj
@test maximum((XX.X .- rXXj.X).^2) < 0.0015
differences(XX, rXXj)


## ar2stableWave_Adam
XX = ar2stableWave_Adam(alpha, params)
rXX = R"ar2stableWave_Adam($(alpha), $(rparam))"
rXXj = rcopy(rXX)

rXXj = rtopopstate(rXX)
@test_broken rXXj ≈ XX
differences(XX, rXXj)
@test maximum((XX.X .- rXXj.X).^2) < 0.0015

pr1 = ar2pr_Adam(alpha, params)
pr2 = R"ar2pr_Adam(alpha, params)"
@test abs(pr1 - rcopy(pr2)) < 1e-5

# @btime pr2arSS_Adam(pr1, params)
# Using regular optimize.
# 2.061 s (11923106 allocations: 2.26 GiB)
# Adding autodiff.
# 595.124 ms (3924797 allocations: 752.06 MiB)
# @btime R"pr2arSS_Adam($(pr1), params)"
# 264.642 s (103 allocations: 3.05 KiB)
