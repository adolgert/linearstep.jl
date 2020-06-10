using Distributions: Beta

function wXageF0(ages)
    1.5 .* ages ./ (2 .+ ages)
end


function bandSparse(N, k::Integer)
    m = zeros(N, N)
    if k > 0
        for i in 1:(N-k)
            m[i, i + k] = 1.0
        end
    else
        for i in 1:(N + k)
            m[i - k, i] = 1.0
        end
    end
    m
end


function bandSparse(N, ks::AbstractArray)
    m = zeros(N, N)
    for k in ks
        if k > 0
            for i in 1:(N-k)
                m[i, i + k] = 1.0
            end
        else
            for i in 1:(N + k)
                m[i - k, i] = 1.0
            end
        end
    end
    m
end


function shiftO(N)
  bandSparse(N, 1)
end


function stageO(N, p)
  bandSparse(N, 0) * (1 - p) .+ bandSparse(N, 1) * p
end


function bdiag(A, B)
  calO = zeros(size(A)[1] + size(B)[1], size(A)[2] + size(B)[2])
  calO[1:size(A)[1], 1:size(A)[2]] .= A
  calO[(size(A)[1] + 1):end, (size(A)[2] + 1):end] .= B
  calO
end


function twoBlockO(A, B, p)
  nA = size(A)[1]
  calO = bdiag(A, B)
  calO[nA, nA + 1] = p
  calO
end

function ages_1066(unit = 365)
  x1 = collect(1:18) * 10
  x2 = maximum(x1) .+ collect(1:18) * (30 + 1 / 1.8)
  x3 = maximum(x2) .+ collect(1:18) * 365
  x4 = maximum(x3) .+ collect(1:12) * 365 * 5
  vcat(x1, x2, x3, x4) / unit
end

function calO_1066()
  A = shiftO(18)
  B = stageO(18, 1 / 3)
  C = stageO(18, 1 / 36.5)
  D = stageO(12, 1 / 36.5 / 5)
  calO = twoBlockO(A, B, 1 / 3)
  calO = twoBlockO(calO, C, 1 / 36.5)
  calO = twoBlockO(calO, D, 1 / 36.5 / 5)
  return(calO)
end

function calP_SoI(N, tau, nu, r)
  s = r .* (1 .- nu)
  v = r .* nu
  # Remain in state with probability s, go to next state probability v
  PSoI = bandSparse(N, 0) .* s .+ bandSparse(N, -1) .* v
  PSoI = vcat(transpose(1 .- r), PSoI) # rbind
  PSoI = hcat(zeros(10), PSoI) # cbind
  # Uninfected remain uninfected.
  PSoI[1, 1] = 1
  PSoI
end


function test_calP_SoI()
  # Both are length 9
  nu = [0.50000000, 0.16666667, 0.14285714, 0.12500000, 0.11111111, 0.11111111,
    0.11111111, 0.08333333, 0.00000000]
  r = [1.0000000, 0.9900498, 0.9607894, 0.9607894, 0.9607894, 0.9607894,
    0.9607894, 0.9607894, 0.9753099]
  calP_SoI(9, 10, nu, r)
end


function SuSblock_Adam(N)
  SuS = bandSparse(N, 0:(N - 1))
  vcat(transpose(zeros(N + 1)), hcat(ones(N), SuS))
end


function ageXimmune(ageInDays)
    ageInDays .> (8 * 365)
end


function xiF_inner(x; p_xi = 2, N = 9, C = 10, D = 0)
    dist = Beta((1 - exp(-p_xi * x)) * C, exp(-p_xi * x) * C + D)
    pdf(dist, x)
end


# I can't tell whether V is an array or a single float when it gets used.
# If it's an array, use pdf.() above, and it will work out.
function xiF_h(pit, V ; p_xi1 = 5, p_xi2 = 2.5, N = 9, eps = 0.001)
  v = 1 - exp(-(pit[2] + p_xi2 * pit[3]) / V)
  xi1 = xiF_inner(v, p_xi = p_xi1)
  xi = xi1 + eps
  xi / sum(xi)
end


function make_parameters_adam()
    tau = 10
    N = 9
    nuVals = [1/2, 1/6, 1/7, 1/8, 1/9, 1/9, 1/9, 1/12, 0]
    rVals = vcat([0, 1/1000], ones(6) / 250, [1/400])
    calDtrue_SoI =vcat([0], ones(9))
    calDlm_SoI = vcat([0, 0.99, .95], 1 .- collect(3:7).^3 / 8.6^3, [.05, .001])
    calDimm_SoI = vcat([0, 0.95, .9, .6, .3], zeros(5))
    calK_SoI = vcat([0, 0], collect(8:-1:1) / 10)
    calBsevere_SoI = vcat([0, .3], zeros(8))
    calBfever_SoI = vcat([0, 1, .7], exp.(-collect(1:6)), [0])
    rho = 0.2
    etaXstage = vcat([0, 0.99, .95], 1 .- collect(3:7).^3 / 8.6^3, [.05, .001]).^2
    rhoXstage = [0, 1, 1, 1, .8, .8, .6, .6, .4, .4]
    ageForImmune = 8
    delta = 1 / 730
    mu1 = 0.001
    mu2 = 0.01
    mu3 = 0.001
    wF = wXageF0
    agesF = ages_1066
    calOF = calO_1066
    # Probability of taking drugs, given _any_ state. It's a background rate.
    Dd = 1 - exp(-delta / tau)
    r = exp.(-rVals * tau)

    PSoI = calP_SoI(N, tau, nuVals, r)
    SuS = SuSblock_Adam(N)
    ages = agesF()
    aXi = (ages .> 8)
    demog = vcat(diff(ages), [15])
    demog = demog / sum(demog)
    demog29 = (ages .>= 2) .& (ages .< 10)
    demog29 = demog29 / sum(demog29)
    U5 = copy(demog)
    O5 = copy(demog)
    O5[ages .< 5] .= 0
    U5[ages .>= 5] .= 0
    calO = calOF()
    w = wF(ages)
    calBfever = vcat(calBfever_SoI, calBfever_SoI, [0, 0])
    trt = rho .* rhoXstage .* etaXstage
    calBtrtdmal = vcat(trt .* calBfever_SoI, trt .* calBfever_SoI, [0, 0])
    calBuntrmal = vcat((1 .- trt) .* calBfever_SoI, (1 .- trt) .* calBfever_SoI, [0, 0])
    calBsevere = vcat(calBsevere_SoI, calBsevere_SoI, [0, 0])
    calBtrtdsev = vcat(trt .* calBsevere_SoI, trt .* calBsevere_SoI, [0, 0])
    calBuntrsev = vcat((1 .- trt) .* calBsevere_SoI, (1 .- trt) .* calBsevere_SoI, [0, 0])
    calK = vcat(calK_SoI, calK_SoI, [0, 0])
    calK[12] = calK[3]
    Dict{Symbol,Any}(
      :tau => tau,
      :N => N,
      :calDtrue => vcat(calDtrue_SoI, collect(1:10), [0, 0]),
      :calDlm => vcat(calDlm_SoI, calDlm_SoI, [0, 0]),
      :calDimm => vcat(calDimm_SoI, calDimm_SoI, [0, 0]),
      :calBfever => calBfever,
      :calBtrtdmal => calBtrtdmal,
      :calBuntrmal => calBuntrmal,
      :calBsevere => calBsevere,
      :calBtrtdsev => calBtrtdsev,
      :calBuntrsev => calBuntrsev,
      :calK => calK,
      :nu => nuVals,
      :r => r,
      :Dd => Dd,
      :rho => rho,
      :etaXstage => etaXstage,
      :rhoXstage => rhoXstage,
      :PSoI => PSoI,
      :SuS => SuS,
      :wF => wF,
      :w => w,
      :ages => ages,
      :ageXimmune => ageXimmune,
      :aXi => aXi,
      :demog => demog,
      :demog29 => demog29,
      :U5 => U5,
      :O5 => O5,
      :calO => calO,
      :mu1 => mu1,
      :mu2 => mu2,
      :mu3 => mu3,
      :xiF => xiF_h
    )
end


### Now we make the Adam model from those parameters

mutable struct CohortState
    X::Array{Float64,2}
    pit::Array{Float64,1}
    V::Float64
    ageInDays::Float64
    alpha1::Float64
    alpha2::Float64
end


function emptyX_Adam()
    X = zeros(22, 1)
    X[1] = 1
    CohortState(X, zeros(3), 0, 0, 0, 0)
end


mutable struct PopulationState
    X::Array{Float64,2}
    pit::Array{Float64,2}
    V::Float64
    ageInDays::Array{Float64,1}
    alpha1::Float64
    alpha2::Float64
end


function emptyXX_Adam(; L = 66)
    X = zeros(22, L)
    X[1, :] .= 1
    pit = zeros(3, L)
    ageInDays = zeros(L)
    PopulationState(X, pit, 0, ageInDays, 0, 0)
end


"""
Acts on members of the cohort state.
"""
function calP_Adam(alpha, pit, V, params)
    SuS = params[:SuS]
    PoI = params[:PoI]
    alpha = minimum(alpha, 1)
    xi = params[:xiF](pit, V)
    E1 = SuS * vcat(0, xi)
    e1 = sum(E1; dims=(1,))
    e1[e1 == 0] .= 1
    hatxi = E1 * Diagonal(1 ./ e1)
    Deta = minimum.(rho * etaXstage .* rhoXstage, 1)
    vecD = Deta .+ (1 .- Deta) .* Dd
    notTreat = Diagonal(1 .- vecD)
    treat = vcat(vecD, vecD)
    B1 = PSoI * notTreat
    B2 = hatxi * notTreat
    chi = vcat([1], cumsum(xi))
    calP = hcat(
        vcat((1 .- alpha .* chi) .* B1, (1 -alpha) .* B2),
        vcat(alpha * chi .* B1, alpha * B2),
        treat,
        0 * treat
    )
    calP = vcat(calP, zeros(dims(calP)[1], 2))
    calP[22, 21] = 1
    calP[1, 22] = 1
    calP
end


function PxX_Adam!(vecX::CohortState, alpha, params)
    vecX.ageInDays += params[:tau]
    w = params[:wF](ageInDays / 365)
    Xt = calP_Adam(w * vecX.alpha, vecX.pit, vecX.V, params) * vecX.X

    vecX.V = alpha + vecX.V * (1 - params[:mu1])
    vecX.pit[1] = vecX.pit[1] * (1 - params[:mu1]) + vecX.alpha2
    vecX.pit[2] = vecX.pit[2] * (1 - params[:mu2]) + params[:calDimm] * vecX.X
    vecX.pit[3] = vecX.pit[3] * (1 - params[:mu3]) + vecX.alpha2 * params[:ageXimmune](vecX.ageInDays)
    vecX.alpha2 = vecX.alpha1
    vecX.alpha1 = w * alpha
    vecX.X = Xt
end


function cohortXX_Adam(alpha, params; mx = 2620)
    # Grow a cohort so you can assign its stages to a population.
    vecX = emptyX_Adam()
    XX = emptyXX_Adam(mx)
    vecX.V = alpha / params[:mu1]
    for i in 1:mx
        XX.X[:, i] .= vecX.X
        XX.pit[:, i] = vecX.pit
        XX.ageInDays[i] = vecX.ageInDays
        PxX_Adam!(vecX, alpha, params)
    end
    XX.alpha2 = vecX.alpha2
    XX.alpha1 = vecX.alpha1
    XX.V = vecX.V
    XX
end


function cohort2ages_Adam(cXX::PopulationState, params)
    aa = params[:tau] * collect(1:size(cXX.X)[2]) / 365
    XX = emptyXX_Adam()
    XX.V = cXX.V
    XX.ageInDays = zeros(66)
    ix = (aa .< params[:ages][1])
    if length(ix) == 1
        XX.X[:, 1] .= cXX.X[:, ix]
        XX.pit[:, 1] = cXX.pit[:, ix]
        XX.ageInDays[1] = cXX.ageInDays[ix]
    else  # length(ix) > 1
        XX.X[:, 1] .= sum(cXX.X[:, ix], dims = 1) / length(ix)
        XX.pit[:, 1] .= sum(cXX.pit[:, ix], dims = 1) / length(ix)
        XX.ageInDays[1] = mean(cXX.ageInDays[ix])
    end

    for i in 2:length(ages)
        ix = (aa .> ages[i - 1]) & (aa .< ages[i])
        if length(ix) == 1
            XX.X[:, i] = cXX.X[:, ix]
            XX.pit[:, i] = cXX.pit[:, ix]
            XX.ageInDays[i] = cXX.ageInDays[ix]
        else
            XX.X[:, i] = sum(cXX.X[:, ix], dims = 1) / length(ix)
            XX.pit[:, i] = sum(cXX.pit[:, ix], dims = 1) / length(ix)
            XX.ageInDays[i] = cXX.ageInDays[ix]
        end
    end
    XX
end


function PtimesX_Adam!(XX::PopulationState, alpha, params)
    ageInDays = (XX.ageInDays + params[:tau]) * params[:calO]
    Xt = copy(XX.X)
    for a in 1:66
        Xt[:, a] .= calP_Adam(
            alpha * params[:w][a],
            XX.pit[:, a],
            XX.V,
            params
            ) * XX.X[:, a]
    end
    Xt *= params[:calO]
    Xt[1, 1] = 1

    XX.V = alpha + XX.V * (1 - params[:mu1])
    XX.pit[1] = XX.pit[1] * (1 - params[:mu1]) + alpha
    XX.pit[2] = XX.pit[2] * (1 - params[:mu2]) + params[:calDimm] * XX.X
    XX.pit[3] = XX.pit[3] * (1 - params[:mu3]) + XX.alpha2 * params[:aXi]

    XX.pit *= params[:calO]
    XX.ageInDays *= params[:calO]
    XX.alpha2 = XX.alpha1
    XX.alpha1 = alpha
    XX.X = Xt
end


function XX2pr29_Adam(XX, params)
    params[:calDlm] * XX.X * params[:demog29]
end


function ar2stableWave_Adam(alpha, params; tol = 0.01)
    cXX = cohortXX_Adam(alpha, params, 2920)
    XX = cohort2ages_Adam(cXX, params)
    for i in 1:100
        PtimesX_Adam!(XX, alpha, params)
    end
    x = XX2pr29_Adam(XX, params)
    df = 1
    while df > tol
        XXl = copy(XX)
        PtimesX_Adam!(XX, alpha, params)
        x = XX2pr29_Adam(XX, params)
        df = sum(abs(XXl$X - XX$X))
    end
    XX
end


function ar2pr_Adam(alpha, params)
    XX = ar2stableWave_Adam(alpha[1], params)
    XX2pr29_Adam(XX, params)
end


function test_ar2pr_Adam()
    params = make_parameters_adam()
    ar2pr_Adam(0.3, params)
end

"""
    pr2arSS_Adam(x, params)

Given a PfPR value and model parameters, return the attack rate
that would produce this PfPR for a population. This is an optimization
step.
"""
function pr2arSS_Adam(x, params)
    ## ar2pr_Adam
end


function pr2ar_Adam(xx, params; alpha0 = nothing, XX0 = nothing, tol = 0.01)
    if alpha0 === nothing
        alpha = pr2arSS_Adam(xx[1], params)
    else
        alpha = alpha0
    end
end
