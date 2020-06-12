

function test_emptyX_Adam()
    X = emptyX_Adam()
    @assert size(X.X) == (22, 1)
    @assert size(X.pit) == (3,)
    @assert isapprox(X.V, 0)
    @assert isapprox(X.ageInDays, 0)
    @assert isapprox(X.alpha1, 0)
    @assert isapprox(X.alpha2, 0)
end

test_emptyX_Adam()

params = make_parameters_adam()

# Test at the ar2pr level.
pr_low = ar2pr_Adam(0.1, params)
@test pr_low < 0.5
pr_high = ar2pr_Adam(0.9, params)
@test pr_high > 0.5
