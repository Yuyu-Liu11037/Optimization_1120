using JuMP, CPLEX, Distributions, LinearAlgebra

function solvePure01WithoutCuts(n::Int, epsilon::Float64)
    mu = rand(Uniform(0,100),1,n)
    a = rand(Uniform(0,100),1,n)
    b = 0.5*sum(a)
    sigma = Vector{Float64}(undef, n)
    for i in 1:n
        sigma[i] = rand(Uniform(0,mu[i]))
    end
    Omega = sqrt((1-epsilon)/epsilon)

    model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_MIP_Strategy_HeuristicEffort" => 0))
    MOI.set(model, MOI.NumberOfThreads(), 1)

    @variable(model, x[1:n], Bin)
    @variable(model, c >= 0)
    @objective(model, Max, sum(mu[i] * x[i] for i in (1:n)') - Omega * c)
    @constraint(model, dot(a,x) <= b)
    @constraint(model, c^2 >= sum(sigma[i]^2 * x[i]^2 for i in (1:n)'))

    optimize!(model)
end