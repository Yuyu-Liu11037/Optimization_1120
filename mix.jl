using JuMP, CPLEX, Distributions, LinearAlgebra

function solveMixed01WithoutCuts(m::Int, epsilon::Float64)
    n = 100
    mu = rand(Uniform(0,100),1,n+m)
    a = rand(Uniform(0,100),1,n+m)
    b = 0.5*sum(a)
    sigma = Vector{Float64}(undef, n+m)
    for i in 1:n+m
        sigma[i] = rand(Uniform(0,mu[i]))
    end

    model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_MIP_Strategy_HeuristicEffort" => 0 , "CPXPARAM_Preprocessing_Presolve" => 0 , "CPX_PARAM_MIPSEARCH" => 1 , "CPX_PARAM_MIRCUTS" => -1 , "CPX_PARAM_FRACCUTS" => -1 , "CPX_PARAM_FLOWCOVERS" => -1))

    @variable(model, x[1:n] , Bin)
    @variable(model, y[1:m] >= 0)
    @variable(model, c >= 0)

    @objective(model, Max, sum(mu[i] * x[i] for i in (1:n)') + sum(mu[i] * y[i] for i in (1:m)') - c)

    @constraint(model, dot(a , [x; y]) <= b)
    @constraint(model, c^2 <= ((1-epsilon)/epsilon) * ((sum(sigma[i]^2 * x[i] for i in (1:n)')) + (sum(sigma[i+n]^2 * y[i] for i in (1:m)'))))
    @constraint(model, y .<= 1)

    JuMP.optimize!(model)
end
    