using JuMP, CPLEX, Distributions, LinearAlgebra

function solvePure01WithoutCuts(n::Int, epsilon::Float64)
    mu = rand(Uniform(0,100),1,n)
    a = rand(Uniform(0,100),1,n)
    b = 0.5*sum(a)
    sigma = Vector{Float64}(undef, n)
    for i in 1:n
        sigma[i] = rand(Uniform(0,mu[i]))
    end

    model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_MIP_Strategy_HeuristicEffort" => 0 , "CPXPARAM_Preprocessing_Presolve" => 0 , "CPX_PARAM_MIPSEARCH" => 1 , "CPX_PARAM_MIRCUTS" => -1 , "CPX_PARAM_FRACCUTS" => -1 , "CPX_PARAM_FLOWCOVERS" => -1))
    
    @variable(model, x[1:n], Bin)
    @variable(model, c >= 0)

    @objective(model, Max, sum(mu[i] * x[i] for i in (1:n)') - c)

    @constraint(model, dot(a,x) <= b)
    @constraint(model, c^2 == ((1-epsilon)/epsilon) * sum(sigma[i]^2 * x[i] for i in (1:n)'))

    JuMP.optimize!(model)
    print(objective_value)
end
    