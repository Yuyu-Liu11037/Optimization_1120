using JuMP, CPLEX, Distributions, LinearAlgebra

model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_MIP_Strategy_HeuristicEffort" => 0))
MOI.set(model, MOI.NumberOfThreads(), 1)
function separate(n::Int64 , x_val , c_val)
    
end

function solvePure01WithCuts(n::Int, epsilon::Float64)
    mu = rand(Uniform(0,100),1,n)
    a = rand(Uniform(0,100),1,n)
    b = 0.5*sum(a)
    sigma = Vector{Float64}(undef, n)
    for i in 1:n
        sigma[i] = rand(Uniform(0,mu[i]))
    end

    @variable(model, x[1:n], Bin)
    @variable(model, c)
    @objective(model, Max, sum(mu[i] * x[i] for i in (1:n)') - c)
    @constraint(model, dot(a,x) <= b)
    @constraint(model, c^2 >= ((1-epsilon)/epsilon) * sum(sigma[i]^2 * x[i]^2 for i in (1:n)'))

    function myCallBackFunction(cb_data::CPLEX.CallbackContext, context_id::Core.Int64)
        if context_id != CPX_CALLBACKCONTEXT_CANDIDATE
            return
        end
        CPLEX.load_callback_variable_primal(cb_data , context_id)
        x_val = callback_value(cb_data , x)
        c_val = callback_value(cb_data , c)
        
        expr = separate()
        if expr !== nothing
            MOI.submit(model, MOI.UserCut(cb_data), @build_constraint(expr <= 0))
        end
    end
    MOI.set(model, CPLEX.CallbackFunction(), myCallBackFunction)
    optimize!(model)
end