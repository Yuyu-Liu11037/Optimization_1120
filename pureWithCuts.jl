using JuMP, CPLEX, Distributions, LinearAlgebra, MathOptInterface


function separate(n::Int , x_val::Vector{Float64})
    function h(x::AbstractArray)
        return sum(mu .* x') - Omega * norm(sigma .* x)
    end
    #compute t = value of objective function at current node
    t = h(x_val)
    #sort entries in x_val, get an array R
    R = sortperm(x_val, rev = true)
    #compute pi
    pi = Vector{Float64}(undef, n)
    x = zeros(n)
    function indicator(R::Vector{Int64}, i::Int, x::Vector{Float64})
        for j in 1:i
            x[R[j]] = 1
        end
        return x
    end
    pi[R[1]] = h(indicator(R , 1 , x))
    for i in 2:n
        pi[R[i]] = h(indicator(R , i , x)) - h(indicator(R , i-1 , x))
    end
    #if pi' * x > t then return pi' * x - t
    #else return
    if pi' * x_val > t
        return pi' * x_val - t
    else
        return 
    end
end

model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_MIP_Strategy_HeuristicEffort" => 0))
MOI.set(model, MOI.NumberOfThreads(), 1)

@variable(model, x[1:n], Bin)
@variable(model, c >= 0)
@objective(model, Max, sum(mu[i] * x[i] for i in (1:n)') - Omega * c)
@constraint(model, dot(a,x) <= b)
@constraint(model, c^2 >= sum(sigma[i]^2 * x[i]^2 for i in (1:n)'))

function myCallBackFunction(cb_data::CPLEX.CallbackContext , context_id::Clong)
    if context_id != CPX_CALLBACKCONTEXT_RELAXATION
        return
    end
    
    CPLEX.load_callback_variable_primal(cb_data, context_id)

    data_p = Ref{Cint}()
    ret = CPXcallbackgetinfoint(cb_data , CPXCALLBACKINFO_NODECOUNT , data_p)
    if ret != 0
        @warn "error retrieving node_count"
    end
    node_count = data_p[]
    if node_count != 0
        return
    end

    x_val = zeros(Float64 , n)
    for i in 1:n
        x_val[i] = callback_value(cb_data,x[i])
    end
    expr = separate(n, x_val)
    if expr !== nothing
        MOI.submit(model, MOI.UserCut(cb_data), @build_constraint(expr <= 0))
    end
end

MOI.set(model, CPLEX.CallbackFunction(), myCallBackFunction)
optimize!(model)