using JuMP, CPLEX, Distributions, LinearAlgebra, MathOptInterface

function solvePure01WithCuts(n::Int, epsilon::Float64)
    mu = rand(Uniform(0,100),1,n)
    a = rand(Uniform(0,100),1,n)
    b = 0.5*sum(a)
    sigma = Vector{Float64}(undef, n)
    for i in 1:n
        sigma[i] = rand(Uniform(0,mu[i]))
    end
    Omega = sqrt((1-epsilon)/epsilon)

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

    #generate cuts at the root node
    model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_MIP_Strategy_HeuristicEffort" => 0 , "CPXPARAM_MIP_Limits_Nodes" => 0))
    MOI.set(model, MOI.NumberOfThreads(), 1)

    @variable(model, x[1:n], Bin)
    @variable(model, c >= 0)
    @objective(model, Max, sum(mu[i] * x[i] for i in (1:n)') - Omega * c)
    @constraint(model, dot(a,x) <= b)
    @constraint(model, c^2 >= sum(sigma[i]^2 * x[i]^2 for i in (1:n)'))

    #cb_calls = Int32[]
    function myCallBackFunction(cb_data::CPLEX.CallbackContext , context_id::Clong)
        #=push!(cb_calls , context_id)
        if context_id != CPX_CALLBACKCONTEXT_RELAXATION
            return
        end
        CPLEX.load_callback_variable_primal(cb_data,context_id)
        data_p = Ref{Cint}()
        status = CPXgetcallbacknodeinfo(ENV , cb_data.context_id , context_id , 0 , CPXCALLBACKINFO_NODECOUNT , data_p)
        node_count = data_p[]
        if node_count != 0
            return
        end
        =#
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

    set_optimizer_attributes(model , "CPXPARAM_MIP_Limits_Nodes" => 9223372036800000000)
    MOI.set(model, CPLEX.CallbackFunction(), nothing)
    optimize!(model)
end