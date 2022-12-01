using JuMP, CPLEX, Distributions, LinearAlgebra

function solveMixed01WithCuts(m::Int, epsilon::Float64)
    n = 100
    mu = rand(Uniform(0,100),1,n+m)
    a = rand(Uniform(0,100),1,n+m)
    b = 0.5*sum(a)
    sigma = Vector{Float64}(undef, n+m)
    for i in 1:n+m
        sigma[i] = rand(Uniform(0,mu[i]))
    end
    Omega = sqrt((1-epsilon)/epsilon)

    function separate(m::Int , xy_val::Vector{Float64} , T::Vector{Int})
        T_indicator = collect(1:m)
        for i in 1:m
            if T[i] == 0
                T_indicator = 0
            else 
                T_indicator = 1
            end
        end

        function h(x::AbstractArray)
            return sum(mu .* x') - Omega * norm(sigma .* x)
        end
        t = h(xy_val)

        new_xy_val = Vector{Float64}(undef, 100+m)
        for i in 1:100
            new_xy_val[i] = xy_val[i]
        end
        for i in 101:m
            new_xy_val[i] = T_indicator[i-100]
        end
        R = sortperm(new_xy_val, rev = true)
        pi = Vector{Float64}(undef, 100+m)
        x = zeros(100+m)
        function indicator(R::Vector{Int64}, i::Int, x::Vector{Float64})
            for j in 1:i
                x[R[j]] = 1
            end
            return x
        end
        pi[R[1]] = h(indicator(R , 1 , x))
        for i in 2:(100+m)
            pi[R[i]] = h(indicator(R , i , x)) - h(indicator(R , i-1 , x))
        end

        if pi' * xy_val > t
            return pi' * xy_val - t
        else
            return 
        end
    end

    model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_MIP_Strategy_HeuristicEffort" => 0 , "CPXPARAM_MIP_Limits_Nodes" => 0))
    MOI.set(model, MOI.NumberOfThreads(), 1)

    @variable(model, x[1:n] , Bin)
    @variable(model, y[1:m] >= 0)
    @variable(model, c >= 0)
    @objective(model, Max, sum(mu[i] * x[i] for i in (1:n)') + sum(mu[i] * y[i] for i in (1:m)') - Omega * c)
    @constraint(model, dot(a , [x; y]) <= b)
    @constraint(model, c^2 >= sum(sigma[i]^2 * x[i]^2 for i in (1:n)') + sum(sigma[i+n]^2 * y[i]^2 for i in (1:m)'))
    @constraint(model, y .<= 1)

    function myCallBackFunction(cb_data::CPLEX.CallbackContext, context_id::Clong)
        xy_val = zeros(Float64 , 100 + m)  
        for i in 1:100
            xy_val[i] = callback_value(cb_data, x[i])
        end
        for i in 101:(100+m)
            xy_val[i] = callback_value(cb_data, y[i-100])
        end

        #add heuristic to pick T
        T = collect(1:m)
        for i in 1:m
            if xy_val[i+100] <= 0.33
                T[i] = 0
            elseif xy_val[i+100] >= 0.66
                T[i] = T[i]
            else 
                r = rand((0,1))
                if r == 0
                    T[i] = 0
                else
                    T[i] = T[i]
                end
            end
        end

        expr = separate(m, xy_val , T)
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