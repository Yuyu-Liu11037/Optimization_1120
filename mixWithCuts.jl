using JuMP, CPLEX, Distributions, LinearAlgebra

function solveMixed01WithCuts(m::Int, epsilon::Float64)
    n = 100
    mu = rand(Uniform(0,100),1,n+m)
    a = rand(Uniform(0,100),1,n+m)
    b = 0.5*sum(a)
    sigma = Vector{Float64}(undef, n+m)
    for i in 1:(n+m)
        sigma[i] = rand(Uniform(0,mu[i]))
    end
    Omega = sqrt((1-epsilon)/epsilon)

    model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_MIP_Strategy_HeuristicEffort" => 0))
    MOI.set(model, MOI.NumberOfThreads(), 1)
    @variable(model, 0 <= x[1:n] <= 1)
    @variable(model, 0 <= y[1:m] <= 1)
    @variable(model, c >= 0)
    @objective(model, Max, sum(mu[i] * x[i] for i in (1:n)') + sum(mu[i] * y[i] for i in (1:m)') - Omega * c)
    @constraint(model, dot(a , [x; y]) <= b)
    @constraint(model, c^2 >= sum(sigma[i]^2 * x[i]^2 for i in (1:n)') + sum(sigma[i+n]^2 * y[i]^2 for i in (1:m)'))
    @constraint(model, y .<= 1)
    function separate(m::Int , xy_val::Vector{Float64} , T::Vector{Int})
        T_indicator = collect(1:m)
        for i in 1:m
            if T[i] == 0
                T_indicator = 0
            else 
                T_indicator = 1
            end
        end
 
        function f_T(x::AbstractArray)
            return sum(mu[i] * x[i] for i in 1:100) + sum(mu[i+100] for i in T_indicator) - Omega * sqrt(sum(sigma[i] * x[i] for i in 1:100) + sum(sigma[i+100] for i in T_indicator))
        end
        t = sum(mu[i] * xy_val[i] for i in 1:(100+m)) - Omega * norm(sqrt(sigma[i]) * xy_val[i] for i in 1:(100+m))

        x_val = Vector{Float64}(undef, 100)
        for i in 1:100
            x_val[i] = xy_val[i]
        end
    
        R = sortperm(x_val, rev = true)
        pi = Vector{Float64}(undef, 100)
        x = zeros(100)
        function indicator(R::Vector{Int64}, i::Int, x::Vector{Float64})
            for j in 1:i
                x[R[j]] = 1
            end
            return x
        end
        pi[R[1]] = f_T(indicator(R , 1 , x))
        for i in 2:100
            pi[R[i]] = f_T(indicator(R , i , x)) - f_T(indicator(R , i - 1 , x))
        end

        temp = Vector{Float64}(undef, 101)
        for i in 1:(100)
            temp[i] = pi[i]
        end
        temp[101] = t
        T_nonzero = filter(!iszero,T)
        if pi' * x_val + sum(mu[i+100] * xy_val[i+100] for i in 1:m) + norm(sqrt(sigma[i+100]) .* xy_val[i+100] for i in T_nonzero) > t
            return temp
        else
            return 
        end
    end
    xy_0 = Vector{Float64}(undef, 100+m)
    while true
        optimize!(model)
        for i in 1:(100)
            xy_0[i] = JuMP.value(x[i])
        end
        for i in 101:(100+m)
            xy_0[i] = JuMP.value(y[i-100])
        end
        T = collect(1:m)
        for i in 1:m
            if xy_0[i+100] <= 0.33
                T[i] = 0
            elseif xy_0[i+100] >= 0.66
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

        if (separate(m, xy_0 , T) === nothing)
            break
        else
            T_nonzero = filter(!iszero,T)
            MOI.add_constraints(model, x , sum(separate(m, xy_0 , T)[i] * x[i] for i in 0:100) + sum(mu[i+100] * y[i] for i in 1:m) + norm(sqrt(sigma[i+100]) .* y[i] for i in T_nonzero) <= separate(m, xy_0 , T)[101])
            continue
        end
    end




    model2 = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_MIP_Strategy_HeuristicEffort" => 0))
    MOI.set(model2, MOI.NumberOfThreads(), 1)
    @variable(model2, x[1:n] , Bin)
    @variable(model2, 0 <= y[1:m] <= 1)
    @variable(model2, c >= 0)
    for i in 1:(100)
        set_start_value(x[i],xy_0[i])
    end
    for i in 101:(100+m)
        set_start_value(y[i-100],xy_0[i])
    end
    @objective(model2, Max, sum(mu[i] * x[i] for i in (1:n)') + sum(mu[i] * y[i] for i in (1:m)') - Omega * c)
    @constraint(model2, dot(a , [x; y]) <= b)
    @constraint(model2, c^2 >= sum(sigma[i]^2 * x[i]^2 for i in (1:n)') + sum(sigma[i+n]^2 * y[i]^2 for i in (1:m)'))
    @constraint(model2, y .<= 1)
    optimize!(model2)
end