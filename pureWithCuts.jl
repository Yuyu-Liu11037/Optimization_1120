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

    model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_MIP_Strategy_HeuristicEffort" => 0))
    MOI.set(model, MOI.NumberOfThreads(), 1)
    @variable(model, 0 <= x[1:n] <= 1)
    @variable(model, c >= 0)
    @objective(model, Max, sum(mu[i] * x[i] for i in (1:n)') - Omega * c)
    @constraint(model, dot(a,x) <= b)
    @constraint(model, c^2 >= sum(sigma[i]^2 * x[i]^2 for i in (1:n)'))
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
        temp = Vector{Float64}(undef, n+1)
        for i in 1:n
            temp[i] = pi[i]
        end
        temp[n+1] = t
        if pi' * x_val > t
            return temp
        else
            return 
        end
    end
    x_0 = Vector{Float64}(undef, n)
    while true
        optimize!(model)
        for i in 1:n
            x_0[i] = JuMP.value(x[i])
        end
        if (separate(n, x_0) === nothing)
            break
        else
            MOI.add_constraint(model , x , sum(x[i] * separate(n, x_0)[i] for i in (1:n)') <= separate(n, x_0)[1+n])
            continue
        end
    end
    

    model2 = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_MIP_Strategy_HeuristicEffort" => 0))
    MOI.set(model2, MOI.NumberOfThreads(), 1)
    @variable(model2, x[1:n] , Bin)
    @variable(model2, c >= 0)
    for i in 1:n
        set_start_value(x[i],x_0[i])
    end
    @objective(model2, Max, sum(mu[i] * x[i] for i in (1:n)') - Omega * c)
    @constraint(model2, dot(a,x) <= b)
    @constraint(model2, c^2 >= sum(sigma[i]^2 * x[i]^2 for i in (1:n)'))
    JuMP.optimize!(model2)
end