using JuMP,Distributions,LinearAlgebra

function parameterGeneration(n::Int,epsilon::Float64)
    mu = rand(Uniform(0,100),1,n)
    a = rand(Uniform(0,100),1,n)
    b = 0.5*sum(a)
    sigma = Vector{Float64}(undef, n)
    for i in 1:n
        sigma[i] = rand(Uniform(0,mu[i]))
    end
    Omega = sqrt((1-epsilon)/epsilon)

    print("mu=",mu)
    println()
    print("a=",a)
    println()
    print("b=",b)
    println()
    print("sigma=",sigma)
    println()
    print("Omega=",Omega)
end