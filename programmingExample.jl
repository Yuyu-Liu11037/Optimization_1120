function solveULS(path::String; solver=CplexSolver(), valid::Bool = true) 
    T, c, h, K, d = readULS(path)
    m = Model(solver = solver)
    @defVar(m, z[1:T], Bin)
    @defVar(m, q[i = 1:T] >= 0)
    @defVar(m, s[1:T] >= 0)
 
    @setObjective(m, Min, sum{c[t]*q[t] + K[t]*z[t] + h[t]*s[t], t=1:T})
    @addConstraint(m, activation[t = 1:T], q[t] <= sum(d[t:T])*z[t])
    @addConstraint(m, balance[t = 1:T], (t > 1 ? s[t-1] : 0) + q[t] == s[t] + d[t])
 
    separationtime = 0.
    separated = 0
    called = 0
    function lSgenerator(cb)
       tt = time()
       called += 1
       z_val = getValue(z)
       q_val = getValue(q)
       expr = separate(T, d, z_val, q_val, z, q)
       if expr != nothing
          @addUserCut(cb, expr >= 0)
          separated += 1
       end
       separationtime += time()-tt
    end 
 
    if valid
       addCutCallback(m, lSgenerator)
    end
 
    status = solve(m)
    println("Objective value: ", getObjectiveValue(m))
    println("Separation time: $separationtime seconds")
    println("Separated: $separations")
    status
 end