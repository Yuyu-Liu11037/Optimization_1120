using JuMP, CPLEX, Distributions, LinearAlgebra, MathOptInterface

n = 25
mu=[45.69250547342361 75.42385196287698 99.29292649272264 90.12534613807001 26.17284108881387 37.376518954546654 6.679093834213035 20.789717613594117 33.36257056179638 57.15630806045373 97.2117285060904 47.29468363862369 22.392413570844138 20.75640243451883 11.189092604419837 12.261574245224649 46.372574667266065 95.66966310532837 49.99524583912126 48.85535458859736 79.66098491982329 3.762350182946328 60.14339056474354 22.742629392575488 28.690046800538948]
a=[6.947049091686441 34.10597780627774 44.869004669879054 29.860941023080446 60.757472049822326 27.37529499615642 80.9566724381156 8.659148232605684 89.96085644428622 27.190903313825043 50.335534246019286 69.52035018393681 39.396935970853875 40.23794523495874 61.34605733580023 72.85292597978584 53.727823947058695 23.06738181710487 63.485099641386135 42.77497173301993 9.466426876441735 79.09046294067822 93.03621463054364 5.765907772190015 44.48302939883511]
b=579.635193887174
sigma=[30.72637049000794, 34.51003151310712, 5.2563796385928985, 6.019544279350812, 24.22232555767528, 5.46904178227805, 4.623200025597562, 3.438798638458221, 32.69160599008206, 30.938286091797252, 54.97218631365624, 43.62260159645204, 8.12032690914649, 19.003759872706425, 0.31716425369448525, 11.882028217012724, 7.8043550908184365, 15.914865139050015, 45.720867773599814, 20.805877874194547, 58.525512919672494, 1.8791724928999844, 9.225389858848864, 8.606539222558844, 20.45761004837807]
Omega=9.9498743710662

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

@variable(model, x1, Bin)
@variable(model, x2, Bin)
@variable(model, x3, Bin)
@variable(model, x4, Bin)
@variable(model, x5, Bin)
@variable(model, x6, Bin)
@variable(model, x7, Bin)
@variable(model, x8, Bin)
@variable(model, x9, Bin)
@variable(model, x10, Bin)
@variable(model, x11, Bin)
@variable(model, x12, Bin)
@variable(model, x13, Bin)
@variable(model, x14, Bin)
@variable(model, x15, Bin)
@variable(model, x16, Bin)
@variable(model, x17, Bin)
@variable(model, x18, Bin)
@variable(model, x19, Bin)
@variable(model, x20, Bin)
@variable(model, x21, Bin)
@variable(model, x22, Bin)
@variable(model, x23, Bin)
@variable(model, x24, Bin)
@variable(model, x25, Bin)

@variable(model, c >= 0)
@objective(model,
           Max, 
           mu[1]*x1 + mu[2]*x2 + mu[3]*x3 + mu[4]*x4 + mu[5]*x5 + mu[6]*x6 + mu[7]*x7 + mu[8]*x8 + mu[9]*x9 + mu[10]*x10 + mu[11]*x11 + mu[12]*x12 + mu[13]*x13 + mu[14]*x14 + mu[15]*x15 + mu[16]*x16 + mu[17]*x17 + mu[18]*x18 + mu[19]*x19 + mu[20]*x20 + mu[21]*x21 + mu[22]*x22 + mu[23]*x23 + mu[24]*x24 + mu[25]*x25 - Omega * c)
@constraint(model, a[1]*x1 + a[2]*x2 + a[3]*x3 + a[4]*x4 + a[5]*x5 + a[6]*x6 + a[7]*x7 + a[8]*x8 + a[9]*x9 + a[10]*x10 + a[11]*x11 + a[12]*x12 + a[13]*x13 + a[14]*x14 + a[15]*x15 + a[16]*x16 + a[17]*x17 + a[18]*x18 + a[19]*x19 + a[20]*x20 + a[21]*x21 + a[22]*x22 + a[23]*x23 + a[24]*x24 + a[25]*x25 <= b)
@constraint(model, c^2 >= sigma[1]^2*x1^2 + sigma[2]^2*x2^2 + sigma[3]^2*x3^2 + sigma[4]^2*x4^2 + sigma[5]^2*x5^2 + sigma[6]^2*x6^2 + sigma[7]^2*x7^2 + sigma[8]^2*x8^2 + sigma[9]^2*x9^2 + sigma[10]^2*x10^2 + sigma[11]^2*x11^2 + sigma[12]^2*x12^2 + sigma[13]^2*x13^2 + sigma[14]^2*x14^2 + sigma[15]^2*x15^2 + sigma[16]^2*x16^2 + sigma[17]^2*x17^2 + sigma[18]^2*x18^2 + sigma[19]^2*x19^2 + sigma[20]^2*x20^2 + sigma[21]^2*x21^2 + sigma[22]^2*x22^2 + sigma[23]^2*x23^2 + sigma[24]^2*x24^2 + sigma[25]^2*x25^2)

function myCallBackFunction(cb_data::CPLEX.CallbackContext , context_id::Clong)
    if context_id != CPX_CALLBACKCONTEXT_CANDIDATE
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
    x_val[1] = callback_value(cb_data , x1)
    x_val[2] = callback_value(cb_data , x2)
    x_val[3] = callback_value(cb_data , x3)
    x_val[4] = callback_value(cb_data , x4)
    x_val[5] = callback_value(cb_data , x5)
    x_val[6] = callback_value(cb_data , x6)
    x_val[7] = callback_value(cb_data , x7)
    x_val[8] = callback_value(cb_data , x8)
    x_val[9] = callback_value(cb_data , x9)
    x_val[10] = callback_value(cb_data , x10)
    x_val[11] = callback_value(cb_data , x11)
    x_val[12] = callback_value(cb_data , x12)
    x_val[13] = callback_value(cb_data , x13)
    x_val[14] = callback_value(cb_data , x14)
    x_val[15] = callback_value(cb_data , x15)
    x_val[16] = callback_value(cb_data , x16)
    x_val[17] = callback_value(cb_data , x17)
    x_val[18] = callback_value(cb_data , x18)
    x_val[19] = callback_value(cb_data , x19)
    x_val[20] = callback_value(cb_data , x20)
    x_val[21] = callback_value(cb_data , x21)
    x_val[22] = callback_value(cb_data , x22)
    x_val[23] = callback_value(cb_data , x23)
    x_val[24] = callback_value(cb_data , x24)
    x_val[25] = callback_value(cb_data , x25)

    expr = separate(n, x_val)
    if expr !== nothing
        MOI.submit(model, MOI.UserCut(cb_data), @build_constraint(expr <= 0))
    end
end

MOI.set(model, CPLEX.CallbackFunction(), myCallBackFunction)
optimize!(model)

println(JuMP.objective_value(model))