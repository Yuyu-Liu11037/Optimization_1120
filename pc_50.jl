using JuMP, CPLEX, Distributions, LinearAlgebra, MathOptInterface

n = 50
mu=[96.34129658828792 60.312045515854365 93.50583519404879 42.91728829922955 43.20953153344804 93.26665821157356 0.27997918402923316 9.751440312740256 81.69414715309865 9.268382769549188 42.93657105389784 20.649012686945845 72.58548950619951 70.28177767218128 31.482704014586858 88.88762770057438 30.683840718861376 53.68859244632045 67.8781633262335 6.849332993813673 89.24451803235304 21.118623825912564 9.238604000691996 26.56051480107966 69.78740045647841 84.65025717824867 3.536481488483978 80.85293986011737 96.1577830845788 72.58290519016434 7.700539712137355 65.92772649937382 63.57835148085537 46.40181628150889 34.344140795704405 22.639910080694538 48.42841169190737 54.29124733479712 97.26534209336705 45.94021192719716 96.19966477274386 34.49120532278205 65.26640686789818 95.44071425993934 57.88142844501417 56.65819005014415 16.521670823661093 97.25369528697803 27.616478393718847 26.835920795386947]
a=[51.922574994551 38.18726682775318 41.151964283393404 41.44120656250955 23.781504255390296 47.1754425308205 92.34628803310977 75.39993378056694 89.73130902406172 97.93676248330708 81.1549760803715 77.62188865485277 80.96828976119146 49.91953293721539 26.523143381417025 4.963558327685391 37.61182760100852 87.44032277519547 42.865077849298316 97.14738444788813 22.312729761338325 73.27922364440639 89.40610020137237 76.47149000263875 42.15582790703727 3.4447277590617276 49.78447534649786 97.42161673777751 27.306688914939247 23.892967746775383 94.96413577277566 76.95460962993195 16.747256888270023 59.4055481455407 68.41939159730967 44.68563515494473 41.11570221922091 62.70923139630623 47.906343921000214 48.284270991612225 81.70693677330078 68.30831138182188 54.300840226956446 58.97524496450901 47.894819003411385 67.68600900016457 82.54953107342907 66.87515576133349 20.066606750607228 2.7154765277263238]
b=1401.5185798968023
sigma=[12.672369370113515, 2.5492554820143707, 41.31754919675693, 41.439943138675645, 6.409054722666487, 66.28055721029926, 0.1992493489691159, 8.693472407255896, 22.021324155433025, 6.500442592974195, 6.638839754298063, 7.809980917533933, 25.98735827136848, 65.60601359107483, 24.560124373663957, 36.42248160284546, 20.95335699092807, 17.050934476270847, 42.66757723551531, 6.345868814477874, 81.45290641258694, 11.847031901664554, 4.153800237387081, 25.482995628087917, 32.4131173085033, 3.7653005449187944, 3.2966876040214697, 34.48581361157358, 40.227454777101805, 4.205249894416882, 6.251672459833214, 5.242562073872441, 30.662762325585007, 15.277391724683032, 8.078852436272513, 13.635426593382489, 4.020571375293106, 50.88082513501987, 53.56926883418995, 24.082049559780977, 76.51330461354907, 15.118891901098927, 52.3608616308796, 10.315853258361921, 40.91499917364862, 8.318861073606232, 7.889387170571406, 31.150835177349066, 18.297891994755417, 1.8654472701844027]
Omega=4.358898943540673

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
@variable(model, x26, Bin)
@variable(model, x27, Bin)
@variable(model, x28, Bin)
@variable(model, x29, Bin)
@variable(model, x30, Bin)
@variable(model, x31, Bin)
@variable(model, x32, Bin)
@variable(model, x33, Bin)
@variable(model, x34, Bin)
@variable(model, x35, Bin)
@variable(model, x36, Bin)
@variable(model, x37, Bin)
@variable(model, x38, Bin)
@variable(model, x39, Bin)
@variable(model, x40, Bin)
@variable(model, x41, Bin)
@variable(model, x42, Bin)
@variable(model, x43, Bin)
@variable(model, x44, Bin)
@variable(model, x45, Bin)
@variable(model, x46, Bin)
@variable(model, x47, Bin)
@variable(model, x48, Bin)
@variable(model, x49, Bin)
@variable(model, x50, Bin)

@variable(model, c >= 0)
@objective(model,
           Max, 
           mu[1]*x1 + mu[2]*x2 + mu[3]*x3 + mu[4]*x4 + mu[5]*x5 + mu[6]*x6 + mu[7]*x7 + mu[8]*x8 + mu[9]*x9 + mu[10]*x10 + mu[11]*x11 + mu[12]*x12 + mu[13]*x13 + mu[14]*x14 + mu[15]*x15 + mu[16]*x16 + mu[17]*x17 + mu[18]*x18 + mu[19]*x19 + mu[20]*x20 + mu[21]*x21 + mu[22]*x22 + mu[23]*x23 + mu[24]*x24 + mu[25]*x25 + mu[26]*x26 + mu[27]*x27 + mu[28]*x28 + mu[29]*x29 + mu[30]*x30 + mu[31]*x31 + mu[32]*x32 + mu[33]*x33 + mu[34]*x34 + mu[35]*x35 + mu[36]*x36 + mu[37]*x37 + mu[38]*x38 + mu[39]*x39 + mu[40]*x40 + mu[41]*x41 + mu[42]*x42 + mu[43]*x43 + mu[44]*x44 + mu[45]*x45 + mu[46]*x46 + mu[47]*x47 + mu[48]*x48 + mu[49]*x49 + mu[50]*x50 - Omega * c)
@constraint(model, a[1]*x1 + a[2]*x2 + a[3]*x3 + a[4]*x4 + a[5]*x5 + a[6]*x6 + a[7]*x7 + a[8]*x8 + a[9]*x9 + a[10]*x10 + a[11]*x11 + a[12]*x12 + a[13]*x13 + a[14]*x14 + a[15]*x15 + a[16]*x16 + a[17]*x17 + a[18]*x18 + a[19]*x19 + a[20]*x20 + a[21]*x21 + a[22]*x22 + a[23]*x23 + a[24]*x24 + a[25]*x25 +  + a[26]*x26 + a[27]*x27 + a[28]*x28 + a[29]*x29 + a[30]*x30 + a[31]*x31 + a[32]*x32 + a[33]*x33 + a[34]*x34 + a[35]*x35 + a[36]*x36 + a[37]*x37 + a[38]*x38 + a[39]*x39 + a[40]*x40 + a[41]*x41 + a[42]*x42 + a[43]*x43 + a[44]*x44 + a[45]*x45 + a[46]*x46 + a[47]*x47 + a[48]*x48 + a[49]*x49 + a[50]*x50 <= b)
@constraint(model, c^2 >= sigma[1]^2*x1^2 + sigma[2]^2*x2^2 + sigma[3]^2*x3^2 + sigma[4]^2*x4^2 + sigma[5]^2*x5^2 + sigma[6]^2*x6^2 + sigma[7]^2*x7^2 + sigma[8]^2*x8^2 + sigma[9]^2*x9^2 + sigma[10]^2*x10^2 + sigma[11]^2*x11^2 + sigma[12]^2*x12^2 + sigma[13]^2*x13^2 + sigma[14]^2*x14^2 + sigma[15]^2*x15^2 + sigma[16]^2*x16^2 + sigma[17]^2*x17^2 + sigma[18]^2*x18^2 + sigma[19]^2*x19^2 + sigma[20]^2*x20^2 + sigma[21]^2*x21^2 + sigma[22]^2*x22^2 + sigma[23]^2*x23^2 + sigma[24]^2*x24^2 + sigma[25]^2*x25^2 +  + sigma[26]^2*x26^2 + sigma[27]^2*x27^2 + sigma[28]^2*x28^2 + sigma[29]^2*x29^2 + sigma[30]^2*x30^2 + sigma[31]^2*x31^2 + sigma[32]^2*x32^2 + sigma[33]^2*x33^2 + sigma[34]^2*x34^2 + sigma[35]^2*x35^2 + sigma[36]^2*x36^2 + sigma[37]^2*x37^2 + sigma[38]^2*x38^2 + sigma[39]^2*x39^2 + sigma[40]^2*x40^2 + sigma[41]^2*x41^2 + sigma[42]^2*x42^2 + sigma[43]^2*x43^2 + sigma[44]^2*x44^2 + sigma[45]^2*x45^2 + sigma[46]^2*x46^2 + sigma[47]^2*x47^2 + sigma[48]^2*x48^2 + sigma[49]^2*x49^2 + sigma[50]^2*x50^2)

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
    x_val[26] = callback_value(cb_data , x26)
    x_val[27] = callback_value(cb_data , x27)
    x_val[28] = callback_value(cb_data , x28)
    x_val[29] = callback_value(cb_data , x29)
    x_val[30] = callback_value(cb_data , x30)
    x_val[31] = callback_value(cb_data , x31)
    x_val[32] = callback_value(cb_data , x32)
    x_val[33] = callback_value(cb_data , x33)
    x_val[34] = callback_value(cb_data , x34)
    x_val[35] = callback_value(cb_data , x35)
    x_val[36] = callback_value(cb_data , x36)
    x_val[37] = callback_value(cb_data , x37)
    x_val[38] = callback_value(cb_data , x38)
    x_val[39] = callback_value(cb_data , x39)
    x_val[40] = callback_value(cb_data , x40)
    x_val[41] = callback_value(cb_data , x41)
    x_val[42] = callback_value(cb_data , x42)
    x_val[43] = callback_value(cb_data , x43)
    x_val[44] = callback_value(cb_data , x44)
    x_val[45] = callback_value(cb_data , x45)
    x_val[46] = callback_value(cb_data , x46)
    x_val[47] = callback_value(cb_data , x47)
    x_val[48] = callback_value(cb_data , x48)
    x_val[49] = callback_value(cb_data , x49)
    x_val[50] = callback_value(cb_data , x50)

    expr = separate(n, x_val)
    if expr !== nothing
        MOI.submit(model, MOI.UserCut(cb_data), @build_constraint(expr <= 0))
    end
end

MOI.set(model, CPLEX.CallbackFunction(), myCallBackFunction)
optimize!(model)

println(JuMP.objective_value(model))