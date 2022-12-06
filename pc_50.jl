using JuMP, CPLEX, Distributions, LinearAlgebra, MathOptInterface

n = 50
#epsilon=0.05
#= mu=[7.463304985802621 60.290059289592776 45.67293173139571 98.49885938160668 69.39524250490182 64.062591855494 79.6213683235748 98.1257246876236 21.651952126650897 0.8591877179772833 79.82803212858943 54.2023316893392 9.638806142205736 63.23739296260147 61.37918636204854 8.995452912463497 47.002480277352376 46.55827985843328 5.608496717934708 22.232738515066906 68.65353878587011 11.54327327500717 3.0288934992744987 31.859787094007373 4.103034350735058 19.060684873119193 20.842193314936587 70.82889061914855 5.930404023102531 70.33723578403887 92.70480037816084 43.80069471409232 89.80424937107553 67.5895586983341 85.56182224019341 36.791823208986344 57.843832134897625 49.11140682628141 78.63717984846285 5.450550826700118 9.027943615513767 29.04113407053711 7.118392353087866 81.97470636236234 42.8243926940049 12.033878794407183 89.04948126752998 3.1678304789231726 33.26977506285919 49.187749513656]
a=[36.615542184416626 86.98789676002653 42.37073717944699 97.92442930342519 52.687486605779775 60.92916108308779 36.118577529551 11.462288299850876 70.53333425348845 94.30480505238683 81.79806771443955 55.65780374951244 53.85664664755755 63.86278625553741 23.11706272992805 80.89845642770047 62.1689859392732 59.76963968389952 56.82897842624839 0.8187115486856911 27.692881176264482 65.07517572931917 37.43487413833135 63.849682210210524 39.20801146646233 32.197718428818455 26.189806067398468 48.31374651132585 12.213549714034977 30.8666446840683 83.19754998794606 69.20456052718994 92.1883004939269 62.92280695422304 42.83096760917152 34.048612450447614 83.88911057391927 9.034390009935 32.67607108685743 40.03133431266889 69.70332499700089 29.641210897855263 53.62174580572873 69.28396184397087 57.744762658510076 62.574608942853935 73.58959386816512 52.57747781232366 81.58558559662765 22.655957454287044]
b=1317.3777106920427
sigma=[6.6499072993591986, 2.1212229331237764, 30.32879659241038, 0.7332659385111749, 27.84643161814835, 0.9614929250780121, 15.14793864590009, 53.90061881492716, 19.14953734299916, 0.6229787656299595, 35.93616849359539, 14.537636258099136, 5.403576479374238, 40.36776250413815, 23.031110496259615, 7.283329483357539, 19.73575581324186, 10.655700867406459, 3.0955414079692365, 7.197017032209098, 32.120894509911594, 7.4178477927683, 1.1931793582423273, 13.919110574806048, 2.6823854984663, 10.542783751803132, 2.823365939106079, 16.80905941519938, 1.3743682465427518, 63.57788567346229, 58.99550521385476, 25.43817051058899, 50.839588739628674, 38.7176210983329, 59.37722329288694, 31.225638340635353, 19.157137017356874, 13.404329398784597, 4.201863005353369, 0.277696202220773, 3.1765581375544447, 15.185666412491615, 1.536110481437705, 24.84092559705601, 29.451705281223308, 4.831760292600087, 27.951541922987005, 3.0568754022102, 14.940200812662924, 42.77286337098473]
Omega=4.358898943540673 =#

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