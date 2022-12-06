using JuMP, CPLEX, Distributions, LinearAlgebra

n = 100
m = 5
mu=[87.07931051754268 1.885674719768382 91.19989721737635 34.56176127646072 19.55643810168001 58.65196451393668 91.29085333409998 9.672821513223296 37.01604969821951 51.9785280595715 57.24601692936437 34.06979949532624 84.78670438876435 53.63336435367351 87.1252884111481 47.7440166422045 4.8144636930481655 35.90833569864085 10.607478591094521 24.482281968054608 39.150269417699256 83.46304394815066 18.04538884457022 37.90577113522119 73.68614782513113 57.472140743030586 74.06816139205236 95.91707486212266 29.42078369266782 48.136821720971135 10.119072733212231 89.37519880333582 98.54889574645983 66.11210735972625 54.406869918522005 42.79208793235407 15.969022422541768 5.338317392293524 18.849653844473934 31.38991176376863 49.156939322486004 40.563719646662086 32.370266047771764 24.464445241306898 10.989911705578226 2.370896091690955 84.60939765176511 73.55072300821539 99.57752627275413 37.80502524766962 14.046678905296728 66.39443915219906 58.6089978355361 95.50947541458453 26.43768489534546 94.97637476861263 0.9814915745096187 35.821464662781366 60.042203122584645 73.45146238267289 28.277789882000448 15.984345260302057 63.72515274684951 18.559269380449006 42.10088854705586 54.24665057338661 98.4752189280723 23.844311427168986 92.87787118493993 78.97408366351517 0.7754629974755289 10.57231407278081 38.12367815064226 65.38169002802587 18.088913470228675 42.34382847045393 75.13887535728213 41.47786947879314 31.35005682789552 42.00460068187169 99.92321689434036 51.40854263053922 55.129084019797524 45.59520681620868 23.376708439284112 95.73326070864644 83.62370396889465 89.98797281983074 54.947434615144054 69.996883288631 70.02669141544472 25.966246875682643 28.482902986030656 41.29703056738162 72.74520009136593 38.754799321633485 86.06581636892467 58.781751266729245 81.50962545108715 2.144025974637176 59.85037496130696 92.30061687813411 58.62415650703797 64.55162886100568 83.4396968437113]
a=[50.1119995902865 15.132667567399782 4.010341790959949 32.68443851752396 70.6821013782562 43.17163329490044 3.1676351196712083 12.101242739575447 36.071974614903176 3.4100185059406374 17.984417984663693 18.85928418981191 42.850917221929144 67.29722306463239 87.99165111677935 85.73668559430288 30.695428926244528 26.45050291060882 34.5465045129677 61.16699525008232 65.63130888557198 32.01100318901019 99.69339896134137 20.942969651804045 0.33125265786618785 1.8945177833471827 35.237056527735554 14.85574437141397 8.985945618786351 88.77332057302942 52.813201148679134 53.90800079769785 37.031212907207845 29.627203907190214 32.3294904998633 9.415673336022202 82.88215653643168 78.75028698991962 82.82800999900903 1.099200600528638 42.1066703906676 53.98537993716045 7.129257165506786 24.657584978778512 42.51030416829833 66.8460227104298 38.95347271438274 28.68611979966944 87.21334722340482 82.5616988468781 5.775576761939383 6.1905913171843014 85.29258942976988 39.02545112486811 63.76884418421361 37.79100508636964 35.30635436314201 52.06516512188751 97.30137032496232 35.867149121784756 11.931463489635286 25.567764038922924 16.270785914615125 69.19749559449953 56.41431004362608 37.93355063164524 51.156455286986066 69.67871538333354 99.4983306969822 91.26225094190947 19.03261146491827 20.791154754077457 28.3436330216989 13.15546497671044 74.25481101359705 47.900992966449586 96.1235656780894 47.982333105045036 80.89996031847342 10.250810087718376 21.482312107660583 73.5358668673205 9.277578787749496 36.78592262749652 6.9070938089432214 91.53182391496698 40.54804118110981 74.6061914361274 82.65047384455951 0.9500928815269982 72.92133795950625 36.06339114258933 36.371973846836994 2.752702940669305 12.15114134602321 92.85816278331049 80.27135937311417 72.59258812436958 70.9675749964615 13.678653411569897 42.35402433811489 69.95222589674293 1.674166243519959 53.63681746335867 94.40196795079719]
b=2332.3842591442863
sigma=[0.02466635472119416, 0.7751272885646003, 42.849743854825086, 17.816704585644484, 14.127785121902534, 4.89084558984959, 64.68860941200283, 2.3015576962048345, 24.893746590440525, 0.40968592338509924, 13.867135936830303, 23.15183914957121, 27.97788839771546, 47.151613564271806, 41.62609274654936, 28.109797457691876, 2.3488092433007783, 33.40960552782956, 7.915751468945418, 19.120007204514472, 31.512450978627346, 27.10902675513179, 0.32366127311918796, 9.812264764299744, 15.557986305905514, 17.804259852905417, 5.035260795922837, 92.74062276709624, 0.4127910230824477, 15.862126327262938, 7.069641743917185, 72.00774239483813, 33.75378991424091, 60.62539083210399, 48.53588063503427, 12.296646623134666, 12.881642625143478, 2.1466551266443203, 0.12650841245669106, 15.249901637793867, 28.211568285168866, 39.65743649115649, 8.121117829749945, 21.783800446729252, 10.25125237383051, 1.6020245707421281, 69.10655947185933, 71.79728971116913, 19.532749753056454, 37.17207458497488, 9.52857383870254, 0.19369689455366876, 2.872388596868621, 44.832023238431525, 19.781248697468715, 67.21599917247914, 0.5842695665256902, 2.0993026571714077, 37.16570310033916, 64.49874025268134, 7.585170506867685, 13.300095935638875, 4.104093152778556, 17.070474938568857, 11.39055692299904, 14.172507135054728, 6.3902701487598685, 8.683621325796723, 58.29525566360567, 56.27164723595722, 0.17481417821528877, 6.257897789002774, 31.952593103457765, 16.142118234891786, 14.005077559097815, 34.72636783627575, 72.0337316913522, 35.2276486595028, 14.277740350645663, 11.307082863794752, 75.48387126710834, 22.139389035496137, 18.4680923940578, 18.39304309900751, 1.46113700744102, 55.21410430375834, 38.09949745226924, 73.29391201354154, 48.81662291583136, 64.4925115006793, 14.014113219349364, 4.020720421711663, 21.74692067583448, 0.04003753412752601, 65.62281059590936, 24.30373522291825, 58.329154940689804, 5.6585113650707815, 3.964061807165238, 1.7774777444304746, 54.41855767607286, 60.981426131318365, 5.278750666201518, 54.12528735344736, 69.73673278698027]

    #epsilon=0.1
#Omega=3.0
    #epsilon=0.05
Omega=4.358898943540673
    #epsilon=0.04
#Omega=4.898979485566356
    #epsilon=0.03
#Omega=5.686240703077327
    #epsilon=0.02
#Omega=7.0
    #epsilon=0.01
#Omega=9.9498743710662

model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_MIP_Strategy_HeuristicEffort" => 0 , "CPXPARAM_MIP_Strategy_Search" => 1))
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
@variable(model, x51, Bin)
@variable(model, x52, Bin)
@variable(model, x53, Bin)
@variable(model, x54, Bin)
@variable(model, x55, Bin)
@variable(model, x56, Bin)
@variable(model, x57, Bin)
@variable(model, x58, Bin)
@variable(model, x59, Bin)
@variable(model, x60, Bin)
@variable(model, x61, Bin)
@variable(model, x62, Bin)
@variable(model, x63, Bin)
@variable(model, x64, Bin)
@variable(model, x65, Bin)
@variable(model, x66, Bin)
@variable(model, x67, Bin)
@variable(model, x68, Bin)
@variable(model, x69, Bin)
@variable(model, x70, Bin)
@variable(model, x71, Bin)
@variable(model, x72, Bin)
@variable(model, x73, Bin)
@variable(model, x74, Bin)
@variable(model, x75, Bin)
@variable(model, x76, Bin)
@variable(model, x77, Bin)
@variable(model, x78, Bin)
@variable(model, x79, Bin)
@variable(model, x80, Bin)
@variable(model, x81, Bin)
@variable(model, x82, Bin)
@variable(model, x83, Bin)
@variable(model, x84, Bin)
@variable(model, x85, Bin)
@variable(model, x86, Bin)
@variable(model, x87, Bin)
@variable(model, x88, Bin)
@variable(model, x89, Bin)
@variable(model, x90, Bin)
@variable(model, x91, Bin)
@variable(model, x92, Bin)
@variable(model, x93, Bin)
@variable(model, x94, Bin)
@variable(model, x95, Bin)
@variable(model, x96, Bin)
@variable(model, x97, Bin)
@variable(model, x98, Bin)
@variable(model, x99, Bin)
@variable(model, x100, Bin)
@variable(model, 0 <= y1 <= 1)
@variable(model, 0 <= y2 <= 1)
@variable(model, 0 <= y3 <= 1)
@variable(model, 0 <= y4 <= 1)
@variable(model, 0 <= y5 <= 1)
@variable(model, c >= 0)

@objective(model, Max, mu[1]*x1 + mu[2]*x2 + mu[3]*x3 + mu[4]*x4 + mu[5]*x5 + mu[6]*x6 + mu[7]*x7 + mu[8]*x8 + mu[9]*x9 + mu[10]*x10 + mu[11]*x11 + mu[12]*x12 + mu[13]*x13 + mu[14]*x14 + mu[15]*x15 + mu[16]*x16 + mu[17]*x17 + mu[18]*x18 + mu[19]*x19 + mu[20]*x20 + mu[21]*x21 + mu[22]*x22 + mu[23]*x23 + mu[24]*x24 + mu[25]*x25 + mu[26]*x26 + mu[27]*x27 + mu[28]*x28 + mu[29]*x29 + mu[30]*x30 + mu[31]*x31 + mu[32]*x32 + mu[33]*x33 + mu[34]*x34 + mu[35]*x35 + mu[36]*x36 + mu[37]*x37 + mu[38]*x38 + mu[39]*x39 + mu[40]*x40 + mu[41]*x41 + mu[42]*x42 + mu[43]*x43 + mu[44]*x44 + mu[45]*x45 + mu[46]*x46 + mu[47]*x47 + mu[48]*x48 + mu[49]*x49 + mu[50]*x50 + 
mu[51]*x51 + mu[52]*x52 + mu[53]*x53 + mu[54]*x54 + mu[55]*x55 + mu[56]*x56 + mu[57]*x57 + mu[58]*x58 + mu[59]*x59 + mu[60]*x60 + mu[61]*x61 + mu[62]*x62 + mu[63]*x63 + mu[64]*x64 + mu[65]*x65 + mu[66]*x66 + mu[67]*x67 + mu[68]*x68 + mu[69]*x69 + mu[70]*x70 + mu[71]*x71 + mu[72]*x72 + mu[73]*x73 + mu[74]*x74 + mu[75]*x75 + mu[76]*x76 + mu[77]*x77 + mu[78]*x78 + mu[79]*x79 + mu[80]*x80 + mu[81]*x81 + mu[82]*x82 + mu[83]*x83 + mu[84]*x84 + mu[85]*x85 + mu[86]*x86 + mu[87]*x87 + mu[88]*x88 + mu[89]*x89 + mu[90]*x90 + mu[91]*x91 + mu[92]*x92 + mu[93]*x93 + mu[94]*x94 + mu[95]*x95 + mu[96]*x96 + mu[97]*x97 + mu[98]*x98 + mu[99]*x99 + mu[100]*x100 + mu[101]*y1 + mu[102]*y2 + mu[103]*y3 + mu[104]*y4 + mu[105]*y5 - Omega * c)
@constraint(model, a[1]*x1 + a[2]*x2 + a[3]*x3 + a[4]*x4 + a[5]*x5 + a[6]*x6 + a[7]*x7 + a[8]*x8 + a[9]*x9 + a[10]*x10 + a[11]*x11 + a[12]*x12 + a[13]*x13 + a[14]*x14 + a[15]*x15 + a[16]*x16 + a[17]*x17 + a[18]*x18 + a[19]*x19 + a[20]*x20 + a[21]*x21 + a[22]*x22 + a[23]*x23 + a[24]*x24 + a[25]*x25 +  + a[26]*x26 + a[27]*x27 + a[28]*x28 + a[29]*x29 + a[30]*x30 + a[31]*x31 + a[32]*x32 + a[33]*x33 + a[34]*x34 + a[35]*x35 + a[36]*x36 + a[37]*x37 + a[38]*x38 + a[39]*x39 + a[40]*x40 + a[41]*x41 + a[42]*x42 + a[43]*x43 + a[44]*x44 + a[45]*x45 + a[46]*x46 + a[47]*x47 + a[48]*x48 + a[49]*x49 + a[50]*x50  + 
a[51]*x51 + a[52]*x52 + a[53]*x53 + a[54]*x54 + a[55]*x55 + a[56]*x56 + a[57]*x57 + a[58]*x58 + a[59]*x59 + a[60]*x60 + a[61]*x61 + a[62]*x62 + a[63]*x63 + a[64]*x64 + a[65]*x65 + a[66]*x66 + a[67]*x67 + a[68]*x68 + a[69]*x69 + a[70]*x70 + a[71]*x71 + a[72]*x72 + a[73]*x73 + a[74]*x74 + a[75]*x75 + a[76]*x76 + a[77]*x77 + a[78]*x78 + a[79]*x79 + a[80]*x80 + a[81]*x81 + a[82]*x82 + a[83]*x83 + a[84]*x84 + a[85]*x85 + a[86]*x86 + a[87]*x87 + a[88]*x88 + a[89]*x89 + a[90]*x90 + a[91]*x91 + a[92]*x92 + a[93]*x93 + a[94]*x94 + a[95]*x95 + a[96]*x96 + a[97]*x97 + a[98]*x98 + a[99]*x99 + a[100]*x100 + a[101]*y1 + a[102]*y2 + a[103]*y3 + a[104]*y4 + a[105]*y5 <= b)
@constraint(model, c^2 >= sigma[1]^2*x1^2 + sigma[2]^2*x2^2 + sigma[3]^2*x3^2 + sigma[4]^2*x4^2 + sigma[5]^2*x5^2 + sigma[6]^2*x6^2 + sigma[7]^2*x7^2 + sigma[8]^2*x8^2 + sigma[9]^2*x9^2 + sigma[10]^2*x10^2 + sigma[11]^2*x11^2 + sigma[12]^2*x12^2 + sigma[13]^2*x13^2 + sigma[14]^2*x14^2 + sigma[15]^2*x15^2 + sigma[16]^2*x16^2 + sigma[17]^2*x17^2 + sigma[18]^2*x18^2 + sigma[19]^2*x19^2 + sigma[20]^2*x20^2 + sigma[21]^2*x21^2 + sigma[22]^2*x22^2 + sigma[23]^2*x23^2 + sigma[24]^2*x24^2 + sigma[25]^2*x25^2 +  + sigma[26]^2*x26^2 + sigma[27]^2*x27^2 + sigma[28]^2*x28^2 + sigma[29]^2*x29^2 + sigma[30]^2*x30^2 + sigma[31]^2*x31^2 + sigma[32]^2*x32^2 + sigma[33]^2*x33^2 + sigma[34]^2*x34^2 + sigma[35]^2*x35^2 + sigma[36]^2*x36^2 + sigma[37]^2*x37^2 + sigma[38]^2*x38^2 + sigma[39]^2*x39^2 + sigma[40]^2*x40^2 + sigma[41]^2*x41^2 + sigma[42]^2*x42^2 + sigma[43]^2*x43^2 + sigma[44]^2*x44^2 + sigma[45]^2*x45^2 + sigma[46]^2*x46^2 + sigma[47]^2*x47^2 + sigma[48]^2*x48^2 + sigma[49]^2*x49^2 + sigma[50]^2*x50^2  + 
sigma[51]^2*x51^2 + sigma[52]^2*x52^2 + sigma[53]^2*x53^2 + sigma[54]^2*x54^2 + sigma[55]^2*x55^2 + sigma[56]^2*x56^2 + sigma[57]^2*x57^2 + sigma[58]^2*x58^2 + sigma[59]^2*x59^2 + sigma[60]^2*x60^2 + sigma[61]^2*x61^2 + sigma[62]^2*x62^2 + sigma[63]^2*x63^2 + sigma[64]^2*x64^2 + sigma[65]^2*x65^2 + sigma[66]^2*x66^2 + sigma[67]^2*x67^2 + sigma[68]^2*x68^2 + sigma[69]^2*x69^2 + sigma[70]^2*x70^2 + sigma[71]^2*x71^2 + sigma[72]^2*x72^2 + sigma[73]^2*x73^2 + sigma[74]^2*x74^2 + sigma[75]^2*x75^2 + sigma[76]^2*x76^2 + sigma[77]^2*x77^2 + sigma[78]^2*x78^2 + sigma[79]^2*x79^2 + sigma[80]^2*x80^2 + sigma[81]^2*x81^2 + sigma[82]^2*x82^2 + sigma[83]^2*x83^2 + sigma[84]^2*x84^2 + sigma[85]^2*x85^2 + sigma[86]^2*x86^2 + sigma[87]^2*x87^2 + sigma[88]^2*x88^2 + sigma[89]^2*x89^2 + sigma[90]^2*x90^2 + sigma[91]^2*x91^2 + sigma[92]^2*x92^2 + sigma[93]^2*x93^2 + sigma[94]^2*x94^2 + sigma[95]^2*x95^2 + sigma[96]^2*x96^2 + sigma[97]^2*x97^2 + sigma[98]^2*x98^2 + sigma[99]^2*x99^2 + sigma[100]^2*x100^2 + sigma[101]^2*y1^2 + sigma[102]^2*y2^2 + sigma[103]^2*y3^2 + sigma[104]^2*y4^2 + sigma[105]^2*y5^2)
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
    for i in 1:100
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

function myCallBackFunction(cb_data::CPLEX.CallbackContext, context_id::Clong)
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

    x_val = zeros(Float64 , 105)
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
    x_val[51] = callback_value(cb_data , x51)
    x_val[52] = callback_value(cb_data , x52)
    x_val[53] = callback_value(cb_data , x53)
    x_val[54] = callback_value(cb_data , x54)
    x_val[55] = callback_value(cb_data , x55)
    x_val[56] = callback_value(cb_data , x56)
    x_val[57] = callback_value(cb_data , x57)
    x_val[58] = callback_value(cb_data , x58)
    x_val[59] = callback_value(cb_data , x59)
    x_val[60] = callback_value(cb_data , x60)
    x_val[61] = callback_value(cb_data , x61)
    x_val[62] = callback_value(cb_data , x62)
    x_val[63] = callback_value(cb_data , x63)
    x_val[64] = callback_value(cb_data , x64)
    x_val[65] = callback_value(cb_data , x65)
    x_val[66] = callback_value(cb_data , x66)
    x_val[67] = callback_value(cb_data , x67)
    x_val[68] = callback_value(cb_data , x68)
    x_val[69] = callback_value(cb_data , x69)
    x_val[70] = callback_value(cb_data , x70)
    x_val[71] = callback_value(cb_data , x71)
    x_val[72] = callback_value(cb_data , x72)
    x_val[73] = callback_value(cb_data , x73)
    x_val[74] = callback_value(cb_data , x74)
    x_val[75] = callback_value(cb_data , x75)
    x_val[76] = callback_value(cb_data , x76)
    x_val[77] = callback_value(cb_data , x77)
    x_val[78] = callback_value(cb_data , x78)
    x_val[79] = callback_value(cb_data , x79)
    x_val[80] = callback_value(cb_data , x80)
    x_val[81] = callback_value(cb_data , x81)
    x_val[82] = callback_value(cb_data , x82)
    x_val[83] = callback_value(cb_data , x83)
    x_val[84] = callback_value(cb_data , x84)
    x_val[85] = callback_value(cb_data , x85)
    x_val[86] = callback_value(cb_data , x86)
    x_val[87] = callback_value(cb_data , x87)
    x_val[88] = callback_value(cb_data , x88)
    x_val[89] = callback_value(cb_data , x89)
    x_val[90] = callback_value(cb_data , x90)
    x_val[91] = callback_value(cb_data , x91)
    x_val[92] = callback_value(cb_data , x92)
    x_val[93] = callback_value(cb_data , x93)
    x_val[94] = callback_value(cb_data , x94)
    x_val[95] = callback_value(cb_data , x95)
    x_val[96] = callback_value(cb_data , x96)
    x_val[97] = callback_value(cb_data , x97)
    x_val[98] = callback_value(cb_data , x98)
    x_val[99] = callback_value(cb_data , x99)
    x_val[100] = callback_value(cb_data , x100)
    x_val[101] = callback_value(cb_data , y1)
    x_val[102] = callback_value(cb_data , y2)
    x_val[103] = callback_value(cb_data , y3)
    x_val[104] = callback_value(cb_data , y4)
    x_val[105] = callback_value(cb_data , y5)

    #add heuristic to pick T
    T = collect(1:m)
    for i in 1:m
        if x_val[i+100] <= 0.33
            T[i] = 0
        elseif x_val[i+100] >= 0.66
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
    T_nonzero = filter(!iszero,T)
    temp = separate(m, x_val , T)
    y = [y1 , y2 , y3 , y4 , y5]
    if temp !== nothing
        MOI.submit(model, MOI.UserCut(cb_data), @build_constraint(temp[1]*x1 + temp[2]*x2 + temp[3]*x3 + temp[4]*x4 + temp[5]*x5 + temp[6]*x6 + temp[7]*x7 + temp[8]*x8 + temp[9]*x9 + temp[10]*x10 + temp[11]*x11 + temp[12]*x12 + temp[13]*x13 + temp[14]*x14 + temp[15]*x15 + temp[16]*x16 + temp[17]*x17 + temp[18]*x18 + temp[19]*x19 + temp[20]*x20 + temp[21]*x21 + temp[22]*x22 + temp[23]*x23 + temp[24]*x24 + temp[25]*x25 + temp[26]*x26 + temp[27]*x27 + temp[28]*x28 + temp[29]*x29 + temp[30]*x30 + temp[31]*x31 + temp[32]*x32 + temp[33]*x33 + temp[34]*x34 + temp[35]*x35 + temp[36]*x36 + temp[37]*x37 + temp[38]*x38 + temp[39]*x39 + temp[40]*x40 + temp[41]*x41 + temp[42]*x42 + temp[43]*x43 + temp[44]*x44 + temp[45]*x45 + temp[46]*x46 + temp[47]*x47 + temp[48]*x48 + temp[49]*x49 + temp[50]*x50 + temp[51]*x51 + temp[52]*x52 + temp[53]*x53 + temp[54]*x54 + temp[55]*x55 + temp[56]*x56 + temp[57]*x57 + temp[58]*x58 + temp[59]*x59 + temp[60]*x60 + temp[61]*x61 + temp[62]*x62 + temp[63]*x63 + temp[64]*x64 + temp[65]*x65 + temp[66]*x66 + temp[67]*x67 + temp[68]*x68 + temp[69]*x69 + temp[70]*x70 + temp[71]*x71 + temp[72]*x72 + temp[73]*x73 + temp[74]*x74 + temp[75]*x75 + temp[76]*x76 + temp[77]*x77 + temp[78]*x78 + temp[79]*x79 + temp[80]*x80 + temp[81]*x81 + temp[82]*x82 + temp[83]*x83 + temp[84]*x84 + temp[85]*x85 + temp[86]*x86 + temp[87]*x87 + temp[88]*x88 + temp[89]*x89 + temp[90]*x90 + temp[91]*x91 + temp[92]*x92 + temp[93]*x93 + temp[94]*x94 + temp[95]*x95 + temp[96]*x96 + temp[97]*x97 + temp[98]*x98 + temp[99]*x99 + temp[100]*x100 + mu[101]*y1 + mu[102]*y2 + mu[103]*y3 + mu[104]*y4 + mu[105]*y5 + sqrt(sum(sigma[i+100]^2 * y[i]^2 for i in T_nonzero)) <= temp[101]))
    end
end

MOI.set(model, CPLEX.CallbackFunction(), myCallBackFunction)
optimize!(model)
print(objective_value(model))