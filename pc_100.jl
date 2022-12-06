using JuMP, CPLEX, Distributions, LinearAlgebra, MathOptInterface

n = 100
#epsilon=0.01的数据，对比明显不要改了
mu=[85.38883522672094 94.2396924080396 39.61270556699572 54.6013934388746 68.28883125610416 34.48870318289588 93.58589954919688 35.55869606242121 68.56587027614968 87.52314053334783 5.9841589159068125 51.718804512879316 35.90394419491454 50.48239441724578 33.85832106134878 65.27805578914754 26.725820685017098 51.27616379024096 54.27339719790966 93.3226787794927 16.740839047095747 32.260941863069924 2.3108563823444617 69.18225492230961 75.92195472360629 64.15306206473818 9.99361025574881 12.973999723208218 64.50345483771767 11.102922855206355 42.445369229000505 64.33210469797459 34.68710659467792 93.17125857171771 48.724697866775315 21.993396538201072 35.93987599656897 11.177192444589734 36.02514642748953 37.23743019134499 7.9330904829360005 99.60762256212112 25.43373150233734 2.436322557102666 30.881821702121858 27.444693578657574 42.628492759334016 76.10689761270586 35.15415141681768 75.41464403720933 36.650759550884885 67.17589236034786 87.1960238234621 21.91846242557427 64.65560102011314 25.234061980463217 97.58486370854935 69.51205016615462 17.647072965554 80.87123430682928 21.117527471810437 29.988290680622576 60.26849276947296 33.07769127374852 66.7638719227768 10.60403884724911 6.049823799143783 79.61370584565894 25.65137424583057 4.564184589922071 20.117769768986193 49.323709820458774 42.28855962425453 99.80672059171832 33.98733888813593 69.63224253919257 81.30824757713226 54.72933675886213 55.43426220661488 2.711178507050016 10.073270780147137 16.29161773093942 8.572883054270742 57.14350082389049 72.88012982017263 52.91889174528273 68.48038796246423 58.431323252070854 27.919306313166004 85.23335529193278 39.076119620681226 74.3229752375296 72.34867200028408 99.26513149157806 37.87618465579773 88.49457360820895 32.64061473173634 82.31020455265721 97.84365071785004 78.64729306151685]
a=[90.3497447607081 59.423576806470145 79.40796750053683 74.53332196908976 3.199473000265718 2.415053457775851 40.619066491302355 48.96296697278807 67.14912793319272 73.9123620185495 21.789435385622134 33.841236139316386 2.203352447562912 24.15943664187893 94.84816234275794 23.810171972989615 17.21851447737601 44.695003687430955 2.1879497040968032 11.66307589959712 30.733947688405628 24.377798369528204 79.84324911510868 69.17944918519767 42.327925633528615 66.89110985633805 56.965005943881785 16.998677459068666 87.29306317561351 97.07628250936801 9.565280569625434 89.19784455651231 47.66111791624584 31.71433555273586 66.71099345410356 50.228695409415955 56.46590756192092 17.019561986516585 10.093980395394286 75.22188675575948 16.447361423172058 16.696222433187767 47.88829154495097 50.88330158301564 40.713284846580954 21.261530610766876 88.34021437521614 41.06311075529092 81.9310683930407 22.819526931846035 9.33234159404207 26.557580117582134 84.21795644239243 0.07089017473038695 53.090779104845446 85.86363523699337 5.681805879480217 24.653062393052373 30.559562672562834 12.267058852625834 94.75883690856286 4.075489424890055 14.561140780175618 43.722343197736635 28.120249695999078 87.62748549867602 83.963090298869 12.324273007391728 29.734915229110015 19.977794830578045 43.391277660925084 48.535850151735524 1.8719581210814629 48.2013350221861 45.331573606819155 96.85595505574143 95.15240186594065 47.86461525270349 44.67353796792042 46.20949172565477 92.34781949253959 10.501755919453604 19.585850054526922 99.46754887197076 60.734326733636166 94.90100563042411 95.71007818907616 80.90342552613711 93.97266718479825 69.11581754354462 39.99150970606548 56.87140730360507 25.371222214201595 38.844921220561865 8.533557406325142 30.476276223225305 80.8625168399543 83.76247350542361 76.51489919826021 94.08132899846977]
b=2395.8988605669247
sigma=[28.675621550456025, 50.31129804193094, 16.460660451524994, 39.49424514230379, 67.80878348551765, 21.61729434946403, 26.69837368490155, 27.446538418550166, 26.88469658375087, 25.46151515110819, 1.1299013856859363, 45.53542185269795, 33.78374208583661, 49.94152590308894, 23.60648003688312, 14.186233406010407, 15.354900005253965, 27.649615178447768, 49.454219451366285, 26.53067705894005, 12.669211220726005, 29.069740557607332, 0.38647750015376453, 30.378054081456344, 58.085000272911934, 13.757422057566865, 2.405261721486074, 4.929580789236528, 35.28405483249155, 8.240385537157044, 33.17017882671827, 50.32724869852045, 26.853015690110542, 38.94226098240429, 16.937487198319005, 9.53180552029533, 13.5523027693026, 8.519876027955391, 24.780411278855194, 4.929520693382756, 6.746764028475192, 14.916118242762723, 1.6219615558748979, 0.07521461574121294, 20.9364453215617, 3.8255105987028752, 13.991003150482266, 8.899014963970467, 34.6307935013495, 68.00440661439686, 19.268283865317848, 29.74930152587676, 83.2706433162123, 18.14639501938429, 8.231274360685072, 5.028469540246054, 59.02155189592037, 63.76714500279677, 2.931311886222003, 59.489361699778875, 7.665098785066166, 6.050849071530158, 22.46806283588477, 24.545942651619306, 64.006044140932, 7.0404717255684695, 3.4750127029652607, 77.56177474468615, 19.582857975157665, 3.3579833420736995, 0.5334636423291759, 1.0053174540052199, 22.89361896670923, 14.084631472079268, 15.2222037066149, 42.88207608950224, 37.772465857499355, 52.20303836133931, 50.870268090561396, 1.800986564363125, 7.264310760882659, 6.594332158041742, 7.108009246372545, 25.314295200942674, 16.196641436600395, 1.5875215791312058, 35.464097769507724, 54.22551923594623, 1.1194952022039837, 14.483182441245566, 21.56266575505877, 27.307866750507554, 26.428818754613186, 98.76395030954846, 11.149710632752518, 29.902324663087544, 7.064256878040111, 32.555584726066115, 58.45118517170635, 51.5104305148142]
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
    temp = Vector{Float64}(undef, 101)
    for i in 1:100
        temp[i] = pi[i]
    end
    temp[101] = t
    if pi' * x_val > t
        return temp
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
@variable(model, c >= 0)
@objective(model,
           Max, 
           mu[1]*x1 + mu[2]*x2 + mu[3]*x3 + mu[4]*x4 + mu[5]*x5 + mu[6]*x6 + mu[7]*x7 + mu[8]*x8 + mu[9]*x9 + mu[10]*x10 + mu[11]*x11 + mu[12]*x12 + mu[13]*x13 + mu[14]*x14 + mu[15]*x15 + mu[16]*x16 + mu[17]*x17 + mu[18]*x18 + mu[19]*x19 + mu[20]*x20 + mu[21]*x21 + mu[22]*x22 + mu[23]*x23 + mu[24]*x24 + mu[25]*x25 + mu[26]*x26 + mu[27]*x27 + mu[28]*x28 + mu[29]*x29 + mu[30]*x30 + mu[31]*x31 + mu[32]*x32 + mu[33]*x33 + mu[34]*x34 + mu[35]*x35 + mu[36]*x36 + mu[37]*x37 + mu[38]*x38 + mu[39]*x39 + mu[40]*x40 + mu[41]*x41 + mu[42]*x42 + mu[43]*x43 + mu[44]*x44 + mu[45]*x45 + mu[46]*x46 + mu[47]*x47 + mu[48]*x48 + mu[49]*x49 + mu[50]*x50 + 
           mu[51]*x51 + mu[52]*x52 + mu[53]*x53 + mu[54]*x54 + mu[55]*x55 + mu[56]*x56 + mu[57]*x57 + mu[58]*x58 + mu[59]*x59 + mu[60]*x60 + mu[61]*x61 + mu[62]*x62 + mu[63]*x63 + mu[64]*x64 + mu[65]*x65 + mu[66]*x66 + mu[67]*x67 + mu[68]*x68 + mu[69]*x69 + mu[70]*x70 + mu[71]*x71 + mu[72]*x72 + mu[73]*x73 + mu[74]*x74 + mu[75]*x75 + mu[76]*x76 + mu[77]*x77 + mu[78]*x78 + mu[79]*x79 + mu[80]*x80 + mu[81]*x81 + mu[82]*x82 + mu[83]*x83 + mu[84]*x84 + mu[85]*x85 + mu[86]*x86 + mu[87]*x87 + mu[88]*x88 + mu[89]*x89 + mu[90]*x90 + mu[91]*x91 + mu[92]*x92 + mu[93]*x93 + mu[94]*x94 + mu[95]*x95 + mu[96]*x96 + mu[97]*x97 + mu[98]*x98 + mu[99]*x99 + mu[100]*x100 - Omega * c)
@constraint(model, a[1]*x1 + a[2]*x2 + a[3]*x3 + a[4]*x4 + a[5]*x5 + a[6]*x6 + a[7]*x7 + a[8]*x8 + a[9]*x9 + a[10]*x10 + a[11]*x11 + a[12]*x12 + a[13]*x13 + a[14]*x14 + a[15]*x15 + a[16]*x16 + a[17]*x17 + a[18]*x18 + a[19]*x19 + a[20]*x20 + a[21]*x21 + a[22]*x22 + a[23]*x23 + a[24]*x24 + a[25]*x25 +  + a[26]*x26 + a[27]*x27 + a[28]*x28 + a[29]*x29 + a[30]*x30 + a[31]*x31 + a[32]*x32 + a[33]*x33 + a[34]*x34 + a[35]*x35 + a[36]*x36 + a[37]*x37 + a[38]*x38 + a[39]*x39 + a[40]*x40 + a[41]*x41 + a[42]*x42 + a[43]*x43 + a[44]*x44 + a[45]*x45 + a[46]*x46 + a[47]*x47 + a[48]*x48 + a[49]*x49 + a[50]*x50  + 
                   a[51]*x51 + a[52]*x52 + a[53]*x53 + a[54]*x54 + a[55]*x55 + a[56]*x56 + a[57]*x57 + a[58]*x58 + a[59]*x59 + a[60]*x60 + a[61]*x61 + a[62]*x62 + a[63]*x63 + a[64]*x64 + a[65]*x65 + a[66]*x66 + a[67]*x67 + a[68]*x68 + a[69]*x69 + a[70]*x70 + a[71]*x71 + a[72]*x72 + a[73]*x73 + a[74]*x74 + a[75]*x75 + a[76]*x76 + a[77]*x77 + a[78]*x78 + a[79]*x79 + a[80]*x80 + a[81]*x81 + a[82]*x82 + a[83]*x83 + a[84]*x84 + a[85]*x85 + a[86]*x86 + a[87]*x87 + a[88]*x88 + a[89]*x89 + a[90]*x90 + a[91]*x91 + a[92]*x92 + a[93]*x93 + a[94]*x94 + a[95]*x95 + a[96]*x96 + a[97]*x97 + a[98]*x98 + a[99]*x99 + a[100]*x100 <= b)
@constraint(model, c^2 >= sigma[1]^2*x1^2 + sigma[2]^2*x2^2 + sigma[3]^2*x3^2 + sigma[4]^2*x4^2 + sigma[5]^2*x5^2 + sigma[6]^2*x6^2 + sigma[7]^2*x7^2 + sigma[8]^2*x8^2 + sigma[9]^2*x9^2 + sigma[10]^2*x10^2 + sigma[11]^2*x11^2 + sigma[12]^2*x12^2 + sigma[13]^2*x13^2 + sigma[14]^2*x14^2 + sigma[15]^2*x15^2 + sigma[16]^2*x16^2 + sigma[17]^2*x17^2 + sigma[18]^2*x18^2 + sigma[19]^2*x19^2 + sigma[20]^2*x20^2 + sigma[21]^2*x21^2 + sigma[22]^2*x22^2 + sigma[23]^2*x23^2 + sigma[24]^2*x24^2 + sigma[25]^2*x25^2 +  + sigma[26]^2*x26^2 + sigma[27]^2*x27^2 + sigma[28]^2*x28^2 + sigma[29]^2*x29^2 + sigma[30]^2*x30^2 + sigma[31]^2*x31^2 + sigma[32]^2*x32^2 + sigma[33]^2*x33^2 + sigma[34]^2*x34^2 + sigma[35]^2*x35^2 + sigma[36]^2*x36^2 + sigma[37]^2*x37^2 + sigma[38]^2*x38^2 + sigma[39]^2*x39^2 + sigma[40]^2*x40^2 + sigma[41]^2*x41^2 + sigma[42]^2*x42^2 + sigma[43]^2*x43^2 + sigma[44]^2*x44^2 + sigma[45]^2*x45^2 + sigma[46]^2*x46^2 + sigma[47]^2*x47^2 + sigma[48]^2*x48^2 + sigma[49]^2*x49^2 + sigma[50]^2*x50^2  + 
                          sigma[51]^2*x51^2 + sigma[52]^2*x52^2 + sigma[53]^2*x53^2 + sigma[54]^2*x54^2 + sigma[55]^2*x55^2 + sigma[56]^2*x56^2 + sigma[57]^2*x57^2 + sigma[58]^2*x58^2 + sigma[59]^2*x59^2 + sigma[60]^2*x60^2 + sigma[61]^2*x61^2 + sigma[62]^2*x62^2 + sigma[63]^2*x63^2 + sigma[64]^2*x64^2 + sigma[65]^2*x65^2 + sigma[66]^2*x66^2 + sigma[67]^2*x67^2 + sigma[68]^2*x68^2 + sigma[69]^2*x69^2 + sigma[70]^2*x70^2 + sigma[71]^2*x71^2 + sigma[72]^2*x72^2 + sigma[73]^2*x73^2 + sigma[74]^2*x74^2 + sigma[75]^2*x75^2 + sigma[76]^2*x76^2 + sigma[77]^2*x77^2 + sigma[78]^2*x78^2 + sigma[79]^2*x79^2 + sigma[80]^2*x80^2 + sigma[81]^2*x81^2 + sigma[82]^2*x82^2 + sigma[83]^2*x83^2 + sigma[84]^2*x84^2 + sigma[85]^2*x85^2 + sigma[86]^2*x86^2 + sigma[87]^2*x87^2 + sigma[88]^2*x88^2 + sigma[89]^2*x89^2 + sigma[90]^2*x90^2 + sigma[91]^2*x91^2 + sigma[92]^2*x92^2 + sigma[93]^2*x93^2 + sigma[94]^2*x94^2 + sigma[95]^2*x95^2 + sigma[96]^2*x96^2 + sigma[97]^2*x97^2 + sigma[98]^2*x98^2 + sigma[99]^2*x99^2 + sigma[100]^2*x100^2)

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

    temp = separate(n, x_val)
    if temp !== nothing
        MOI.submit(model, MOI.UserCut(cb_data), @build_constraint(sum(temp[i] * x[i] for i in 1:n) <= temp[101]))
    end
end

MOI.set(model, CPLEX.CallbackFunction(), myCallBackFunction)

optimize!(model)
println(JuMP.objective_value(model))