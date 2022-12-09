using JuMP, CPLEX, Distributions, LinearAlgebra

n = 100
m = 10
#very good performance
mu=[83.09205559925192 20.14619899972744 15.999629728643072 74.34302435905094 19.567161513089314 99.85540596818825 96.6440218441089 96.31961971426314 83.37652823329626 82.02729415029515 8.903874540847035 20.925819191681573 6.304410662313997 74.90562557743715 70.72583948355134 17.453086412302532 65.88388320113746 5.972765068172203 64.68543982185871 44.80283341138047 34.87729338248653 35.02251687726814 93.45812005654672 27.122040173306928 50.498851141542566 13.115604692068672 73.07470544362421 35.32800073086907 88.49987811233585 87.63732472655525 74.33892389844056 97.08296994338119 31.626479725839097 5.895782113484971 94.35440172223036 90.0477596370582 27.838634811907237 32.055720368507544 58.884780223672706 5.572236747236092 95.3231969554463 0.36033938763886386 67.86754059836126 61.48118408612654 65.4891218863798 45.277719681545584 18.370435634958472 39.27489575692486 88.17801796908009 97.51501185547912 1.321361883503125 24.897205622808826 99.78435721066735 33.44746826840613 22.0428309410474 27.474989330943355 36.87401494672176 20.687112436824005 39.81970138694256 29.523368757143476 93.3894888049713 58.60581582258704 50.76650287353488 42.026053422368705 93.015774681532 88.46990441521788 51.02834476617918 57.323904340357444 37.80292162259882 65.12445826732997 86.36030227466966 44.10318585594328 45.99792483549599 93.11170222747779 49.36298079977862 43.88662289703411 11.615299698473613 22.293283862618484 43.245464556881416 66.61949701244122 69.66439931731495 25.971307786051977 70.73724693800216 81.0340293036835 41.36947291847196 19.83023173319366 5.805855573509433 35.852998620039166 48.02455832567887 89.50628950888478 72.8017235325686 42.90512691236535 0.8164367466408362 24.79461383839362 78.34552866064826 87.35390998489856 70.94704286690596 38.707633295878296 77.92996505927609 80.32710656752296 40.381864941817454 91.78211465450971 31.629569430471992 83.71078648587856 49.40343077086842 24.668104403914615 8.5792081621351 82.16829822095202 59.1687719854446 50.961936057247904]
a=[46.985179084429504 11.986693083851085 20.844800081548843 62.71927637203092 50.852666721578096 23.362006208210275 34.69196177212323 92.10115208393836 85.41692325865122 37.0109670321764 38.29124731833461 52.14470166460568 81.38020491129943 4.433325381460474 16.431059271621262 75.24169345306187 65.10892553879209 40.629915945394615 90.89825110017959 75.99444212066467 31.081175651401615 10.60285648134397 81.39078275355139 74.25719222033506 80.93563066793683 3.5431904672471703 10.311811136245918 53.42618263485529 59.228863203781266 42.983764941972225 99.75355203300033 18.386652387617865 52.444963286740546 66.48541409840257 1.017305695659243 58.958426833644474 72.98767228084392 1.249688804740945 50.391150146651384 97.73184169060376 71.6721062977154 60.912433816574286 38.811270830686105 70.79205272937597 63.39024171948316 66.8590229929298 96.5043253605501 18.353414378003187 97.0225938166638 26.36458247943968 22.33567375498302 76.89185154973858 9.262922751933289 17.51533546971992 27.768953207963055 15.161941513670573 56.49604904748625 32.683236436480655 44.87802626150458 4.331572310577303 81.67010428909435 21.631219075779317 51.580018546675156 90.07385787956682 14.289983237150494 32.047113216455934 55.26815353478587 56.60004743676596 25.768150067690854 56.08376113768946 39.55768515007414 59.761066980655706 4.657241951412317 86.36795215482611 49.51550435859182 98.80194437056153 58.21790515231176 90.29035590967844 83.5992438065169 3.0018326138659845 97.87913668153786 21.41101558250762 33.53744708682226 89.99720318295248 9.088217102558904 13.854249084761971 69.3972008354932 12.720726958189355 56.12025424696982 7.722192320864451 51.46986293717538 93.29907605425188 9.648327068940322 58.73024055214094 16.18008918425875 0.660297678440791 10.196328953193712 3.56680711097479 26.538005752448115 24.149834681490812 30.866649915606992 54.45120839574938 77.3747688891801 97.86014421844106 27.11820889529418 49.16861580505739 96.92492630998416 30.03076286186953 14.327228616958909 3.905868320924011]
b=2569.3395763505973
sigma=[68.57162880843815, 6.517491323921962, 10.442341255508296, 29.45092174346329, 7.972373169104754, 67.15976926119781, 59.81211826724379, 20.022005858238693, 63.16670304484237, 5.3667675250573526, 0.8540102694891571, 18.385623681300807, 4.3113097107734815, 47.730271960893795, 15.50463186019617, 17.043763307996482, 6.045956870668331, 5.387433488397819, 48.39080717360607, 12.721073439941199, 18.78667648096671, 27.147247770296765, 80.71753954007247, 18.04914740150733, 11.015080636427607, 6.252948397885532, 8.141714385279345, 18.21147697296982, 73.2793373887968, 61.04449560273146, 40.28377434730351, 30.257091133277118, 22.91922521371884, 1.229987286246495, 85.72282650332987, 38.45458896543819, 10.969680936624066, 2.726511405948786, 58.00085589154239, 0.9700441175908574, 64.47853199973902, 0.340978384580995, 52.337042734201496, 48.6062784081993, 32.821741723961225, 18.352904888680026, 6.707540289325013, 22.310995300478282, 42.42604393037605, 52.62203878387636, 0.3758911725137508, 0.2994475978600505, 11.899468459261307, 11.996947262879575, 9.76145553267894, 8.551629788868023, 17.825201129327322, 14.945342987839341, 37.2151587757661, 11.890976667493375, 3.6036960411931376, 9.461829119062006, 18.47907257618378, 35.757164613548326, 14.133946680071027, 55.44745029384055, 45.14307169042162, 40.03347282317019, 9.408962761024066, 22.286819517701332, 29.07125274992057, 17.747926617497598, 28.66225030471637, 9.917111549684735, 25.53586296718013, 37.043660138412775, 4.078473707185989, 3.1276024731135954, 15.548292212052854, 34.03673360518685, 62.305644959154264, 9.541637294624858, 55.37930072496348, 12.06957476435424, 28.225892588184372, 15.079271593815871, 0.1373867553217474, 23.171285819658397, 47.70859871116184, 52.20372501973269, 32.53965824689922, 12.88187581306314, 0.6814420869355657, 24.516205454883995, 65.54796381676933, 78.8849320777254, 4.726235571351146, 9.795094743313513, 24.769111197196334, 62.063560647818804, 2.0297885193023895, 76.28152750743929, 10.150653770970873, 42.051399255675584, 29.878460712051965, 20.973199739518638, 3.4920200673934354, 25.723546377337232, 29.507554145396256, 3.985343452946898]

Omega=4.358898943540673


model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_MIP_Strategy_HeuristicEffort" => 0 , "CPXPARAM_MIP_Strategy_Search" => 1, "CPX_PARAM_CUTSFACTOR" => 1 , "CPXPARAM_Preprocessing_Linear" => 0 , "CPXPARAM_MIP_Interval" => 1))
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
@variable(model, 0 <= y6 <= 1)
@variable(model, 0 <= y8 <= 1)
@variable(model, 0 <= y9 <= 1)
@variable(model, 0 <= y7 <= 1)
@variable(model, 0 <= y10 <= 1)
@variable(model, c >= 0)

@objective(model, Max, mu[1]*x1 + mu[2]*x2 + mu[3]*x3 + mu[4]*x4 + mu[5]*x5 + mu[6]*x6 + mu[7]*x7 + mu[8]*x8 + mu[9]*x9 + mu[10]*x10 + mu[11]*x11 + mu[12]*x12 + mu[13]*x13 + mu[14]*x14 + mu[15]*x15 + mu[16]*x16 + mu[17]*x17 + mu[18]*x18 + mu[19]*x19 + mu[20]*x20 + mu[21]*x21 + mu[22]*x22 + mu[23]*x23 + mu[24]*x24 + mu[25]*x25 + mu[26]*x26 + mu[27]*x27 + mu[28]*x28 + mu[29]*x29 + mu[30]*x30 + mu[31]*x31 + mu[32]*x32 + mu[33]*x33 + mu[34]*x34 + mu[35]*x35 + mu[36]*x36 + mu[37]*x37 + mu[38]*x38 + mu[39]*x39 + mu[40]*x40 + mu[41]*x41 + mu[42]*x42 + mu[43]*x43 + mu[44]*x44 + mu[45]*x45 + mu[46]*x46 + mu[47]*x47 + mu[48]*x48 + mu[49]*x49 + mu[50]*x50 + 
mu[51]*x51 + mu[52]*x52 + mu[53]*x53 + mu[54]*x54 + mu[55]*x55 + mu[56]*x56 + mu[57]*x57 + mu[58]*x58 + mu[59]*x59 + mu[60]*x60 + mu[61]*x61 + mu[62]*x62 + mu[63]*x63 + mu[64]*x64 + mu[65]*x65 + mu[66]*x66 + mu[67]*x67 + mu[68]*x68 + mu[69]*x69 + mu[70]*x70 + mu[71]*x71 + mu[72]*x72 + mu[73]*x73 + mu[74]*x74 + mu[75]*x75 + mu[76]*x76 + mu[77]*x77 + mu[78]*x78 + mu[79]*x79 + mu[80]*x80 + mu[81]*x81 + mu[82]*x82 + mu[83]*x83 + mu[84]*x84 + mu[85]*x85 + mu[86]*x86 + mu[87]*x87 + mu[88]*x88 + mu[89]*x89 + mu[90]*x90 + mu[91]*x91 + mu[92]*x92 + mu[93]*x93 + mu[94]*x94 + mu[95]*x95 + mu[96]*x96 + mu[97]*x97 + mu[98]*x98 + mu[99]*x99 + mu[100]*x100 + mu[101]*y1 + mu[102]*y2 + mu[103]*y3 + mu[104]*y4 + mu[105]*y5 + mu[106]*y6 + mu[107]*y7 + mu[108]*y8 + mu[109]*y9 + mu[110]*y10 - Omega * c)
@constraint(model, a[1]*x1 + a[2]*x2 + a[3]*x3 + a[4]*x4 + a[5]*x5 + a[6]*x6 + a[7]*x7 + a[8]*x8 + a[9]*x9 + a[10]*x10 + a[11]*x11 + a[12]*x12 + a[13]*x13 + a[14]*x14 + a[15]*x15 + a[16]*x16 + a[17]*x17 + a[18]*x18 + a[19]*x19 + a[20]*x20 + a[21]*x21 + a[22]*x22 + a[23]*x23 + a[24]*x24 + a[25]*x25 +  + a[26]*x26 + a[27]*x27 + a[28]*x28 + a[29]*x29 + a[30]*x30 + a[31]*x31 + a[32]*x32 + a[33]*x33 + a[34]*x34 + a[35]*x35 + a[36]*x36 + a[37]*x37 + a[38]*x38 + a[39]*x39 + a[40]*x40 + a[41]*x41 + a[42]*x42 + a[43]*x43 + a[44]*x44 + a[45]*x45 + a[46]*x46 + a[47]*x47 + a[48]*x48 + a[49]*x49 + a[50]*x50  + 
a[51]*x51 + a[52]*x52 + a[53]*x53 + a[54]*x54 + a[55]*x55 + a[56]*x56 + a[57]*x57 + a[58]*x58 + a[59]*x59 + a[60]*x60 + a[61]*x61 + a[62]*x62 + a[63]*x63 + a[64]*x64 + a[65]*x65 + a[66]*x66 + a[67]*x67 + a[68]*x68 + a[69]*x69 + a[70]*x70 + a[71]*x71 + a[72]*x72 + a[73]*x73 + a[74]*x74 + a[75]*x75 + a[76]*x76 + a[77]*x77 + a[78]*x78 + a[79]*x79 + a[80]*x80 + a[81]*x81 + a[82]*x82 + a[83]*x83 + a[84]*x84 + a[85]*x85 + a[86]*x86 + a[87]*x87 + a[88]*x88 + a[89]*x89 + a[90]*x90 + a[91]*x91 + a[92]*x92 + a[93]*x93 + a[94]*x94 + a[95]*x95 + a[96]*x96 + a[97]*x97 + a[98]*x98 + a[99]*x99 + a[100]*x100 + a[101]*y1 + a[102]*y2 + a[103]*y3 + a[104]*y4 + a[105]*y5  + a[106]*y6 + a[107]*y7 + a[108]*y8 + a[109]*y9 + a[110]*y10 <= b)
@constraint(model, c^2 >= sigma[1]^2*x1^2 + sigma[2]^2*x2^2 + sigma[3]^2*x3^2 + sigma[4]^2*x4^2 + sigma[5]^2*x5^2 + sigma[6]^2*x6^2 + sigma[7]^2*x7^2 + sigma[8]^2*x8^2 + sigma[9]^2*x9^2 + sigma[10]^2*x10^2 + sigma[11]^2*x11^2 + sigma[12]^2*x12^2 + sigma[13]^2*x13^2 + sigma[14]^2*x14^2 + sigma[15]^2*x15^2 + sigma[16]^2*x16^2 + sigma[17]^2*x17^2 + sigma[18]^2*x18^2 + sigma[19]^2*x19^2 + sigma[20]^2*x20^2 + sigma[21]^2*x21^2 + sigma[22]^2*x22^2 + sigma[23]^2*x23^2 + sigma[24]^2*x24^2 + sigma[25]^2*x25^2 +  + sigma[26]^2*x26^2 + sigma[27]^2*x27^2 + sigma[28]^2*x28^2 + sigma[29]^2*x29^2 + sigma[30]^2*x30^2 + sigma[31]^2*x31^2 + sigma[32]^2*x32^2 + sigma[33]^2*x33^2 + sigma[34]^2*x34^2 + sigma[35]^2*x35^2 + sigma[36]^2*x36^2 + sigma[37]^2*x37^2 + sigma[38]^2*x38^2 + sigma[39]^2*x39^2 + sigma[40]^2*x40^2 + sigma[41]^2*x41^2 + sigma[42]^2*x42^2 + sigma[43]^2*x43^2 + sigma[44]^2*x44^2 + sigma[45]^2*x45^2 + sigma[46]^2*x46^2 + sigma[47]^2*x47^2 + sigma[48]^2*x48^2 + sigma[49]^2*x49^2 + sigma[50]^2*x50^2  + 
sigma[51]^2*x51^2 + sigma[52]^2*x52^2 + sigma[53]^2*x53^2 + sigma[54]^2*x54^2 + sigma[55]^2*x55^2 + sigma[56]^2*x56^2 + sigma[57]^2*x57^2 + sigma[58]^2*x58^2 + sigma[59]^2*x59^2 + sigma[60]^2*x60^2 + sigma[61]^2*x61^2 + sigma[62]^2*x62^2 + sigma[63]^2*x63^2 + sigma[64]^2*x64^2 + sigma[65]^2*x65^2 + sigma[66]^2*x66^2 + sigma[67]^2*x67^2 + sigma[68]^2*x68^2 + sigma[69]^2*x69^2 + sigma[70]^2*x70^2 + sigma[71]^2*x71^2 + sigma[72]^2*x72^2 + sigma[73]^2*x73^2 + sigma[74]^2*x74^2 + sigma[75]^2*x75^2 + sigma[76]^2*x76^2 + sigma[77]^2*x77^2 + sigma[78]^2*x78^2 + sigma[79]^2*x79^2 + sigma[80]^2*x80^2 + sigma[81]^2*x81^2 + sigma[82]^2*x82^2 + sigma[83]^2*x83^2 + sigma[84]^2*x84^2 + sigma[85]^2*x85^2 + sigma[86]^2*x86^2 + sigma[87]^2*x87^2 + sigma[88]^2*x88^2 + sigma[89]^2*x89^2 + sigma[90]^2*x90^2 + sigma[91]^2*x91^2 + sigma[92]^2*x92^2 + sigma[93]^2*x93^2 + sigma[94]^2*x94^2 + sigma[95]^2*x95^2 + sigma[96]^2*x96^2 + sigma[97]^2*x97^2 + sigma[98]^2*x98^2 + sigma[99]^2*x99^2 + sigma[100]^2*x100^2 + sigma[101]^2*y1^2 + sigma[102]^2*y2^2 + sigma[103]^2*y3^2 + sigma[104]^2*y4^2 + sigma[105]^2*y5^2 + sigma[106]^2*y6^2 + sigma[107]^2*y7^2 + sigma[108]^2*y8^2 + sigma[109]^2*y9^2 + sigma[110]^2*y10^2)

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
    t = sum(mu[i] * xy_val[i] for i in 1:(100+m)) - Omega * sqrt(sigma[i]^2 * xy_val[i]^2 for i in 1:(100+m))

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
    x_val[106] = callback_value(cb_data , y6)
    x_val[107] = callback_value(cb_data , y7)
    x_val[108] = callback_value(cb_data , y8)
    x_val[109] = callback_value(cb_data , y9)
    x_val[110] = callback_value(cb_data , y10)

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
    y = [y1 , y2 , y3 , y4 , y5 , y6 , y7, y8 , y9 , y10]
    temp = separate(m, x_val , T)
    if temp !== nothing
        MOI.submit(model, MOI.UserCut(cb_data), @build_constraint(temp[1]*x1 + temp[2]*x2 + temp[3]*x3 + temp[4]*x4 + temp[5]*x5 + temp[6]*x6 + temp[7]*x7 + temp[8]*x8 + temp[9]*x9 + temp[10]*x10 + temp[11]*x11 + temp[12]*x12 + temp[13]*x13 + temp[14]*x14 + temp[15]*x15 + temp[16]*x16 + temp[17]*x17 + temp[18]*x18 + temp[19]*x19 + temp[20]*x20 + temp[21]*x21 + temp[22]*x22 + temp[23]*x23 + temp[24]*x24 + temp[25]*x25 + temp[26]*x26 + temp[27]*x27 + temp[28]*x28 + temp[29]*x29 + temp[30]*x30 + temp[31]*x31 + temp[32]*x32 + temp[33]*x33 + temp[34]*x34 + temp[35]*x35 + temp[36]*x36 + temp[37]*x37 + temp[38]*x38 + temp[39]*x39 + temp[40]*x40 + temp[41]*x41 + temp[42]*x42 + temp[43]*x43 + temp[44]*x44 + temp[45]*x45 + temp[46]*x46 + temp[47]*x47 + temp[48]*x48 + temp[49]*x49 + temp[50]*x50 + temp[51]*x51 + temp[52]*x52 + temp[53]*x53 + temp[54]*x54 + temp[55]*x55 + temp[56]*x56 + temp[57]*x57 + temp[58]*x58 + temp[59]*x59 + temp[60]*x60 + temp[61]*x61 + temp[62]*x62 + temp[63]*x63 + temp[64]*x64 + temp[65]*x65 + temp[66]*x66 + temp[67]*x67 + temp[68]*x68 + temp[69]*x69 + temp[70]*x70 + temp[71]*x71 + temp[72]*x72 + temp[73]*x73 + temp[74]*x74 + temp[75]*x75 + temp[76]*x76 + temp[77]*x77 + temp[78]*x78 + temp[79]*x79 + temp[80]*x80 + temp[81]*x81 + temp[82]*x82 + temp[83]*x83 + temp[84]*x84 + temp[85]*x85 + temp[86]*x86 + temp[87]*x87 + temp[88]*x88 + temp[89]*x89 + temp[90]*x90 + temp[91]*x91 + temp[92]*x92 + temp[93]*x93 + temp[94]*x94 + temp[95]*x95 + temp[96]*x96 + temp[97]*x97 + temp[98]*x98 + temp[99]*x99 + temp[100]*x100 + mu[101]*y1 + mu[102]*y2 + mu[103]*y3 + mu[104]*y4 + mu[105]*y5 + mu[106]*y6 + mu[107]*y7 + mu[108]*y8 + mu[109]*y9 + mu[110]*y10 + sqrt(sum(sigma[i+100]^2 * y[i]^2 for i in T_nonzero)) <= temp[101]))
    end
end

MOI.set(model, CPLEX.CallbackFunction(), myCallBackFunction)

optimize!(model)
print(objective_value(model))