using JuMP, CPLEX, Distributions, LinearAlgebra

n = 100
m = 10
mu=[16.03533670006101 14.139550484878537 37.371738381869555 71.54628061907586 43.787804038906096 5.654976750342833 63.59792145247775 25.198138187764975 61.67065084711767 40.87910582850293 34.94742785284175 7.757810511240926 61.39276478434812 62.700789593168125 5.089639641585741 1.1356791679302924 54.900233120200056 62.528337193400795 65.29824723418604 57.78407331734182 15.393761609753598 26.118493008963927 84.86488038813052 58.25866947834921 69.8854682992065 32.41244530204567 73.12333910872053 52.844617032246234 93.70411920248849 66.31074557134757 72.19920627295429 18.189453757134423 2.854988224807997 70.49521610285015 78.17542266284269 51.77452781848429 20.351010588569395 93.96336374127421 43.463415468493814 55.98957304917259 43.245655296042486 41.498834334803924 55.67052931496548 93.36482772042486 24.00841048006528 91.9085601090717 25.806934192194376 81.98868183775751 94.14392235460848 30.43064274636039 73.9314680996564 21.15417962244558 79.80330556046597 5.720648111319681 65.58024380896745 3.655097585109246 98.12327215621892 12.152163683851247 59.66725025193498 26.849724981705503 47.03690421344814 78.80728070995022 49.62464786277886 96.39072279985949 90.57546453613115 31.87901577550719 73.08337782237925 17.47450433037815 95.1688320438786 22.472209696774293 85.35118875427302 13.91706943011558 94.00745327804077 10.261789823039912 92.48956823446548 18.6849381674681 93.34286125235293 42.381371227487264 98.82952988060855 44.061596735792776 60.96758449288008 87.09408598211581 35.62604163498278 51.07719278662631 64.41367143971519 31.448463785730308 28.024347124190708 61.178097456802746 6.695440598497315 80.0773741437086 83.02173304550344 12.703432350168875 73.43141155183217 78.73778975401459 14.846335590841042 25.382056142257127 18.352657414605844 95.9553491070656 50.00321639116314 22.12133448292314 90.55164962005593 68.02087241944226 44.60976624895121 57.17653204021396 78.32014812502621 70.80483326244857 99.89778613173488 42.04490951266252 12.925736200260873 7.477603623960838]
a=[18.691979144206506 25.863634398442027 9.185317242027812 38.63205160326088 83.86698862241177 56.071112227520715 47.87888946280615 95.44869891792678 17.96880962194085 87.73991179264968 94.15833289845553 67.57750288698078 61.73475084710758 52.74662388113691 88.22727568498627 62.355975865315095 65.38693195254403 83.03585284568285 40.9462568374806 51.64666699470787 54.12631288000578 97.66813848333311 61.733324593583696 18.53058524558039 99.80107047323757 65.78212083965127 3.5124307018194356 43.01370439170758 1.4530842075272665 61.778901532369304 73.40061470647657 60.936534961759946 33.533053341146 21.278466121134333 41.425375618437656 63.7360866107718 31.579861805661658 94.33411232097707 89.7680807920979 0.18936967053089226 88.25516513054038 2.5600299646583724 36.9217331953393 76.70009246433986 37.3071099403219 3.945509893197452 43.81815764048584 0.8427921442840214 37.31959666612544 11.737089023997116 93.31552207764062 14.261651553040355 80.43604544737357 12.29609412575009 12.218349740046108 29.65175203351854 45.195616218370546 36.599362228681386 47.79363318785007 38.01144815428524 68.38256145188923 25.053678047739513 26.77079961701867 57.58497315120744 63.063158484418395 44.22919374519642 33.885452401958716 58.93702946991207 1.5151273922614839 16.756153884367862 76.41419219631051 20.14363517454535 76.74171555677896 9.066322058536214 35.592864511372724 46.12005320789648 36.016015740726026 21.956098235738075 32.97003370711753 9.821146290134708 11.429207430940625 59.44874217975814 44.1618865274997 86.76355336962868 57.682113736553276 37.160024832494564 19.37754068099994 20.85493069230777 54.650252807407874 48.515414753184935 66.3323737171165 4.447610038051952 67.12914820567512 31.072175073788443 18.111301703135783 82.04589983297811 74.63005561234834 51.40077608058515 4.1424836523326185 42.07847087008842 3.0729392632106167 48.398155385484685 89.07837571119967 19.824511876962436 3.0094037353716185 2.738096826788905 35.36836791269206 48.24917927620385 59.029461013521015 89.6681368276242]
b=2478.4111689031533
sigma=[1.3286984932220296, 3.4215718050786057, 19.9754256244379, 11.507227563402942, 39.95912587161312, 2.3074469018806516, 44.03990436988234, 20.635028987823098, 28.722333351175024, 4.696937168193808, 34.06259166298545, 4.854686694636684, 59.931962916978385, 60.81884201667623, 3.6865357571327535, 0.2459098374108936, 37.326706110572445, 32.911889365105736, 3.0419097252596643, 21.95339851777996, 14.251641822370246, 6.736784837720443, 36.26312904789307, 5.612279240021581, 68.49074501265393, 18.056988062905255, 12.112078036003016, 19.387522833919142, 7.4546572571385425, 44.310451468512376, 26.86946709777274, 6.327327490165894, 2.4876634248860348, 66.98653681764856, 47.30347202031208, 16.436903902744934, 1.6432113532057422, 71.25612917901717, 2.855170326721882, 12.818708071423357, 23.430068420442865, 37.56551000562247, 51.14686441621404, 69.12379531506541, 7.2304237663084, 8.529984990130274, 16.849977982090806, 73.32613110092804, 21.612071515505566, 18.44490716680761, 53.11681087987432, 6.9411690423643595, 51.31279814377119, 4.92974096136975, 15.46315784814143, 3.156072252989054, 55.46178987200934, 1.5847153019577946, 10.24284772789808, 0.9907787658036507, 13.475076584025272, 19.57193935476399, 37.796559865116755, 30.128883620464745, 11.238694367059235, 8.376985839150842, 29.289837598294387, 15.283074405989987, 59.481090784575684, 16.538732288708417, 66.41879753915057, 1.997797972987544, 91.82956717875243, 7.9143776542894155, 22.336667679544497, 11.522066254686866, 89.14228099102198, 10.878263268627757, 94.75997348994129, 11.189523956245539, 48.15535830573565, 71.31146309957919, 8.150726347334865, 43.2506332464661, 62.05308044655904, 5.4192183944986105, 8.694040206750282, 59.45543936861825, 2.838338593462094, 75.10863931466918, 10.9139448771451, 5.257587547959743, 70.87113016505243, 24.267146881935634, 5.502803497517955, 11.316972074721422, 9.199560214641751, 70.28010107862997, 8.842735871826568, 17.64051481902161, 78.74700114441629, 36.13704249789792, 30.34566829131177, 9.260742020842642, 25.397269659165016, 42.21701875969531, 41.5097566201296, 12.392023286215583, 6.6974248513991235, 2.4858791172850077]

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