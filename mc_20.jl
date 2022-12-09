using JuMP, CPLEX, Distributions, LinearAlgebra

nUCUTS = 0
n = 100
m = 20
mu=[94.76383715123572 33.23043588583412 54.69886833778017 1.8102648236075103 17.602943494470747 42.45315984540249 23.113349674450102 12.916735092914887 17.149123104125806 49.395438804069016 40.16291028139308 3.5886604414051004 24.947956307433827 96.38705728106778 37.78085168201988 30.63070469447916 49.98929908240013 15.301570275737475 91.93836126341604 6.026907088580547 67.95194309873422 29.291073816269318 49.51441900854909 25.79215096786661 10.050937432624728 36.14127180796677 44.292653678525554 59.28275169405437 51.15092572942713 90.32564261722386 30.050898329301866 21.760020161009706 25.119856524276674 16.319016961473697 46.80474714658348 39.036321464716764 91.98134633786792 13.326075430306295 84.13341347567341 34.08189911365661 78.9948830410087 49.335693145967106 23.967383417847365 37.048067203771716 1.8478568693371167 25.848547374348275 24.67764163210543 23.537839908535872 10.415111621842254 33.778874754487255 34.817002475555235 32.85099602779192 0.08423756301806362 79.71743325402288 20.808076621964688 41.638048508340354 73.74708941344763 28.46249575947023 99.7375360630426 66.10054719498287 95.80180100696208 93.56841238357038 36.190991921427 66.19491151001235 21.323481165134773 69.06512898956333 90.44088822255122 69.13364035829028 14.197341549220289 18.87376174783163 56.277604036332654 44.99192412978337 27.205542755825473 59.43104533069422 10.870245332969796 52.36328270456953 3.149734913975133 75.87955447957503 76.81626983344115 55.00741638970007 90.98512218715312 51.83459143791589 99.70512923440975 68.57461856290332 9.844507754933684 93.08496004719751 35.00053415758347 36.376579702165614 77.49443824205632 51.04604264824164 20.75247476363079 31.997627298886233 25.330634571128485 47.944851361166606 62.625798319030125 99.90292830227483 12.49113156474122 88.28169798455725 40.271729295326594 59.074070423602365 34.17481276746178 26.43287719626337 84.96759039419904 79.93803723681295 2.9172092168409325 70.06403562926315 95.03808495140946 2.655899168099729 23.649972663832408 15.457793913379925 9.128049041730568 61.129795637616354 92.91188704842016 10.001138342454164 62.531642589561166 54.38265920570422 20.76210060667475 86.80239449997691 67.79750130710144 87.46015476594042]
a=[31.386036632177582 63.655737265939784 66.87389398998167 32.196926402158866 18.173987691634153 85.36144865442336 35.934921569035616 94.93576419353865 92.52776891335209 91.57106323400865 10.54730972247626 71.9889859188689 78.77740404351742 65.63319505703159 68.95700869565854 52.42942720426276 62.332815421845076 58.67382045817789 70.48969436142724 59.36413209242903 44.81422309697602 12.494426001099178 67.61258669457018 80.07761067787081 59.557452041930226 21.8729813071164 81.54233502653967 92.00092550648799 49.531205757979095 89.33133956017768 6.133659207015207 76.02427279187576 90.7124286989904 97.11806426038837 83.33893860261405 52.97127044814896 48.786822286859234 85.1953214978682 35.27716228661228 41.273688523901576 19.753498695348526 28.79722020644604 13.859560910144086 56.151009310516166 27.321842716109114 15.00989347104702 57.01122761048059 18.201575444874273 64.28385560047288 88.04636964478519 35.82173397776005 40.18990281768321 97.49135619260645 4.544679393411655 59.77713489285333 97.81922465208447 98.42061759399895 53.66802620955484 71.58722962637331 43.3282755700884 50.04105582040669 52.08404463914233 5.623388582855338 42.25259729157985 90.27774351519957 4.737951034619581 56.57443350792201 11.951310233041589 79.06906795293851 56.138523444751975 15.857992388018161 87.03605544905288 57.538936402463214 75.50956328599175 12.566545927978723 57.93588160398587 83.53179269184994 86.77654546902642 16.92409745013801 59.98317257598185 44.54534208490363 10.7188315766503 62.31531353815703 9.519434706661434 31.767538976860322 55.1421266632529 50.021063356368686 92.92586701274165 25.881482514676602 4.224983673214988 86.18752951711527 0.36279767203369895 93.44203862627195 2.9622589488566 45.74700037595723 17.2583822251679 93.46216013370497 37.2064515694272 37.35768484618477 98.11695367803746 66.72787107487453 16.31032260338693 55.073225663645964 12.519867875857194 69.69731534203927 6.073541346942612 94.33795093416633 68.76854328223659 75.76722336091663 43.10317934252696 92.76941310188715 52.28772985431157 41.881340354133265 68.4776205746832 64.01255771437532 71.91331833510245 88.21577414259227 82.98149195249131 13.575376250147453 97.34322983502862]
b=3250.9875621190704
sigma=[10.383265368775083, 23.377038847608222, 45.14101545980379, 1.292723491784178, 3.802377760474864, 31.198735958490932, 14.658549272552142, 0.757730334724595, 7.934684604536421, 38.52307077858826, 30.229035170555033, 1.8819547889262582, 14.72166784032876, 23.760369959316645, 27.430100300268048, 4.297765911311384, 45.829675835915815, 1.14941524182568, 32.91571122624648, 5.442963307495941, 50.846494414459066, 4.619263566054056, 29.81216902520306, 25.618896682865522, 6.745812037073139, 19.07205983899779, 1.0699002397162047, 27.280076726708337, 47.06272865037031, 56.19697314553977, 20.053403952574822, 10.35426786877869, 4.697602710029764, 1.0652562986119833, 18.20223275631279, 21.503559886361995, 7.255392659063088, 3.170496694408836, 20.75162645223102, 24.346523223144853, 70.59632948643457, 17.98352467058198, 20.977518822151925, 7.751758144475114, 1.23490165798387, 1.7445432148033369, 21.436052413265624, 0.5184398513814551, 2.0958358892281566, 30.08703666289813, 24.721538404997737, 25.897635571798425, 0.08000001662502806, 55.52648930897068, 9.348883704729502, 11.925441261984563, 43.612496177864344, 9.817651042789185, 12.39916848023758, 52.04466918664284, 40.935059764397664, 63.77764099430104, 14.57327554628545, 45.33437888769901, 5.389498718303169, 62.19387908983415, 19.06077127336902, 15.88541445814202, 10.771899154568711, 18.828828462096443, 41.43211769290071, 27.19722911934656, 23.603506180083674, 27.06062071546473, 6.246058154409388, 4.299031053412334, 1.3557258255683988, 20.637783298893066, 76.63709525493215, 16.258355815018717, 61.97029839514889, 40.59869291249107, 91.91415483215879, 20.823553386739807, 8.26129493656327, 63.750527201566136, 5.104001234235243, 29.174784395189178, 1.3417793449236055, 50.974493895103365, 6.807961342859707, 21.599816047898727, 10.12951407952531, 27.64831437983299, 16.37910580648235, 73.20381415109163, 6.24852210602977, 48.78268993788264, 28.794147191067733, 5.144708273090075, 12.985448409147967, 19.752165545092332, 13.009539031785483, 15.915561670340757, 0.6729811896528397, 41.590148741517496, 91.5955880930206, 1.5597335136577113, 2.336459789921828, 5.7697963006131925, 2.5212688865016797, 39.322814076925496, 51.65730924862792, 1.6572034438018388, 42.54913001682749, 8.404572981383584, 9.363210722317454, 2.835334823352845, 17.766094164925843, 4.183501159739565]
Omega=9.9498743710662






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
@variable(model, 0 <= y7 <= 1)
@variable(model, 0 <= y8 <= 1)
@variable(model, 0 <= y9 <= 1)
@variable(model, 0 <= y10 <= 1)
@variable(model, 0 <= y11 <= 1)
@variable(model, 0 <= y12 <= 1)
@variable(model, 0 <= y13 <= 1)
@variable(model, 0 <= y14 <= 1)
@variable(model, 0 <= y15 <= 1)
@variable(model, 0 <= y16 <= 1)
@variable(model, 0 <= y17 <= 1)
@variable(model, 0 <= y18 <= 1)
@variable(model, 0 <= y19 <= 1)
@variable(model, 0 <= y20 <= 1)
@variable(model, c >= 0)

@objective(model, Max, mu[1]*x1 + mu[2]*x2 + mu[3]*x3 + mu[4]*x4 + mu[5]*x5 + mu[6]*x6 + mu[7]*x7 + mu[8]*x8 + mu[9]*x9 + mu[10]*x10 + mu[11]*x11 + mu[12]*x12 + mu[13]*x13 + mu[14]*x14 + mu[15]*x15 + mu[16]*x16 + mu[17]*x17 + mu[18]*x18 + mu[19]*x19 + mu[20]*x20 + mu[21]*x21 + mu[22]*x22 + mu[23]*x23 + mu[24]*x24 + mu[25]*x25 + mu[26]*x26 + mu[27]*x27 + mu[28]*x28 + mu[29]*x29 + mu[30]*x30 + mu[31]*x31 + mu[32]*x32 + mu[33]*x33 + mu[34]*x34 + mu[35]*x35 + mu[36]*x36 + mu[37]*x37 + mu[38]*x38 + mu[39]*x39 + mu[40]*x40 + mu[41]*x41 + mu[42]*x42 + mu[43]*x43 + mu[44]*x44 + mu[45]*x45 + mu[46]*x46 + mu[47]*x47 + mu[48]*x48 + mu[49]*x49 + mu[50]*x50 + 
mu[51]*x51 + mu[52]*x52 + mu[53]*x53 + mu[54]*x54 + mu[55]*x55 + mu[56]*x56 + mu[57]*x57 + mu[58]*x58 + mu[59]*x59 + mu[60]*x60 + mu[61]*x61 + mu[62]*x62 + mu[63]*x63 + mu[64]*x64 + mu[65]*x65 + mu[66]*x66 + mu[67]*x67 + mu[68]*x68 + mu[69]*x69 + mu[70]*x70 + mu[71]*x71 + mu[72]*x72 + mu[73]*x73 + mu[74]*x74 + mu[75]*x75 + mu[76]*x76 + mu[77]*x77 + mu[78]*x78 + mu[79]*x79 + mu[80]*x80 + mu[81]*x81 + mu[82]*x82 + mu[83]*x83 + mu[84]*x84 + mu[85]*x85 + mu[86]*x86 + mu[87]*x87 + mu[88]*x88 + mu[89]*x89 + mu[90]*x90 + mu[91]*x91 + mu[92]*x92 + mu[93]*x93 + mu[94]*x94 + mu[95]*x95 + mu[96]*x96 + mu[97]*x97 + mu[98]*x98 + mu[99]*x99 + mu[100]*x100 + mu[101]*y1 + mu[102]*y2 + mu[103]*y3 + mu[104]*y4 + mu[105]*y5 + mu[106]*y6 + mu[107]*y7 + mu[108]*y8 + mu[109]*y9 + mu[110]*y10 + mu[111]*y11 + mu[112]*y12 + mu[113]*y13 + mu[114]*y14 + mu[115]*y15 + mu[116]*y16 + mu[117]*y17 + mu[118]*y18 + mu[119]*y19 + mu[120]*y20 - Omega * c)

@constraint(model, a[1]*x1 + a[2]*x2 + a[3]*x3 + a[4]*x4 + a[5]*x5 + a[6]*x6 + a[7]*x7 + a[8]*x8 + a[9]*x9 + a[10]*x10 + a[11]*x11 + a[12]*x12 + a[13]*x13 + a[14]*x14 + a[15]*x15 + a[16]*x16 + a[17]*x17 + a[18]*x18 + a[19]*x19 + a[20]*x20 + a[21]*x21 + a[22]*x22 + a[23]*x23 + a[24]*x24 + a[25]*x25 +  + a[26]*x26 + a[27]*x27 + a[28]*x28 + a[29]*x29 + a[30]*x30 + a[31]*x31 + a[32]*x32 + a[33]*x33 + a[34]*x34 + a[35]*x35 + a[36]*x36 + a[37]*x37 + a[38]*x38 + a[39]*x39 + a[40]*x40 + a[41]*x41 + a[42]*x42 + a[43]*x43 + a[44]*x44 + a[45]*x45 + a[46]*x46 + a[47]*x47 + a[48]*x48 + a[49]*x49 + a[50]*x50  + 
a[51]*x51 + a[52]*x52 + a[53]*x53 + a[54]*x54 + a[55]*x55 + a[56]*x56 + a[57]*x57 + a[58]*x58 + a[59]*x59 + a[60]*x60 + a[61]*x61 + a[62]*x62 + a[63]*x63 + a[64]*x64 + a[65]*x65 + a[66]*x66 + a[67]*x67 + a[68]*x68 + a[69]*x69 + a[70]*x70 + a[71]*x71 + a[72]*x72 + a[73]*x73 + a[74]*x74 + a[75]*x75 + a[76]*x76 + a[77]*x77 + a[78]*x78 + a[79]*x79 + a[80]*x80 + a[81]*x81 + a[82]*x82 + a[83]*x83 + a[84]*x84 + a[85]*x85 + a[86]*x86 + a[87]*x87 + a[88]*x88 + a[89]*x89 + a[90]*x90 + a[91]*x91 + a[92]*x92 + a[93]*x93 + a[94]*x94 + a[95]*x95 + a[96]*x96 + a[97]*x97 + a[98]*x98 + a[99]*x99 + a[100]*x100 + a[101]*y1 + a[102]*y2 + a[103]*y3 + a[104]*y4 + a[105]*y5  + a[106]*y6 + a[107]*y7 + a[108]*y8 + a[109]*y9 + a[110]*y10 + a[111]*y11 + a[112]*y12 + a[113]*y13 + a[114]*y14 + a[115]*y15 + a[116]*y16 + a[117]*y17 + a[118]*y18 + a[119]*y19 + a[120]*y20 <= b)

@constraint(model, c^2 >= sigma[1]^2*x1^2 + sigma[2]^2*x2^2 + sigma[3]^2*x3^2 + sigma[4]^2*x4^2 + sigma[5]^2*x5^2 + sigma[6]^2*x6^2 + sigma[7]^2*x7^2 + sigma[8]^2*x8^2 + sigma[9]^2*x9^2 + sigma[10]^2*x10^2 + sigma[11]^2*x11^2 + sigma[12]^2*x12^2 + sigma[13]^2*x13^2 + sigma[14]^2*x14^2 + sigma[15]^2*x15^2 + sigma[16]^2*x16^2 + sigma[17]^2*x17^2 + sigma[18]^2*x18^2 + sigma[19]^2*x19^2 + sigma[20]^2*x20^2 + sigma[21]^2*x21^2 + sigma[22]^2*x22^2 + sigma[23]^2*x23^2 + sigma[24]^2*x24^2 + sigma[25]^2*x25^2 +  + sigma[26]^2*x26^2 + sigma[27]^2*x27^2 + sigma[28]^2*x28^2 + sigma[29]^2*x29^2 + sigma[30]^2*x30^2 + sigma[31]^2*x31^2 + sigma[32]^2*x32^2 + sigma[33]^2*x33^2 + sigma[34]^2*x34^2 + sigma[35]^2*x35^2 + sigma[36]^2*x36^2 + sigma[37]^2*x37^2 + sigma[38]^2*x38^2 + sigma[39]^2*x39^2 + sigma[40]^2*x40^2 + sigma[41]^2*x41^2 + sigma[42]^2*x42^2 + sigma[43]^2*x43^2 + sigma[44]^2*x44^2 + sigma[45]^2*x45^2 + sigma[46]^2*x46^2 + sigma[47]^2*x47^2 + sigma[48]^2*x48^2 + sigma[49]^2*x49^2 + sigma[50]^2*x50^2  + 
sigma[51]^2*x51^2 + sigma[52]^2*x52^2 + sigma[53]^2*x53^2 + sigma[54]^2*x54^2 + sigma[55]^2*x55^2 + sigma[56]^2*x56^2 + sigma[57]^2*x57^2 + sigma[58]^2*x58^2 + sigma[59]^2*x59^2 + sigma[60]^2*x60^2 + sigma[61]^2*x61^2 + sigma[62]^2*x62^2 + sigma[63]^2*x63^2 + sigma[64]^2*x64^2 + sigma[65]^2*x65^2 + sigma[66]^2*x66^2 + sigma[67]^2*x67^2 + sigma[68]^2*x68^2 + sigma[69]^2*x69^2 + sigma[70]^2*x70^2 + sigma[71]^2*x71^2 + sigma[72]^2*x72^2 + sigma[73]^2*x73^2 + sigma[74]^2*x74^2 + sigma[75]^2*x75^2 + sigma[76]^2*x76^2 + sigma[77]^2*x77^2 + sigma[78]^2*x78^2 + sigma[79]^2*x79^2 + sigma[80]^2*x80^2 + sigma[81]^2*x81^2 + sigma[82]^2*x82^2 + sigma[83]^2*x83^2 + sigma[84]^2*x84^2 + sigma[85]^2*x85^2 + sigma[86]^2*x86^2 + sigma[87]^2*x87^2 + sigma[88]^2*x88^2 + sigma[89]^2*x89^2 + sigma[90]^2*x90^2 + sigma[91]^2*x91^2 + sigma[92]^2*x92^2 + sigma[93]^2*x93^2 + sigma[94]^2*x94^2 + sigma[95]^2*x95^2 + sigma[96]^2*x96^2 + sigma[97]^2*x97^2 + sigma[98]^2*x98^2 + sigma[99]^2*x99^2 + sigma[100]^2*x100^2 + sigma[101]^2*y1^2 + sigma[102]^2*y2^2 + sigma[103]^2*y3^2 + sigma[104]^2*y4^2 + sigma[105]^2*y5^2 + sigma[106]^2*y6^2 + sigma[107]^2*y7^2 + sigma[108]^2*y8^2 + sigma[109]^2*y9^2 + sigma[110]^2*y10^2 + sigma[111]^2*y11^2 + sigma[112]^2*y12^2 + sigma[113]^2*y13^2 + sigma[114]^2*y14^2 + sigma[115]^2*y15^2 + sigma[116]^2*y16^2 + sigma[117]^2*y17^2 + sigma[118]^2*y18^2 + sigma[119]^2*y19^2 + sigma[120]^2*y20^2)

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
    if pi' * x_val + sum(mu[i+100] * xy_val[i+100] for i in 1:m) + norm(sigma[i+100] .* xy_val[i+100] for i in T_nonzero) > t
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

    x_val = zeros(Float64 , 120)
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
    x_val[111] = callback_value(cb_data , y11)
    x_val[112] = callback_value(cb_data , y12)
    x_val[113] = callback_value(cb_data , y13)
    x_val[114] = callback_value(cb_data , y14)
    x_val[115] = callback_value(cb_data , y15)
    x_val[116] = callback_value(cb_data , y16)
    x_val[117] = callback_value(cb_data , y17)
    x_val[118] = callback_value(cb_data , y18)
    x_val[119] = callback_value(cb_data , y19)
    x_val[120] = callback_value(cb_data , y20)

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
    y = [y1 , y2 , y3 , y4 , y5 , y6 , y7, y8 , y9 , y10 , y11 , y12 , y13 , y14 , y15 , y16 , y17 , y18 , y19 , y20]
    temp = separate(m, x_val , T)
    if temp !== nothing
        MOI.submit(model, MOI.UserCut(cb_data), @build_constraint(temp[1]*x1 + temp[2]*x2 + temp[3]*x3 + temp[4]*x4 + temp[5]*x5 + temp[6]*x6 + temp[7]*x7 + temp[8]*x8 + temp[9]*x9 + temp[10]*x10 + temp[11]*x11 + temp[12]*x12 + temp[13]*x13 + temp[14]*x14 + temp[15]*x15 + temp[16]*x16 + temp[17]*x17 + temp[18]*x18 + temp[19]*x19 + temp[20]*x20 + temp[21]*x21 + temp[22]*x22 + temp[23]*x23 + temp[24]*x24 + temp[25]*x25 + temp[26]*x26 + temp[27]*x27 + temp[28]*x28 + temp[29]*x29 + temp[30]*x30 + temp[31]*x31 + temp[32]*x32 + temp[33]*x33 + temp[34]*x34 + temp[35]*x35 + temp[36]*x36 + temp[37]*x37 + temp[38]*x38 + temp[39]*x39 + temp[40]*x40 + temp[41]*x41 + temp[42]*x42 + temp[43]*x43 + temp[44]*x44 + temp[45]*x45 + temp[46]*x46 + temp[47]*x47 + temp[48]*x48 + temp[49]*x49 + temp[50]*x50 + temp[51]*x51 + temp[52]*x52 + temp[53]*x53 + temp[54]*x54 + temp[55]*x55 + temp[56]*x56 + temp[57]*x57 + temp[58]*x58 + temp[59]*x59 + temp[60]*x60 + temp[61]*x61 + temp[62]*x62 + temp[63]*x63 + temp[64]*x64 + temp[65]*x65 + temp[66]*x66 + temp[67]*x67 + temp[68]*x68 + temp[69]*x69 + temp[70]*x70 + temp[71]*x71 + temp[72]*x72 + temp[73]*x73 + temp[74]*x74 + temp[75]*x75 + temp[76]*x76 + temp[77]*x77 + temp[78]*x78 + temp[79]*x79 + temp[80]*x80 + temp[81]*x81 + temp[82]*x82 + temp[83]*x83 + temp[84]*x84 + temp[85]*x85 + temp[86]*x86 + temp[87]*x87 + temp[88]*x88 + temp[89]*x89 + temp[90]*x90 + temp[91]*x91 + temp[92]*x92 + temp[93]*x93 + temp[94]*x94 + temp[95]*x95 + temp[96]*x96 + temp[97]*x97 + temp[98]*x98 + temp[99]*x99 + temp[100]*x100 + mu[101]*y1 + mu[102]*y2 + mu[103]*y3 + mu[104]*y4 + mu[105]*y5 + mu[106]*y6 + mu[107]*y7 + mu[108]*y8 + mu[109]*y9 + mu[110]*y10 + mu[111]*y11 + mu[112]*y12 + mu[113]*y13 + mu[114]*y14 + mu[115]*y15 + mu[116]*y16 + mu[117]*y17 + mu[118]*y18 + mu[119]*y19 + mu[120]*y20 + sqrt(sum(sigma[i+100]^2 * y[i]^2 for i in T_nonzero)) <= temp[101]))
        nUCUTS += 1
    end
end

MOI.set(model, CPLEX.CallbackFunction(), myCallBackFunction)
optimize!(model)
println(objective_value(model))
println(nUCUTS)