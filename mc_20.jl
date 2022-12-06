using JuMP, CPLEX, Distributions, LinearAlgebra

n = 100
m = 20
#可以用的数据集
mu=[0.38339242343766733 91.41151494606311 2.6100060623990062 66.2528283710079 83.45011505243087 86.1896140635096 62.423894054469876 13.465418255755491 30.811261033023097 90.12709934161293 50.48804330663251 40.31826941196013 73.15044460372143 5.703770682177234 32.1171106248045 60.62832659205634 55.22591550273705 70.54856921613501 9.886718802720818 76.40489723957596 11.889802150452688 40.76815240933514 46.90322337257111 16.43590333014511 33.167206703379115 7.574835036794047 4.878641495780334 28.53787499181094 74.25960904626876 56.74794334136212 78.68838914601439 58.040632713829545 87.53648015612268 54.84278102865251 78.54807062545251 21.486248151286613 5.770970124965036 35.77094215181986 76.22991831711973 73.48900980945083 15.109699344414796 77.8523776642454 2.657907057761688 67.2977113957565 43.20436256276907 17.902399940336633 29.98919552230479 10.616518109135498 42.01426379403447 67.15263269460534 24.244656645091645 12.169876669840784 29.43902509707782 27.387749756574276 39.63193098885844 2.0010015613370125 18.807271896696243 73.15630321185377 42.35021196774594 63.08888326443388 83.75814427510781 37.8028597542391 88.32581330401564 38.527450831897504 48.57205076034871 77.03463800422435 97.84448293593772 74.81980719272536 65.72881832310964 43.05080282531981 52.70290979097877 19.781850603700224 54.259167981347765 63.875489378852365 78.54179692697075 11.195378168415127 18.141759457896576 35.00307772058577 25.682622261430666 17.708142004598983 0.25889151955443834 2.378514298224954 13.565015968639981 22.36601019554475 70.06807869307397 8.249652405440589 46.36379113298166 44.44776874277523 24.01865223792865 17.614582528197797 0.039437798839037086 28.62718078871168 24.190522773070043 90.16479951698992 76.45045252384807 78.90930511891156 54.8320057641688 81.73655120667543 63.59894338265166 73.20696956709067 2.480913076082414 38.30102771327657 5.628533529923841 67.52457218502272 85.23206388285224 73.84787163965728 3.6594759742133998 94.9652236710631 39.96828797432194 98.54618126966443 86.4979846270759 53.878244150208964 70.0725505226638 53.776828334844474 12.90738395001323 76.3208208262586 37.22169818761615 10.24089814714494 34.49679513754168 2.4022363260779445]
a=[42.30323681437248 92.95235235874566 90.74349447541437 55.74552802933338 19.238442203941975 8.35599544601383 58.366507174776316 21.46300190596123 77.51275616425296 18.583705736296718 93.63130057842683 0.39410422205657936 18.962868898985352 88.77425896004742 28.048732373633467 6.976939803958471 49.4954837951703 70.5303078436172 22.72777944424402 62.38222294448552 72.04650994352014 17.69864318269383 97.8673134469352 91.70819126489998 42.21662937749704 42.11262983253575 11.114671014819166 75.32037166993982 77.41527300801319 58.62614694221319 89.31647133158033 14.588429997705143 36.61628042035324 24.56425724914212 47.87346735130451 71.83719610575871 95.15568505153814 42.99007079440417 14.774736546089917 73.00459702008466 46.31984504094305 13.090452779152994 48.23558189184604 58.40586270285989 30.192996213468028 61.93723961534562 27.825071851879677 62.479342487018684 56.53287048934672 53.189613984419196 76.81017663678065 15.894552360502988 2.3449005930240796 21.959857515133052 79.72315215301794 0.0023039805049029916 63.87638148884774 37.34387165185331 66.16175510052868 31.500960543503396 93.1066422233411 56.643453927071285 37.756807585126516 78.55300191087358 25.89604497040329 88.77989429001603 78.08563622990738 61.816984301896305 31.01452888564954 45.291801216152194 36.25374641602923 8.068643099650197 91.95130363038564 97.30348365125973 69.1167573965387 60.23851676534665 46.59527074165677 19.54409046253611 86.56555653614758 21.026221298966806 41.056265559719854 84.75836733955524 78.2174190403974 22.04405782470732 11.555862058675736 1.0420041482975484 33.052995743567315 37.6815697582331 93.73125201473277 44.55349142205559 69.42792896358138 89.4230691495286 18.63942264331473 69.03884370666232 92.50741806022981 89.158480362625 25.424703422569205 10.970797747518402 62.764121674632264 74.95118780235155 28.59170394043661 1.964225592791169 77.68524622856147 40.08839588802334 95.61196627138148 18.210911369294756 96.19681942202459 76.80356733635506 62.310141884734804 9.029629677612295 86.35438609879064 0.7799014423154782 29.78571704787064 28.83710219382546 59.71009248049892 28.761681693944453 88.850856803749 98.41728921422232 93.66182992962074 89.92987680181895]
b=3086.5261825522607
sigma=[0.16276843599396432, 80.71277767786226, 0.023542324019695287, 13.22080455172537, 31.308979786041423, 52.60283614463664, 4.873312479011793, 10.10003089043063, 13.273413327451463, 50.96334384910643, 37.22457337252345, 2.927952132814765, 4.936254574407271, 1.3478808714836812, 7.458266055567498, 5.667327109415715, 17.998317013045334, 67.41955995962236, 4.2932434626908, 0.15555335482431135, 9.24703082423249, 35.31467264063869, 46.2517157917034, 13.031417572115107, 10.412232722279132, 1.8037037824319733, 4.055820371495161, 23.19205722453624, 49.906071126255654, 48.13673833096207, 62.85042601798859, 28.97070580955364, 5.428503625022692, 38.26067006525232, 28.921855411721218, 19.514229526712356, 1.1424205478444123, 33.40288134121236, 52.352974565945836, 22.636749853834253, 12.771723363560623, 58.277316394613884, 0.07911568905175613, 44.23714814330599, 33.53319220634923, 16.747874824627495, 18.0848740407201, 4.031307424446334, 40.50967487501204, 16.34436696990247, 18.688958852504626, 10.919519507819059, 11.541485462291789, 23.1624761752117, 34.259033846557045, 1.1919920022383925, 13.603711650082483, 43.37234585995339, 1.5934183584738737, 44.649486936906364, 67.94958493843518, 16.379801299115034, 38.86477946960343, 28.280636228722553, 2.1790141574564688, 45.67074607220288, 16.907102321162736, 71.12081918652686, 21.26446909131254, 22.2837715520172, 47.92122871189987, 7.822468983705742, 1.3263392910905603, 0.9071303673134026, 1.5618947516536625, 5.796124586734006, 1.5899370747222723, 16.6008378573599, 16.54809410418781, 1.6324490546729407, 0.09781720959479281, 1.1001590675569233, 7.9284402364464, 19.867474036593684, 43.87771348226338, 8.206681768382948, 46.04256308303038, 10.86463798494377, 23.506808256728085, 0.5405423286507993, 0.0008824638101219297, 27.700319186471948, 8.174136893641403, 2.090764316312222, 59.44835827587633, 56.26736989739023, 10.626559969955228, 1.298609174889638, 37.46713956027814, 24.601280855295126, 0.7645594322977879, 27.029326059140246, 1.0738807409272966, 18.818697841677064, 9.607574484351716, 29.720300292723877, 3.628874379535872, 19.77747575721166, 9.07223680590321, 70.53090775591592, 8.852852941667557, 26.212102731042208, 55.58012594121697, 32.04330018263353, 7.023082542601025, 42.303491373698336, 28.5385637295543, 4.483557318921641, 10.55280515743832, 1.6036538151961723]
    #epsilon=0.1
Omega=3.0
    #epsilon=0.05
#Omega=4.358898943540673
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
    end
end

MOI.set(model, CPLEX.CallbackFunction(), myCallBackFunction)
optimize!(model)
print(objective_value(model))