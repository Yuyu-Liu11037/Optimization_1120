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
@variable(model, x[1:n] , Bin)
@variable(model, 0<= y[1:m] <= 1)
@variable(model, c >= 0)
@objective(model, Max, sum(mu[i] * x[i] for i in (1:n)') + sum(mu[100+i] * y[i] for i in (1:m)') - Omega * c)
@constraint(model, sum(a[i] * x[i] for i in (1:n)') + sum(a[100+i] * y[i] for i in (1:m)') <= b)
@constraint(model, c^2 >= sum(sigma[i]^2 * x[i]^2 for i in (1:n)') + sum(sigma[i+n]^2 * y[i]^2 for i in (1:m)'))

optimize!(model)
print(objective_value(model))