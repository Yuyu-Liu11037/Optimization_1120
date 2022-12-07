using JuMP, CPLEX, Distributions, LinearAlgebra

n = 25
mu=[45.69250547342361 75.42385196287698 99.29292649272264 90.12534613807001 26.17284108881387 37.376518954546654 6.679093834213035 20.789717613594117 33.36257056179638 57.15630806045373 97.2117285060904 47.29468363862369 22.392413570844138 20.75640243451883 11.189092604419837 12.261574245224649 46.372574667266065 95.66966310532837 49.99524583912126 48.85535458859736 79.66098491982329 3.762350182946328 60.14339056474354 22.742629392575488 28.690046800538948]
a=[6.947049091686441 34.10597780627774 44.869004669879054 29.860941023080446 60.757472049822326 27.37529499615642 80.9566724381156 8.659148232605684 89.96085644428622 27.190903313825043 50.335534246019286 69.52035018393681 39.396935970853875 40.23794523495874 61.34605733580023 72.85292597978584 53.727823947058695 23.06738181710487 63.485099641386135 42.77497173301993 9.466426876441735 79.09046294067822 93.03621463054364 5.765907772190015 44.48302939883511]
b=579.635193887174
sigma=[30.72637049000794, 34.51003151310712, 5.2563796385928985, 6.019544279350812, 24.22232555767528, 5.46904178227805, 4.623200025597562, 3.438798638458221, 32.69160599008206, 30.938286091797252, 54.97218631365624, 43.62260159645204, 8.12032690914649, 19.003759872706425, 0.31716425369448525, 11.882028217012724, 7.8043550908184365, 15.914865139050015, 45.720867773599814, 20.805877874194547, 58.525512919672494, 1.8791724928999844, 9.225389858848864, 8.606539222558844, 20.45761004837807]
Omega=9.9498743710662

model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_MIP_Strategy_HeuristicEffort" => 0 , "CPXPARAM_Preprocessing_Presolve" => 0))
MOI.set(model, MOI.NumberOfThreads(), 1)

@variable(model, x[1:n], Bin)
@variable(model, c >= 0)
@objective(model, Max, sum(mu[i] * x[i] for i in (1:n)') - Omega * c)
@constraint(model, dot(a,x) <= b)
@constraint(model, c^2 >= sum(sigma[i]^2 * x[i]^2 for i in (1:n)'))

optimize!(model)

println(JuMP.objective_value(model))