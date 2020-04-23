using JuMP, Cbc, Juniper, Ipopt,Random
model = read_from_file("/HDD/Workspace/CLT/mps/data/R100701005_2_cplex.mps")
nl_solver= optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
mip_solver = optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0 , "seconds" => 100, "feasibilityPump" => "on","maxSolutions" => 1)
optimizer = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>nl_solver, "mip_solver"=>mip_solver)
set_optimizer(model, optimizer)

#optimize!(model)


const J = Juniper
#=
function init_juniper_problem!(jp::JuniperProblem, model::MOI.AbstractOptimizer)
     - model의 타입이 맞지 않음
=#
typeof(backend(model))
J.init_juniper_problem!(J.JuniperProblem(), backend(model)) # 타입이 맞지 않아 error

#JuMP.objective_bound(backend(model))

#=
코드 이해를 위한 간단한 함수 분석
-JuMP 소스코드까지 확인해야 코드 수정이 가능할 듯 합니다.
=#
#패키지의 'objective_function'은 JUMP로 이루어져있다. 
objective_function(model)
#---------------------------------------------------------------
FunType=JuMP.objective_function_type(model)
MOIFunType = moi_function_type(FunType)
func=MOI.get(backend(model), MOI.ObjectiveFunction{MOIFunType}())
jump_function(model, func )
#---------------------------------------------------------------
objective_function_type(model)
