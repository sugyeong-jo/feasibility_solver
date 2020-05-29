#= 
현재 필요한 문제 형태는 다음 코드에서 A와 같습니다.
A는 다음과 같이 구성됩니다.
- type: sparse matrix
- 형태 : [variable index, constraint index, coefficent]

A 구축을 좀 더 쉽게 핸들링할 수 있게 다음 dictionary가 필요할 것 같습니다.
- Variable name - variable index 
- constraint name - constraint index

=#
#####
println("================패키지 부르기=================")
using Clp
using ProgressBars
using MathOptInterface
using JuMP, Cbc
const MOI = MathOptInterface
using Cbc
using MathProgBase
using GLPKMathProgInterface
using Random
import JuMP
using MathOptInterface
#######
println("================문제 생성=================")
####################################
########  MPS problem read
####################################

m = MathProgBase.LinearQuadraticModel(GLPKSolverMIP())
path = "/HDD/Workspace/CLT/mps/processing/CPLEX_file/4_R170725003_veri_cplex.mps"
MathProgBase.loadproblem!(m,path)

c = MathProgBase.getobj(m)
A = MathProgBase.getconstrmatrix(m) # 이 형태가 필요합니다.
nrow,ncol = size(A)
xlb = MathProgBase.getvarLB(m)
xub = MathProgBase.getvarUB(m)
l = MathProgBase.getconstrLB(m)
u = MathProgBase.getconstrUB(m)
t = MathProgBase.getvartype(m)


#########################################################################
#=
다음'model'로 이어 A를 만들어야 합니다.
=#

model = read_from_file(path)
optimizer = Cbc.Optimizer()
set_optimizer(model, optimizer)

# 여기에서부터 A를 만들어야 하는데.. 정말 어렵네요... 저는 cbc.jl 패키지를 참고하고 있습니다.
# https://github.com/JuliaOpt/Cbc.jl/blob/master/src/MOI_wrapper/MOI_wrapper.jl

MOI.Utilities.IndexMap()
var_org=all_variables(model)
less_than_constraints = all_constraints(m, GenericAffExpr{Float64,VariableRef}, MOI.LessThan{Float64})

#...


