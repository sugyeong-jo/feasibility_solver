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
#using Pkg
#Pkg.clone("https://github.com/odow/LPWriter.jl")
#using LPWriter
include("FeasSolFinder_func_ver.1.0.jl")


#read_from_file("/HDD/Workspace/CLT/mps/processing/CPLEX_file/4_R170725003_veri_cplex.mps")
#bounds_to_set()
#LPWriter.read("/HDD/Workspace/CLT/mps/processing/CPLEX_file/4_R170725003_veri_cplex.mps")
println("================문제 생성=================")
####################################
########  MPS problem read
####################################

m = MathProgBase.LinearQuadraticModel(GLPKSolverMIP())
print("The file name is: ")
#filename = readline()::String
#filepath = string("/HDD/Workspace/CLT/mps/processing/CPLEX_file/",filename)
filepath = string("/HDD/Workspace/CLT/mps/processing/CPLEX_file/4_R170725003_veri_cplex.mps")

MathProgBase.loadproblem!(m,"/HDD/Workspace/CLT/mps/processing/CPLEX_file/4_R170725003_veri_cplex.mps")
#MathProgBase.loadproblem!(m,"/HDD/Workspace/CLT/FP/data/rococoB10-011000.mps")

c = MathProgBase.getobj(m)
A = MathProgBase.getconstrmatrix(m)
nrow,ncol = size(A)
xlb = MathProgBase.getvarLB(m)
xub = MathProgBase.getvarUB(m)
l = MathProgBase.getconstrLB(m)
u = MathProgBase.getconstrUB(m)
t = MathProgBase.getvartype(m)


####################################
########  MPS problem to model
####################################
m=Model(optimizer_with_attributes(Cbc.Optimizer,  "maxSolutions" => 1, "logLevel "=>1, "seconds" => 600,"allowableGap "=>70))
dest = Cbc.Optimizer()
index_x = 1: ncol
const_x = 1: nrow


# variable 만들기
@variable(m, x[i in index_x], lower_bound = xlb[i], upper_bound=xub[i])
var=all_variables(m)

dict_xlb = Dict()
dict_xub = Dict()
for i in index_x
    dict_xlb[var[i]] = lower_bound(var[i])
    dict_xub[var[i]] = upper_bound(var[i])
end


# type 만들어 주기
# integer 조건 type dictionary 만들어주기
dict_IB_type= Dict() # integer type (binary and integer) dictionary
dict_bin_type = Dict() # binary dictionary
dict_int_type = Dict() # integer dictionary
dict_cont_type = Dict() # cont dictionary
dict_var_index = Dict()
for i in 1:length(var)
    dict_var_index[var[i]] = i
end

for i in tqdm(index_x)
    if t[i] != :Cont
        if (xlb[i]==0 && xub[i]==1) == true
            dict_IB_type[var[i]] = :Bin
        else
            dict_IB_type[var[i]] = :Int
        end
    else
        continue
    end
end


for i in tqdm(index_x)
    if t[i] == :Cont
        dict_cont_type[var[i]] = t[i]                       
    else
        if (xlb[i]==0 && xub[i]==1) == true
            dict_bin_type[var[i]] = :Bin
        else
            dict_int_type[var[i]] = :Int
        end
    end
end


"""
# objective function 만들기
@objective(m, Min, sum( c[i]*var[i] for i in index_x) )
"""


###################################
# JuMP로 MPS 바로 적용해보기!
# constraint 핸들링이 힘듦....



####################################
#test
#filename = "R100701005_2_cplex.mps"

model = read_from_file("/HDD/Workspace/CLT/mps/processing/CPLEX_file/4_R170725003_veri_cplex.mps",  format = MOI.FORMAT_MOF)
all_variables(model)
MOI.ConstraintName()
MOI.get(model,MOI.ListOfConstraintIndices{MOI.SingleVariable, MOI.Integer}())
optimizer = Cbc.Optimizer()
set_optimizer(model, optimizer)
MOI.read!("/HDD/Workspace/CLT/mps/processing/CPLEX_file/4_R170725003_veri_cplex.mps")
cbc_dest = Cbc.Optimizer()
dest = Cbc.Optimizer
src =  model
typeof(m)
println(model)
using Clp
Clp_loadProblem(model)
model = MathProgBase.LinearQuadraticModel(ClpSolver())
model = MathProgBase.LinearQuadraticModel(GLPKSolverMIP())
MathProgBase.loadproblem!(model,"/HDD/Workspace/CLT/mps/processing/CPLEX_file/testMPS_small.mps")
MathProgBase.loadproblem!(model,"/HDD/Workspace/CLT/mps/processing/CPLEX_file/4_R170725003_veri_cplex.mps")
MOI.Utilities.IndexMap()
c = MathProgBase.getobj(model)
A = MathProgBase.getconstrmatrix(model)
nrow,ncol = size(A)
xlb = MathProgBase.getvarLB(model)
xub = MathProgBase.getvarUB(model)
l = MathProgBase.getconstrLB(model)
u = MathProgBase.getconstrUB(model)
t = MathProgBase.getvartype(model)
getObjectiveValue(model)
copy(model)
column_value(mapping, index::MOI.VariableIndex)

#####################################
m = read_from_file("/HDD/Workspace/CLT/mps/processing/CPLEX_file/4_R170725003_veri_cplex.mps")
optimizer=optimizer_with_attributes(Cbc.Optimizer,  "maxSolutions" => 1, "logLevel "=>1, "seconds" => 300,"allowableGap "=>70)
set_optimizer(m, optimizer)
var_org = all_variables(m)

model=MathOptInterface.FileFormats.Model(filename = "/HDD/Workspace/CLT/mps/processing/CPLEX_file/4_R170725003_veri_cplex.mps")
MathOptInterface.Utilities.get_bounds(m, {lb},var_org)

dict_var_index_to_org = Dict()
dict_var_org_to_index = Dict()
for x in tqdm(1:ncol)
    dict_var_org_to_index[var_org[x]] = x
end
for x in tqdm(1:ncol)
    dict_var_index_to_org[x] = var_org[x]
end


dict_var_fix = Dict()
dict_bin_type_org = Dict()
dict_int_type_org = Dict()
dict_IB_type_org = Dict()
for x in tqdm(keys(var))
    if is_fixed(x) == true
        dict_var_fix[x] = fix_value(x)
    end
    if is_binary(x) == true
        dict_bin_type_org[x] = :Bin
        dict_IB_type_org[x] = :Bin

    end
    if is_integer(x) == true
        dict_int_type_org[x] = :Int
        dict_IB_type_org[x] = :Int
    end

end

typeof(var_org[1])
keys(dict_var_fix)
dict_var_fix[var_org[3]]
less_than_constraints = all_constraints(m, GenericAffExpr{Float64,VariableRef}, MOI.LessThan{Float64})
less_than_constraints = all_constraints(m, GenericAffExpr{Float64,VariableRef}, MOI.LessThan{Float64})
less_than_constraints = all_constraints(m, GenericAffExpr{Float64,VariableRef}, MOI.LessThan{Float64})
less_than_constraints[1]

function JuMPContainer(::Type{T}, args...) where T
    return JuMPArray(Array{T}(undef, length.(args)...), args...)
end
function JuMPContainer(::Type{T}, args::UnitRange...) where T
    if all(first.(args) .== 1)
        return Array{T}(undef, length.(args)...)
    else
        return JuMPArray(Array{T}(undef, length.(args)...), args...)
    end
end
JuMPContainer(less_than_constraints[1],1)
list_of_constraint_types(m)[1][1]

typeof(less_than_constraints[1]{MathOpt Interface.ScalarAffineFunction{Float64}})

typeof(constraint_object(less_than_constraints[1]).func)
[constraint_object(less_than_constraints[1]).func]
name(less_than_constraints[1])



collect(linear_terms(constraint_object(less_than_constraints[1]).func))[4][2]


length(linear_terms(constraint_object(less_than_constraints[1]).func))
constraint_object(less_than_constraints[1]).func[1]
for (a,b,c) in less_than_constraints[1]
    println(a,b,c)
end

index(less_than_constraints[100])
# unfix the variable
for i in tqdm(index_x)
    #println(var_org[i]) 
    try 
        unfix(var_org[i])
    catch
        continue
    end
end


# bound initial 
for i in tqdm(var_org)
    try
        JuMP.set_lower_bound(i,dict_xlb[dict_var_org_index[i]])
        JuMP.set_upper_bound(i,dict_xub[dict_var_org_index[i]])
        
    catch
        continue
    end    
end

for i in tqdm(keys(dict_bin_type_org))
    JuMP.set_lower_bound(i,0)
    JuMP.set_upper_bound(i,1)    
end


#unset
for i in tqdm(keys(dict_bin_type_org))
    if dict_IB_type_org[i]==:Bin
        JuMP.unset_binary(i)
    elseif dict_IB_type_org[i] == :Int
        JuMP.unset_integer(i)
    else
        continue
    end
end

dict_xlb_org = Dict()
dict_xub_org = Dict()
for i in tqdm(var_org)
    if has_lower_bound(i) == true
        dict_xlb_org[i] = lower_bound(i)
    else 
        dict_xlb_org[i] = -Inf
    end
    if has_upper_bound(i) == true
        dict_xub_org[i] = upper_bound(i)
    else 
        dict_xub_org[i] = +Inf
    end
end



###################################
"""
#constraint 만들기
for i in tqdm(1:nrow) 
    const_s=A[i,:,:]#constraint_specific
    const_s_n=1:const_s.colptr[2]-1
    const_s_coef=const_s.nzval[const_s_n,]#특정 constraint에 해당하는 coefficient
    const_s_var_index=const_s.rowval[const_s_n]#특정 constraint에 해당하는 variable index
    if length(const_s_var_index)==0
        continue
    else
        @constraint(m,l[i] <= sum((x[const_s_var_index[j]]*const_s_coef[j]) for j = 1:length(const_s_coef) ) <= u[i])
    end
end

all_constraints(m, GenericAffExpr{Float64,VariableRef}, MOI.Interval{Float64})
orig_m=copy(m)
#m = copy(orig_m)


# bound initial 
for i in tqdm(var)
    JuMP.set_lower_bound(i,dict_xlb[i])
    JuMP.set_upper_bound(i,dict_xub[i])    
end

for i in tqdm(keys(dict_bin_type))
    JuMP.set_lower_bound(i,0)
    JuMP.set_upper_bound(i,1)    
end
"""

############################################################
# 흠.. ? 이렇게만 해도 되네?
println("================Running!!!!!=================")
@time begin
    Initial_bound()
    try 
        UnSet_Type()
    catch  end

    optimize!(m)
    termination_status(m)==MOI.OPTIMAL
    solution_k, dict_xlb_k,dict_xub_k, dict_LP_int, solution_k_1=LP_solve(m)
    dict_infeasible_index, n_infeasible_var=Infeasible_Check(solution_k)
    Update(dict_LP_int)
    CP(dict_infeasible_index)
    optimize!(m)
    termination_status(m)==MOI.OPTIMAL    
    solution_k_check, dict_xlb_k_check,dict_xub_k_check,dict_LP_int_check, solution_k_1_check=LP_solve(m)
    dict_infeasible_index_check, n_infeasible_var_check=Infeasible_Check(solution_k_check)
    Reupdate(dict_xlb, dict_xub_k_check)
    Set_Type()
    solution_k, dict_xlb_k,dict_xub_k, dict_LP_int, solution_k_1=LP_solve(m,dict_xlb, dict_xub_k_check)
    termination_status(m)!=MOI.OPTIMAL

end
#############################################################################
# result
#############################################################################
sol_check_bound = Array{Any}(undef,0)
sol_check_type = Array{Any}(undef,0)
for x in keys(solution_k)
    push!(sol_check_bound,dict_xlb[x]<=solution_k[x]<=dict_xub[x])
end

for x in keys(dict_bin_type)
    if (solution_k[x] ==0 && solution_k[x] ==1) == true
        push!(sol_check_type, (x,solution_k[x], :Bin))
    else
        continue
    end
end

for x in keys(dict_int_type)
    if isapprox(abs(solution_k[x] - round(solution_k[x])), 0, atol=1e-8) ==false  
        push!(sol_check_type, (x,solution_k[x], :Int))
    else
        continue
    end
end


println("Is it optimal?: ", termination_status(m)!=MOI.OPTIMAL)
println("The the number of unsatisfied variable is:")
println("    - The the number of unsatisfied bound is:", length(findall(x->false, sol_check_bound)))
println("    - The the number of unsatisfied type is:", length(findall(x->false, sol_check_type)))