
##constraint propagation 적용

using ProgressBars
using MathOptInterface
using JuMP, Cbc
const MOI = MathOptInterface
using Cbc
using MathProgBase
using GLPKMathProgInterface
####################################
########  MPS problem read
####################################
m = MathProgBase.LinearQuadraticModel(GLPKSolverMIP())
MathProgBase.loadproblem!(m,"/HDD/Workspace/CLT/mps/data/R100701005_2_cplex.mps")
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
#m = Model(Cbc.Optimizer)
m=Model(optimizer_with_attributes(Cbc.Optimizer,  "maxSolutions" => 1, "logLevel "=>1, "seconds" => 300,"allowableGap "=>70))

index_x = 1: ncol
const_x = 1: nrow

# variable 만들기
@variable(m, x[i in index_x], lower_bound = xlb[i], upper_bound=xub[i])
var=all_variables(m)
xlb_dict = Dict()
xub_dict = Dict()
for i in index_x
    xlb_dict[var[i]] = xlb[i]
    xub_dict[var[i]] = xub[i]
end


# objective function 만들기
@objective(m, Min, sum( c[i]*var[i] for i in index_x) )



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

# type 만들어 주기
# integer 조건 type dictionary 만들어주기
type_dict = Dict() # integer type (binary and integer) dictionary
bin_type_dict = Dict() # binary dictionary
int_type_dict = Dict() # integer dictionary
cont_type_dict = Dict() # cont dictionary
var_index_dict = Dict()
for i in 1:length(var)
    var_index_dict[var[i]] = i
end

for i in tqdm(index_x)
    if t[i] != :Cont
        if (xlb[i]==0 && xub[i]==1) == true
            type_dict[var[i]] = :Bin
        else
            type_dict[var[i]] = :Int
        end
    else
        continue
    end
end


for i in tqdm(index_x)
    if t[i] == :Cont
        cont_type_dict[var[i]] = t[i]                       
    else
        if (xlb[i]==0 && xub[i]==1) == true
            bin_type_dict[var[i]] = :Bin
        else
            int_type_dict[var[i]] = :Int
        end
    end
end

for i in tqdm(var)
    JuMP.set_lower_bound(i,xlb_dict[i])
    JuMP.set_upper_bound(i,xub_dict[i])    
end


for i in tqdm(keys(type_dict))
    if type_dict[i]==:Bin
        #JuMP.set_binary(i)
        JuMP.set_lower_bound(i,0)
        JuMP.set_upper_bound(i,1)    
    elseif type_dict[i] == :Int
        #JuMP.set_integer(i)
    else
        continue
    end
end

#########################################
## LP relaxation 타이트 시키기!
#########################################
#FP
# solution_k = LP솔루션 이후 최초 해
# solution_k_1 = rounding 한 후 해
# solution_k_2 = 몇 개를 rounding 반대로 해 준 해

# step 1 LP relaxation
JuMP.optimize!(m)
termination_status(m)==MOI.OPTIMAL


solution_k = Dict() # LP솔루션 저장
for i in index_x
    solution_k[var[i]] = JuMP.value(var[i])
end

# infeasible 한 index dictionary
dict_infeasible_index = Dict()
for i in keys(bin_type_dict)
    if solution_k[i] ==  float(0)
        continue
    elseif solution_k[i] ==  float(1)
        continue
    else 
        dict_infeasible_index[i] = :Bin
    end
end
for i in keys(int_type_dict)
    if solution_k[i]!=  round(solution_k[i]) 
        dict_infeasible_index[i] = :Int
    end
end
n_infeasible_var= length(keys(dict_infeasible_index)) #n_infeasible_var: infeasible한  개수


# LP로 구해진 integer set
dict_LP_int = Dict()
solution_k_1 = Dict()
for i in keys(type_dict)
    if solution_k[i]!=round(solution_k[i])
        solution_k_1[i] = round(solution_k[i])
    else dict_LP_int[i] = solution_k[i]
    end
end
# 다시 한 번 확인해보기

# upper bound는 느슨하게 업데이트
for i in keys(dict_LP_int) #8437개
    if lower_bound(i)<dict_LP_int[i]
        JuMP.set_lower_bound(i,dict_LP_int[i])
    end
    if upper_bound(i)>dict_LP_int[i]+0.01
        JuMP.set_upper_bound(i,dict_LP_int[i]+0.01)
    end
end


JuMP.optimize!(m)
termination_status(m)==MOI.OPTIMAL



##############################################
#feasible check
##############################################
infeasible_index = Dict()
for i in keys(bin_type_dict)
    if solution_k[i] ==  float(0)
        continue
    elseif solution_k[i] ==  float(1)
        continue
    else 
        infeasible_index[i] = :Bin
    end
end

for i in keys(int_type_dict)
    if solution_k[i]!=  round(solution_k[i]) 
        infeasible_index[i] = :Int
    end
end

################################################
#  constraint propagation
################################################
#=
각 variable 최대 최소 구하기
1. 해당 variable이 속해 있는 모든 제약조건 찾음
2. 제약조건으로 constraint programming 하기
=#
# 특정 VARIABLE속한 모든 제약조건 찾기
#var_s=A[:,var_index_dict[collect(keys(infeasible_index))[i]],:] # 특정 variable이 속한 모든 constraint

const_list =  Array{Int}(undef,0)
for i in tqdm(keys(infeasible_index))
    #print(i)
    var_s=A[:,var_index_dict[i],:] # 특정 variable이 속한 모든 constraint
    var_s_n=1:var_s.colptr[2]-1 # 속해 있는 constraint 개수
    for j in var_s.rowval[var_s_n]
        #println(j)
        push!(const_list, j)
    end
end

infeasible_var_const=unique(const_list)
sort!(infeasible_var_const)

"""
for i in 1:length(infeasible_var_const)
    print(infeasible_var_const[i],"|")
end
"""


for constraint_index in tqdm(infeasible_var_const)
    #print("const:",constraint_index,"|")
    const_s=A[constraint_index,:,:]
    const_s_n=1:const_s.colptr[2]-1 # 속해 있는 variable 개수
    const_s_coef=const_s.nzval[const_s_n,]#특정 constraint에 해당하는 coefficient
    const_s_var_index=const_s.rowval[const_s_n]#특정 constraint에 해당하는 variable index
    L_min=  Array{Float64}(undef,0)
    for j in const_s_n
        if const_s_coef[j] >0
            #print("variable iteration:",j,"|")
            push!(L_min,const_s_coef[j]*xlb[const_s_var_index[j]])
        else
            push!(L_min,const_s_coef[j]*xub[const_s_var_index[j]])
        end
    end
        L_min=sum(L_min)

    for j in const_s_n
        if const_s_coef[j] >0
            u_new = xlb[const_s_var_index[j]]+(u[constraint_index]-L_min)/const_s_coef[j]
            if upper_bound(var[const_s_var_index[j]]) > u_new
                println(u_new)
                set_upper_bound(var[const_s_var_index[j]],u_new)
            end
        else
            u_new = xub[const_s_var_index[j]] + (u[constraint_index]-L_min)/const_s_coef[j]
            if lower_bound(var[const_s_var_index[j]]) < u_new
                println(u_new)
                set_lower_bound(var[const_s_var_index[j]],u_new)
    
            end

        end
    end
end

optimize!(m)
termination_status(m)==MOI.OPTIMAL

##################################################################################
# FP
##################################################################################
# solution_k = LP솔루션 이후 최초 해
# solution_k_1 = rounding 한 후 해
# solution_k_2 = 몇 개를 rounding 반대로 해 준 해

# step 1 LP relaxation
JuMP.optimize!(m)
termination_status(m)==MOI.OPTIMAL


solution_k = Dict() # LP솔루션 저장
for i in index_x
    solution_k[var[i]] = JuMP.value(var[i])
end

# infeasible 한 index dictionary
dict_infeasible_index = Dict()
for i in keys(bin_type_dict)
    if solution_k[i] ==  float(0)
        continue
    elseif solution_k[i] ==  float(1)
        continue
    else 
        dict_infeasible_index[i] = :Bin
    end
end
for i in keys(int_type_dict)
    if solution_k[i]!=  round(solution_k[i]) 
        dict_infeasible_index[i] = :Int
    end
end
n_infeasible_var= length(keys(dict_infeasible_index)) #n_infeasible_var: infeasible한  개수


# LP로 구해진 integer set
dict_LP_int = Dict()
solution_k_1 = Dict()
for i in keys(type_dict)
    if solution_k[i]!=round(solution_k[i])
        solution_k_1[i] = round(solution_k[i])
    else dict_LP_int[i] = solution_k[i]
    end
end
# 다시 한 번 확인해보기

# upper bound는 느슨하게 업데이트
for i in keys(dict_LP_int) #8437개
    if lower_bound(i)<dict_LP_int[i]
        JuMP.set_lower_bound(i,dict_LP_int[i])
    end
    #if upper_bound(i)>dict_LP_int[i]+0.01
    #    JuMP.set_upper_bound(i,dict_LP_int[i]+0.01)
    #end

end

JuMP.optimize!(m)
termination_status(m)==MOI.OPTIMAL



# step 2/3 rounding
# only for integer
solution_k_1 = Dict()
for i in keys(bin_type_dict)
    if solution_k[i]!=round(solution_k[i])
        solution_k_1[i] = round(solution_k[i])
    end
end
for i in keys(int_type_dict)
    if solution_k[i]!=round(solution_k[i])
        solution_k_1[i] = round(solution_k[i])
    end
end


# step 4 projection
score = Array{Float64}(undef,0)
#score_list = Dict()
score_list = Any[]

for i in keys(solution_k_1)
    if solution_k_1[i] == upper_bound(i)
        val = upper_bound(i)-solution_k[i]
        push!(score, abs(val) )
        #score_list[i] =val
        push!(score_list, (i, abs(val)))
    elseif solution_k_1[i] == lower_bound(i)
        val = solution_k[i]-lower_bound(i)
        push!(score, abs(val))
        #score_list[i] = val
        push!(score_list, (i, abs(val)))
    else
        val = abs(solution_k[i]-solution_k_1[i])
        push!(score,val)
        #score_list[i] = val
        push!(score_list, (i, abs(val)))
    end
end
top_score_list=sort!(score_list, by = x -> abs(x[2]), rev = false)
dist2 = sum(score)


T = 10

solution_k_2 = Dict()
#copy(solution_k)
for i in 1:T
    x_idx=top_score_list[i][1]
    if (solution_k[x_idx]-solution_k_1[x_idx])<0
        solution_k_2[x_idx] = solution_k_1[x_idx]-1
    else
        solution_k_2[x_idx] = solution_k_1[x_idx]+1
    end
end

# step 8 update
for i in keys(solution_k_2)
    if lower_bound(i)<=solution_k_2[i]<=upper_bound(i)
        JuMP.set_lower_bound(i,solution_k_2[i])
    end
    
end

optimize!(m) # 새로운 LP
termination_status(m)==MOI.OPTIMAL


#####################################################################################################
