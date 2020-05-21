
##constraint propagation 적용

using ProgressBars
using MathOptInterface
using JuMP, Cbc
const MOI = MathOptInterface
using Cbc
using MathProgBase
using GLPKMathProgInterface
using Random
####################################
########  MPS problem read
####################################
m = MathProgBase.LinearQuadraticModel(GLPKSolverMIP())
MathProgBase.loadproblem!(m,"/HDD/Workspace/CLT/mps/processing/CPLEX_file/R100701005_2_cplex.mps")
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
import JuMP
using MathOptInterface

#m = Model(Cbc.Optimizer)
m=Model(optimizer_with_attributes(Cbc.Optimizer,  "maxSolutions" => 1, "logLevel "=>1, "seconds" => 300,"allowableGap "=>70))


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


# bound initial 
for i in tqdm(var)
    JuMP.set_lower_bound(i,dict_xlb[i])
    JuMP.set_upper_bound(i,dict_xub[i])    
end

for i in tqdm(keys(dict_bin_type))
    JuMP.set_lower_bound(i,0)
    JuMP.set_upper_bound(i,1)    
end


function Set_Type()
    for i in tqdm(keys(dict_IB_type))
        if dict_IB_type[i]==:Bin
            JuMP.set_binary(i)
        elseif dict_IB_type[i] == :Int
            JuMP.set_integer(i)
        else
            continue
        end
    end
end

function UnSet_Type()
    for i in tqdm(keys(dict_IB_type))
        if dict_IB_type[i]==:Bin
            JuMP.unset_binary(i)
        elseif dict_IB_type[i] == :Int
            JuMP.unset_integer(i)
        else
            continue
        end
    end
end


##############################################
function Initial_bound()
    global var
    for i in tqdm(var)
        JuMP.set_lower_bound(i,dict_xlb[i])
        JuMP.set_upper_bound(i,dict_xub[i])    
    end
end

# upper bound는 느슨하게 업데이트
function Update(dict_LP_int)
    for i in keys(dict_LP_int) #8437개
        if lower_bound(i)<=dict_LP_int[i]<=upper_bound(i)
            JuMP.set_lower_bound(i,dict_LP_int[i])
        end
    end
end


function LPUpdate(dict_LP_int)
    for i in keys(dict_LP_int) #8437개
        if lower_bound(i)<=dict_LP_int[i]<=upper_bound(i)-0.01
            JuMP.set_lower_bound(i,dict_LP_int[i])
            JuMP.set_upper_bound(i,dict_LP_int[i]+0.01)
        end
    end
end


# bound 재조정
function Reupdate(dict_xlb_k, dict_xub_k)
    Initial_bound()
    for i in keys(dict_xlb_k) 
        set_lower_bound(i,dict_xlb_k[i])
    end
    for i in keys(dict_xub_k) 
        set_upper_bound(i,dict_xub_k[i])
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

function LP_solve(m,dict_xlb_k,dict_xub_k)
    global dict_IB_type
    Reupdate(dict_xlb_k, dict_xub_k)
    optimize!(m)        
    solution_k = Dict() # LP솔루션 저장
    for i in index_x
        solution_k[var[i]] = JuMP.value(var[i])
    end   
    dict_LP_int = Dict()     # LP로 구해진 integer set
    solution_k_1 = Dict()    # LP로 구해진 integer set을 round 한 solution set
    for i in keys(dict_IB_type)
        if solution_k[i]!= round(solution_k[i])
            solution_k_1[i] = round(solution_k[i])
        else dict_LP_int[i] = solution_k[i]
        end
    end
    dict_xlb_k = Dict() # lower bound dictionary
    for i in var
        dict_xlb_k[i]=lower_bound(i)
    end
    dict_xub_k = Dict() # lower bound dictionary
    for i in var
        dict_xub_k[i]=upper_bound(i)
    end
    return solution_k, dict_xlb_k, dict_xub_k, dict_LP_int, solution_k_1
end


function LP_solve(m)
    optimize!(m)        
    solution_k = Dict() # LP솔루션 저장
    for i in index_x
        solution_k[var[i]] = JuMP.value(var[i])
    end   
    dict_LP_int = Dict()     # LP로 구해진 integer set
    solution_k_1 = Dict()    # LP로 구해진 integer set을 round 한 solution set
    for i in keys(dict_IB_type)
        if solution_k[i]!=round(solution_k[i])
            solution_k_1[i] = round(solution_k[i])
        else dict_LP_int[i] = solution_k[i]
        end
    end
    dict_xlb_k = Dict() # lower bound dictionary
    for i in var
        dict_xlb_k[i]=lower_bound(i)
    end
    dict_xub_k = Dict() # lower bound dictionary
    for i in var
        dict_xub_k[i]=upper_bound(i)
    end
    return solution_k, dict_xlb_k, dict_xub_k, dict_LP_int, solution_k_1
end


##############################################
#feasible check
##############################################
# infeasible 한 index dictionary
function Infeasible_Check(solution_k)
    dict_infeasible_index = Dict()
    for i in keys(dict_bin_type)
        if solution_k[i] ==  float(0)
            continue
        elseif solution_k[i] ==  float(1)
            continue
        else 
            dict_infeasible_index[i] = :Bin
        end
    end
    for i in keys(dict_int_type)
        if solution_k[i] !=  round(solution_k[i]) 
            dict_infeasible_index[i] = :Int
        end
    end
    n_infeasible_var= length(keys(dict_infeasible_index)) #n_infeasible_var: infeasible한  개수
    return dict_infeasible_index, n_infeasible_var
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

function CP(dict_infeasible_index)
    const_list =  Array{Int}(undef,0)
    println("The all constraints including infeasible variables")
    for i in tqdm(keys(dict_infeasible_index))
        #print(i)
        var_s=A[:,dict_var_index[i],:] # 특정 variable이 속한 모든 constraint
        var_s_n=1:var_s.colptr[2]-1 # 속해 있는 constraint 개수
        for j in var_s.rowval[var_s_n]
            #println(j)
            push!(const_list, j)
        end
    end

    infeasible_var_const=unique(const_list)
    #sort!(infeasible_var_const)

    println("Constraint Propagarion")
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
            if const_s_coef[j] >= 0
                u_new = xlb[const_s_var_index[j]] + (u[constraint_index]-L_min)/const_s_coef[j]
                if upper_bound(var[const_s_var_index[j]]) > u_new
                    set_upper_bound(var[const_s_var_index[j]],u_new)
                end
            else
                u_new = xub[const_s_var_index[j]] + (u[constraint_index]-L_min)/const_s_coef[j]
                if lower_bound(var[const_s_var_index[j]]) < u_new
                    set_lower_bound(var[const_s_var_index[j]],u_new)
        
                end
            end
        end
     end
end


function Projection(solution_k_1)
    # step 4 projection
    # score_list: feasible LP의 score list!
    # score_list_1: round 후 score
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
    top_score_list=sort!(score_list, by = x -> abs(x[2]), rev = true)
    dist = sum(score)
    return dist, top_score_list
end

function Round_change(dict_x_indx,solution_k,solution_k_1)
    solution_k_2 = Dict()
    for x_idx in dict_x_indx
        if (solution_k[x_idx]-solution_k_1[x_idx])<0
            solution_k_2[x_idx] = solution_k_1[x_idx]-1
        else
            solution_k_2[x_idx] = solution_k_1[x_idx]+1
        end
    end
    return solution_k_2
end



    
###########################################################################################
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

# LP strengthen
while true
    global solution_k, dict_xlb_k, dict_xub_k, dict_LP_int, solution_k_1, solution_k_check,dict_xlb_k_check,dict_xub_k_check, dict_LP_int_check, solution_k_1_check
    global dict_infeasible_index, n_infeasible_var, dict_infeasible_index_check, n_infeasible_var_check
    if n_infeasible_var>n_infeasible_var_check
        solution_k, dict_xlb_k,dict_xub_k,dict_LP_int, solution_k_1 = solution_k_check, dict_xlb_k_check,dict_xub_k_check,dict_LP_int_check, solution_k_1_check
        dict_infeasible_index, n_infeasible_var = dict_infeasible_index_check, n_infeasible_var_check
        Update(dict_LP_int)
        CP(dict_infeasible_index)
        optimize!(m)
        if termination_status(m)!=MOI.OPTIMAL
            Reupdate(dict_xlb_k, dict_xub_k)
            solution_k_check, dict_xlb_k_check,dict_xub_k_check,dict_LP_int_check, solution_k_1_check=LP_solve(m)
            dict_infeasible_index, n_infeasible_var=Infeasible_Check(solution_k_check)
            break
        end    
        solution_k_check, dict_xlb_k_check, dict_xub_k_check,dict_LP_int_check, solution_k_1_check=LP_solve(m)
        dict_infeasible_index_check, n_infeasible_var_check=Infeasible_Check(solution_k_check)
        global solution_k, dict_xlb_k, dict_xub_k, dict_LP_int, solution_k_1, solution_k_check,dict_xlb_k_check,dict_xub_k_check, dict_LP_int_check, solution_k_1_check
        global dict_infeasible_index, n_infeasible_var, dict_infeasible_index_check, n_infeasible_var_check
        continue
    else
        break
    end
end

#####
#LP 최적 해 저장!
solution_k_Best,dict_xlb_k_Best,dict_xub_k_Best, dict_LP_int_Best, solution_k_1_Best = solution_k,dict_xlb_k,dict_xub_k,  dict_LP_int, solution_k_1
dict_infeasible_index_Best, n_infeasible_var_Best=dict_infeasible_index, n_infeasible_var
#

solution_k,dict_xlb_k,dict_xub_k,  dict_LP_int, solution_k_1 = solution_k_Best,dict_xlb_k_Best,dict_xub_k_Best, dict_LP_int_Best, solution_k_1_Best
dict_infeasible_index, n_infeasible_var = dict_infeasible_index_Best, n_infeasible_var_Best
####
# dist Top T 업데이트
N=5
T = 5
for i in 1:N
    try
        global m
        global solution_k, dict_xlb_k, dict_xub_k, dict_LP_int, solution_k_1
        global solution_k_check,dict_xlb_k_check,dict_xub_k_check, dict_LP_int_check, solution_k_1_check
        global solution_k_Best, dict_xlb_k_Best, dict_xub_k_Best, dict_LP_int_Best, solution_k_1_Best
        global dict_infeasible_index, n_infeasible_var
        global dict_infeasible_index_check, n_infeasible_var_check
        global dict_infeasible_index_Best, n_infeasible_var_Best
        global T    

        solution_k, dict_xlb_k,dict_xub_k, dict_LP_int, solution_k_1=LP_solve(m,dict_xlb_k_Best,dict_xub_k_Best)
        dict_infeasible_index, n_infeasible_var=Infeasible_Check(solution_k)
        dist, top_score_list=Projection(solution_k_1)
    
        
        x_indx = Array{Any}(undef,0)
        for i in 1:T
            push!(x_indx,top_score_list[i][1]) 
        end
        dict_round_select = Dict() 
        for i in x_indx
            dict_round_select[i]=solution_k_1[i]
        end
        Update(dict_round_select)
    #    CP(dict_infeasible_index)
        
        solution_k_check, dict_xlb_k_check, dict_xub_k_check,dict_LP_int_check, solution_k_1_check=LP_solve(m)
        dict_infeasible_index_check, n_infeasible_var_check=Infeasible_Check(solution_k_check)
        dist, top_score_list=Projection(solution_k_1_check)
        global T, dist  
        println("===================================================================================================================")
        println("   ")
        println("[T, N]: [" ,T,",", N, "]|n_infeasible_var: ", n_infeasible_var_check, "|dist: ", dist)
        println("   ")
        println("===================================================================================================================")

        #####
        #LP 최적 해 저장!
        if n_infeasible_var_check < n_infeasible_var_Best

            solution_k_Best,dict_xlb_k_Best,dict_xub_k_Best, dict_LP_int_Best, solution_k_1_Best = solution_k_check,dict_xlb_k_check,dict_xub_k_check,  dict_LP_int_check, solution_k_1_check
            dict_infeasible_index_Best, n_infeasible_var_Best=dict_infeasible_index_check, n_infeasible_var_check
            dist_Best = dist
            #

            solution_k,dict_xlb_k,dict_xub_k,  dict_LP_int, solution_k_1 = solution_k_Best,dict_xlb_k_Best,dict_xub_k_Best, dict_LP_int_Best, solution_k_1_Best
            dict_infeasible_index, n_infeasible_var = dict_infeasible_index_Best, n_infeasible_var_Best
        end
        global solution_k, dict_xlb_k, dict_xub_k, dict_LP_int, solution_k_1
        global solution_k_check,dict_xlb_k_check,dict_xub_k_check, dict_LP_int_check, solution_k_1_check
        global solution_k_Best, dict_xlb_k_Best, dict_xub_k_Best, dict_LP_int_Best, solution_k_1_Best
        global dict_infeasible_index, n_infeasible_var
        global dict_infeasible_index_check, n_infeasible_var_check
        global dict_infeasible_index_Best, n_infeasible_var_Best
        global T, dist    
        println("===================================================================================================================")
        println("   ")
        println("[T, N]: [" ,T,",", N, "]|n_infeasible_var_Best: ", n_infeasible_var_Best, "|dist: ", dist_Best)
        println("   ")
        println("===================================================================================================================")
        continue
    catch
            continue
    
    end
    
end

################################################
# Random    
################################################
N=5
T = 2
for i in 1:N
    try

        global m
        global solution_k, dict_xlb_k, dict_xub_k, dict_LP_int, solution_k_1
        global solution_k_check,dict_xlb_k_check,dict_xub_k_check, dict_LP_int_check, solution_k_1_check
        global solution_k_Best, dict_xlb_k_Best, dict_xub_k_Best, dict_LP_int_Best, solution_k_1_Best
        global dict_infeasible_index, n_infeasible_var
        global dict_infeasible_index_check, n_infeasible_var_check
        global dict_infeasible_index_Best, n_infeasible_var_Best
        global T    

        solution_k, dict_xlb_k,dict_xub_k, dict_LP_int, solution_k_1=LP_solve(m,dict_xlb_k_Best,dict_xub)
        dict_infeasible_index, n_infeasible_var=Infeasible_Check(solution_k)
        dist, top_score_list=Projection(solution_k_1)
        #CP(dict_infeasible_index)
        dist_Best = dist
    
        
        x_indx=rand(keys(dict_infeasible_index), T )
        solution_k_2 = Round_change(x_indx,solution_k_1,solution_k_1)
        Update(solution_k_2)       
        solution_k_check, dict_xlb_k_check, dict_xub_k_check,dict_LP_int_check, solution_k_1_check=LP_solve(m)
        dict_infeasible_index_check, n_infeasible_var_check=Infeasible_Check(solution_k_check)
        dist, top_score_list=Projection(solution_k_1_check)
        println("===================================================================================================================")
        println("   ")
        println("[T, N]: [" ,T,",", N, "]|n_infeasible_var: ", n_infeasible_var_check, "|dist: ", dist)
        println("   ")
        println("===================================================================================================================")

        #####
        #LP 최적 해 저장!
        if n_infeasible_var_check < n_infeasible_var_Best

            solution_k_Best,dict_xlb_k_Best,dict_xub_k_Best, dict_LP_int_Best, solution_k_1_Best = solution_k_check,dict_xlb_k_check,dict_xub_k_check,  dict_LP_int_check, solution_k_1_check
            dict_infeasible_index_Best, n_infeasible_var_Best=dict_infeasible_index_check, n_infeasible_var_check
            dist_Best = dist
            #

            solution_k,dict_xlb_k,dict_xub_k,  dict_LP_int, solution_k_1 = solution_k_Best,dict_xlb_k_Best,dict_xub_k_Best, dict_LP_int_Best, solution_k_1_Best
            dict_infeasible_index, n_infeasible_var = dict_infeasible_index_Best, n_infeasible_var_Best
            CP(dict_infeasible_index)
        end
        global solution_k, dict_xlb_k, dict_xub_k, dict_LP_int, solution_k_1
        global solution_k_check,dict_xlb_k_check,dict_xub_k_check, dict_LP_int_check, solution_k_1_check
        global solution_k_Best, dict_xlb_k_Best, dict_xub_k_Best, dict_LP_int_Best, solution_k_1_Best
        global dict_infeasible_index, n_infeasible_var
        global dict_infeasible_index_check, n_infeasible_var_check
        global dict_infeasible_index_Best, n_infeasible_var_Best
        global T    
        println("===================================================================================================================")
        println("   ")
        println("[T, N]: [" ,T,",", N, "]|n_infeasible_var_Best: ", n_infeasible_var_Best, "|dist: ", dist_Best)
        println("   ")
        println("===================================================================================================================")
        continue
    catch
        continue
    end
    

end
        

Reupdate(dict_xlb, dict_xub_k_check)
Set_Type()
optimize!(m)
termination_status(m)!=MOI.OPTIMAL


@time begin
    
end 


#UnSet_Type()
############################################################
# 흠.. ? 이렇게만 해도 되네?
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
    optimize!(m)
    termination_status(m)!=MOI.OPTIMAL
end
##################################################


@time begin
    Set_Type()
    optimize!(m)
    termination_status(m)==MOI.OPTIMAL
end