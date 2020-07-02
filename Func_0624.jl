using ProgressBars
using MathOptInterface
using JuMP, Cbc
#const MOI = MathOptInterface
using Cbc
using MathProgBase
using GLPKMathProgInterface
using Random


function Set_Type()
    for i in tqdm([k for (k,v) in var if v==:Bin])
        JuMP.set_binary(i)
    end
    for i in tqdm([k for (k,v) in var if v==:Int])
        JuMP.set_integer(i)
    end
end

function UnSet_Type()
    for i in tqdm([k for (k,v) in var if v==:Bin])
        JuMP.unset_binary(i)
    end
    for i in tqdm([k for (k,v) in var if v==:Int])
        JuMP.unset_integer(i)
    end
end

function Initial_bound()
    for i in tqdm(keys(var))
        JuMP.set_lower_bound(i,var_lb[i])
        JuMP.set_upper_bound(i,var_ub[i])
    end
end

function Initial_bound_toy()
    for i in tqdm(keys(var))
        JuMP.set_lower_bound(i,var_lb[i])
        JuMP.set_upper_bound(i,var_ub[i])
    end
end

################################
function Update_UB(dict_LP_int)
    println("-----upper bound update")
    for i in keys(dict_LP_int) #8437개
        try 
            if lower_bound(i) < dict_LP_int[i] < upper_bound(i)
                JuMP.set_upper_bound(i,dict_LP_int[i])
            end
        catch end
    end
end

function Update_LB(dict_LP_int)
    println("-----lower bound update")
    for i in keys(dict_LP_int) #8437개
        try 
            if lower_bound(i) < dict_LP_int[i] < upper_bound(i)
                JuMP.set_lower_bound(i,dict_LP_int[i])
            end
        catch end
    end
end


function Update_UL(dict_LP_int)
    println("-----lower/upper bound update")
    for i in keys(dict_LP_int) #8437개
        if lower_bound(i) <= dict_LP_int[i] <= upper_bound(i)
            JuMP.set_lower_bound(i,dict_LP_int[i])
            JuMP.set_upper_bound(i,dict_LP_int[i])
        end
    end
end


# bound 재조정
"""
function Reupdate(dict_xlb_k, dict_xub_k)
    Initial_bound()
    for i in keys(dict_xlb_k) 
        if dict_xlb_k[i] > var_lb[i]
            set_lower_bound(i,dict_xlb_k[i])
        else continue 
        end
    end
    for i in keys(dict_xub_k) 
        if dict_xub_k[i] < var_ub[i]
            set_upper_bound(i,dict_xub_k[i])
        else continue 
        end
    end
end
"""
function Reupdate(dict_xlb_k, dict_xub_k)
    println("-----bound re-update")
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

function LP_solve(m)
    println("-----LP solver")
    optimize!(m)        
    solution_k = Dict() # LP솔루션 저장
    for i in all_variables(m)
        solution_k[i] = JuMP.value(i)
    end   
    dict_LP_int = Dict()     # LP로 구해진 integer set
    solution_k_1 = Dict()    # LP로 구해진 integer set을 round 한 solution set
    for i in [k for (k,v) in var if v==:Bin]
        if solution_k[i] ==  float(0)  ##
            dict_LP_int[i] = solution_k[i]            
        elseif solution_k[i] ==  float(1)
            dict_LP_int[i] = solution_k[i]
        else 
            solution_k_1[i] = round(solution_k[i])
        end
    end

    var_set = [k for (k,v) in var if v==:Int]
    for i in var_set
        if solution_k[i]!=round(solution_k[i])
            solution_k_1[i] = round(solution_k[i])
        else dict_LP_int[i] = solution_k[i]
        end
    end

    global var
 
    dict_xlb_k = Dict() # lower bound dictionary
    dict_xub_k = Dict() # lower bound dictionary
    for i in tqdm(keys(var))
        dict_xlb_k[i] = lower_bound(i)
        dict_xub_k[i] = upper_bound(i)
    end
    return solution_k, dict_xlb_k, dict_xub_k, dict_LP_int, solution_k_1
end


function LP_solve(m,dict_xlb_k,dict_xub_k)
    Reupdate(dict_xlb_k, dict_xub_k)
    println("-----LP solver")
    optimize!(m)        
    solution_k = Dict() # LP솔루션 저장
    for i in keys(var)
        solution_k[i] = JuMP.value(i)
    end   
    println("Step 2| if x* is integer, return(x*) --> dict_LP_int")
    println("Step 3| let x~:=[x*] (rounding of x*) --> solution_k_1")
    dict_LP_int = Dict()     # LP로 구해진 integer set
    solution_k_1 = Dict()    # LP로 구해진 integer set을 round 한 solution set
    for i in [k for (k,v) in var if v==:Bin]
        if solution_k[i] ==  float(0)
            dict_LP_int[i] = solution_k[i]            
        elseif solution_k[i] ==  float(1)
            dict_LP_int[i] = solution_k[i]
        else 
            solution_k_1[i] = round(solution_k[i])
        end
    end

    var_set = [k for (k,v) in var if v==:Int]
    for i in var_set
        if solution_k[i]!=round(solution_k[i])
            solution_k_1[i] = round(solution_k[i])
        else dict_LP_int[i] = solution_k[i]
        end
    end
 
    dict_xlb_k = Dict() # lower bound dictionary
    dict_xub_k = Dict() # lower bound dictionary
    for i in tqdm(keys(var))
        dict_xlb_k[i] = lower_bound(i)
        dict_xub_k[i] = upper_bound(i)
    end
    return solution_k, dict_xlb_k, dict_xub_k, dict_LP_int, solution_k_1
end

##############################################
#feasible check
##############################################
# infeasible 한 index dictionary
function Infeasible_Check(solution_k)
    dict_infeasible_index = Dict()
    for i in [k for (k,v) in var if v==:Int]
        if solution_k[i] !=  round(solution_k[i]) 
            dict_infeasible_index[i] = :Int
        end
    end
    for i in [k for (k,v) in var if v==:Bin]
        if solution_k[i] !=  round(solution_k[i]) 
            dict_infeasible_index[i] = :Bin
        end
    end

    n_infeasible_var= length(keys(dict_infeasible_index)) #n_infeasible_var: infeasible한  개수
    println("Number of infeasible variable: $n_infeasible_var.")
    return dict_infeasible_index, n_infeasible_var
end

################################################
#  Projection
################################################
# step 4 projection
# score_list: feasible LP의 score list!
# score_list_1: round 후 score
function Projection(solution_k,solution_k_1)
    println("Step 11| Projection")
    score = Array{Float64}(undef,0)
    score_list = Any[]
    for i in keys(solution_k_1)
        if solution_k_1[i] == upper_bound(i)
            val = upper_bound(i)-solution_k[i]
            push!(score, abs(val) )
            push!(score_list, (i, abs(val)))

        elseif solution_k_1[i] == lower_bound(i)
            val = solution_k[i]-lower_bound(i)
            push!(score, abs(val))
            push!(score_list, (i, abs(val)))
            
        else
            val = abs(solution_k[i]-solution_k_1[i])
            push!(score,val)
            push!(score_list, (i, abs(val)))
        end
    end
    dist = sum(score)

    top_score_list=sort!(score_list, by = x -> abs(x[2]), rev = true)
    top_idx = Array{Any}(undef,0)
    for i in 1:length(top_score_list)
        push!(top_idx,top_score_list[i][1]) 
    end

    println("The dist is $dist")
    return top_score_list,top_idx, dist
end

function Round_change(dict_x_idx,solution_k,solution_k_1)
    println("Step 12| move the T components x~ with largest score")
    solution_k_2 = Dict()
    
    for x_idx in dict_x_idx
        if var[x_idx] == :Bin
            if (solution_k[x_idx]-solution_k_1[x_idx])<0
                solution_k_2[x_idx] = solution_k_1[x_idx]-1
            else
                solution_k_2[x_idx] = solution_k_1[x_idx]+1
            end
        else
            solution_k_2[x_idx] = rand(solution_k_1[x_idx]:upper_bound(x_idx))            
        end
    end
    return solution_k_2
end

function Round_change_pert(dict_x_idx,solution_k,solution_k_1)
    println("Step ##| purtubation")
    solution_k_2 = Dict()
    w = rand(Float64, 1)[1]
    if w>0.5
        t_w = 2*w*(1-w)
    else
        t_w = 1-2*w*(1-w)
    end

    for x_idx in dict_x_idx
        if (solution_k[x_idx]-solution_k_1[x_idx])<0
            solution_k_2[x_idx] = round(solution_k_1[x_idx]-t_w)
        else
            solution_k_2[x_idx] = round(solution_k_1[x_idx]+t_w)
        end
    end
    return solution_k_2
end

function dict_idx_sol(x_idx, solution_k_1)
    dict_round_select = Dict() 
    for i in x_idx
        dict_round_select[i]=solution_k_1[i]
    end
    return dict_round_select
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
    for i in keys(dict_infeasible_index)
        var_s=A[:,var_idx[i],:] # 특정 variable이 속한 모든 constraint
        var_s_n=1:var_s.colptr[2]-1 # 속해 있는 constraint 개수
        for j in var_s.rowval[var_s_n]
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
            if const_s_coef[j] >= 0
                #print("variable iteration:",j,"|")
                push!(L_min,const_s_coef[j]*var_lb[idx_var[const_s_var_index[j]]])
            else
                push!(L_min,const_s_coef[j]*var_ub[idx_var[const_s_var_index[j]]])
            end
        end
        L_min=sum(L_min)

        for j in const_s_n
            if const_s_coef[j] >= 0
                u_new = var_lb[idx_var[const_s_var_index[j]]] + (u[constraint_index]-L_min)/const_s_coef[j]
                if upper_bound(idx_var[const_s_var_index[j]]) >  u_new
                    set_upper_bound(idx_var[const_s_var_index[j]],u_new)
                end
            else
                #u_new = var_ub[idx_var[const_s_var_index[j]]] + (u[constraint_index]-L_min)/const_s_coef[j]
                #if lower_bound(idx_var[const_s_var_index[j]]) < u_new
                #    set_lower_bound(idx_var[const_s_var_index[j]],u_new)        
                #end
                continue
            end
        end
     end
end


#=
Input: dict_infeasible_index, infeasible index로 variable의 이름이 인풋으로 들어갑니다.

[part 1] 
 : input으로 들어온 모든 variable이 속한 constraint을 모두 sort
  - infeasible variable index 중 하나를 선택
  -> 이 variable index가 속한 모든 constraint를 구함
  -> 모든 constraint의 variable sorting
  
  - 이 과정을 모든 infeasible variable index 하나하나마다 하여 unique함

[part 2]
 : infeasiebl variable들 하나하나를 bound strengthening 
 - paper에서 (1), (2)를 구현

Output: bound에 맞게 lower bound/ upper bound 수정됨


이 부분을 구현할 예정입니다.
##constraint propagation
Pn = 
Pf = 
Q = 
R = 

=#

function BoundStrength(dict_infeasible_index)
    const_list =  Array{Int}(undef,0)
    println("The all constraints including infeasible variables")
    for i in keys(dict_infeasible_index)
        var_s=A[:,var_idx[i],:] # 특정 variable이 속한 모든 constraint
        var_s_n=1:var_s.colptr[2]-1 # 속해 있는 constraint 개수
        for j in var_s.rowval[var_s_n]
            push!(const_list, j)            
        end
    end

    infeasible_var_const=unique(const_list) #dict infeasible index가 속한 모든 constraint에 속한 variable set
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
            if const_s_coef[j] >= 0
                #print("variable iteration:",j,"|")
                push!(L_min,const_s_coef[j]*lb[const_s_var_index[j]])
            else
                push!(L_min,const_s_coef[j]*ub[const_s_var_index[j]])
            end
        end
        L_min=sum(L_min)

        for j in const_s_n
            if const_s_coef[j] >= 0
                u_new = lb[const_s_var_index[j]] + (u[constraint_index]-L_min)/const_s_coef[j]
                if upper_bound(idx_var[const_s_var_index[j]]) >  u_new
                    set_upper_bound(idx_var[const_s_var_index[j]],u_new)
                end
            else
                u_new = var_ub[idx_var[const_s_var_index[j]]] + (u[constraint_index]-L_min)/const_s_coef[j]
                if lower_bound(idx_var[const_s_var_index[j]]) < u_new
                    set_lower_bound(idx_var[const_s_var_index[j]],u_new)        
                end
                continue
            end
        end
     end
end

function BoundStrength(constraint_index)
    global A

    #println("bound strengthening")
    const_s=A[constraint_index,:,:]
    n_var = const_s.colptr[2]-1
    const_s_n=1:n_var # 속해 있는 variable 개수
    #println("The number of strengthening variables: $n_var ")
    const_s_coef=const_s.nzval[const_s_n,]#특정 constraint에 해당하는 coefficient
    const_s_var_index=const_s.rowval[const_s_n]#특정 constraint에 해당하는 variable index
    L_min=  Array{Float64}(undef,0)
    for j in const_s_n
        if const_s_coef[j] >= 0
            #print("variable iteration:",j,"|")
            push!(L_min,const_s_coef[j]*lb[const_s_var_index[j]])
        else
            push!(L_min,const_s_coef[j]*ub[const_s_var_index[j]])
        end
    end
    L_min=sum(L_min)

    for j in const_s_n
        if const_s_coef[j] >= 0
            u_new = var_lb[idx_var[const_s_var_index[j]]] + (u[constraint_index]-L_min)/const_s_coef[j]
            if upper_bound(idx_var[const_s_var_index[j]]) >  u_new
                set_upper_bound(idx_var[const_s_var_index[j]],u_new)
            end
        else
            u_new = var_ub[idx_var[const_s_var_index[j]]] + (u[constraint_index]-L_min)/const_s_coef[j]
            if lower_bound(idx_var[const_s_var_index[j]]) < u_new
                set_lower_bound(idx_var[const_s_var_index[j]],u_new)        
            end
            continue
        end
    end
end

##########################################################################################################################

##########################################################################################################################
function Round_FP2()
    idx_int_bin = deepcopy(vcat([k for (k,v) in var if v==:Int],  [k for (k,v) in var if v==:Bin]))
    
end
"""
function Propagate(D)
    Q =
    R =
    while length(Q) != 0
        p = Pop(Q)
        D = 
    end
end
"""
###################################################################################################################################################
function Result(path, solution_k, nIter)    
    sol_check_bound = Array{Any}(undef,0)
    sol_check_type = Array{Any}(undef,0)
    for x in keys(solution_k)
        push!(sol_check_bound,var_lb[x]<=solution_k[x]<=var_ub[x])
    end

    for x in [k for (k,v) in var if v==:Bin]
        if isapprox(solution_k[x], 0, atol=1e-3) ==true
            continue 
        elseif isapprox(solution_k[x], 1, atol=1e-3) ==true
            continue
        else
            println(solution_k[x])
            push!(sol_check_type, (x,solution_k[x], :Bin))
        end
    end

    for x in [k for (k,v) in var if v==:Int]
        if isapprox(abs(solution_k[x] - round(solution_k[x])), 0, atol=1e-3) ==false  
            push!(sol_check_type, (x,solution_k[x], :Int))
        else
            continue
        end
    end

    check_optimal = "FP 1.0"
    check_bound = length(findall(x->false, sol_check_bound))
    check_type = length(sol_check_type)
    if check_type != 0
        obj_value = "null"
    else
        obj_value = objective_value(m)
    end
    t_total = t_problemLaod+t_run
    println("The problem load time is : $t_problemLaod")
    println("The running time is : $t_run s")
    println("The total time is : $t_total")
    println("Is it optimal?: $check_optimal")
    println("The FP Iteration is : $nIter")
    println("The objective value is : $obj_value" )
    println("The the number of unsatisfied variable is:")
    println("    - The the number of unsatisfied bound is: $check_bound" )
    println("    - The the number of unsatisfied type is: $check_type")
    df=DataFrame(
        Name = [filename],
        Total = [t_total],
        Problem = [t_problemLaod],
        Run = [t_run],
        Optimal = [check_optimal],
        Check_Bound = [check_bound],
        check_type = [check_type],
        obj_value = [obj_value],
        nIter = [nIter]
    )

    CSV.write(path, df) 

end


function check_result()
    println("===================================================================================================================")
    println("   ")
    println("[T, nIter, pertub]: [$T, $nIter, $KK]| n_infeasible_var: $n_infeasible_var, |dist: $dist")
    println("   ")
    println("===================================================================================================================")   
end


####
# Perturbation
####

function Perturbation(sol0, sol_tilde0, TT, itg_idx0)
    dictionary1 = Dict(zip(itg_idx0, (sol0-sol_tilde0)[itg_idx0]))
    sortedDict = sort(collect(dictionary1), by =x->abs(x[2]), rev = true);
    println((sol0-sol_tilde0)[15565])
    println(minimum((sol0-sol_tilde0)[itg_idx0]))
    println(sortedDict[1:3])

    for j in 1:TT
        topTT = sortedDict[j]
        if topTT[2,]>0
            sol_tilde0[topTT[1,]] += 1
        else
            sol_tilde0[topTT[1,]] += -1
        end
    end
    return sol_tilde0
end

#########################################################################################
using MathOptInterface
using JuMP, Cbc
#const MOI = MathOptInterface
using SparseArrays
using ArgParse
using QPSReader

##########################
s = ArgParseSettings()
s.description = "Find Feasible Solution"
s.commands_are_required = true
s.version = "1.0"
s.add_version = true

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "filename"
            help = "a positional argument"
            required = true
        "Sec"
            help = "a solving time"
            arg_type = Int
            required = true
    end

    return parse_args(s)
end


function MPS_read_full(filepath, Sec)
    m = read_from_file(filepath)
    optimizer=optimizer_with_attributes(Cbc.Optimizer,  "maxSolutions" => 1, "logLevel "=>1, "seconds" => Sec ,"allowableGap "=>70)
    set_optimizer(m, optimizer)

    # constraint
    list_of_constraint_types(m)
    con = Dict()
    con_rhs = Dict()
    idx_con = Dict()
    con_idx = Dict()

    println("Equal To constraints load")
    constraint = all_constraints(m, GenericAffExpr{Float64,VariableRef}, MathOptInterface.EqualTo{Float64})
    for i in (1:length(constraint))
        con_name = name(constraint[i])
        con_rhs[con_name] = constraint_object(constraint[i]).set.value
        con[con_name]= :Equal
        idx_con[constraint[i].index.value] = con_name
        con_idx[con_name] = constraint[i].index.value
    end

    println("Greater Than constraints load")
    constraint = all_constraints(m, GenericAffExpr{Float64,VariableRef}, MathOptInterface.GreaterThan{Float64})
    for i in (1:length(constraint))
        con_name = name(constraint[i])
        con_rhs[con_name] = constraint_object(constraint[i]).set.lower
        con[con_name]= :Greater
        idx_con[constraint[i].index.value] = con_name
        con_idx[con_name] = constraint[i].index.value
    end

    println("Less Than constraints load")
    constraint = all_constraints(m, GenericAffExpr{Float64,VariableRef}, MathOptInterface.LessThan{Float64})
    for i in (1:length(constraint))
        con_name = name(constraint[i])
        con_rhs[con_name] = constraint_object(constraint[i]).set.upper
        con[con_name]= :Less
        idx_con[constraint[i].index.value] = con_name
        con_idx[con_name] = constraint[i].index.value
    end


    # variable
    var = Dict()
    var_idx = Dict()
    idx_var = Dict()
    var_lb = Dict()
    var_ub = Dict()
    var_ref = Dict()
    println("All variable load")
    idx = [1:1:length(all_variables(m));]
    var_name = all_variables(m)
    for i in 1:length(all_variables(m)) 
        var_=var_name[i]
        var[var_] = :Con
        var_idx[var_] = idx[i] 
        idx_var[idx[i]] = var_
        var_lb[var_] = -Inf
        var_ub[var_] = Inf
        var_ref[var_] = :except
    end




    println("Interval variables load")
    variable = all_constraints(m, VariableRef, MathOptInterface.Interval{Float64})
    for i in (1:length(variable))
        var_name = constraint_object(variable[i]).func
        var_lb[var_name] = constraint_object(variable[i]).set.lower
        var_ub[var_name] = constraint_object(variable[i]).set.upper
        var_ref[var_name] = :Interval
    end

    println("Equal To variables load")
    variable = all_constraints(m, VariableRef, MathOptInterface.EqualTo{Float64})
    for i in (1:length(variable))
        var_name = constraint_object(variable[i]).func
        var_lb[var_name] = constraint_object(variable[i]).set.value
        var_ub[var_name] = constraint_object(variable[i]).set.value
        var_ref[var_name] = :EqualTo
    end

    println("Integer variables load")
    variable = all_constraints(m, VariableRef, MathOptInterface.Integer)
    for i in (1:length(variable))
        var_name = constraint_object(variable[i]).func
        var[var_name] = :Int
    end

    println("Binary variables load")
    variable = all_constraints(m, VariableRef, MathOptInterface.ZeroOne)
    for i in (1:length(variable))
        var_name = constraint_object(variable[i]).func
        var[var_name] = :Bin
        var_lb[var_name] = 0
        var_ub[var_name] = 1
    end

    println("Constraint load")
    # Sparse matrix
    I_I = Int[]
    J = Int[]
    V = Float64[]
    u = Dict()
    l = Dict()


    con_set = [k for (k,v) in con if v==:Less]
    for i in (con_set)
        con_term =collect(linear_terms(constraint_object(constraint_by_name(m, i )).func))
        u[con_idx[i]] = con_rhs[i]
        l[con_idx[i]] = -Inf
        for j in 1:length(con_term)
            push!(I_I, con_idx[i])
            push!(J, var_idx[con_term[j][2]])
            push!(V, con_term[j][1])
        end
    end

    con_set = [k for (k,v) in con if v==:Greater]
    for i in (con_set)
        con_term =collect(linear_terms(constraint_object(constraint_by_name(m, i )).func))
        u[con_idx[i]] = -(con_rhs[i])
        l[con_idx[i]] = -Inf

        for j in 1:length(con_term)
            push!(I_I, con_idx[i])
            push!(J, var_idx[con_term[j][2]])
            push!(V, -(con_term[j][1]))
        end
    end
    con_set = [k for (k,v) in con if v==:Equal]
    for i in (con_set)
        con_term =collect(linear_terms(constraint_object(constraint_by_name(m, i )).func))
        u[con_idx[i]] = con_rhs[i]
        l[con_idx[i]] = con_rhs[i]
        for j in 1:length(con_term)
            push!(I_I, con_idx[i])
            push!(J, var_idx[con_term[j][2]])
            push!(V, con_term[j][1])
        end
    end

    A = sparse(I_I,J,V)

    global m, var, var_lb, var_ub
    #Interval, EqualTo를 lower bound/ upper bound로 할당
    variable = all_constraints(m, VariableRef, MathOptInterface.EqualTo{Float64})
    for i in variable
        delete(m,i)
    end
    variable = all_constraints(m, VariableRef, MathOptInterface.Interval{Float64})
    for i in variable
        delete(m,i)
    end
    Initial_bound()


    x = all_variables(m)
    lb = lower_bound.(x)
    ub = upper_bound.(x)



end


function MPS_read_FP(filepath, Sec)  
    m = Model()
    println("================MPS load=================")
    println("The file name is: $filename ")
    println("The limit time of solver is: $Sec s") 
    println("=========================================") 
    m = read_from_file(filepath)
    optimizer=optimizer_with_attributes(Cbc.Optimizer,  "maxSolutions" => 1, "logLevel "=>1, "seconds" => Sec ,"allowableGap "=>70)
    set_optimizer(m, optimizer)

    qps = readqps(filepath)
    var_idx_qps = qps.varindices


    # constraint
    con_idx = qps.conindices
    idx_con = Dict()
    u = qps.ucon
    l = qps.lcon
    
    for con_name in keys(con_idx)
        idx_con[con_idx[con_name]] = con_name
    end
    A = sparse(qps.arows, qps.acols, qps.avals, qps.ncon, qps.nvar)


    #########################################
    # variable
    var = Dict()
    var_idx = Dict()
    idx_var = Dict()
    var_lb = Dict()
    var_ub = Dict()
    qps_var = Dict()

    for var_name in (all_variables(m))
        var[var_name] = :Con
        var_idx[var_name] = var_idx_qps[string(var_name)]
        idx_var[var_idx_qps[string(var_name)]] = var_name
        qps_var[string(var_name)] = var_name
    end

    variable = all_constraints(m, VariableRef, MathOptInterface.Integer)
    for i in (1:length(variable))
        var_name = constraint_object(variable[i]).func
        var[var_name] = :Int
    end
    variable = all_constraints(m, VariableRef, MathOptInterface.ZeroOne)
    for i in (1:length(variable))
        var_name = constraint_object(variable[i]).func
        var[var_name] = :Bin
    end
    var_lb = Dict()
    var_ub = Dict()
    for i in 1:qps.nvar
        var_name=qps.varnames[i]
        var_lb[qps_var[var_name]] = qps.lvar[i]
        var_ub[qps_var[var_name]] = qps.uvar[i]
    end
    ##########################################
    for var_name in (all_variables(m))
        var[var_name] = :Con
        var_idx[var_name] = var_idx_qps[string(var_name)]
        idx_var[var_idx_qps[string(var_name)]] = var_name
        qps_var[string(var_name)] = var_name
    end

    variable = all_constraints(m, VariableRef, MathOptInterface.Integer)
    for i in (1:length(variable))
        var_name = constraint_object(variable[i]).func
        var[var_name] = :Int
    end
    variable = all_constraints(m, VariableRef, MathOptInterface.ZeroOne)
    for i in (1:length(variable))
        var_name = constraint_object(variable[i]).func
        var[var_name] = :Bin
    end

    itg = Any[]
    for i in qps.varnames
        if var[qps_var[i]]==:Con
            push!(itg, false)
        else 
            push!(itg, true)
        end
    end

    
    ## 이 부분이 추가되었습니다. (2020.06.23)
    # variable 순서대로 VariableRef 형식으로 나오는 array
    var_name = Any[]
    for i in qps.varnames
        push!(var_name, qps_var[i])
    end

    global m, var, var_lb, var_ub
    #Interval, EqualTo를 lower bound/ upper bound로 할당
    variable = all_constraints(m, VariableRef, MathOptInterface.EqualTo{Float64})
    for i in variable
        delete(m,i)
    end
    variable = all_constraints(m, VariableRef, MathOptInterface.Interval{Float64})
    for i in variable
        delete(m,i)
    end
    Initial_bound()

    return lb, ub, itg, l, u, A, m,var_name, con_idx, idx_con, var_idx_qps, var, var_idx, idx_var, var_lb, var_ub, qps_var, qps
end

function MPS_read(filepath, Sec)  
    m = Model()
    println("================MPS load=================")
    println("The file name is: $filename ")
    println("The limit time of solver is: $Sec s") 
    println("=========================================") 
    m = read_from_file(filepath)
    optimizer=optimizer_with_attributes(Cbc.Optimizer,  "maxSolutions" => 1, "logLevel "=>1, "seconds" => Sec ,"allowableGap "=>70)
    set_optimizer(m, optimizer)

    qps = readqps(filepath)
    var_idx_qps = qps.varindices


    # constraint
    #con_idx = qps.conindices
    #idx_con = Dict()
    #u = qps.ucon
    #l = qps.lcon
    
    #for con_name in keys(con_idx)
    #    idx_con[con_idx[con_name]] = con_name
    #end
    #A = sparse(qps.arows, qps.acols, qps.avals, qps.ncon, qps.nvar)


    #########################################
    # variable
    var = Dict()
    var_idx = Dict()
    idx_var = Dict()
    var_lb = Dict()
    var_ub = Dict()
    qps_var = Dict()


    for var_name in (all_variables(m))
        var[var_name] = var_name
        var_idx[var_name] = var_idx_qps[string(var_name)]
        idx_var[var_idx_qps[string(var_name)]] = var_name
        qps_var[string(var_name)] = var_name
    end
    """
    variable = all_constraints(m, VariableRef, MathOptInterface.Integer)
    for i in (1:length(variable))
        var_name = constraint_object(variable[i]).func
        var[var_name] = :Int
    end
    variable = all_constraints(m, VariableRef, MathOptInterface.ZeroOne)
    for i in (1:length(variable))
        var_name = constraint_object(variable[i]).func
        var[var_name] = :Bin
    end

    itg = Any[]
    for i in qps.varnames
        if var[qps_var[i]]==:Con
            push!(itg, false)
        else 
            push!(itg, true)
        end
    end

    """

    var_lb = Dict()
    var_ub = Dict()
    for i in 1:qps.nvar
        var_name=qps.varnames[i]
        var_lb[qps_var[var_name]] = qps.lvar[i]
        var_ub[qps_var[var_name]] = qps.uvar[i]
    end
    ##########################################
    
    ## 이 부분이 추가되었습니다. (2020.06.23)
    # variable 순서대로 VariableRef 형식으로 나오는 array
    #var_name = Any[]
    #for i in qps.varnames
    #    push!(var_name, qps_var[i])
    #end

    global m, var, var_lb, var_ub
    #Interval, EqualTo를 lower bound/ upper bound로 할당
    variable = all_constraints(m, VariableRef, MathOptInterface.EqualTo{Float64})
    for i in variable
        delete(m,i)
    end
    variable = all_constraints(m, VariableRef, MathOptInterface.Interval{Float64})
    for i in variable
        delete(m,i)
    end

    Initial_bound()
    return m
end