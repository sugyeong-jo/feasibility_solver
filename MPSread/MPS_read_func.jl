using ProgressBars
using MathOptInterface
using JuMP, Cbc
const MOI = MathOptInterface
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

################################
# upper bound는 느슨하게 업데이트
function Update(dict_LP_int)
    for i in keys(dict_LP_int) #8437개
        try 
            if lower_bound(i)+0.01 < dict_LP_int[i] < upper_bound(i)-0.01
                JuMP.set_upper_bound(i,dict_LP_int[i])
            end
        catch end
    end
end


function LPUpdate(dict_LP_int)
    for i in keys(dict_LP_int) #8437개
        if lower_bound(i)+0.01 < dict_LP_int[i] < upper_bound(i)-0.01
            JuMP.set_lower_bound(i,dict_LP_int[i])
            JuMP.set_upper_bound(i,dict_LP_int[i])
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

function LP_solve(m)
    optimize!(m)        
    solution_k = Dict() # LP솔루션 저장
    for i in keys(var)
        solution_k[i] = JuMP.value(i)
    end   
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


function LP_solve(m,dict_xlb_k,dict_xub_k)
    Reupdate(dict_xlb_k, dict_xub_k)
    optimize!(m)        
    solution_k = Dict() # LP솔루션 저장
    for i in keys(var)
        solution_k[i] = JuMP.value(i)
    end   
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
        var_s=A[:,var_idx[i],:] # 특정 variable이 속한 모든 constraint
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
                if upper_bound(idx_var[const_s_var_index[j]]) >  u_new+0.01
                    set_upper_bound(idx_var[const_s_var_index[j]],u_new)
                end
            else
                u_new = var_ub[idx_var[const_s_var_index[j]]] + (u[constraint_index]-L_min)/const_s_coef[j]
                if lower_bound(idx_var[const_s_var_index[j]]) < u_new-0.01
                    set_lower_bound(idx_var[const_s_var_index[j]],u_new)
        
                end
            end
        end
     end
end


