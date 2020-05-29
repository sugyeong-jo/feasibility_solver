using ProgressBars
using MathOptInterface
using JuMP, Cbc
const MOI = MathOptInterface
using Cbc
using MathProgBase
using GLPKMathProgInterface
using Random


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



    
