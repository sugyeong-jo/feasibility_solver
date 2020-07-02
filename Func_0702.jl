using ProgressBars
using MathOptInterface
using JuMP, Cbc
#const MOI = MathOptInterface
using Cbc
using MathProgBase
using GLPKMathProgInterface
using Random

using SparseArrays
using ArgParse
using QPSReader

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

    variable = all_constraints(m, VariableRef, MathOptInterface.EqualTo{Float64})
    for i in variable
        delete(m,i)
    end
    variable = all_constraints(m, VariableRef, MathOptInterface.Interval{Float64})
    for i in variable
        delete(m,i)
    end

    Initial_bound()


    qps = readqps(filepath)
    var_qps = Dict()
    for i in keys(var)
        var_qps[i] = string(i)
    end

    Coeff = Dict()
    for i in 1:length(var)
        Coeff[qps.varnames[i]] = qps.c[i]
    end

    c = Array{Float64}(undef,0)
    for i in 1:length(var)
        push!(c,Coeff[var_qps[idx_var[i]]])
    end


    x = all_variables(m)
    ub = upper_bound.(x)
    lb = lower_bound.(x)
    itg_idx = is_binary.(x)+is_integer.(x)


    return m, con_idx, idx_con, A, l, u, var, var_idx, idx_var, var_lb, var_ub, c,itg_idx

end


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
 
    dict_xlb_k = Dict() # lower bound dictionary
    dict_xub_k = Dict() # lower bound dictionary
    for i in tqdm(keys(var))
        dict_xlb_k[i] = lower_bound(i)
        dict_xub_k[i] = upper_bound(i)
    end
    return solution_k, dict_xlb_k, dict_xub_k, dict_LP_int, solution_k_1
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
            u_new = var_ub[idx_var[const_s_var_index[j]]] + (u[constraint_index]-L_min)/const_s_coef[j]
            if lower_bound(idx_var[const_s_var_index[j]]) < u_new
                set_lower_bound(idx_var[const_s_var_index[j]],u_new)        
            end
            continue
        end
    end
end


