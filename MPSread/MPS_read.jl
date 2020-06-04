println("================패키지 부르기=================")
using ProgressBars
using MathOptInterface
using JuMP, Cbc
const MOI = MathOptInterface
using SparseArrays

#using MathProgBase
#using GLPKMathProgInterface
#using Random
#using Pkg
#Pkg.clone("https://github.com/odow/LPWriter.jl")
#using LPWriter
include("/HDD/sugyeong/Julia/feasibility_solver/MPSread/MPS_read_func.jl")
##########################
println("================MPS load=================")
val, t_problemLaod, bytes, gctime, memallocs = @timed begin
    m = read_from_file("/HDD/Workspace/CLT/mps/processing/CPLEX_file/R100701001_2_cplex.mps")
#    m = read_from_file("/HDD/Workspace/CLT/mps/processing/CPLEX_file/1_R170718004_veri_cplex.mps")
#    m = read_from_file("/HDD/Workspace/CLT/mps/processing/CPLEX_file/4_R170725003_veri_cplex.mps")

    optimizer=optimizer_with_attributes(Cbc.Optimizer,  "maxSolutions" => 1, "logLevel "=>1, "seconds" => 500,"allowableGap "=>70)
    set_optimizer(m, optimizer)

    # constraint
    list_of_constraint_types(m)
    con = Dict()
    con_rhs = Dict()
    idx_con = Dict()
    con_idx = Dict()

    constraint = all_constraints(m, GenericAffExpr{Float64,VariableRef}, MathOptInterface.EqualTo{Float64})
    for i in tqdm(1:length(constraint))
        con_name = name(constraint[i])
        con_rhs[con_name] = constraint_object(constraint[i]).set.value
        con[con_name]= :Equal
        idx_con[constraint[i].index.value] = con_name
        con_idx[con_name] = constraint[i].index.value
    end
    constraint = all_constraints(m, GenericAffExpr{Float64,VariableRef}, MathOptInterface.GreaterThan{Float64})
    for i in tqdm(1:length(constraint))
        con_name = name(constraint[i])
        con_rhs[con_name] = constraint_object(constraint[i]).set.lower
        con[con_name]= :Greater
        idx_con[constraint[i].index.value] = con_name
        con_idx[con_name] = constraint[i].index.value
    end
    constraint = all_constraints(m, GenericAffExpr{Float64,VariableRef}, MathOptInterface.LessThan{Float64})
    for i in tqdm(1:length(constraint))
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
    for var_name in tqdm(all_variables(m))
        var[var_name] = :Con
        var_idx[var_name] = var_name.index.value
        idx_var[var_name.index.value] = var_name
        var_lb[var_name] = -Inf
        var_ub[var_name] = Inf
        var_ref[var_name] = :except
    end

    variable = all_constraints(m, VariableRef, MathOptInterface.Interval{Float64})
    for i in tqdm(1:length(variable))
        var_name = constraint_object(variable[i]).func
        var_lb[var_name] = constraint_object(variable[i]).set.lower
        var_ub[var_name] = constraint_object(variable[i]).set.upper
        var_ref[var_name] = :Interval
    end


    variable = all_constraints(m, VariableRef, MathOptInterface.EqualTo{Float64})
    for i in tqdm(1:length(variable))
        var_name = constraint_object(variable[i]).func
        var_lb[var_name] = constraint_object(variable[i]).set.value
        var_ub[var_name] = constraint_object(variable[i]).set.value
        var_ref[var_name] = :EqualTo
    end

    variable = all_constraints(m, VariableRef, MathOptInterface.Integer)
    for i in tqdm(1:length(variable))
        var_name = constraint_object(variable[i]).func
        var[var_name] = :Int
    end

    variable = all_constraints(m, VariableRef, MathOptInterface.ZeroOne)
    for i in tqdm(1:length(variable))
        var_name = constraint_object(variable[i]).func
        var[var_name] = :Bin
        var_lb[var_name] = 0
        var_ub[var_name] = 1
    end


    # Sparse matrix
    I = Int[]
    J = Int[]
    V = Float64[]
    u = Dict()
    l = Dict()


    con_set = [k for (k,v) in con if v==:Less]
    for i in tqdm(con_set)
        con_term =collect(linear_terms(constraint_object(constraint_by_name(m, i )).func))
        u[con_idx[i]] = con_rhs[i]
        l[con_idx[i]] = -Inf
        for j in 1:length(con_term)
            push!(I, con_idx[i])
            push!(J, var_idx[con_term[j][2]])
            push!(V, con_term[j][1])
        end
    end

    con_set = [k for (k,v) in con if v==:Greater]
    for i in tqdm(con_set)
        con_term =collect(linear_terms(constraint_object(constraint_by_name(m, i )).func))
        u[con_idx[i]] = -(con_rhs[i])
        l[con_idx[i]] = -Inf

        for j in 1:length(con_term)
            push!(I, con_idx[i])
            push!(J, var_idx[con_term[j][2]])
            push!(V, -(con_term[j][1]))
        end
    end
    con_set = [k for (k,v) in con if v==:Equal]
    for i in tqdm(con_set)
        con_term =collect(linear_terms(constraint_object(constraint_by_name(m, i )).func))
        u[con_idx[i]] = con_rhs[i]
        l[con_idx[i]] = con_rhs[i]
        for j in 1:length(con_term)
            push!(I, con_idx[i])
            push!(J, var_idx[con_term[j][2]])
            push!(V, con_term[j][1])
        end
    end

    A = sparse(I,J,V)
end

###############################################################################################################################################################
println("================Running!!!!!=================")
val, t_run, bytes, gctime, memallocs = @timed begin
    variable = all_constraints(m, VariableRef, MathOptInterface.EqualTo{Float64})
    for i in variable
        delete(m,i)
    end
    variable = all_constraints(m, VariableRef, MathOptInterface.Interval{Float64})
    for i in variable
        delete(m,i)
    end

    println("================FP & CP=================")
    Initial_bound()
    try 
        UnSet_Type()
    catch  end

    optimize!(m)
    termination_status(m)==MOI.OPTIMAL
    solution_k, dict_xlb_k,dict_xub_k, dict_LP_int, solution_k_1=LP_solve(m)
    dict_infeasible_index, n_infeasible_var=Infeasible_Check(solution_k)
    Update(dict_LP_int)
    solution_k_check, dict_xlb_k_check,dict_xub_k_check,dict_LP_int_check, solution_k_1_check=LP_solve(m)
    dict_infeasible_index, n_infeasible_var=Infeasible_Check(solution_k)
    LPUpdate(dict_LP_int_check)
    solution_k, dict_xlb_k,dict_xub_k, dict_LP_int, solution_k_1=LP_solve(m)
    dict_infeasible_index, n_infeasible_var=Infeasible_Check(solution_k)
    CP(dict_infeasible_index)
    solution_k_check, dict_xlb_k_check,dict_xub_k_check,dict_LP_int_check, solution_k_1_check=LP_solve(m)
    dict_infeasible_index, n_infeasible_var=Infeasible_Check(solution_k)
    #Reupdate(dict_xlb_k_check, dict_xub_k_check)
    println("================MILP CBC=================")
    Set_Type()
    solution_k, dict_xlb_k,dict_xub_k, dict_LP_int, solution_k_1=LP_solve(m,var_lb, dict_xub_k_check)
end

###############################################################################################################################################################
#############################################################################
# result
#############################################################################
sol_check_bound = Array{Any}(undef,0)
sol_check_type = Array{Any}(undef,0)
for x in keys(solution_k)
    push!(sol_check_bound,var_lb[x]<=solution_k[x]<=var_ub[x])
end

for x in [k for (k,v) in var if v==:Bin]
    if isapprox(abs(solution_k[x] - round(solution_k[x])), 0, atol=1e-8) ==false  
        push!(sol_check_type, (x,solution_k[x], :Int))
    else
        continue
    end
end


for x in [k for (k,v) in var if v==:Int]
    if isapprox(abs(solution_k[x] - round(solution_k[x])), 0, atol=1e-8) ==false  
        push!(sol_check_type, (x,solution_k[x], :Int))
    else
        continue
    end
end

check_optimal = termination_status(m)!=MOI.OPTIMAL
check_bound = length(findall(x->false, sol_check_bound))
check_type = length(findall(x->false, sol_check_type))
t_total = t_problemLaod+t_run
println("The problem load time is :",t_problemLaod)
println("The running time is : ", t_run,"s")
println("The total time is : ", t_total)
println("Is it optimal?: ", check_optimal)
println("The the number of unsatisfied variable is:")
println("    - The the number of unsatisfied bound is:",check_bound )
println("    - The the number of unsatisfied type is:", check_type)

df=DataFrame(
    Name = [filename],
    Total = [t_total],
    Problem = [t_problemLaod],
    Run = [t_run],
    Optimal = [check_optimal],
    Check_Bound = [check_bound],
    check_type = [check_type]
)




"""
solution_k[collect(keys(dict_infeasible_index))[1]]
for i in keys(dict_infeasible_index)
    println(solution_k[i])
end

for x in [k for (k,v) in var if v==:Bin]
    if (solution_k[x] ==0 && solution_k[x] ==1) == true
        push!(sol_check_type, (x,solution_k[x], :Bin))
    else
        continue
    end
end
"""
