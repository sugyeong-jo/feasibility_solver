println("================패키지 부르기=================")
using ProgressBars
using MathOptInterface
using JuMP, Cbc
const MOI = MathOptInterface
using SparseArrays
using DataFrames
using CSV
using ArgParse
using QPSReader

include("SG_func_2.jl")
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

println("================MPS load=================")

parsed_args = parse_commandline()
filename=string(parsed_args["filename"])
Sec=parsed_args["Sec"]


#filename = string("R100701001_2_cplex.mps")
#Sec = 2000

#=
      File read
=#

println("The file name is: $filename ")
println("The limit time of solver is: $Sec s")
println("The max solution is: $MaxSol ")

filepath = string("/HDD/Workspace/CLT/mps/processing/CPLEX_file/",filename)
val, t_problemLaod, bytes, gctime, memallocs = @timed begin
    qps = readqps(filepath)
    qps.varnames
    qps.varindices
    qps.connames
    qps.contypes
    qps.conindices
    qps.uvar
    qps.lvar
    qps.ucon
    qps.lcon
    A = sparse(qps.arows, qps.acols, qps.avals, qps.ncon, qps.nvar)
end

println("========================================")
val, t_problemLaod, bytes, gctime, memallocs = @timed begin
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
    for var_name in (all_variables(m))
        var[var_name] = :Con
        var_idx[var_name] = var_name.index.value
        idx_var[var_name.index.value] = var_name
        var_lb[var_name] = -Inf
        var_ub[var_name] = Inf
        var_ref[var_name] = :except
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
    I = Int[]
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
            push!(I, con_idx[i])
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
            push!(I, con_idx[i])
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
            push!(I, con_idx[i])
            push!(J, var_idx[con_term[j][2]])
            push!(V, con_term[j][1])
        end
    end

    A = sparse(I,J,V)
end

###############################################################################################################################################################
#=
Solving!
=#
path = string("/home/sugyeong/HDD/sugyeong/Julia/feasibility_solver/result/",filename,".csv")
println("================Running!!!!!=================")
function SG()
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
    solution_k, dict_xlb_k,dict_xub_k, dict_LP_int, solution_k_1=LP_solve(m)
    dict_infeasible_index, n_infeasible_var=Infeasible_Check(solution_k)
    Update(dict_LP_int) 
    solution_k, dict_xlb_k,dict_xub_k, dict_LP_int, solution_k_1=LP_solve(m)
    dict_infeasible_index, n_infeasible_var=Infeasible_Check(solution_k)
    Update(dict_LP_int)    
    CP(dict_infeasible_index)
    solution_k_check, dict_xlb_k_check,dict_xub_k_check,dict_LP_int_check, solution_k_1_check=LP_solve(m)
    dict_infeasible_index_check, n_infeasible_var_check=Infeasible_Check(solution_k_check) 
 
    println("================MILP CBC=================")
    Set_Type()
    try
        solution_k, dict_xlb_k,dict_xub_k, dict_LP_int, solution_k_1=LP_solve(m, var_lb, dict_xub_k_check) 
        return solution_k

    catch  ex 
        println(ex)
        println("________________________Update Lower bound________________________")
        Initial_bound()
        try 
            UnSet_Type()
        catch  end
        
        optimize!(m)
        termination_status(m)==MOI.OPTIMAL
        solution_k, dict_xlb_k,dict_xub_k, dict_LP_int, solution_k_1=LP_solve(m)
        dict_infeasible_index, n_infeasible_var=Infeasible_Check(solution_k)
        Update_LB(dict_LP_int)
        CP(dict_infeasible_index)
        solution_k_check, dict_xlb_k_check,dict_xub_k_check,dict_LP_int_check, solution_k_1_check=LP_solve(m)
        dict_infeasible_index_check, n_infeasible_var_check=Infeasible_Check(solution_k_check) 
     
        println("================MILP CBC=================")
        Set_Type()
        solution_k, dict_xlb_k,dict_xub_k, dict_LP_int, solution_k_1=LP_solve(m, var_lb, dict_xub_k_check) 
        return solution_k

    end
end
val, t_run, bytes, gctime, memallocs = @timed begin
    solution_k = SG()

end

###############################################################################################################################################################
#############################################################################
# result
Result(path, solution_k)

#############################################################################
###########################################
##########################################