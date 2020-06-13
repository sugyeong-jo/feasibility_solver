println("================패키지 부르기=================")
using ProgressBars
using MathOptInterface
using JuMP, Cbc
const MOI = MathOptInterface
using SparseArrays
using DataFrames
using CSV
using ArgParse


include("SG_fin_func.jl")
include("utility.jl")

println("================MPS load=================")
parsed_args = parse_commandline()
filename=string(parsed_args["filename"])
Sec=parsed_args["Sec"]

#filename = string("R100701006_2_cplex.mps")
#Sec = 2000

println("The file name is: $filename ")
println("The limit time of solver is: $Sec s")

filepath = string("/HDD/Workspace/CLT/mps/processing/CPLEX_file/",filename)

println("========================================")
val, t_problemLaod, bytes, gctime, memallocs = @timed begin
    m, con_idx, idx_con, A, l, u, var_idx_qps, var, var_idx, idx_var, var_lb, var_ub, qps_var = MPS_read(filepath, Sec)
end


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
    solution_k, dict_xlb_k,dict_xub_k, dict_LP_int, solution_k_1=LP_solve(m)
    dict_infeasible_index, n_infeasible_var=Infeasible_Check(solution_k)
    Update_UB(dict_LP_int)
    CP(dict_infeasible_index)
    dict_xlb_k = Dict() # lower bound dictionary
    dict_xub_k = Dict() # lower bound dictionary
    for i in keys(var)
        dict_xlb_k[i] = lower_bound(i)
        dict_xub_k[i] = upper_bound(i)
    end
    dict_xub_k_check = dict_xub_k
    dict_xlb_k_check = dict_xlb_k
    
    println("================MILP CBC=================")
    Set_Type()
    try
        solution_k, dict_xlb_k,dict_xub_k, dict_LP_int, solution_k_1=LP_solve(m, var_lb, dict_xub_k_check)
        return solution_k
    catch  ex
        println(ex)
        println("________________________Update Lower bound________________________")
        solution_k, dict_xlb_k,dict_xub_k, dict_LP_int, solution_k_1=LP_solve(m, dict_xlb_k_check, var_ub)
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
