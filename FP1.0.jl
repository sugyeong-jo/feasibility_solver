println("================패키지 부르기=================")
using ProgressBars
using MathOptInterface
using JuMP, Cbc
const MOI = MathOptInterface
using SparseArrays
using DataFrames
using CSV
using ArgParse


include("FP1.0_func.jl")
include("utility.jl")

parsed_args = parse_commandline()
filename=string(parsed_args["filename"])
Sec=parsed_args["Sec"]

#filename = string("R100701005_2_cplex.mps")
#filename = string("atlanta-ip.mps")
#filename = string("b1c1s1.mps")

#Sec = 2000

#filepath = string("/HDD/Workspace/CLT/FP/data/",filename)
#filepath = string("/HDD/Workspace/CLT/mps/processing/CPLEX_file/",filename)

val, t_problemLaod, bytes, gctime, memallocs = @timed begin
    m, con_idx, idx_con, A, l, u, var_idx_qps, var, var_idx, idx_var, var_lb, var_ub, qps_var = MPS_read(filepath, Sec)
end

path = string("/home/sugyeong/HDD/sugyeong/Julia/feasibility_solver/result/FP1_$filename.csv")
println("================Running!!!!!=================")
function FP1(T, maxIter, P)
    variable = all_constraints(m, VariableRef, MathOptInterface.EqualTo{Float64})
    for i in variable
        delete(m,i)
    end
    variable = all_constraints(m, VariableRef, MathOptInterface.Interval{Float64})
    for i in variable
        delete(m,i)
    end
    println("Step 1| Initialized as a minimum-cost solution of the LP relaxation")
    Initial_bound()
    try
        UnSet_Type()
    catch  end
    println("----------------Pumping Cycle--------------------")
    solution_k, dict_xlb_k,dict_xub_k, dict_LP_int, solution_k_1=LP_solve(m)
    dict_infeasible_index, n_infeasible_var=Infeasible_Check(solution_k)
    Update_UL(dict_LP_int)
    top_score_list,top_idx ,dist = Projection(solution_k,solution_k_1)
    
    
    nIter = 0
    KK = Int[]
    while dist>0 && nIter <maxIter
        nIter = nIter+1
        try
            solution_k, dict_xlb_k, dict_xub_k,dict_LP_int, solution_k_1=LP_solve(m)
            dict_infeasible_index, n_infeasible_var=Infeasible_Check(solution_k)
            top_score_list,top_idx ,dist =Projection(solution_k,solution_k_1)

            if dist > 0
                KK =Int[]
                purturb = sum(KK)
                TT = Int(round(rand((T/2,3*T/2))))
                if TT>=length(top_idx)
                    TT = length(top_idx)
                end
                dict_x_idx = top_idx[1:TT]
                Update_UL(Round_change(dict_x_idx,solution_k,solution_k_1))
                
                println("===================================================================================================================")
                println("   ")
                println("[T, nIter, purturb]: [$T, $nIter, $purturb]| n_infeasible_var: $n_infeasible_var, |dist: $dist")
                println("   ")
                println("===================================================================================================================")   
            
            else
                println("===================================================================================================================")
                println("   ")
                println("[T, nIter, purturb]: [$T, $nIter, $purturb]| n_infeasible_var: $n_infeasible_var, |dist: $dist")
                println("   ")
                println("===================================================================================================================")   
                return solution_k, nIter  
                break
            end
        catch
            purturb = sum(push!(KK, 1))
            if purturb < P
                TT = Int(round(rand((T/2,3*T/2))))
                if TT>=length(top_idx)
                    TT = length(top_idx)
                end
                dict_x_idx = top_idx[1:TT]
                Update_UL(Round_change_pert(dict_x_idx,solution_k,solution_k_1))
                println("===================================================================================================================")
                println("   ")
                println("[T, nIter, perturb]: [$T, $nIter, $purturb]| n_infeasible_var: $n_infeasible_var, |dist: $dist")
                println("   ")
                println("===================================================================================================================")   
                continue
            else
                TT = Int(round(rand((T/2,3*T/2))))
                if TT>=length(top_idx)
                    TT = length(top_idx)
                end
                dict_x_idx = rand(top_idx, TT)
                Update_UL(Round_change_pert(dict_x_idx,solution_k,solution_k_1))
                Update_UL(Round_change_pert(dict_x_idx,solution_k,solution_k_1))
                println("===================================================================================================================")
                println("  Perturb ")
                println("[T, nIter, perturb]: [$T, $nIter, $purturb]| n_infeasible_var: $n_infeasible_var, |dist: $dist")
                println("   ")
                println("===================================================================================================================")   
                continue

            end
        end


    end    
    return solution_k, nIter 
    println("end")

end

val, t_run, bytes, gctime, memallocs = @timed begin
    solution_k, nIter = FP1(1000, 60, 30)
end

###############################################################################################################################################################
#############################################################################
# result
Result(path, solution_k, nIter)
