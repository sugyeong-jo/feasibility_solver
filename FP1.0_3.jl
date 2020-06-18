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

# for CMD
parsed_args = parse_commandline()
filename=string(parsed_args["filename"])
Sec=parsed_args["Sec"]

# for Julia
#filename = string("R100701005_2_cplex.mps")
#filename = string("atlanta-ip.mps")
#filename = string("blp-ic98.mps")
#filename = string("blp-ar98.mps")
#filename = string("b1c1s1.mps")
#filename = string("1_R170718004_veri_cplex.mps")


#Sec = 2000

filepath = string("/HDD/Workspace/CLT/FP/data/",filename)
#filepath = string("/HDD/Workspace/CLT/mps/processing/CPLEX_file/",filename)

val, t_problemLaod, bytes, gctime, memallocs = @timed begin
    m, con_idx, idx_con, A, l, u, var_idx_qps, var, var_idx, idx_var, var_lb, var_ub, qps_var = MPS_read(filepath, Sec)
end

dist_LP=copy(m)
objective_function(dist_LP,AffExpr)
dist_LP= Model()
variable_list = string(collect(keys(solution_k_1)))
for i in keys(solution_k_1)
    @variable(dist_LP, keys(solution_k_1))
end

@objective(dist_LP, Min, sum( x[i] for i in 1:length(solution_k_1) ))

collect(keys(solution_k_1))[1]
#[k for (k,v) in var if v==:Int]

path = string("/home/sugyeong/HDD/sugyeong/Julia/feasibility_solver/result/FP1_$filename.csv")
println("================Running!!!!!=================")
function FP1(T, maxIter, P) #P == kk, 70
    maxIter_bin = Int(round(maxIter*0.5))
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
    dict_infeasible_index, n_infeasible_var=Infeasible_Check(solution_k) #
    #Update_UL(dict_LP_int)
    #Update_UL(solution_k_1)
    top_score_list,top_idx ,dist = Projection(solution_k,solution_k_1)
    top_idx_bin = Any[] 
    top_idx_int = Any[]
    for i in top_idx
        if var[i] == :Bin
            push!(top_idx_bin, i)
        else
            push!(top_idx_int, i)
        end
    end    

    nIter = 0
    KK = Int[]
    result = "Infeas"
    # Binary Step
    n_Bin = length(top_idx_bin)
    while n_Bin >0 && nIter < maxIter
        nIter = nIter+1
        #Binary
        if nIter <= maxIter_bin
            try 
                purturb = sum(KK)
                solution_k, dict_xlb_k, dict_xub_k,dict_LP_int, solution_k_1=LP_solve(m)
                dict_infeasible_index, n_infeasible_var=Infeasible_Check(solution_k)
                top_score_list,top_idx ,dist =Projection(solution_k,solution_k_1)
                #####
                #if dist < dist_Best
                solution_k_Best, dict_xlb_k_Best, dict_xub_k_Best,dict_LP_int_Best, solution_k_1_Best = solution_k, dict_xlb_k, dict_xub_k,dict_LP_int, solution_k_1
                global solution_k_Best, dict_xlb_k_Best, dict_xub_k_Best,dict_LP_int_Best, solution_k_1_Best
                KK = Int[]
                top_idx_bin = Any[]
                top_idx_int = Any[]
                for i in top_idx
                    if var[i] == :Bin
                        push!(top_idx_bin, i)
                    else
                        push!(top_idx_int, i)
                    end
                end

                n_Bin = length(top_idx_bin)   
                if n_Bin == 0
                    break
                else
                    if n_Bin >0 &&   nIter <maxIter
                        TT = Int(round(rand((T/2,3*T/2))))
                        if TT>=n_Bin
                            TT = length(top_idx_bin)
                        end
                        dict_x_idx = top_idx_bin[1:TT]
                        #Update_UL(dict_LP_int)
                        #Update_UL(solution_k_1)
                        Update_UL(Round_change(dict_x_idx,solution_k,solution_k_1))
                        solution_k, dict_xlb_k, dict_xub_k,dict_LP_int, solution_k_1=LP_solve(m)


                        println("===================================================================================================================")
                        println("   ")
                        println("[T, nIter, purturb]: [$T, $nIter, $purturb]| n_infeasible_var: $n_infeasible_var, |dist: $dist")
                        println("   ")
                        println("===================================================================================================================")
                        result = "InFeas"   
                    else
                        println("===================================================================================================================")
                        println("Feas!")
                        println("[T, nIter, purturb]: [$T, $nIter, $purturb]| n_infeasible_var: $n_infeasible_var, |dist: $dist")
                        println("   ")
                        println("===================================================================================================================")   
                        result = "Feas"
                        break
                    end
                end
            catch
                purturb = sum(push!(KK, 1))
                if purturb < P
                    TT = Int(round(rand((T/2,3*T/2))))
                    if TT>=length(top_idx_bin)
                        T = length(top_idx_bin)
                        TT = Int(round(rand((T/2,3*T/2))))
                    end
                    dict_x_idx = top_idx_bin[1:TT]
                    Update_UL(Round_change_pert(dict_x_idx,solution_k,solution_k_1))
                    println("===================================================================================================================")
                    println("   ")
                    println("[T, nIter, perturb]: [$T, $nIter, $purturb]| n_infeasible_var: $n_infeasible_var, |dist: $dist")
                    println("   ")
                    println("===================================================================================================================")   
                    result = "Infeas"
                    continue
                else
                    TT = Int(round(rand((T/2,3*T/2))))
                    if TT>=length(top_idx_bin)
                        T = length(top_idx_bin)
                        TT = Int(round(rand((T/2,3*T/2))))
                    end
                    dict_x_idx = rand(top_idx_bin, TT)
                    Update_UL(Round_change_pert(dict_x_idx,solution_k,solution_k_1))
                    Update_UL(Round_change_pert(dict_x_idx,solution_k,solution_k_1))
                    println("===================================================================================================================")
                    println("  Perturb ")
                    println("[T, nIter, perturb]: [$T, $nIter, $purturb]| n_infeasible_var: $n_infeasible_var, |dist: $dist")
                    println("   ")
                    println("===================================================================================================================")   
                    result = "Infeas"
                    continue
                end
            end    
        else
            #Integer
            println(nIter)
            if nIter == (maxIter_bin+1)
                println("reset!!!")
                solution_k, dict_xlb_k, dict_xub_k,dict_LP_int, solution_k_1=LP_solve(m, dict_xlb_k_Best, dict_xub_k_Best)
                Reupdate(dict_xlb_k_Best, dict_xub_k_Best)
                KK = Int[]
            end

            try 
                println("Integer part")
                purturb = sum(KK)
                solution_k, dict_xlb_k, dict_xub_k,dict_LP_int, solution_k_1=LP_solve(m)
                dict_infeasible_index, n_infeasible_var=Infeasible_Check(solution_k)
                top_score_list,top_idx ,dist =Projection(solution_k,solution_k_1)
                solution_k_Best, dict_xlb_k_Best, dict_xub_k_Best,dict_LP_int_Best, solution_k_1_Best = solution_k, dict_xlb_k, dict_xub_k,dict_LP_int, solution_k_1
                global solution_k_Best, dict_xlb_k_Best, dict_xub_k_Best,dict_LP_int_Best, solution_k_1_Best
                Reupdate(dict_xlb_k_Best, dict_xub_k_Best)
                KK = Int[]
                top_idx_bin = Any[]
                top_idx_int = Any[]
                for i in top_idx
                    if var[i] == :Bin
                        push!(top_idx_bin, i)
                    else
                        push!(top_idx_int, i)
                    end
                end
                n_Int = length(top_idx_int)
                println(n_Int)
                if n_Int == 0
                    break
                else

                    if n_Int >0  &&   nIter <maxIter
                        TT = Int(round(rand((T/2,3*T/2))))
                        if TT>=n_Int
                            T = length(top_idx_int)
                            TT = Int(round(rand((T/2,3*T/2))))
                        end
                        dict_x_idx = top_idx_int[1:TT]
                        Update_UL(Round_change(dict_x_idx,solution_k,solution_k_1))
                        println("===================================================================================================================")
                        println(" Int  ")
                        println("[T, nIter, purturb]: [$T, $nIter, $purturb]| n_infeasible_var: $n_infeasible_var, |dist: $dist")
                        println("   ")
                        println("===================================================================================================================")   
                        result = "Infeas"
                        continue
        
                    else
                        println("===================================================================================================================")
                        println(" InT_feas  ")
                        println("[T, nIter, purturb]: [$T, $nIter, $purturb]| n_infeasible_var: $n_infeasible_var, |dist: $dist")
                        println("   ")
                        println("===================================================================================================================")   
                        result = "Feas"
        
                        break
                    end
                end
            catch
                print("====================")
                
                if n_Int == 0
                    break
                else
                    purturb = sum(push!(KK, 1))
                    if purturb < P
                        TT = Int(round(rand((T/2,3*T/2))))
                        if TT>=length(top_idx_int)
                            T = length(top_idx_int)
                            TT = Int(round(rand((T/2,3*T/2))))
                        end
                        dict_x_idx = top_idx_int[1:TT]
                        Update_UL(Round_change_pert(dict_x_idx,solution_k,solution_k_1))
                        println("===================================================================================================================")
                        println(" Int_2  ")
                        println("[T, nIter, perturb]: [$T, $nIter, $purturb]| n_infeasible_var: $n_infeasible_var, |dist: $dist")
                        println("   ")
                        println("===================================================================================================================")   
                        result = "Infeas"
        
                        continue
                    else
                        TT = Int(round(rand((T/2,3*T/2))))
                        if TT>=length(top_idx_int)
                            T = length(top_idx_int)
                            TT = Int(round(rand((T/2,3*T/2))))
                        end
                        dict_x_idx = rand(top_idx_int, TT)
                        Update_UL(Round_change_pert(dict_x_idx,solution_k,solution_k_1))
                        println("===================================================================================================================")
                        println("  Int_Perturb ")
                        println("[T, nIter, perturb]: [$T, $nIter, $purturb]| n_infeasible_var: $n_infeasible_var, |dist: $dist")
                        println("   ")
                        println("===================================================================================================================")   
                        result = "Infeas"
                        continue
                    end
                end
            end
        end
    end
    println("===================================================================================================================")
    println("  End ")
    println("[T, nIter, perturb]: [$T, $nIter]| n_infeasible_var: $n_infeasible_var, |dist: $dist")
    println("   ")
    println("===================================================================================================================")   
    return solution_k, nIter, result
    println("end")
end

val, t_run, bytes, gctime, memallocs = @timed begin
    T = 20
    maxIter = 30
    P = 10
    NIter = 0
    NmaxIter = 5
    result = "Infeas"
    while result == "Infeas" &&  NIter < NmaxIter
        global  NIter, result, maxIter, T, P        
        NIter = NIter+1
        solution_k, nIter, result = FP1(T, maxIter, P)
        println("Restart! $NIter")
        global  NIter, result, nIter        
    end

end

###############################################################################################################################################################
#############################################################################
# result
Result(path, solution_k, nIter*NIter)

