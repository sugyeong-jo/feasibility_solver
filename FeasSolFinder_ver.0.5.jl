println("================패키지 부르기=================")
using ProgressBars
using MathOptInterface
using JuMP, Cbc
const MOI = MathOptInterface
using Cbc
using MathProgBase
using GLPKMathProgInterface
using Random
import JuMP
using MathOptInterface
using DataFrames
using CSV
using ArgParse

include("FeasSolFinder_func_ver.0.5.jl")

println("================문제 생성=================")
####################################
########  MPS problem read
####################################
s = ArgParseSettings()
s.description = "Finde Feasible Solution"
s.commands_are_required = true
s.version = "1.0"
s.add_version = true

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "filename"
            help = "a positional argument"
            required = true
    end

    return parse_args(s)
end


parsed_args = parse_commandline()
filename=string(parsed_args["filename"])
println("The file name is: $filename ")


###########################################
#filename = "R100701005_2_cplex.mps"
filepath = string("/HDD/Workspace/CLT/mps/processing/CPLEX_file/",filename)
m = MathProgBase.LinearQuadraticModel(GLPKSolverMIP())
MathProgBase.loadproblem!(m,filepath)
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
val, t_problemLaod, bytes, gctime, memallocs = @timed begin
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
        
end



############################################################
# 흠.. ? 이렇게만 해도 되네?
println("================Running!!!!!=================")
val, t_run, bytes, gctime, memallocs = @timed begin
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
            

    Reupdate(dict_xlb, dict_xub_k_check)
    Set_Type()
    optimize!(m)
    termination_status(m)!=MOI.OPTIMAL

end
#############################################################################
# result
#############################################################################
sol_check_bound = Array{Any}(undef,0)
sol_check_type = Array{Any}(undef,0)
for x in keys(solution_k)
    push!(sol_check_bound,dict_xlb[x]<=solution_k[x]<=dict_xub[x])
end

for x in keys(dict_bin_type)
    if (solution_k[x] ==0 && solution_k[x] ==1) == true
        push!(sol_check_type, (x,solution_k[x], :Bin))
    else
        continue
    end
end

for x in keys(dict_int_type)
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

path = string("/home/sugyeong/HDD/sugyeong/Julia/feasibility_solver/result/",filename,".csv")
CSV.write(path, df) 