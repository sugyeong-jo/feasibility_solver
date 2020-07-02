println("================패키지 부르기=================")
using ProgressBars
using MathOptInterface
using JuMP, Cbc
const MOI = MathOptInterface
using SparseArrays
using DataFrames
using CSV
using ArgParse
using JuMP  # Need to say it whenever we use JuMP
using Cbc # Loading the GLPK module for using its solver
using MathOptInterface
#using Plots
using LinearAlgebra


include("Func_0625.jl")

filename = string("atlanta-ip.mps")
filename = string("b1c1s1.mps")
#filename = string("timtab1.mps")
#filepath = string("/Users/sungilkim/Dropbox/Sungil_project/CyberLogitec/연장과제/data/",filename)
#filepath = string("C:\\Users\\sungi\\Dropbox\\Sungil_project\\CyberLogitec\\연장과제\\data\\",filename)


#filename = string("b1c1s1.mps")
Sec = 20
filepath = string("/HDD/Workspace/CLT/FP/data/",filename)
val, t_problemLaod, bytes, gctime, memallocs = @timed begin
    m, con_idx, idx_con, A, l, u, var, var_idx, idx_var, var_lb, var_ub,coeff, itg_idx = MPS_read_full(filepath, Sec)
end


"""

#BoundStrength(1281)
Initial_bound()
try
    UnSet_Type()
catch
end

solution_k, dict_xlb_k, dict_xub_k, dict_LP_int, solution_k_1 = LP_solve(m)

fixed_var_element = rand(collect(keys(dict_LP_int)),1)[1]

filter((k,v) -> k == fixed_var_element,solution_k) # fixed variable 있는 solution dictionary

Update_UL(filter((k,v) -> k == fixed_var_element,solution_k))


##############################################
# 이밑을 반복하면 fiexed되는 variable을 구할 수 있습니다.
for i in tqdm(1:length(keys(idx_con))) #시간이 다소 소요됩니다.
    BoundStrength(i)
end

fixed_var = Any[]
non_fixed_var = Any[]
for i in tqdm([k for (k,v) in var if v==:Bin])
    if isapprox(abs(lower_bound(i)-upper_bound(i)),0,  atol=1e-3)
        push!(fixed_var,i)
    else
        push!(non_fixed_var, i)
    end
end


# fixed 된 variable을 제외한 다른 LP로 나온 variable 중 integer로 나온 값을 하나 fixed 시키면서 업데이트하면 늘어날 것 같습니다.
# 하지만 논문의 방법이 아니므로 논문을 더 읽어보려 합니다. 



"""
##########################################################################


function distance_fucntion(s1,s2)
    return sum(abs.(s1-s2)) #or norm(s1-s2,1)
end

function objective_set(sol_tilde,ub,lb,itg_idx)
    d_idx = []
    u_idx = []
    l_idx = []
    #c=zeros(length(sol_tilde))
    for i in 1:length(sol_tilde)
        if itg_idx[i] == 1
            if ub[i] == sol_tilde[i]
                append!(u_idx,i)
            elseif lb[i] == sol_tilde[i]
                append!(l_idx,i)
            else
                append!(d_idx,i)
            end
        end
    end
    return d_idx, u_idx, l_idx
end

function Rounding(sol,itg_idx)
    sol_tilde = copy(sol)
    sol_tilde[itg_idx.==1] = round.(sol[itg_idx.==1],RoundNearest)
    return sol_tilde
end


"""
x = all_variables(m)
ub = upper_bound.(x)
lb = lower_bound.(x)
global itg_idx = is_binary.(x)+is_integer.(x)
"""
Initial_bound()
try
    UnSet_Type()
catch
end

optimize!(m)
termination_status(m)

sol = Any[]
for i in all_variables(m)
    push!(sol,JuMP.value(i))
end   

global sol
#global sol = JuMP.value.(x)
global sol_tilde = Rounding(sol,itg_idx)

dist = distance_fucntion(sol,sol_tilde)
#plot(sol-sol_tilde, marker = :circle, legend = false, grid = false)

global nIter=0
maxIter = 20
#global alpha = 1

while dist>0 && nIter<maxIter

    global nIter += 1
    global sol
    global sol_tilde

    ## model copy
    m2 = copy(m)
    set_optimizer(m2, Cbc.Optimizer)
    set_optimizer_attribute(m2, MOI.Silent(), true)


    x = all_variables(m2)
    ub = upper_bound.(x)
    lb = lower_bound.(x)
    #itg_idx = is_binary.(x)+is_integer.(x)

    d_idx, u_idx, l_idx = objective_set(sol_tilde,ub,lb,itg_idx)

    #plot(is_binary.(x),marker=:circle)
    #plot!(is_integer.(x),marker=:circle)

    alpha = rand(Float64, 1)[1]
    global alpha
    dist = distance_fucntion(sol,sol_tilde)

    @variable(m2, d[d_idx]>=0)
    #@objective(m2,Min,c'*x+ones(length(d))'*d)
    #@objective(m2,Min,(1-alpha)/norm(dist,1)*(sum(ub[u_idx]-x[u_idx])+sum(x[l_idx]-lb[l_idx])+sum(d))+(alpha/norm(coeff,1))*sum(coeff'*x))
    @objective(m2,Min,sum(ub[u_idx]-x[u_idx])+sum(x[l_idx]-lb[l_idx])+sum(d))

    @constraint(m2, con1[i = d_idx], x[i]- sol_tilde[i] <= d[i])
    @constraint(m2, con2[i = d_idx], sol_tilde[i] - x[i] <= d[i])

    optimize!(m2)

    #all_constraints(m2, VariableRef, MathOptInterface.LessThan{Float64})
    #constraint_by_name(m2, "con1")

    #objective_function(m2)
    sol = JuMP.value.(x) #keys(var))
    #sol = JuMP.value.(keys(var))
    #sold = JuMP.value.(d)
    #objective_value(m2)
    #c'*sol1+ones(length(sold))'*sold

    global dist = objective_value(m2)#distance_fucntion(sol,sol_tilde)


    if dist>0#isapprox(dist, 0, atol = 1e-3) == false
        global itg_idx

        if sum(round.(sol[itg_idx.==1],RoundNearest)-sol_tilde[itg_idx.==1])>0
            #update
            sol_tilde = Rounding(sol,itg_idx)

            #dist = distance_fucntion(sol,sol_tilde)
        else
            dictionary1 = Dict(zip(Vector(1:length(itg_idx))[itg_idx.==1],(sol-sol_tilde)[itg_idx.==1]))
            sortedDic = sort(collect(dictionary1), by=x->abs(x[2]),rev=true);

            T = 10
            TT = floor(Int, rand(T/2:3T/2,1)[1])

            for j in 1:TT
                topTT = sortedDic[j]
                if topTT[2,]>0
                    sol_tilde[topTT[1,]] += 1
                else
                    sol_tilde[topTT[1,]] += -1
                end
            end
        end

    end

    println("dist = ", dist)#dist)
end


# Printing the optimal solutions obtained
println("Feasible MIP Solutions:")
#println("x = ", sol)
println("dist = ", dist)


x = all_variables(m)
solution_k = Dict()
for i in 1:length(var)
    solution_k[x[i]] = sol[i]
end
sol_check_bound = Array{Any}(undef,0)
sol_check_type = Array{Any}(undef,0)

for i in 1:length(var)    
    push!(sol_check_bound,var_lb[x[i]]<=sol[i]<=var_ub[x[i]])
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
#t_total = t_problemLaod+t_run
#println("The problem load time is : $t_problemLaod")
#println("The running time is : $t_run s")
#println("The total time is : $t_total")
println("Is it optimal?: $check_optimal")
println("The FP Iteration is : $nIter")
println("The objective value is : $obj_value" )
println("The the number of unsatisfied variable is:")
println("    - The the number of unsatisfied bound is: $check_bound" )
println("    - The the number of unsatisfied type is: $check_type")
