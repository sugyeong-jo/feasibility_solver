using JuMP  # Need to say it whenever we use JuMP
using Cbc # Loading the GLPK module for using its solver  
using MathOptInterface
#using Plots
using LinearAlgebra
include("utility_FP.jl")
include("Func_0622.jl")



function Rounding(sol,itg_idx)
    sol_tilde = copy(sol)
    sol_tilde[itg_idx] = round.(sol[itg_idx],RoundNearest)
    return sol_tilde
end


function distance_fucntion(s1,s2)
    return sum(abs.(s1-s2)) #or norm(s1-s2,1)
end

function objective_set(sol_tilde,ub,lb,itg_idx)
    d_idx = []
    u_idx = []
    l_idx = []
    c=zeros(length(sol_tilde))
    for i in 1:length(sol_tilde)
        if i in itg_idx
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

function Perturb(sol0, sol_tilde0,TT,itg_idx0)
    dictionary1 = Dict(zip(itg_idx0,(sol0-sol_tilde0)[itg_idx0]))
    sortedDic = sort(collect(dictionary1), by=x->abs(x[2]),rev=true);
    println((sol0-sol_tilde0)[1565])
    println(minimum((sol0-sol_tilde0)[itg_idx0]))
    println(sortedDic[1:3])

    for j in 1:TT
        #println(j)
        topTT = sortedDic[j]

        if topTT[2,]>0
            sol_tilde0[topTT[1,]] += 1
            #println(sol_tilde0[topTT[1,]])
        else
            sol_tilde0[topTT[1,]] += -1
            #println(sol_tilde0[topTT[1,]])
        end
    end
    #output = copy(sol_tilde0)
    #sol_tilde0 = 0
    return sol_tilde0
end


##

Sec = 2000


filename = string("b1c1s1.mps")
filename = string("atlanta-ip.mps")
filepath = string("/HDD/Workspace/CLT/FP/data/",filename)

val, t_problemLaod, bytes, gctime, memallocs = @timed begin
    lb, ub, itg, l, u, A, m, var_name, con_idx, idx_con, var_idx_qps, var, var_idx, idx_var, var_lb, var_ub, qps_var, qps = MPS_read_FP(filepath, Sec)
end

#optimize!(m)
println("relaxed objective value is = ",  JuMP.objective_value(m)) #JuMP.value.(keys(var)), " and ",


plot(ub, marker = :circle, legend = false, grid = false)
plot!(lb, marker = :circle, legend = false, grid = false)


intvar = all_constraints(m, VariableRef, MathOptInterface.Integer)
binvar = all_constraints(m, VariableRef, MathOptInterface.ZeroOne)

int_bin = deepcopy(vcat(intvar,binvar))
delete(m,vcat(intvar,binvar))



optimize!(m)
println("relaxed objective value is = ",  JuMP.objective_value(m)) #JuMP.value.(keys(var)), " and ",


#y = collect(values(var))
#ind = findall(x -> x == :Int, y)
#value.(keys(var))
#intlropt = value.(collect(keys(var))[ind])
#vv = all_variables(m)


global itg_idx = findall(xx->xx==true, itg)
global var_name

x= var_name #all_variables(m_org)


global sol = JuMP.value.(x)
global sol_tilde = Rounding(sol,itg_idx)

dist = distance_fucntion(sol,sol_tilde)
plot(sol-sol_tilde, marker = :circle, legend = false, grid = false)

global nIter=0
maxIter = 100
#global alpha = 1

while dist>0 && nIter<maxIter

    global nIter += 1
    global sol
    global sol_tilde
    global alpha
    ## model copy
    m2 = copy(m)
    set_optimizer(m2, Cbc.Optimizer)
    set_optimizer_attribute(m2, MOI.Silent(), true)

    d_idx, u_idx, l_idx = objective_set(sol_tilde,ub,lb,itg_idx)

    x = all_variables(m)
    x = var_name


    @variable(m2, d[d_idx]>=0)
    #@objective(m2,Min,c'*x+ones(length(d))'*d)
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

        if sum(round.(sol[itg_idx],RoundNearest)-sol_tilde[itg_idx])>0
            #update
            sol_tilde = Rounding(sol,itg_idx)

            #dist = distance_fucntion(sol,sol_tilde)
        else
            dictionary1 = Dict(zip(itg_idx,(sol-sol_tilde)[itg_idx]))
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
println("x = ", sol)
println("dist = ", dist)
