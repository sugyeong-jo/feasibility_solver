using JuMP  # Need to say it whenever we use JuMP
using Cbc # Loading the GLPK module for using its solver
using MathOptInterface
using Plots
using LinearAlgebra
######################################


include("Func_0624.jl")

#filename = string("atlanta-ip.mps")
#filename = string("timtab1.mps")
#filepath = string("/Users/sungilkim/Dropbox/Sungil_project/CyberLogitec/연장과제/data/",filename)
#filepath = string("C:\\Users\\sungi\\Dropbox\\Sungil_project\\CyberLogitec\\연장과제\\data\\",filename)


filename = string("b1c1s1.mps")
filepath = string("/HDD/Workspace/CLT/FP/data/",filename)


##

function Rounding(sol,itg_idx)
    sol_tilde = copy(sol)
    sol_tilde[itg_idx.==1] = round.(sol[itg_idx.==1],RoundNearest)
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



val, t_problemLaod, bytes, gctime, memallocs = @timed begin
   m = MPS_read(filepath, Sec)
end

#optimize!(m)
#println("relaxed objective value is = ",  JuMP.objective_value(m)) #JuMP.value.(keys(var)), " and ",


plot(ub, marker = :circle, legend = false, grid = false)
plot!(lb, marker = :circle, legend = false, grid = false)


intvar = all_constraints(m, VariableRef, MathOptInterface.Integer)
binvar = all_constraints(m, VariableRef, MathOptInterface.ZeroOne)

int_bin = deepcopy(vcat(intvar,binvar))


x = all_variables(m)
ub = upper_bound.(x)
lb = lower_bound.(x)
global itg_idx = is_binary.(x)+is_integer.(x)



delete(m,vcat(intvar,binvar))



optimize!(m)
println("relaxed objective value is = ",  JuMP.objective_value(m)) #JuMP.value.(keys(var)), " and ",


#y = collect(values(var))
#ind = findall(x -> x == :Int, y)
#value.(keys(var))
#intlropt = value.(collect(keys(var))[ind])
#vv = all_variables(m)


#global itg_idx = findall(xx->xx==true, itg)


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


    x = all_variables(m2)
    ub = upper_bound.(x)
    lb = lower_bound.(x)
    #itg_idx = is_binary.(x)+is_integer.(x)

    d_idx, u_idx, l_idx = objective_set(sol_tilde,ub,lb,itg_idx)

    #plot(is_binary.(x),marker=:circle)
    #plot!(is_integer.(x),marker=:circle)

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
println("x = ", sol)
println("dist = ", dist)
