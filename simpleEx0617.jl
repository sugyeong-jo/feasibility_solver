using JuMP  # Need to say it whenever we use JuMP
using Cbc # Loading the GLPK module for using its solver
using MathOptInterface


########################################

ub = [10,10,1,20, 30]
lb = [0,0,0,-Inf, 0]
itg = [false,true,true, true, true]

m_org = Model()
set_optimizer(m_org, Cbc.Optimizer)
set_optimizer_attribute(m_org, MOI.Silent(), true)
@variable(m_org, x[i=1:5], upper_bound=ub[i],lower_bound=lb[i],integer=itg[i])
@objective(m_org, Max, [1,2,5, 10, -1]'*x)
@constraint(m_org, constraint1, [1,1,3,2,-5]'*x <= -5)
@constraint(m_org, constraint2, [1,-2,-7,5,1]'*x <= 10)
@constraint(m_org, constraint3, [2,5,8,5,1]'*x >= 4)
optimize!(m_org)

println("true MIP solution x = ", JuMP.value.(x), " and ", JuMP.objective_value(m_org))


include("utility_FP.jl")

#filename = string("atlanta-ip.mps")
filename = string("b1c1s1.mps")
filepath = string("/HDD/Workspace/CLT/FP/data/",filename)

Sec = 2000
#=   
lb/ub = variable lower/upper bound
itg = integer type (t/f)
l/u = constraint lower/upper bound
A = constraint matrix
=# 
val, t_problemLaod, bytes, gctime, memallocs = @timed begin
    lb, ub, itg, l, u, A, m_org, con_idx, idx_con, var_idx_qps, var, var_idx, idx_var, var_lb, var_ub, qps_var, qps = MPS_read_FP(filepath, Sec) 
end


intvar = all_constraints(m_org, VariableRef, MathOptInterface.Integer)
delete(m_org,intvar)

intvar = all_constraints(m_org, VariableRef, MathOptInterface.ZeroOne)
delete(m_org,intvar)

optimize!(m_org)
println("relaxed solution is x = ", JuMP.value.(x), " and ", JuMP.objective_value(m_org))


#y = collect(values(var))
#ind = findall(x -> x == :Int, y)
#value.(keys(var))
#intlropt = value.(collect(keys(var))[ind])
#vv = all_variables(m)


sol = Any[]
for i in keys(var)
    push!(sol,JuMP.value(i))
end   
global sol
global sol_tilde = copy(sol)
for i in 1:length(sol)
    if itg[i] == true
        sol_tilde[i] = round(sol_tilde[i],RoundNearest)
    end
end

#println(sol_tilde)
#println(sol)


function distance_fucntion(s1,s2)
    return sum(abs.(s1-s2))
end

dist = distance_fucntion(sol,sol_tilde)



nIter=0
maxIter = 100
while distance_fucntion(sol,sol_tilde)>0 && nIter<maxIter
    m2 = copy(m_org)

    set_optimizer(m2, Cbc.Optimizer)
    set_optimizer_attribute(m2, MOI.Silent(), true)



    global nIter += 1

    d_idx = []
    c=zeros(length(sol))
    for i in 1:length(sol_tilde)
        if ub[i] == sol_tilde[i]
            c[i] = -1
        elseif lb[i] == sol_tilde[i]
            c[i] = 1
        else
            append!(d_idx,i)
        end
    end


    #####

    
    x= all_variables(m2)

    @variable(m2, d[d_idx]>=0)
    @objective(m2,Min,c'*x+ones(length(d))'*d)

    @constraint(m2, con1[i = d_idx], x[i]- sol_tilde[i] <= d[i])
    @constraint(m2, con2[i = d_idx], sol_tilde[i] - x[i] <= d[i])


    intvar = all_constraints(m2, VariableRef, MathOptInterface.Integer)
    delete(m2,intvar)
    optimize!(m2)

    sol = Any[]
    for i in keys(var)
        push!(sol,JuMP.value(i))
    end   
    global sol

    global dist = distance_fucntion(sol,sol_tilde)

    if isapprox(dist, 0, atol = 1e-3) == false
        itg_idx = findall(x->x==true, itg)
        

        if sum(round.(sol[itg_idx],RoundNearest)-sol_tilde[itg_idx])>0

            #update
            global sol_tilde = copy(sol)
            for i in 1:length(sol)
                if itg[i] == true
                    sol_tilde[i] = round(sol[i],RoundNearest)
                end
            end
            global dist = distance_fucntion(sol,sol_tilde)
        else

            dictionary1 = Dict(zip(itg_idx,(sol-sol_tilde)[itg_idx]))
            TT = 2
            sortedDic = sort(collect(dictionary1), by=x->abs(x[2]),rev=true)

            for j in 1:TT
                topTT = sortedDic[j]
                if topTT[2,]>0
                    sol_tilde[topTT[1,]] += 1
                else sol_tilde[topTT[1,]] += -1
                end
            end
        end
    else
        break
    end


    #println(dist)
    println(dist)
    println("Feasible MIP Solutions:")
    println("x = ", sol)


end


# Printing the optimal solutions obtained
println("Feasible MIP Solutions:")
println("x = ", sol)
