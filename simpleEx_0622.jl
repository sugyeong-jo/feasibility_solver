using JuMP  # Need to say it whenever we use JuMP
using Cbc # Loading the GLPK module for using its solver
using MathOptInterface



include("utility_FP.jl")

filename = string("atlanta-ip.mps")
filepath = string("/HDD/Workspace/CLT/FP/data/",filename)

Sec = 2000
#=   
lb/ub = variable lower/upper bound
itg = integer type (t/f)
l/u = constraint lower/upper bound
A = constraint matrix
=# 
val, t_problemLaod, bytes, gctime, memallocs = @timed begin
    lb, ub, itg, l, u, A, m, con_idx, idx_con, var_idx_qps, var, var_idx, idx_var, var_lb, var_ub, qps_var, qps = MPS_read_FP(filepath, Sec) 
end


ub = [10,10,1]
lb = [0,0,0]
itg = [false,true,false]

m2 = Model()
set_optimizer(m2, Cbc.Optimizer)
@variable(m2, x[i=1:3], upper_bound=ub[i],lower_bound=lb[i],integer=itg[i])
@objective(m2, Max, [1,2,5]'*x)
@constraint(m2, constraint1, [-1,1,3]'*x <= -5)
@constraint(m2, constraint2, [1,3,-7]'*x <= 10)
optimize!(m2)
println("true MIP solution x = ", JuMP.value.(x))

intvar = all_constraints(m2, VariableRef, MathOptInterface.Integer)
delete(m2,intvar)
optimize!(m2)
println("x = ", JuMP.value.(x))


#y = collect(values(var))
#ind = findall(x -> x == :Int, y)
#value.(keys(var))
#intlropt = value.(collect(keys(var))[ind])
#vv = all_variables(m)

global sol = JuMP.value.(x)
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



j=1

while dist>0 && j<10
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
    m2 = Model()
    set_optimizer(m2, Cbc.Optimizer)

    @variable(m2, x[i=1:3], upper_bound=ub[i],lower_bound=lb[i],integer=itg[i])
    @variable(m2, d[d_idx]>=0)

    @objective(m2,Min,c'*x+ones(length(d))'*d)

    @constraint(m2, constraint1, [-1,1,3]'*x <= -5)
    @constraint(m2, constraint2, [1,3,-7]'*x <= 10)

    @constraint(m2, con1[i = d_idx], x[i]- sol_tilde[i] <= d[i])
    @constraint(m2, con2[i = d_idx], sol_tilde[i] - x[i] <= d[i])


    intvar = all_constraints(m2, VariableRef, MathOptInterface.Integer)
    delete(m2,intvar)
    optimize!(m2)


    ########

    global sol = JuMP.value.(x)
    global sol_tilde = copy(sol)
    for i in 1:length(sol)
        if itg[i] == true
            sol_tilde[i] = round(sol_tilde[i],RoundNearest)
        end
    end


    global dist = distance_fucntion(sol,sol_tilde)

    #println(dist)
    #print(j)
    global j += 1

end


# Printing the optimal solutions obtained
println("Feasible MIP Solutions:")
println("x = ", sol)


var

