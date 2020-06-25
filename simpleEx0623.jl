using JuMP  # Need to say it whenever we use JuMP
using Cbc # Loading the GLPK module for using its solver
using MathOptInterface
include("utility_FP.jl")
include("Func_0622.jl")


########################################
"""
# Toy problem  

ub = [10,10,1,20, 30]
lb = [0,0,0,-Inf, 0]
itg = [false,true,true, true, true]


m = Model()
set_optimizer(m, Cbc.Optimizer)
set_optimizer_attribute(m, MOI.Silent(), true)
@variable(m, x[i=1:5], upper_bound=ub[i],lower_bound=lb[i],integer=itg[i])
@objective(m, Max, [1,2,5, 10, -1]'*x)
@constraint(m, constraint1, [1,1,3,2,-5]'*x <= -5)
@constraint(m, constraint2, [1,-2,-7,5,1]'*x <= 10)
@constraint(m, constraint3, [2,5,8,5,1]'*x >= 4)
optimize!(m)

################################################################
## model to dictionary!
################################################################
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

str_var_idx = Dict()
for i in 1:length(all_variables(m))
    str_var_idx[string(all_variables(m)[i])] = all_variables(m)[i]
end

for i in 1:length(lb)
    var_lb[str_var_idx["x[$i]"]] = lb[i] 
    var_ub[str_var_idx["x[$i]"]] = ub[i] 
    var_idx[str_var_idx["x[$i]"]] = i
end


println("true MIP solution x = ", JuMP.value.(x), " and ", JuMP.objective_value(m))

################################################################
## LP relaxation
################################################################
variable = all_constraints(m, VariableRef, MathOptInterface.EqualTo{Float64})
for i in variable
    delete(m,i)
end
variable = all_constraints(m, VariableRef, MathOptInterface.Interval{Float64})
for i in variable
    delete(m,i)
end

Initial_bound()
try
    UnSet_Type()
catch  end


optimize!(m)
println("relaxed solution is x = ", JuMP.value, " and ", JuMP.objective_value(m))


sol = Any[]
solution_k = Dict()
for i in keys(var)
    push!(sol,JuMP.value(i))
    solution_k[i] = JuMP.value(i)
end   

global sol
global sol_tilde = copy(sol)
for i in 1:length(sol)
    if itg[i] == true
        sol_tilde[i] = round(sol_tilde[i],RoundNearest)
    end
end


################################################################
## Bound Strengthening part
dict_infeasible_index, n_infeasible_var=Infeasible_Check(solution_k)
BoundStrength(dict_infeasible_index)
################################################################


"""


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
    lb, ub, itg, l, u, A, m, var_name, con_idx, idx_con, var_idx_qps, var, var_idx, idx_var, var_lb, var_ub, qps_var, qps = MPS_read_FP(filepath, Sec) 
end

################################################################
## LP relaxation
################################################################
# LP relaxation으로 변환 (만약 type지정하고 싶다면 Set_Type() 사용)
try
    UnSet_Type()
catch  end

optimize!(m)
println("relaxed solution is x = ", JuMP.value, " and ", JuMP.objective_value(m))

################################################################
## Initial step
################################################################
sol = Any[]
solution_k = Dict()
for i in keys(var)
    push!(sol,JuMP.value(i))
    solution_k[i] = JuMP.value(i)
end   

global sol
global sol_tilde = copy(sol)
for i in 1:length(sol)
    if itg[i] == true
        sol_tilde[i] = round(sol_tilde[i],RoundNearest)
    end
end


################################################################
## Bound Strengthening part
dict_infeasible_index, n_infeasible_var=Infeasible_Check(solution_k)
BoundStrength(dict_infeasible_index)
#################################################################


function distance_fucntion(s1,s2)
    return sum(abs.(s1-s2))
end

dist = distance_fucntion(sol,sol_tilde)



nIter=0
maxIter = 50
while distance_fucntion(sol,sol_tilde)>0 && nIter<maxIter
    m2 = copy(m)

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
    solution_k = Dict()
    for i in keys(var)
        push!(sol,JuMP.value(i))
        solution_k[i] = JuMP.value(i)
    end 
    dict_infeasible_index, n_infeasible_var=Infeasible_Check(solution_k)
    global sol, solution_k

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
    #println("Feasible MIP Solutions:")
    #println("x = ", sol)


end


# Printing the optimal solutions obtained
println("Feasible MIP Solutions:")
println("x = ", sol)

