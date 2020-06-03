println("================패키지 부르기=================")
using ProgressBars
using MathOptInterface
using JuMP, Cbc
const MOI = MathOptInterface
using SparseArrays

#using MathProgBase
#using GLPKMathProgInterface
#using Random
#using Pkg
#Pkg.clone("https://github.com/odow/LPWriter.jl")
#using LPWriter
include("FeasSolFinder_func_ver.1.0.jl")
##########################

m = read_from_file("/HDD/Workspace/CLT/mps/processing/CPLEX_file/4_R170725003_veri_cplex.mps")
optimizer=optimizer_with_attributes(Cbc.Optimizer,  "maxSolutions" => 1, "logLevel "=>1, "seconds" => 300,"allowableGap "=>70)
set_optimizer(m, optimizer)

# constraint
list_of_constraint_types(m)
con = Dict()
con_rhs = Dict()
idx_con = Dict()
con_idx = Dict()

constraint = all_constraints(m, GenericAffExpr{Float64,VariableRef}, MathOptInterface.EqualTo{Float64})
for i in tqdm(1:length(constraint))
    con_name = name(constraint[i])
    con_rhs[con_name] = constraint_object(constraint[i]).set.value
    con[con_name]= :Equal
    idx_con[constraint[i].index.value] = con_name
    con_idx[con_name] = constraint[i].index.value
end
constraint = all_constraints(m, GenericAffExpr{Float64,VariableRef}, MathOptInterface.GreaterThan{Float64})
for i in tqdm(1:length(constraint))
    con_name = name(constraint[i])
    con_rhs[con_name] = constraint_object(constraint[i]).set.lower
    con[con_name]= :Greater
    idx_con[constraint[i].index.value] = con_name
    con_idx[con_name] = constraint[i].index.value
end
constraint = all_constraints(m, GenericAffExpr{Float64,VariableRef}, MathOptInterface.LessThan{Float64})
for i in tqdm(1:length(constraint))
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

for i in tqdm(all_variables(m))
    var_name = string(i)
    var[i] = :Con
    var_idx[i] = i.index.value
    idx_var[i.index.value] = i
    var_lb[i] = -Inf
    var_ub[i] = Inf
end

variable = all_constraints(m, VariableRef, MathOptInterface.Interval{Float64})
for i in tqdm(1:length(variable))
    var_name = constraint_object(variable[i]).func
    var_lb[var_name] = constraint_object(variable[i]).set.lower
    var_ub[var_name] = constraint_object(variable[i]).set.upper
end

variable = all_constraints(m, VariableRef, MathOptInterface.GreaterThan{Float64})
for i in tqdm(1:length(variable))
    var_name = constraint_object(variable[i]).func
    var_lb[var_name] = constraint_object(variable[i]).set.lower
end
variable = all_constraints(m, VariableRef, MathOptInterface.LessThan{Float64})
for i in tqdm(1:length(variable))
    var_name = constraint_object(variable[i]).func
    var_ub[var_name] = constraint_object(variable[i]).set.upper
end
variable = all_constraints(m, VariableRef, MathOptInterface.Integer)
for i in tqdm(1:length(variable))
    var_name = constraint_object(variable[i]).func
    var[var_name] = :Int
end

variable = all_constraints(m, VariableRef, MathOptInterface.ZeroOne)
for i in tqdm(1:length(variable))
    var_name = constraint_object(variable[i]).func
    var[var_name] = :Bin
    var_lb[var_name] = 0
    var_ub[var_name] = 1
end


# Sparse matrix
I = Int[]
J = Int[]
V = Float64[]
u = Dict()
l = Dict()


con_set = [k for (k,v) in con if v==:Less]
for i in tqdm(con_set)
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
for i in tqdm(con_set)
    u[con_idx[i]] = -(con_rhs[i])
    l[con_idx[i]] = -Inf
    con_term =collect(linear_terms(constraint_object(constraint_by_name(m, i )).func))
    for j in 1:length(con_term)
        push!(I, con_idx[i])
        push!(J, var_idx[con_term[j][2]])
        push!(V, -(con_term[j][1]))
    end
end
con_set = [k for (k,v) in con if v==:Equal]
for i in tqdm(con_set)
    u[con_idx[i]] = con_rhs[i]
    l[con_idx[i]] = con_rhs[i]
    con_term =collect(linear_terms(constraint_object(constraint_by_name(m, i )).func))
    for j in 1:length(con_term)
        push!(I, con_idx[i])
        push!(J, var_idx[con_term[j][2]])
        push!(V, con_term[j][1])
    end
end

A = sparse(I,J,V)

# unfix the variable
for i in tqdm(keys(var))
    #println(var_org[i]) 
    try 
        unfix(i)
    catch
        continue
    end
end

dict_bound_check = Dict()
for i in keys(var)
    try
        if has_lower_bound(i)
            dict_bound_check[i] = var[i]
        end
    catch end
    try
        if has_upper_bound(i)
            dict_bound_check[i] = var[i]
        end
    catch end

end
[k for (k,v) in dict_bound_check if v==:Con]
##############################

Initial_bound()
UnSet_Type()
optimize!(m)

var_ = [k for (k,v) in var if v==:Bin]
has_lower_bound(var_[1])
is_binary(var_[3])

JuMP.set_binary(var_[1])


is_integer(var_[1])




set_fail = Dict()
for i in tqdm(keys(var))
    try
        JuMP.set_lower_bound(i,var_lb[i])
        JuMP.set_upper_bound(i,var_ub[i])
    catch 
        set_fail[i] = var[i]    
    end    
end

set_fail = Dict()
for i in tqdm(keys(var))
    JuMP.set_lower_bound(i,var_lb[i])
    JuMP.set_upper_bound(i,var_ub[i])
end

x = collect(keys(set_fail))[1]

set_upper_bound(x,var_ub[x])
con[idx_con[2843]]



var_bin = [k for (k,v) in var if v==:Bin]
var_bin_array = Array{Any}(undef,0)
for i in 1:length(var_bin)
    push!(var_bin_array,var_ub[var_bin[i]])
end
# collect(keys(con)) : constraint index set


dict_var_idx = Dict()
for i in tqdm(1:length(var))
    dict_var_idx[var[i]] = i
end


constraint_by_name(m, var[1])

con_term =collect(linear_terms(constraint_object(constraint_by_name(m, con_set[1])).func))


var = Dict()

dict_var_fix = Dict()
dict_bin_type = Dict()
dict_int_type = Dict()
dict_IB_type = Dict()
for x in tqdm(collect(var))
    if is_fixed(x) == true
        dict_var_fix[x] = fix_value(x)
    end
    if is_binary(x) == true
        dict_bin_type[x] = :Bin
        dict_IB_type[x] = :Bin

    end
    if is_integer(x) == true
        dict_int_type[x] = :Int
        dict_IB_type[x] = :Int
    end
end

# unfix the variable
for i in tqdm(var)
    try 
        unfix(i)
    catch
        continue
    end
end


# binary bound setting
for i in tqdm(keys(dict_bin_type))
    JuMP.set_lower_bound(i,0)
    JuMP.set_upper_bound(i,1)    
end

dict_xlb = Dict()
dict_xub = Dict()
for i in tqdm(var)
    if has_lower_bound(i) == true
        dict_xlb[i] = lower_bound(i)
    else 
        dict_xlb[i] = -Inf
    end
    if has_upper_bound(i) == true
        dict_xub[i] = upper_bound(i)
    else 
        dict_xub[i] = +Inf
    end
    
end


Initial_bound()
UnSet_Type()
Set_Type()

list_of_constraint_types(m)[1][1]



dict_xub[dict_var_str["Y_SGSIN1701W_NLRTM_DG"]]

i = 1
con=all_constraints(m, GenericAffExpr{Float64,VariableRef}, MathOptInterface.EqualTo{Float64})[i]
con_name = name(con)



con_term= collect(linear_terms(constraint_object(con).func))
con_var_num = length(con_term)
j = 1
cond_var = con_term[j][2]
con_var_coef = con_term[j][1]

I = Int[]
J = Int[]
V = Float64[]



fieldnames(MOI.Interval{Float64})

collect(linear_terms(constraint_object(all_constraints(m, GenericAffExpr{Float64,VariableRef}, MathOptInterface.EqualTo{Float64})).func))


collect(constraint_object(all_constraints(m, VariableRef, MathOptInterface.Interval{Float64})[1]).set)
constraint_object(all_constraints(m, VariableRef, MathOptInterface.Interval{Float64})[1]).set.upper