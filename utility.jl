using MathOptInterface
using JuMP, Cbc
const MOI = MathOptInterface
using SparseArrays
using ArgParse
using QPSReader

##########################
s = ArgParseSettings()
s.description = "Find Feasible Solution"
s.commands_are_required = true
s.version = "1.0"
s.add_version = true

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "filename"
            help = "a positional argument"
            required = true
        "Sec"
            help = "a solving time"
            arg_type = Int
            required = true
    end

    return parse_args(s)
end

"""
# without package

function MPS_read(filepath, Sec)
    m = read_from_file(filepath)
    optimizer=optimizer_with_attributes(Cbc.Optimizer,  "maxSolutions" => 1, "logLevel "=>1, "seconds" => Sec ,"allowableGap "=>70)
    set_optimizer(m, optimizer)

    # constraint
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
end
"""

function MPS_read(filepath, Sec)  
    println("================MPS load=================")
    println("The file name is: $filename ")
    println("The limit time of solver is: $Sec s") 
    println("=========================================") 
    m = read_from_file(filepath)
    optimizer=optimizer_with_attributes(Cbc.Optimizer,  "maxSolutions" => 1, "logLevel "=>1, "seconds" => Sec ,"allowableGap "=>70)
    set_optimizer(m, optimizer)

    qps = readqps(filepath)
    var_idx_qps = qps.varindices

    # constraint
    con_idx = qps.conindices
    idx_con = Dict()
    u = Dict()
    l = Dict()
    for con_name in keys(con_idx)
        idx_con[con_idx[con_name]] = con_name
    end
    A = sparse(qps.arows, qps.acols, qps.avals, qps.ncon, qps.nvar)
    for i in 1:qps.ncon
        con_name=con_idx[qps.connames[i]]
        l[con_name] = qps.lcon[i]
        u[con_name] = qps.ucon[i]
    end




    # variable
    var = Dict()
    var_idx = Dict()
    idx_var = Dict()
    var_lb = Dict()
    var_ub = Dict()
    qps_var = Dict()

    for var_name in (all_variables(m))
        var[var_name] = :Con
        var_idx[var_name] = var_idx_qps[string(var_name)]
        idx_var[var_idx_qps[string(var_name)]] = var_name
        qps_var[string(var_name)] = var_name
    end

    variable = all_constraints(m, VariableRef, MathOptInterface.Integer)
    for i in (1:length(variable))
        var_name = constraint_object(variable[i]).func
        var[var_name] = :Int
    end
    variable = all_constraints(m, VariableRef, MathOptInterface.ZeroOne)
    for i in (1:length(variable))
        var_name = constraint_object(variable[i]).func
        var[var_name] = :Bin
    end
    var_lb = Dict()
    var_ub = Dict()
    for i in 1:qps.nvar
        var_name=qps.varnames[i]
        var_lb[qps_var[var_name]] = qps.lvar[i]
        var_ub[qps_var[var_name]] = qps.uvar[i]
    end
    return m, con_idx, idx_con, A, l, u, var_idx_qps, var, var_idx, idx_var, var_lb, var_ub, qps_var
end
