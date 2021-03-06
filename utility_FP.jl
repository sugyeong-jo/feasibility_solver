using MathOptInterface
using JuMP, Cbc
#const MOI = MathOptInterface
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



function MPS_read_FP(filepath, Sec)  
    m = Model()
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
    u = qps.ucon
    l = qps.lcon
    
    for con_name in keys(con_idx)
        idx_con[con_idx[con_name]] = con_name
    end
    A = sparse(qps.arows, qps.acols, qps.avals, qps.ncon, qps.nvar)


    #########################################
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
    ##########################################
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

    itg = Any[]
    for i in qps.varnames
        if var[qps_var[i]]==:Con
            push!(itg, false)
        else 
            push!(itg, true)
        end
    end

    
    ## 이 부분이 추가되었습니다. (2020.06.23)
    # variable 순서대로 VariableRef 형식으로 나오는 array
    var_name = Any[]
    for i in qps.varnames
        push!(var_name, qps_var[i])
    end

    global m, var, var_lb, var_ub
    #Interval, EqualTo를 lower bound/ upper bound로 할당
    variable = all_constraints(m, VariableRef, MathOptInterface.EqualTo{Float64})
    for i in variable
        delete(m,i)
    end
    variable = all_constraints(m, VariableRef, MathOptInterface.Interval{Float64})
    for i in variable
        delete(m,i)
    end
    Initial_bound()

    return lb, ub, itg, l, u, A, m,var_name, con_idx, idx_con, var_idx_qps, var, var_idx, idx_var, var_lb, var_ub, qps_var, qps
end

function MPS_read(filepath, Sec)  
    m = Model()
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
    #con_idx = qps.conindices
    #idx_con = Dict()
    #u = qps.ucon
    #l = qps.lcon
    
    #for con_name in keys(con_idx)
    #    idx_con[con_idx[con_name]] = con_name
    #end
    #A = sparse(qps.arows, qps.acols, qps.avals, qps.ncon, qps.nvar)


    #########################################
    # variable
    var = Dict()
    var_idx = Dict()
    idx_var = Dict()
    var_lb = Dict()
    var_ub = Dict()
    qps_var = Dict()


    for var_name in (all_variables(m))
        var[var_name] = var_name
        var_idx[var_name] = var_idx_qps[string(var_name)]
        idx_var[var_idx_qps[string(var_name)]] = var_name
        qps_var[string(var_name)] = var_name
    end
    """
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

    itg = Any[]
    for i in qps.varnames
        if var[qps_var[i]]==:Con
            push!(itg, false)
        else 
            push!(itg, true)
        end
    end

    """

    var_lb = Dict()
    var_ub = Dict()
    for i in 1:qps.nvar
        var_name=qps.varnames[i]
        var_lb[qps_var[var_name]] = qps.lvar[i]
        var_ub[qps_var[var_name]] = qps.uvar[i]
    end
    ##########################################
    
    ## 이 부분이 추가되었습니다. (2020.06.23)
    # variable 순서대로 VariableRef 형식으로 나오는 array
    #var_name = Any[]
    #for i in qps.varnames
    #    push!(var_name, qps_var[i])
    #end

    global m, var, var_lb, var_ub
    #Interval, EqualTo를 lower bound/ upper bound로 할당
    variable = all_constraints(m, VariableRef, MathOptInterface.EqualTo{Float64})
    for i in variable
        delete(m,i)
    end
    variable = all_constraints(m, VariableRef, MathOptInterface.Interval{Float64})
    for i in variable
        delete(m,i)
    end

    Initial_bound()
    return m
end