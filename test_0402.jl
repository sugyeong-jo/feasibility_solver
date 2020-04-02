using JuMP
using MathOptInterface
using Cbc
const MOI = MathOptInterface

# functions
const SVF = MOI.SingleVariable
const SAF = MOI.ScalarAffineFunction{Float64}
const SQF = MOI.ScalarQuadraticFunction{Float64}
const VECTOR = MOI.VectorOfVariables

function get_type_dict(obj)
    T = typeof(obj)
    type_dict = Dict{Symbol,Type}()
    for (name, typ) in zip(fieldnames(T), T.types)
        type_dict[name] = typ
    end
    return type_dict
end


mutable struct JuniperProblem 
    nl_solver           :: Any
    nl_solver_options   :: Vector{Pair}
   
    model               :: JuMP.Model

    relaxation_status   :: MOI.TerminationStatusCode
    relaxation_objval   :: Float64
    relaxation_solution :: Vector{Float64}

    status              :: MOI.TerminationStatusCode
    objval              :: Float64
    best_bound          :: Float64

    x                   :: Vector{JuMP.VariableRef}
    primal_start        :: Vector{Real}
    num_constr          :: Int64
    num_nl_constr       :: Int64
    num_q_constr        :: Int64
    num_l_constr        :: Int64
    num_var             :: Int64
    l_var               :: Vector{Float64}
    u_var               :: Vector{Float64}

    has_nl_objective    :: Bool
    nlp_evaluator       :: MOI.AbstractNLPEvaluator

    objective           :: Union{SVF, SAF, SQF, Nothing}

    disc2var_idx        :: Vector{Int64}
    var2disc_idx        :: Vector{Int64}

    var_type            :: Vector{Symbol}
    obj_sense           :: Symbol
    num_disc_var        :: Int64

    solution            :: Vector{Float64}

    soltime             :: Float64
    #options             :: SolverOptions
    #solutions           :: Vector{SolutionObj}
    nsolutions          :: Int64

    mip_solver          :: Any
    mip_solver_options  :: Vector{Pair}

    relaxation_time     :: Float64
    start_time          :: Float64

    # Info  
    nintvars            :: Int64
    nbinvars            :: Int64
    nnodes              :: Int64
    ncuts               :: Int64
    nbranches           :: Int64
    nlevels             :: Int64

    fpump_info          :: Dict{Symbol,Float64}

    # debug 
    debugDict           :: Dict{Symbol,Any}

    JuniperProblem() = new()
end

function init_juniper_problem!(jp::JuniperProblem, model::MOI.AbstractOptimizer)
    num_variables = length(model.variable_info)
    num_linear_le_constraints = length(model.linear_le_constraints)
    num_linear_ge_constraints = length(model.linear_ge_constraints)
    num_linear_eq_constraints = length(model.linear_eq_constraints)
    num_quadratic_le_constraints = length(model.quadratic_le_constraints)
    num_quadratic_ge_constraints = length(model.quadratic_ge_constraints)
    num_quadratic_eq_constraints = length(model.quadratic_eq_constraints)

    jp.status = MOI.OPTIMIZE_NOT_CALLED
    jp.relaxation_status = MOI.OPTIMIZE_NOT_CALLED
    jp.has_nl_objective = model.nlp_data.has_objective
    jp.nlp_evaluator = model.nlp_data.evaluator
    jp.objective = model.objective

    jp.objval = NaN
    jp.best_bound = NaN
    jp.solution = fill(NaN, num_variables)
    jp.nsolutions = 0
    jp.solutions = []
    jp.num_disc_var = 0
    jp.nintvars = 0
    jp.nbinvars = 0
    jp.nnodes = 1 # is set to one for the root node
    jp.ncuts = 0
    jp.nbranches = 0
    jp.nlevels = 1
    jp.relaxation_time = 0.0

    jp.start_time = time()

    jp.nl_solver = model.options.nl_solver
    nl_vec_opts = Vector{Pair}()
    if isa(jp.nl_solver, MOI.OptimizerWithAttributes)
        for arg in model.options.nl_solver.params
            push!(nl_vec_opts, arg)
        end
    end

    jp.nl_solver_options = nl_vec_opts

    if model.options.mip_solver !== nothing
        jp.mip_solver = model.options.mip_solver
        mip_vec_opts = Vector{Pair}()
        if isa(jp.mip_solver, MOI.OptimizerWithAttributes)
            for arg in model.options.mip_solver.params
                push!(mip_vec_opts, arg)
            end
        end
    
        jp.mip_solver_options = mip_vec_opts
    end
    jp.options = model.options    
    if model.sense == MOI.MIN_SENSE 
        jp.obj_sense = :Min
    else
        jp.obj_sense = :Max
    end
    jp.l_var = info_array_of_variables(model.variable_info, :lower_bound)
    jp.u_var = info_array_of_variables(model.variable_info, :upper_bound)
    integer_bool_arr = info_array_of_variables(model.variable_info, :is_integer)
    binary_bool_arr = info_array_of_variables(model.variable_info, :is_binary)
    primal_start_arr = info_array_of_variables(model.variable_info, :start)
    jp.primal_start = primal_start_arr
    jp.nintvars = sum(integer_bool_arr)
    jp.nbinvars = sum(binary_bool_arr)
    jp.num_disc_var = sum(integer_bool_arr)+sum(binary_bool_arr)
    jp.num_var = length(model.variable_info)
    jp.var_type = [:Cont for i in 1:jp.num_var]
    jp.var_type[integer_bool_arr .== true] .= :Int
    jp.var_type[binary_bool_arr .== true] .= :Bin
    jp.disc2var_idx = zeros(jp.num_disc_var)
    jp.var2disc_idx = zeros(jp.num_var)
    int_i = 1
    for i=1:jp.num_var
        if jp.var_type[i] != :Cont
            jp.disc2var_idx[int_i] = i
            jp.var2disc_idx[i] = int_i
            int_i += 1
        end
    end
    jp.num_l_constr = num_linear_le_constraints+num_linear_ge_constraints+num_linear_eq_constraints
    jp.num_q_constr = num_quadratic_le_constraints+num_quadratic_ge_constraints+num_quadratic_eq_constraints
    jp.num_nl_constr = length(model.nlp_data.constraint_bounds)
    jp.num_constr = jp.num_l_constr+jp.num_q_constr+jp.num_nl_constr
end

function info_array_of_variables(variable_info::Vector{VariableInfo}, attr::Symbol)
    len_var_info = length(variable_info)
    type_dict = get_type_dict(variable_info[1])
    result = Array{type_dict[attr], 1}(undef, len_var_info)
    for i = 1:len_var_info
        result[i] = getfield(variable_info[i], attr)
        # if type is binary then set bounds correctly
        if result[i] < 0 && attr == :lower_bound && getfield(variable_info[i], :is_binary)
            result[i] = 0
        end
        if result[i] > 1 && attr == :upper_bound && getfield(variable_info[i], :is_binary)
            result[i] = 1
        end
    end
    return result
end




init_juniper_problem!(JuniperProblem(), model::MOI.AbstractOptimizer)



#cd("/HDD/Workspace/CLT/FP/data")
cd("/HDD/Workspace/CLT/mps/data")
pwd()
println("============================flag! read======================================")

#model = read_from_file("mzzv42z.mps")
model = read_from_file("R100701005_2_cplex.mps")
optimizer = Cbc.Optimizer
set_optimizer(model, optimizer)
#set_optimizer(model, model_manual)
MOIU.attach_optimizer(model)


#------------------------------------------------
#listing all variables
var = all_variables(model)
x = []
for i in var
    push!(x, has_lower_bound(i))
end
sum(x) #0

for i in var
    push!(x, has_upper_bound(i))
end
sum(x) #0
#result> variable마다 UB, LB는 정해져 있지 않음
#---------------------------------------------------


#==
EDA
==#
objective_sense(model)
objective_function(model)
objective_function_type(model)

list_of_constraint_types(model)
less_than_constraints = all_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.LessThan{Float64})
grater_than_constraints = all_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.GreaterThan{Float64})
equal_to_constraints = all_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.EqualTo{Float64})
binary_variable = all_constraints(model, VariableRef, MathOptInterface.ZeroOne)
integer_variable = all_constraints(model, VariableRef, MathOptInterface.Integer)
typeof(integer_variable[1])

con = constraint_object(less_than_constraints[1])
con.func
typeof(con.func)
con.set
name(less_than_constraints[1])

#==
구현 해야 할 것

1. 핸들링하기 쉬운 형태로 만들기
 - constraint에 속해있는 variable들을 각각 선택할 수 있게 만드는 함수
 - constraint의 rhs를 따로 데어내 사전으로 만듦
 
2. mip-generator 
 - variable type을 모두 제거 한 후 LP로 풀기
 - 참고 코드
 function generate_mip(optimizer, m, nlp_sol, tabu_list, start_fpump)
    mip_optimizer = m.mip_solver.optimizer_constructor()
    mip_model = Model(m.mip_solver)
    @variable(mip_model, mx[i = 1:m.num_var], 
        binary = m.var_type[i] == :Bin, 
        integer = m.var_type[i] == :Int
    )

    # only add bounds for non binary variables
    for i=1:m.num_var
        if m.var_type[i] != :Bin
           @constraint(mip_model, m.l_var[i] <= mx[i] <= m.u_var[i])
        end
        
        if m.var_type[i] == :Bin && (m.l_var[i] > 0 || m.u_var[i] < 1)
            # must be 1
            if m.l_var[i] > 0
                @constraint(mip_model, mx[i] == 1)
            else # or 0
                @constraint(mip_model, mx[i] == 0)
            end
        end
    end

    backend = JuMP.backend(mip_model);

    llc = optimizer.linear_le_constraints
    lgc = optimizer.linear_ge_constraints
    lec = optimizer.linear_eq_constraints
    for constr_type in [llc, lgc, lec]
        for constr in constr_type
            MOI.add_constraint(backend, constr[1], constr[2])
        end
    end

    @variable(mip_model, mabsx[i=1:m.num_disc_var] >= 0)
    for i=1:m.num_disc_var
        vi = m.disc2var_idx[i]
        @constraint(mip_model, mabsx[i] >= mx[vi]-nlp_sol[vi])
        @constraint(mip_model, mabsx[i] >= -mx[vi]+nlp_sol[vi])
    end

    # How long is the tabu list
    num_sols = 0
    for i=1:tabu_list.length
        if !isnan(tabu_list.sols[i][1])
            num_sols += 1
        else
            break
        end
    end

    # If there solutions in the tabu list => avoid them
    if num_sols > 0
        @variable(mip_model, z1[j=1:m.num_disc_var,k=1:num_sols], Bin)
        @variable(mip_model, z2[j=1:m.num_disc_var,k=1:num_sols], Bin)
        v = tabu_list.sols
        for k=1:num_sols, j=1:m.num_disc_var
            i = m.disc2var_idx[j]
            lbi = m.l_var[i] > typemin(Int64) ? m.l_var[i] : typemin(Int64)
            ubi = m.u_var[i] < typemax(Int64) ? m.u_var[i] : typemax(Int64)
            @constraint(mip_model, z1[j,k]+z2[j,k] <= 1)
            @constraint(mip_model, (lbi - v[k][j])*z1[j,k]+z2[j,k]+v[k][j] <= mx[i])
            @constraint(mip_model, mx[i] <= v[k][j] - z1[j,k] + (ubi-v[k][j])*z2[j,k])
        end
        for k=1:num_sols
            @constraint(mip_model, sum(z1[j,k]+z2[j,k] for j=1:m.num_disc_var) >= 1)
        end
    end

    @objective(mip_model, Min, sum(mabsx[i] for i=1:m.num_disc_var))

    # Break the mip solver if it takes too long or throw a warning when this option isn't available 
    current_time = time()-start_fpump  
    time_left = m.options.feasibility_pump_time_limit-current_time
    time_left < 0 && (time_left = 1.0)

    # set time limit if supported
    old_time_limit = set_time_limit!(mip_optimizer, time_left)
    
    status, backend = optimize_get_status_backend(mip_model)

    # reset time limit
    set_time_limit!(mip_optimizer, old_time_limit)

    obj_val = NaN
    values = fill(NaN,m.num_var)
    if state_is_optimal(status; allow_almost=m.options.allow_almost_solved)
        obj_val = JuMP.objective_value(mip_model)

        # round mip values
        values = JuMP.value.(mx)
        for i=1:m.num_disc_var
            vi = m.disc2var_idx[i]
            values[vi] = round(values[vi])
        end
    end

  
    return status, values, obj_val
end

3.  MOI 핸들링 정리 필요

==#
