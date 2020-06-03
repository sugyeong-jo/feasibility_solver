import MathOptInterface
import SparseArrays

const MOI = MathOptInterface

# Supported scalar sets
const SCALAR_SETS = Union{
    MOI.GreaterThan{Float64},
    MOI.LessThan{Float64},
    MOI.EqualTo{Float64},
    MOI.Interval{Float64}
}

# Maps Clp's parameters to getter/setter function
const CLP_OPTION_MAP = Dict(
    :PrimalTolerance => (Clp_primalTolerance, Clp_setPrimalTolerance),
    :DualTolerance => (Clp_dualTolerance, Clp_setDualTolerance),
    :DualObjectiveLimit => (Clp_dualObjectiveLimit, Clp_setDualObjectiveLimit),
    :MaximumIterations => (
        # TODO(odow): Clp doesn't follow convention with maximumIterations.
        maximumIterations, Clp_setMaximumIterations
    ),
    :MaximumSeconds => (Clp_maximumSeconds, Clp_setMaximumSeconds),
    :LogLevel => (Clp_logLevel, Clp_setLogLevel),
    :Scaling => (Clp_scalingFlag, Clp_scaling),
    :Perturbation => (Clp_perturbation, Clp_setPerturbation),
    :Algorithm => (Clp_algorithm, Clp_setAlgorithm)
)

const SOLVE_OPTION_MAP = Dict(
   :PresolveType => (ClpSolve_getPresolveType, ClpSolve_setPresolveType),
   :SolveType => (ClpSolve_getSolveType, ClpSolve_setSolveType),
   :InfeasibleReturn => (ClpSolve_infeasibleReturn, ClpSolve_setInfeasibleReturn)
)

mutable struct Optimizer <: MOI.AbstractOptimizer
    inner::Ptr{Cvoid}
    solver_options::Ptr{Cvoid}
    options_set::Set{Symbol}
    optimize_called::Bool
    solve_time::Float64
    # Work-around for upstream bug in Clp:
    maximumSeconds::Float64

    function Optimizer(; kwargs...)
        model = new(
            Clp_newModel(),
            ClpSolve_new(),
            Set{Symbol}(),
            false,
            0.0,
            -1.0,
        )
        for (key, value) in kwargs
            MOI.set(model, MOI.RawParameter(String(key)), value)
        end
        finalizer(model) do m
            Clp_deleteModel(m.inner)
            ClpSolve_delete(m.solver_options)
        end
        return model
    end
end

# ====================
#   empty functions
# ====================

function MOI.is_empty(model::Optimizer)
    # A problem is empty if it has no variable and no linear constraints
    return (Clp_getNumRows(model.inner) == 0) && (Clp_getNumCols(model.inner) == 0)
end

function MOI.empty!(model::Optimizer)
    old_model = model.inner
    model.inner = Clp_newModel()
    model.optimize_called = false
    model.solve_time = 0.0
    # Copy parameters from old model into new model
    for key in model.options_set
        getter, setter = CLP_OPTION_MAP[key]
        setter(model.inner, getter(old_model))
    end
    Clp_deleteModel(old_model)
    # Work-around for maximumSeconds
    Clp_setMaximumSeconds(model.inner, model.maximumSeconds)
    return
end

MOI.get(::Optimizer, ::MOI.SolverName) = "Clp"

# TODO: improve type-stability of `MOI.RawParameter`-related methods.
MOI.supports(::Optimizer, param::MOI.RawParameter) = true

function MOI.set(model::Optimizer, param::MOI.RawParameter, value)
    key = Symbol(param.name)
    if haskey(CLP_OPTION_MAP, key)
        push!(model.options_set, key)
        CLP_OPTION_MAP[key][2](model.inner, value)
    elseif haskey(SOLVE_OPTION_MAP, key)
        SOLVE_OPTION_MAP[key][2](model.solver_options, value)
    else
        throw(MOI.UnsupportedAttribute(param))
    end
    return
end

function MOI.get(model::Optimizer, param::MOI.RawParameter)
    key = Symbol(param.name)
    if haskey(CLP_OPTION_MAP, key)
        return CLP_OPTION_MAP[key][1](model.inner)
    elseif haskey(SOLVE_OPTION_MAP, key)
        return SOLVE_OPTION_MAP[key][1](model.solver_options)
    else
        throw(MOI.UnsupportedAttribute(param))
    end
end

MOI.supports(::Optimizer, ::MOI.Silent) = true
function MOI.set(model::Optimizer, ::MOI.Silent, value::Bool)
    push!(model.options_set, :LogLevel)
    Clp_setLogLevel(model.inner, value ? 0 : 1)
    return
end

function MOI.get(model::Optimizer, ::MOI.Silent)
    return Clp_logLevel(model.inner) == 0
end

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, value)
    push!(model.options_set, :MaximumSeconds)
    value = value === nothing ? -1.0 : value
    Clp_setMaximumSeconds(model.inner, value)
    model.maximumSeconds = value
    return
end

function MOI.get(model::Optimizer, ::MOI.TimeLimitSec)
    # TODO(odow): replace with `Clp_maximumSeconds(model.inner)` when upstream
    # is fixed.
    return model.maximumSeconds
end

MOI.supports(::Optimizer, ::MOI.NumberOfThreads) = false

# ========================================
#   Supported constraints and objectives
# ========================================

function MOI.supports_constraint(
    ::Optimizer, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{<:SCALAR_SETS}
)
    return true
end

function MOI.supports_constraint(
    ::Optimizer, ::Type{MOI.SingleVariable}, ::Type{<:SCALAR_SETS}
)
    return true
end

MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true

function MOI.supports(
    ::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}
)
    return true
end

# =======================
#   `copy_to` function
# =======================

_add_bounds(::Vector{Float64}, ub, i, s::MOI.LessThan{Float64}) = ub[i] = s.upper
_add_bounds(lb, ::Vector{Float64}, i, s::MOI.GreaterThan{Float64}) = lb[i] = s.lower
_add_bounds(lb, ub, i, s::MOI.EqualTo{Float64}) = lb[i], ub[i] = s.value, s.value
_add_bounds(lb, ub, i, s::MOI.Interval{Float64}) = lb[i], ub[i] = s.lower, s.upper

function _extract_bound_data(src, mapping, lb, ub, S)
    for con_index in MOI.get(
        src, MOI.ListOfConstraintIndices{MOI.SingleVariable, S}()
    )
        f = MOI.get(src, MOI.ConstraintFunction(), con_index)
        s = MOI.get(src, MOI.ConstraintSet(), con_index)
        column = mapping.varmap[f.variable].value
        _add_bounds(lb, ub, column, s)
        mapping.conmap[con_index] = MOI.ConstraintIndex{MOI.SingleVariable, S}(column)
    end
end

function _copy_to_columns(dest, src, mapping)
    x_src = MOI.get(src, MOI.ListOfVariableIndices())
    N = Cint(length(x_src))
    for i = 1:N
        mapping.varmap[x_src[i]] = MOI.VariableIndex(i)
    end

    fobj = MOI.get(src, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    c = fill(0.0, N)
    for term in fobj.terms
        i = mapping.varmap[term.variable_index].value
        c[i] += term.coefficient
    end
    # Clp seems to negates the objective offset
    Clp_setObjectiveOffset(dest.inner, -fobj.constant)
    return N, c
end

_bounds(s::MOI.GreaterThan{Float64}) = (s.lower, Inf)
_bounds(s::MOI.LessThan{Float64}) = (-Inf, s.upper)
_bounds(s::MOI.EqualTo{Float64}) = (s.value, s.value)
_bounds(s::MOI.Interval{Float64}) = (s.lower, s.upper)

function _extract_row_data(src, mapping, lb, ub, I, J, V, S)
    row = length(I) == 0 ? 1 : I[end] + 1
    for c_index in MOI.get(
        src, MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64}, S}()
    )
        print(c_index)
        f = MOI.get(src, MOI.ConstraintFunction(), c_index)
        l, u = _bounds(MOI.get(src, MOI.ConstraintSet(), c_index))
        push!(lb, l - f.constant)
        push!(ub, u - f.constant)
        for term in f.terms
            push!(I, row)
            push!(J, Cint(mapping.varmap[term.variable_index].value))
            push!(V, term.coefficient)
        end
        mapping.conmap[c_index] = MOI.ConstraintIndex{
            MOI.ScalarAffineFunction{Float64}, S
        }(length(ub))
        row += 1
    end
    return
end

function MOI.copy_to(
    dest::Optimizer,
    src::MOI.ModelLike;
    copy_names::Bool = false
)
    @assert MOI.is_empty(dest)
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        println(S)
        if !MOI.supports_constraint(dest, F, S)
            throw(MOI.UnsupportedConstraint{F, S}("Clp.Optimizer does not support constraints of type $F-in-$S."))
        end
    end
    fobj_type = MOI.get(src, MOI.ObjectiveFunctionType())
    if !MOI.supports(dest, MOI.ObjectiveFunction{fobj_type}())
        throw(MOI.UnsupportedAttribute(MOI.ObjectiveFunction(fobj_type)))
    end
    mapping = MOI.Utilities.IndexMap()
    N, c = _copy_to_columns(dest, src, mapping)
    cl, cu = fill(-Inf, N), fill(Inf, N)
    rl, ru, I, J, V = Float64[], Float64[], Cint[], Cint[], Float64[]
    for S in (
        MOI.GreaterThan{Float64},
        MOI.LessThan{Float64},
        MOI.EqualTo{Float64},
        MOI.Interval{Float64},
    )
        _extract_bound_data(src, mapping, cl, cu, S)
        _extract_row_data(src, mapping, rl, ru, I, J, V, S)
    end
    M = Cint(length(rl))
    A = SparseArrays.sparse(I, J, V, M, N)
    Clp.Clp_loadProblem(
        dest.inner,
        A.n,
        A.m,
        A.colptr .- Cint(1),
        A.rowval .- Cint(1),
        A.nzval,
        cl,
        cu,
        c,
        rl,
        ru
    )

    sense = MOI.get(src, MOI.ObjectiveSense())
    if sense == MOI.MIN_SENSE
        Clp_setObjSense(dest.inner, 1)
    elseif sense == MOI.MAX_SENSE
        Clp_setObjSense(dest.inner, -1)
    else
        @assert sense == MOI.FEASIBILITY_SENSE
        Clp_setObjSense(dest.inner, 0)
    end
    return mapping
end