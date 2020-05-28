function MOI.copy_to(
    dest::Optimizer,
    src::MOI.ModelLike;
    copy_names::Bool = false
)
    @assert MOI.is_empty(dest)
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
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
    Clp_loadProblem(
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
