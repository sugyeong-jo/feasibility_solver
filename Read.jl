println("******************************************************************")
println("*****                     OOESAlgorithm                      *****")
println("*****           A comprehensive Julia package for            *****")
println("***** Bi-Objective Mixed Integer Linear Programming Problems *****")
println("***** To report any bug, email us at:                        *****")
println("*****    amsierra@mail.usf.edu   &  hcharkhgard@usf.edu      *****")
println("*****                                                        *****")
println("***** To Support us, please cite:                            *****")
println("***** Sierra-Altamiranda, A., Charkhgard, H., 2018. A new    *****")
println("***** exact algorithm to optimize a linear function over     *****")
println("***** the set of efficient solutions for bi-objective mixed  *****")
println("***** integer linear programming. INFORMS Journal on         *****")
println("***** Computing To appear.                                   *****")
println("*****                                                        *****")
println("***** Sierra-Altamiranda, A. and Charkhgard, H., 2018.       *****")
println("***** OOESAlgorithm.jl: A julia package for optimizing a     *****")
println("***** linear function over the set of efficient solutions    *****")
println("***** for bi-objective mixed integer linear programming.     *****")
println("******************************************************************")

# OOES-알고리즘에서 jump로 mps파일의 sense와 instance받아오는 코드가 있어 참고함
# 하지만 type이나, constraint의 coefficient를 얻어오지는 않음 (integer조건이 없음)
# 여기서도 바로 핸들링 하지 않고 instance를 받아와 따로 구축하는 듯 함


function compute_objective_function_value!(solution::OOESolution, instance::MOOInstance)
    @inbounds for i in 1:size(instance.c)[1]
        if length(solution.obj_vals) >= i
            solution.obj_vals[i] = sum(instance.c[i, :] .* solution.vars)
        else
            push!(solution.obj_vals, sum(instance.c[i, :] .* solution.vars))
        end
    end
end

@inbounds function read_a_moo_instance_from_a_mathprogbase_model(model::MathProgBase.AbstractMathProgModel, sense::Vector{Symbol}, temp_val::Bool)
    var_types = MathProgBase.getvartype(model)
    v_lb = MathProgBase.getvarLB(model)
    v_ub = MathProgBase.getvarUB(model)
    for i in 1:length(v_lb)
        if v_lb[i] == 0.0 && v_ub[i] == 1.0 && var_types[i] != :Cont
            var_types[i] = :Bin
        end
    end
    A = MathProgBase.getconstrmatrix(model)
    cons_lb = MathProgBase.getconstrLB(model)
    cons_ub = MathProgBase.getconstrUB(model)
    m, n = size(A)
    c = zeros(length(sense), n)
    c_prime = zeros(length(sense), n)
    if temp_val == false
	    cons_lb = cons_lb[1:end-length(sense)]
	    cons_ub = cons_ub[1:end-length(sense)]
	    c[2:end, :] = A[end-length(sense)+1:end-length(sense)+2, :]
	    c[1, :] = A[end, :]
	    A = A[1:end-length(sense), :]
    else
	    i = 1; j = 0
            while i < length(cons_lb) + 1
		if cons_ub[i] == 0 && cons_lb[i] == 0
			j += 1
 			if i == 1
				cons_lb = cons_lb[2:end]
				cons_ub = cons_ub[2:end]
				for k in 1:n
				   c_prime[j, k] = A[1, k]
				end
				A = A[2:end, :]
			elseif i < length(cons_lb)
				cons_lb = vcat(cons_lb[1:i-1],cons_lb[i+1:end])
				cons_ub = vcat(cons_ub[1:i-1],cons_ub[i+1:end])
				for k in 1:n
				   c_prime[j, k] = A[i, k]
				end
				A = vcat(A[1:i-1, :],A[i+1:end, :])
			else
				cons_lb = cons_lb[1:end-1]
				cons_ub = cons_ub[1:end-1]
				for k in 1:n
				   c_prime[j, k] = A[end, k]
				end
				A = A[1:end-1, :]
			end
		else
			i += 1
		end
	    end
	    for i in 1:length(sense)
		for j in 1:n
			if i<length(sense)
				c[i+1,j] = c_prime[i, j]
			else
				c[1,j] = c_prime[i, j]
			end
		end
	    end
    end
    for i in 1:length(sense)
        if sense[i] == :Max
            c[i, :] = -1.0*c[i, :]
        end
    end
    m, n = size(A)
    for i in 1:m
        if cons_ub[i] != Inf && cons_lb[i] == -Inf
            cons_lb[i] = -1.0*cons_ub[i]
            cons_ub[i] = Inf
            A[i, :] = -1.0*A[i, :]
        end
    end
    sparsity = length(findall(x -> x=="0.0", A))/(m*n) 
    if sparsity >= 0.5
        A = sparse(A)
    end
    instance = MOOInstance(var_types, v_lb, v_ub, c, A, cons_lb, cons_ub)
    instance, sense
end

@inbounds function read_an_instance_from_a_jump_model(model::JuMP.Model, sense::Vector{Symbol}, version::Bool)
    if version == true
	model2 = backend(model)
	model_ = GLPK.Optimizer()
	MathOptInterface.copy_to(model_, model2)
	MathOptInterface.write_to_file(model_, "temp.lp")
	instance, sense = read_an_instance_from_a_lp_or_a_mps_file("temp.lp", sense, temp_val=true)
    else
	writeLP(model, "temp.lp")
	instance, sense = read_an_instance_from_a_lp_or_a_mps_file("temp.lp", sense, temp_val=true)
    end
    rm("temp.lp")
    instance, sense
end

@inbounds function read_an_instance_from_a_jump_model(model::JuMP.Model, sense::Vector{Symbol})
    version = Pkg.installed()["JuMP"]
    if version < v"0.19.0"
	read_an_instance_from_a_jump_model(model, sense, false)
    else
	read_an_instance_from_a_jump_model(model, sense, true)
    end
end

@inbounds function read_an_instance_from_a_lp_or_a_mps_file(filename::String, sense::Vector{Symbol}; temp_val::Bool = false)
    model = MathProgBase.LinearQuadraticModel(GLPKSolverMIP())
    MathProgBase.loadproblem!(model, filename)
    read_a_moo_instance_from_a_mathprogbase_model(model, sense, temp_val)
end

model = MathProgBase.LinearQuadraticModel(GLPKSolverMIP())
MathProgBase.loadproblem!(model,"/HDD/Workspace/CLT/mps/processing/CPLEX_file/R100701002_2_cplex.mps")


model = read_from_file("/HDD/Workspace/CLT/mps/processing/CPLEX_file/R100701002_2_cplex.mps")
optimizer=optimizer_with_attributes(Cbc.Optimizer,  "maxSolutions" => 1, "logLevel "=>1, "seconds" => 300,"allowableGap "=>70)
set_optimizer(model, optimizer)
var_org = all_variables(model)
dict_xlb = Dict()
dict_xub = Dict()
for i in index_x
    dict_xlb[var[i]] = lower_bound(var[i])
    dict_xub[var[i]] = upper_bound(var[i])
end

x = []
for i in var
    push!(x, has_lower_bound(i))
end
sum(x) #0

for i in var
    push!(x, has_upper_bound(i))
end
sum(x) #0
i = var[2]

i=1

