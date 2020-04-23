# FP
using MathOptInterface, JuMP, Cbc
const MOI = MathOptInterface
using Cbc


m = Model(Cbc.Optimizer)

c = [ -1; -2; -5]
A = [-1  1  3;
      1  3 -7]
b = [-5; 10]
lb = [0;0;0]
ub = [10;Inf;Inf]
index_x = 1:3
index_constraints = 1:2
                       
@variable(m, x[i in index_x], lower_bound = lb[i], upper_bound=ub[i])

var=all_variables(m)
t = Dict(var[2]=>:Int, var[3]=>:Bin)

for i in keys(t)
    if t[i]==:Bin
        JuMP.set_binary(i)
        JuMP.set_lower_bound(i,0)
        JuMP.set_upper_bound(i,1)    
    elseif t[i] == :Int
        JuMP.set_integer(i)
    end
end

for i in keys(t)
    if t[i]==:Bin
        JuMP.unset_binary(i)    
    elseif t[i] == :Int
        JuMP.unset_integer(i)
    end
end

@objective(m, Min, sum( c[i]*var[i] for i in index_x) )
@constraint(m, constraint[j in index_constraints],
               sum( A[j,i]*var[i] for i in index_x ) <= b[j] )

print(m)

# step 1 LP relaxation
JuMP.optimize!(m)
### LP솔루션 저장
solution_k = Dict()
for i in index_x
    solution_k[var[i]] = JuMP.value(var[i])
end

# step 2/3 rounding
### only for integer
solution_k_1 = Dict()
for i in keys(t)
    if solution_k[i]!=round(solution_k[i])
        solution_k_1[i] = round(solution_k[i])
    end
end

# step 6 projection
score = Array{Float64}(undef,0)
for i in keys(solution_k_1)
    if solution_k_1[i] == upper_bound(i)
        push!(score, (upper_bound(i)-solution_k[i]))
    elseif solution_k_1[i] == lower_bound(i)
        push!(score, (solution_k[i])-lower_bound(i))
    else
        push!(score,abs(solution_k[i]-solution_k_1[i])) 
    end
end
dist=sum(score)

maxIter = 100
nIter = 0

# step 4 while
while true
    if dist==0 || nIter > maxIter
        break
    else

        # step 1 LP relaxation
        JuMP.optimize!(m)
        ### LP솔루션 저장
        solution_k = Dict()
        for i in index_x
            solution_k[var[i]] = JuMP.value(var[i])
        end

        # step 2/3 rounding
        ### only for integer
        solution_k_1 = Dict()
        for i in keys(t)
            if solution_k[i]!=round(solution_k[i])
                solution_k_1[i] = round(solution_k[i])
            end
        end

        # step 5 
        
        
        # step 6 projection
        score = Array{Float64}(undef,0)
        for i in keys(solution_k_1)
            if solution_k_1[i] == upper_bound(i)
                push!(score, (upper_bound(i)-solution_k[i]))
            elseif solution_k_1[i] == lower_bound(i)
                push!(score, (solution_k[i])-lower_bound(i))
            else
                push!(score,abs(solution_k[i]-solution_k_1[i])) 
            end
        end


        # step 7 update
        for i in keys(solution_k_1)
            if abs(solution_k[i]-solution_k_1[i])<0
                solution_k[i] = solution_k_1[i]-1
            else
                solution_k[i] = solution_k_1[i]+1
            end
        end

        # step 8 update 2
        for i in keys(solution_k_1)
            if lower_bound(i)<solution_k_1[i]
                JuMP.set_lower_bound(i,solution_k_1[i])
            end
        end

        println(solution_k)
        println("solution_update")
        for i in keys(solution_k_1)
            solution_k[i]=solution_k_1[i]
        end
    end
    println("The total iteration is $nIter ")
    global nIter += 1
    global dist = sum(score)  
    global solution_k
    println("This distance is $dist.")
    println("The solution is ",solution_k) 
    

end

JuMP.objective_value(m)


