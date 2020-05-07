# FP
using Pkg
Pkg.add("ProgressBars")
using ProgressBars
using MathOptInterface
using JuMP, Cbc
const MOI = MathOptInterface
using Cbc
using MathProgBase
using GLPKMathProgInterface
####################################
########  MPS problem read
####################################
m = MathProgBase.LinearQuadraticModel(GLPKSolverMIP())
MathProgBase.loadproblem!(m,"/HDD/Workspace/CLT/mps/data/R100701005_2_cplex.mps")
c = MathProgBase.getobj(m)
A = MathProgBase.getconstrmatrix(m)
nrow,ncol = size(A)
xlb = MathProgBase.getvarLB(m)
xub = MathProgBase.getvarUB(m)
l = MathProgBase.getconstrLB(m)
u = MathProgBase.getconstrUB(m)
t = MathProgBase.getvartype(m)

####################################
########  MPS problem to model
####################################
m = Model(Cbc.Optimizer)
index_x = 1: ncol
const_x = 1: nrow

# variable 만들기
@variable(m, x[i in index_x], lower_bound = xlb[i], upper_bound=xub[i])
var=all_variables(m)

# objective function 만들기
@objective(m, Min, sum( c[i]*var[i] for i in index_x) )

#constraint 만들기
for i in tqdm(1:nrow) 
    const_s=A[i,:,:]#constraint_specific
    const_s_n=1:const_s.colptr[2]-1
    const_s_coef=const_s.nzval[const_s_n,]#특정 constraint에 해당하는 coefficient
    const_s_var_index=const_s.rowval[const_s_n]#특정 constraint에 해당하는 variable index
    if length(const_s_var_index)==0
        continue
    else
        @constraint(m,  l[i] <= sum((x[const_s_var_index[j]]*const_s_coef[j]) for j = 1:length(const_s_coef) ) <= u[i])
    end
end


orig_m=copy(m)
#m = copy(orig_m)

# type 만들어 주기
# integer 조건 type dictionary 만들어주기
type_dict = Dict()
for i in tqdm(index_x)
    if t[i] != :Cont
        type_dict[var[i]] = t[i]
    else
        continue
    end
end

for i in tqdm(index_x)
    if t[i]==:Bin
        JuMP.set_binary(var[i])
        JuMP.set_lower_bound(var[i],0)
        JuMP.set_upper_bound(var[i],1)    
    elseif t[i] == :Int
        JuMP.set_integer(var[i])
    end
end

#LP를 위해 type 제거
for i in tqdm(index_x)
    if t[i]==:Bin
        JuMP.unset_binary(var[i])
    elseif t[i] == :Int
        JuMP.unset_integer(var[i])
    end
end



#########################################
## LP relaxation 완전 타이트 시키기!
#########################################

JuMP.optimize!(m)
### LP솔루션 저장
solution_k = Dict()
for i in index_x
    solution_k[var[i]] = JuMP.value(var[i])
end
# step 2/3 rounding
### only for integer
initial_solution_fixed = Dict()
solution_k_1 = Dict()
for i in keys(type_dict)
    if solution_k[i]!=round(solution_k[i])
        solution_k_1[i] = round(solution_k[i])
    else initial_solution_fixed[i] = solution_k[i]
    end
end
# 다시 한 번 확인해보기

# upper bound도 업데이트해주니 infeasible 발생! (알수가 없네 정말..)
for i in keys(initial_solution_fixed) #8437개
    if lower_bound(i)<initial_solution_fixed[i]
        JuMP.set_lower_bound(i,initial_solution_fixed[i])
    end
end

JuMP.optimize!(m) #값은 더 안 좋음

### LP솔루션 저장
solution_k = Dict()
for i in index_x
    solution_k[var[i]] = JuMP.value(var[i])
end

second_solution_fixed = Dict() #8437개 똑같넹~ 하지만 fixed된 값들은 다 다름..!
solution_k_1 = Dict()
for i in keys(type_dict)
    if solution_k[i]!=round(solution_k[i])
        solution_k_1[i] = round(solution_k[i])
    else second_solution_fixed[i] = solution_k[i]
    end
end

global initial_solution_fixed
global second_solution_fixed
while true        
    if keys(initial_solution_fixed) == keys(second_solution_fixed)
        break
    else        
        initial_solution_fixed = second_solution_fixed
        # upper bound도 업데이트해주니 infeasible 발생! (알수가 없네 정말..)
        for i in keys(initial_solution_fixed) #8437개
            if lower_bound(i)<initial_solution_fixed[i]
                JuMP.set_lower_bound(i,initial_solution_fixed[i])
            end
        end

        # step 1 LP relaxation
        JuMP.optimize!(m)
        ### LP솔루션 저장
        solution_k = Dict()
        for i in index_x
            solution_k[var[i]] = JuMP.value(var[i])
        end
        # step 2/3 rounding
        ### only for integer
        second_solution_fixed = Dict()
        solution_k_1 = Dict()
        for i in keys(type_dict)
            if solution_k[i]!=round(solution_k[i])
                solution_k_1[i] = round(solution_k[i])
            else second_solution_fixed[i] = solution_k[i]
            end
        end

        # tigthed the bound!
        # upper bound도 업데이트해주니 infeasible 발생! (알수가 없네 정말..)
        for i in keys(second_solution_fixed) #8437개
            if lower_bound(i)<second_solution_fixed[i]
                JuMP.set_lower_bound(i,second_solution_fixed[i])
            end
        end
        #JuMP.objective_value(m)
        if keys(initial_solution_fixed) == keys(second_solution_fixed)
            break
        else
            global initial_solution_fixed
            global second_solution_fixed
            println("end")
            continue 
        end
    end
end

################### LP bound 타이트 끝! ##########################
JuMP.objective_value(m)
##################################################################

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
        for i in keys(type_dict)
            if solution_k[i]!=round(solution_k[i])
                solution_k_1[i] = round(solution_k[i])
            end
        end

        # step 4 projection
        score = Array{Float64}(undef,0)
        #score_list = Dict()
        score_list = Any[]

        for i in keys(solution_k_1)
            if solution_k_1[i] == upper_bound(i)
                val = upper_bound(i)-solution_k[i]
                push!(score, val )
                #score_list[i] =val
                push!(score_list, (i, val))
            elseif solution_k_1[i] == lower_bound(i)
                val = solution_k[i]-lower_bound(i)
                push!(score, val)
                #score_list[i] = val
                push!(score_list, (i, val))
            else
                val = abs(solution_k[i]-solution_k_1[i])
                push!(score,val)
                #score_list[i] = val
                push!(score_list, (i, val))
            end
        end
        top_score_list=sort!(score_list, by = x -> x[2], rev = true)
        

        #####################################################여기에서 부터 다시 코딩 시작!######################################################
        # step 5 update
        # 가장 큰 score차이나는 top 20
        T = 20

        # solution_k_1: round 한 결과
        # solution_k_2: round 한 것 중 상위 몇개만  원래 솔루션에서 바꾼 결과
        solution_k_2 = copy(solution_k)
        for i in 1:T
            x_idx=top_score_list[i][1]
            if abs(solution_k[x_idx]-solution_k_1[x_idx])<0
                solution_k_2[x_idx] = solution_k_1[x_idx]-1
            else
                solution_k_2[x_idx] = solution_k_1[x_idx]+1
            end
        end

        # step 8 update 2
        for i in keys(solution_k_2)
            if lower_bound(i)<solution_k_2[i]
                JuMP.set_lower_bound(i,solution_k_2[i])
            end
        end



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
        dist2=sum(score)

        m = copy(orig_m)
        optimizer = Cbc.Optimizer
        set_optimizer(m, optimizer)

        if dist2 < dist     
            # step 8 update 2
            for i in keys(solution_k_2)
                if lower_bound(i)<solution_k_2[i]
                    JuMP.set_lower_bound(i,solution_k_2[i])
                end
            end

            println(solution_k)
            println("solution_update")
            for i in keys(solution_k_1)
                solution_k[i]=solution_k_1[i]
                println("The total iteration is $nIter ")
                global nIter += 1
                global dist = sum(score)  
                global solution_k
                println("This distance is $dist.")
                println("The solution is ",solution_k)
        
            end
        else 
            dist = dist2
            global dist
            continue
        end

    end 
    

end

JuMP.objective_value(m)


