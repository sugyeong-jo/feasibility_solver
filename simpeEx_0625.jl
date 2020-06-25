println("================패키지 부르기=================")
using ProgressBars
using MathOptInterface
using JuMP, Cbc
const MOI = MathOptInterface
using SparseArrays
using DataFrames
using CSV
using ArgParse
using JuMP  # Need to say it whenever we use JuMP
using Cbc # Loading the GLPK module for using its solver
using MathOptInterface
#using Plots
using LinearAlgebra


include("Func_0625.jl")

filename = string("atlanta-ip.mps")
#filename = string("timtab1.mps")
#filepath = string("/Users/sungilkim/Dropbox/Sungil_project/CyberLogitec/연장과제/data/",filename)
#filepath = string("C:\\Users\\sungi\\Dropbox\\Sungil_project\\CyberLogitec\\연장과제\\data\\",filename)


#filename = string("b1c1s1.mps")
Sec = 20
filepath = string("/HDD/Workspace/CLT/FP/data/",filename)
val, t_problemLaod, bytes, gctime, memallocs = @timed begin
    m, con_idx, idx_con, A, l, u, var, var_idx, idx_var, var_lb, var_ub = MPS_read_full(filepath, Sec)
end



#BoundStrength(1281)
Initial_bound()
try
    UnSet_Type()
catch
end

solution_k, dict_xlb_k, dict_xub_k, dict_LP_int, solution_k_1 = LP_solve(m)

fixed_var_element = rand(collect(keys(dict_LP_int)),1)[1]

filter((k,v) -> k == fixed_var_element,solution_k) # fixed variable 있는 solution dictionary

Update_UL(filter((k,v) -> k == fixed_var_element,solution_k))


##############################################
# 이밑을 반복하면 fiexed되는 variable을 구할 수 있습니다.
for i in tqdm(1:length(keys(idx_con))) #시간이 다소 소요됩니다.
    BoundStrength(i)
end

fixed_var = Any[]
non_fixed_var = Any[]
for i in tqdm([k for (k,v) in var if v==:Bin])
    if isapprox(abs(lower_bound(i)-upper_bound(i)),0,  atol=1e-3)
        push!(fixed_var,i)
    else
        push!(non_fixed_var, i)
    end
end


# fixed 된 variable을 제외한 다른 LP로 나온 variable 중 integer로 나온 값을 하나 fixed 시키면서 업데이트하면 늘어날 것 같습니다.
# 하지만 논문의 방법이 아니므로 논문을 더 읽어보려 합니다. 