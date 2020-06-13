println("================패키지 부르기=================")
using ProgressBars
using MathOptInterface
using JuMP, Cbc
const MOI = MathOptInterface
using SparseArrays
using DataFrames
using CSV
using ArgParse

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

println("================MPS load=================")

parsed_args = parse_commandline()
filename=string(parsed_args["filename"])
Sec=parsed_args["Sec"]

#filename = string("R100701005_2_cplex.mps")
#Sec = 30

#=
      File read
=#

println("The file name is: $filename ")
println("The limit time of solver is: $Sec s")

filepath = string("/HDD/Workspace/CLT/mps/processing/CPLEX_file/",filename)

println("========================================")
val, t_problemLaod, bytes, gctime, memallocs = @timed begin
    m = read_from_file(filepath)
    optimizer=optimizer_with_attributes(Cbc.Optimizer, "feasibilityPump" => "on", "maxSolutions" => 1, "PASSF"=>10000, "seconds" => Sec ,"allowableGap "=>70)
    set_optimizer(m, optimizer)

end

###############################################################################################################################################################
#=
Solving!
=#
path = string("/home/sugyeong/HDD/sugyeong/Julia/feasibility_solver/cbc_result/",filename,".csv")
println("================Running!!!!!=================")
function SG()
    optimize!(m)
    termination_status(m)==MOI.OPTIMAL
end

val, t_run, bytes, gctime, memallocs = @timed begin
    SG()
end

###############################################################################################################################################################
#############################################################################
# result

check_optimal = termination_status(m)
check_bound = "null"
check_type = "null"
obj_value = objective_value(m)
t_total = t_problemLaod+t_run
println("The problem load time is :",t_problemLaod)
println("The running time is : ", t_run,"s")
println("The total time is : ", t_total)
println("Is it optimal?: ", check_optimal)
println("The objective value is :", obj_value )
println("The the number of unsatisfied variable is:")
println("    - The the number of unsatisfied bound is:",check_bound )
println("    - The the number of unsatisfied type is:", check_type)
df=DataFrame(
    Name = [filename],
    Total = [t_total],
    Problem = [t_problemLaod],
    Run = [t_run],
    Optimal = [check_optimal],
    Check_Bound = [check_bound],
    check_type = [check_type],
    obj_value = [obj_value]
)

CSV.write(path, df) 



#############################################################################
###########################################
##########################################