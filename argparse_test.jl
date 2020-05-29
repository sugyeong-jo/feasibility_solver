using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "filename"
            help = "a positional argument"
            required = true
    end

    return parse_args(s)
end


parsed_args = parse_commandline()
println(parsed_args[filename])

println("Parsed args:")
for (arg,val) in parsed_args
    println("  $arg  =>  $val")
    filename = val
end

println(filename)