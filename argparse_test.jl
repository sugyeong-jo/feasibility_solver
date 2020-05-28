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
filename=string(parsed_args["filename"])



