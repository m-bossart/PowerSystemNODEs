try
    [1 2] * [3 2 1]
catch e
    @error e "Something went wrong" exception = (e, catch_backtrace())
end
