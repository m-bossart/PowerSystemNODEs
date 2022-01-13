using Dates
open(string("timing", rand()), "w") do io
    print(io, now())
    print(io, "\n")
    sleep(10)
    print(io, now())
    print(io, "\n")
    sleep(10)
    print(io, now())
    print(io, "\n")
    sleep(10)
    print(io, now())
    print(io, "\n")
    sleep(10)
end
