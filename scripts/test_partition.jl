using IterTools

as = [1, 2, 3, 4, 5, 6]
bs = [7, 8, 9, 10, 11, 12]

for (a, b) in zip(partition(as, 1), partition(bs, 1))
    @show a
    @show b
end
