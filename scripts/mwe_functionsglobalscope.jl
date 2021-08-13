global a = 1.0
function aa()
    global a = 2.0
end
function bb()
    print(a)
end

print(a)
aa()
bb()
