using DifferentialEquations
using Plots
function parameterized_lorenz!(du,u,p,t)
    du[1] = p[1]*(u[2]-u[1])
    du[2] = u[1]*(p[2]-u[3]) - u[2]
    du[3] = u[1]*u[2] - p[3]*u[3]
end

global u₀ = [1.0,0.0,0.0]
#u0 = [1.0,0.0,0.0]
tspan = (0.0,100.0)
p = [10.0,28.0,8/3]
prob = ODEProblem(parameterized_lorenz!,u0,tspan,p)
sol = solve(prob)
plot(sol,vars=(1,2,3))


function predict(p)
    pnew = vcat(p, 8/3)
    prob2 = remake(prob, p=pnew, u0= u₀)
    sol2 = solve(prob2)
    display(plot(sol2,vars=(1,2,3)))
    return Array(sol2)
end 

predict([10.0,28.0])
predict([10.0,0.0])