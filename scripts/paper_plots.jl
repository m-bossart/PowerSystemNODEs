using LaTeXStrings
plot()

tspan = (0.0, 1.0)
step = 5e-4
tlin= tspan[1]:step:tspan[2]


p1a = plot( tlin, Vm.(tlin), xlabel=L"t \: \mathrm{[s]}", ylabel=L"v_{tvib} [p.u.]", width=1.5)
p1b = plot( tlin, Vθ.(tlin),  xlabel=L"t \: \mathrm{[s]}", ylabel=L"\theta_{tvib} [rad]", width=1.5)
p1 = plot(p1a,p1b, layout=(2,1))



p2a = plot(tsteps, ode_data[1,:], xlabel=L"t \: \mathrm{[s]}", ylabel=L"i^r_{tvib} \: \mathrm{[pu]}", width=1.5)
p2b = plot(tsteps, ode_data[1,:], xlabel=L"t \: \mathrm{[s]}", ylabel=L"i^r_{tvib} \: \mathrm{[pu]}", width=1.5, xscale=:log)
p2c = plot(tsteps, ode_data[2,:], xlabel=L"t \: \mathrm{[s]}", ylabel=L"i^i_{tvib} \: \mathrm{[pu]}", width=1.5)
p2d = plot(tsteps, ode_data[2,:], xlabel=L"t \: \mathrm{[s]}", ylabel=L"i^i_{tvib} \: \mathrm{[pu]}", width=1.5, xscale=:log)
p2 = plot(p2a,p2b,p2c,p2d, layout=(2,2))

png(p1,"V_tvib")
png(p2,"I_truth")

#IN MORNING - INSERT THE BEST COMPRISON MODEL...
data = readdlm(string("figs/","run4_comp.txt"), ',', Float64, '\n')
p3a = plot(data[:,1], data[:,2], xlabel=L"t \: \mathrm{[s]}", ylabel=L"i^r_{tvib} \: \mathrm{[pu]}", width=1.5, label="true solution",xscale=:log)
plot!(p3a, data[:,1], data[:,4], label="UODE model")
plot!(p3a, data[:,1], data[:,10], label="average inverter model")
p3b = plot(data[:,1], data[:,8],  xlabel=L"t \: \mathrm{[s]}", ylabel=L"i^r_{nn} \: \mathrm{[pu]}", width=1.5,xscale=:log)
p3= plot(p3a,p3b,layout=(2,1))
png(p3, "Ir_compare")

p3a = plot(data[:,1], data[:,3], xlabel=L"t \: \mathrm{[s]}", ylabel=L"i^r_{tvib} \: \mathrm{[pu]}", width=1.5, label="true solution",xscale=:log)
plot!(p3a, data[:,1], data[:,5], label="UODE model")
plot!(p3a, data[:,1], data[:,11], label="average inverter model")
p3b = plot(data[:,1], data[:,9],  xlabel=L"t \: \mathrm{[s]}", ylabel=L"i^r_{nn} \: \mathrm{[pu]}", width=1.5,xscale=:log)
p3= plot(p3a,p3b,layout=(2,1))
png(p3, "Ii_compare")


l1 = readdlm(string("run1_loss.txt"), ',', Float64, '\n')
l2 = readdlm(string("run2_loss.txt"), ',', Float64, '\n')
l3 = readdlm(string("run3_loss.txt"), ',', Float64, '\n')
l4 = readdlm(string("run4_loss.txt"), ',', Float64, '\n')
l5 = readdlm(string("run5_loss.txt"), ',', Float64, '\n')
l6 = readdlm(string("run6_loss.txt"), ',', Float64, '\n')
l7 = readdlm(string("run7_loss.txt"), ',', Float64, '\n')

p4= plot(l1, xlabel="iteration", ylabel="loss", width=1.5)
plot!(p4,l2)
plot!(p4,l3)
plot!(p4,l4)
plot!(p4,l5)
png(p4, "loss_compare")