#Script to verify going from fft to the sin/cos form needed for
#The variable periodic source in PSID.
using Plots
using FFTW
v(t) = 0.5 + 3* sin(500*2*pi*t) +1.5*cos(7*2*pi*t)
tspan = (0.0, 2.0)
step = 1e-3
tsteps = tspan[1]:step:tspan[2]
V = v.(tsteps)
p1 = plot(tsteps,V, label = "original")

N = length(V)
fs = (N-1)/(tspan[2]-tspan[1])
freqs = fftfreq(N, fs)
freqs_pos = freqs[freqs .>= 0]
F_V = fft(V)
F_V = F_V[freqs .>= 0]
F_V = F_V/N
F_V[2:end]= F_V[2:end]*2
V_reconstruct = zeros(length(tsteps))
V_reconstruct = V_reconstruct .+ abs(F_V[1])
for (i,f) in enumerate(freqs_pos[2:end])
    global V_reconstruct
    V_reconstruct -= imag(F_V[i+1])* sin.(f .*2 .* pi .* tsteps)
    V_reconstruct += real(F_V[i+1])* cos.(f .*2 .* pi .* tsteps)
    print(f,"\n")
end
p2 = scatter(freqs_pos,abs.(F_V))
plot!(p1,tsteps, V_reconstruct,label="reconstructed")
plot(p1,p2, layout = (2,1))
