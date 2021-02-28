
Ψ = qnmfunction(s=2,l=2,m=2,n=2,a=0.1)
Ψ2 = qnmfunction(s=-1,l=2,m=2,n=2,a=0.1)

Ψ2 = qnmfunction(s=2,l=2,m=-2,n=3,a=0.7)
Ψ2

using Plots
using LaTeXStrings
r = (Ψ.R.r₊+0.2):0.001:10
rʼ = 8.0:0.1:10
plotlyjs()
using Statistics
let lw=1, height = 3
Ψ1 = qnmfunction(s=-1,l=2,m=2,n=2,a=0.1)
Ψ2 = qnmfunction(s=2,l=2,m=-2,n=0,a=0.2)
rʼ = 2.0:0.1:10
wave1 = real.(Ψ1.(rʼ,0.01,0.01,2.0))
wave2 = real.(Ψ2.(rʼ,0.01,0.01,2.0))
wh1 = mean(abs.(wave1))
wh2 = mean(abs.(wave2))
#wh = height*max(wave1height,wave2height)
p = plot(rʼ,wave1./wh1,
        ylim = [-height,height],
        lab = L"{\ }_{%$(Ψ1.s)}\Psi_{%$(Ψ1.l),%$(Ψ1.m)}",
        linewidth = lw)

plot!(p,rʼ,wave2./wh2,
        ylim = [-height,height],
        lab = L"{\ }_{%$(Ψ2.s)}\Psi_{%$(Ψ2.l),%$(Ψ2.m)}",
        linewidth = lw)

end

gr()
Ψm2 = qnmfunction(s=-2,l=2,m=2,n=0,a=0.1)
Ψm1 = qnmfunction(s=-2,l=2,m=1,n=0,a=0.1)
Ψ0 = qnmfunction(s=-2,l=2,m=0,n=0,a=0.1)
Ψp1 = qnmfunction(s=-2,l=2,m=-1,n=0,a=0.1)
Ψp2 = qnmfunction(s=-2,l=2,m=-2,n=0,a=0.1)
ΨArray = [Ψm2,Ψm1,Ψ0,Ψp1,Ψp2]
rʼ = 2.0:0.1:10
lw=1; height = 3;
anim = @animate for t ∈ 0:0.2:3
        N = length(ΨArray)
        waves = [];
        whs = [];
        for i in 1:N
            thewave = real.(ΨArray[i].R.(rʼ)*exp(-ΨArray[i].ω*im*t))
            push!(waves,thewave)
            push!(whs,mean(abs.(thewave)))
        end

    wh = height*maximum(whs)
    p = plot(rʼ,waves[1],
            ylim = [-200,200],
            lab = L"{\ }_{%$(ΨArray[1].s)}R_{%$(ΨArray[1].l),%$(ΨArray[1].m)}",
            linewidth = lw)
        for i in 2:N
            plot!(p,rʼ,waves[i],
                ylim = [-200,200],
                lab = L"{\ }_{%$(ΨArray[i].s)}R_{%$(ΨArray[i].l),%$(ΨArray[i].m)}",
                linewidth = lw)
        end
        xlabel!(p,L"r (\textrm{in units of M})")
        ylabel!(p,L"{\ }_{s}R_{lm}")
        title!(p, "Quasinormal Mode Radial Function time dependence \n Does blow up at infinity")

end
display(anim)
gif(anim, "anim_fps15.gif", fps = 30)

x = 0:0.1:100
ω = 0.4 - im*0.5
anim2 = @animate for t ∈ 0:1:100
        plot(x, map(r-> real(exp(-im*ω*(t-r))),x), ylim = [-1000,1000])
end
gif(anim2, "anim_fps15.gif", fps = 30)


Ψ2(Ψ2.R.r₊+0.000001)

function factorial2(n)
    if n > 20
        return n
    return factorial(n)
end

function SterlingsApprox(n)
    return sqrt(2*pi*n)*(n/e)^n
end

#=
using PyCall
qnm = pyimport("qnm")
grav_220 = qnm.modes_cache(s=-2,l=2,m=2,n=2)
omega, A, C = grav_220(a=0.1)
println(omega)
Ψ1 = qnmfunction(s=-2,l=2,m=2,n=2,a=0.1)
println(Ψ1.ω)
grav_220 = qnm.modes_cache(s=-2,l=2,m=-2,n=2)
omega, A, C = grav_220(a=0.1)
println(omega)
Ψ2 = qnmfunction(s=-2,l=2,m=-2,n=2,a=0.1)
println(Ψ2.ω)
=#
