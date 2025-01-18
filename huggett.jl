using Parameters
# using Interpolations
using LinearInterpolations
using Interpolations
gamma = 2.0
w=5        # wages
beta=1.03^(-1/12)       # discount factor
r=-0.01      # interest rate
phi = -1.0
a_min= -3.0
a_max=18
N = 100
z_vals = [0.5,1]
z_size = 2
Pi = [0.4 0.6; 0.6 0.4]
mu = 0.5
agrid = range(a_min,a_max,length = N)

function time_g()
    global gg, gb
    ag1 = copy(agrid)
    ab1 = copy(agrid)
    ag = copy(agrid)
    ab = copy(agrid)
    metric = 1
    while metric >1e-6
        gg = linear_interpolation(ag1,agrid,extrapolation_bc = Line())
        
        gb = linear_interpolation(ab1,agrid,extrapolation_bc = Line())

        Eg =@. beta*(1+r)*( Pi[1,1]*(agrid*(1+r)+w-max(gg[agrid],phi))^(-gamma)+Pi[2,1]*(agrid*(1+r)+mu*w-max(gb[agrid],phi))^(-gamma) )
        Eb =@. beta*(1+r)*( Pi[1,2]*(agrid*(1+r)+w-max(gg[agrid],phi))^(-gamma)+Pi[1,1]*(agrid*(1+r)+mu*w-max(gb[agrid],phi))^(-gamma) )
        
        ag1 =@. (Eg^(-1/gamma)-w+agrid)/(1+r)
        ab1 =@. (Eb^(-1/gamma)-mu*w+agrid)/(1+r)
        metric = maximum([abs.(ag1-ag) abs.(ab1-ab)])
        
        ag = ag1 
        ab = ab1
    end
    return gg, gb 
end

gg, gb = time_g()
a1g_star, a1b_star = max.(gg(agrid),phi), max.(gb(agrid),phi)

begin    
    using Plots
    plot(agrid,agrid,linestyle = :dash, color="#FF5733",label ="")
    plot!(agrid,a1g_star,color = :blue,xlabel = "current asset",label = "z = 1")
    plot!(agrid,a1b_star,color = :red,ylabel = "next period asset",label = "z = 0.5")
    savefig("policyk.pdf")
end
z_state = reshape(z_vals,1,2)
c_star =@. 5 * z_state + (1 + r) * agrid - [a1b_star a1g_star]
plot(agrid,c_star[:,1],label = "z = 0.5")
plot!(agrid,c_star[:,2],label = "z = 1")

F = linear_interpolation(agrid, agrid)