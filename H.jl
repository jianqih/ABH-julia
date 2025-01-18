# Huggett model version 1 
# partial equilibrium 
using SparseArrays
using Parameters
@with_kw mutable struct model
    gamma = 2
    r = 0.03 # interest rate 
    rho = 0.05 
    z1 = 0.1
    z2 = 0.2
    z = [z1 z2]
    lambda1 = 0.02
    lambda2 = 0.03
    la = [lambda1 lambda2]
    na = 500
    amin = -0.02
    amax = 2
    agrid = range(amin,amax, length = na)
    da = (amax-amin)/(na-1)
    amat = agrid*ones(1,2)
    zmat = ones(na)*z
    maxit = 200
    crit = 1e-6
end


function ss(v0,param)
    @unpack_model param
    v = v0;
    dVf = zeros(na,2);
    dVb = zeros(na,2);
    c = zeros(na,2);

    # INITIAL GUESS
    v0[:,1] = @. (z[1] + r*agrid).^(1-gamma)/(1-gamma)/rho;
    v0[:,2] = @. (z[2] + r*agrid).^(1-gamma)/(1-gamma)/rho;
    for i = 1:maxit
        V = v;
        dVf[1:na-1,:] = (V[2:na,:]-V[1:na-1,:])/da;
        dVf[na,:] .= 0;

        dVb[2:na,:] = (V[2:na,:]-V[1:na-1,:])/da;
        dVb[1,:] = @. (z + r*amin)^(-gamma);

        I_concave = dVb .> dVf;

        cf = dVf.^(-1/gamma);
        sf = @. zmat + r*amat - cf;

        cb = dVb.^(-1/gamma);
        sb = @. zmat + r*amat - cb;

        c0 = @. zmat + r.*amat;
        dV0 = c0.^(-gamma);

        If = sf .> 0 # pos of forward saving
        Ib = sb .< 0 # neg of backward saving
        I0 = @. (1 - If - Ib); # one neg one pos savings

        Ib[na,:] .= 1; # if reach upper bound, only consume can not save more. 
        If[na,:] .= 0; 

        dV_upwind = dVf.*If + dVb.*Ib + dV0.*I0;
        c = dV_upwind.^(-1/gamma)
        u = @. c^(1-gamma)/(1-gamma)
        V_switch = [V[:,2] V[:,1]];  # exchange the value positions
        
        laM = ones(na)*la;
        Vchange = @. u + dV_upwind*(zmat + r*amat - c) + laM*(V_switch - V) - rho*V;

        Delta = @. 0.9*da/max(z2 + r*agrid);

        v_new = @. v + Delta*Vchange;

        if maximum(abs.(Vchange))< crit
            println("success!")
            return v
        else
            v = v_new
        end
    end
    return c,v
end

param = model()
@unpack na,r,agrid,gamma,rho,zmat,amat = param
v0 = zeros(na,2)
c,v = ss(v0,param)
saving = zmat + r*amat - c;
@unpack amin,amax = param
using Plots
plot(agrid,v, xlims = (amin,amax),xlabel = "a",ylabel = "\$V_i(a)\$")
plot(agrid,saving,xlims = (amin,amax))



