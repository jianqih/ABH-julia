## Version 1.0
## Modified from Ben Moll 
using Parameters
using SparseArrays
using LinearAlgebra

@with_kw struct model
    gamma = 2
    rho = 0.05
    d = 0.05
    al = 1/3
    Aprod = 0.1
    z1 = 1
    z2 = 2*z1
    z = [z1 z2]
    lambda1 = 1/3
    lambda2 = 1/3
    lambda = [lambda1,lambda2]
    z_ave = (z1*lambda2 + z2*lambda1)/(lambda1 + lambda2)
    N= 1000
    amin = 0
    amax = 20
    a = collect(range(amin,amax,N))
    da = (amax-amin)/(N-1)

    aa = [a a]
    zz = ones(N)*z
    maxit= 100
    crit = 10^(-6)
    Delta = 1000
end

param =model()
@unpack N,lambda,rho = param
dVf = zeros(N,2)
dVb = zeros(N,2)
c = zeros(N,2)


Aswitch = [-sparse(I,N,N)*lambda[1] sparse(I,N,N)*lambda[1];sparse(I,N,N)*lambda[2] -sparse(I,N,N)*lambda[2]];


function HJB(v0,w,r,param)
    @unpack crit,Delta, Aprod, d, z_ave,amin,amax,gamma,da,z,a,maxit = param 
    v = v0;
    for n = 1:maxit 
        dVf[1:N-1,:] = (v[2:N,:]-v[1:N-1,:])/da; #forward difference 
        dVf[N,:] =@. (w*z + r*amax)^(-gamma);
        dVb[2:N,:] = (v[2:N,:]-v[1:N-1,:])/da;
        dVb[1,:] =@. (w*z + r*amin)^(-gamma); # state constraint boundary condition
        cf = dVf.^(-1/gamma);
        ssf =@. w*zz + r*aa - cf;
        cb = dVb.^(-1/gamma);
        ssb = @. w*zz + r*aa - cb;
        
        c0 =@. w*zz + r*aa;

        If = ssf .> 0;
        Ib = ssb .< 0;
        I0 = 1 .- If .- Ib;

        c = cf.*If + cb.*Ib + c0.*I0;
        u = c.^(1-gamma)./(1-gamma);

        X = -min.(ssb,0)/da;
        Y = -max.(ssf,0)/da + min.(ssb,0)/da;
        Z = max.(ssf,0)/da;

        # updiag = zeros(2N-1);
        # centdiag = zeros(2N-1);
        # lowdiag = zeros(2N-1);

        # for j = 1:2
        #     centdiag[I*(j-1)+1:I*j] .= Y[:,j]
        #     lowdiag[I*(j-1)+1:I*j - 1] .= X[2:end, j]
        #     updiag[I*(j-1)+1:I*j - 1] .= Z[1:end-1, j]
        # end
        A1=spdiagm(0 => Y[:,1], -1 => X[2:N,1], 1 => Z[1:N-1,1])
        A2=spdiagm(0 => Y[:,2], -1 => X[2:N,2], 1 => Z[1:N-1,2])
        A = [A1 spzeros(N,N);spzeros(N,N) A2] + Aswitch

        # if maximum(abs.(sum.(A)))
        #     println("Improper Transition Matrix")
        # end

        B = (1/Delta + rho)*sparse(I,2N,2N) - A;


        b = vec(u+(1/Delta)*v);
        v_new = reshape(B\b, (N,2));


        Vchange = v_new - v;
        v = v_new;
        dist = maximum(abs.(Vchange));
        if dist < crit 
            println("Value function converged, Iteration = $n")
            return A,v,c
        end
    end
end

function KFE(A;param)
    @unpack N,da,a = param
    AT = A';
    b = zeros(2N,1);

    # fix one value to prevent matrix from being singular
    b[1] = 0.1;
    AT[1,1] = 1
    AT[1,2:end] .= 0

    # solve linear system
    gg = AT\b;
    g_sum = gg'*ones(2N,1)*da;
    gg = gg./g_sum;

    g = [gg[1:N] gg[N+1:2N]]
    
    # get aggregate capital
    K =sum(g'*a*da)
    return K,g
end

@unpack zz,aa,gamma,N,Aprod,al,z_ave,d= param
# Fokker-Planck Equation 
excess_demand = 1

crit_S = 10^(-5);

# initial guess
rmax = 0.049;
r = 0.04;
w = 0.05;

r0 = 0.03;
rmin = 0.01;
rmax = 0.99*rho;
for it in 1:1000
    KD = (al*Aprod/(r + d))^(1/(1-al))*z_ave;
    w = (1-al)*Aprod*KD.^al*z_ave^(-al);
    v0 =@. (w*zz + r*aa).^(1-gamma)/(1-gamma)/rho;
    
    A,v,c = HJB(v0,w,r,param)
    KS,g = KFE(A;param)
    excess_demand = KS - KD;
    println("The excess demand = $excess_demand")
    if excess_demand > crit_S 
        println("Excess supply")
        rmax = r;
        r = 0.5*(r+rmin);
    elseif excess_demand < -crit_S
        println("Excess Demand");
        rmin = r;
        r = 0.5*(r+rmax);
    elseif abs(excess_demand)< crit_S
        println("Equilibrium found, interest rate = $r")
        break
    end
end


# Plotting
begin
    using Plots
    amax1 = 5;
    amin1 = amin-0.1;
    savings =@. w*zz + r.*aa - c
    plot(a,savings[:,1],ylim = (-0.03,0.05))
    plot!(a,savings[:,2],xlim = (amin1,amax1),yaxis = "Savings",xaxis = "Wealth")
end

begin
    plot(a,g[:,1])
    plot!(a,g[:,2],xlim = (amin1,amax1),yaxis = "Densities \$g_i(a)\$",xaxis = "Wealth, \$a\$")
end


