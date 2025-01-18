using SparseArrays
function F(k; AL = 0.5, AH = 1.0, α = 1/3, kappa = 2)
    # input: capital k (nk,)
    FL =@. AL*k^α
    FH =@. AH*max(k - kappa,0)^α
    return max(FL, FH)
end

k = range(0, 6,100)
y = F.(k)
begin
    using Plots
    plot(k,y, linewidth = 2)
    y1 =@. 0.5*k^(1/3)
    y2 = @. 1.0*max(k-2, 0)^(1/3)
    plot!(k, y1,linestyle = :dash, linewidth = 2)
    plot!(k,y2,linestyle = :dash, linewidth = 2)
    ylabel!("F(k)")
    xlabel!("k")
end


s = 2
a = 0.3
r = 0.05
d = 0.05
AH = 0.6
AL = 0.4
kappa = 2 #fixed cost

kssH = (a*AH/(r+d))^(1/(1-a)) + kappa
kstar = kappa/(1-(AL/AH)^(1/a))

I = 1000
kmin = 0.001*kssH
kmax = 1.3*kssH
k = range(kmin,kmax,length = I)
dk = (kmax-kmin)/(I-1)

yH =@. AH*max(k - kappa,0)^a
yL =@. AL*k^a
y = max.(yH,yL)

plot(k,y)
plot!(k,yH)
plot!(k,yL)


maxit=1000
crit = 1e-6
Delta = 1000

dVf = zeros(I)
dVb = zeros(I)
c = zeros(I)

v0 =@. -(k^a)^(-1.0)/r
v = v0
for i in 1:maxit
    V = v 
    dVf[1:I-1] = (V[2:I]-V[1:I-1])/dk
    dVf[I] = (y[I] - d*kmax)^(-s)
    #backward difference
    dVb[2:I] = (V[2:I] - V[1:I-1])/dk
    dVb[1] = (y[1] - d*kmin)^(-s) # state constraint, for stability

    cf = dVf.^(-1/s);
    sf =@. y - d*k - cf

    cb = dVf.^(-1/s);
    sb = @. y - d*k - cb

    c0 =@. y - d*k;
    dV0 = c0.^(-s);

    If = cf >0
    Ib = cb <0
    I0 = 1 .- If .- Ib

    dV_upwind = cf*If + cb*Ib + c0*I0
    c = dV_upwind.^(-1/s)

    u =@. c^(1-s)/(1-s)
    X = -Ib*sb/dk
    Y = -If*sf/dk + Ib*sb/dk
    Z = If*sf/dk 

    A = 
    B = (r+1/Delta)*sparse(I,I,I)
    b = (u+V/Delta)
    V = B/b 

    Vchange = V - v 
    v = V 
    if (maximum(abs.(Vchange))<crit)
        println("Value functin converged, Interation = $i")
        break
    end
end