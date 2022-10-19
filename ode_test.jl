# solve the IVP y' = -3y, with initial condition y(0) = 17 on 
# the time domain 0 ≤ t ≤ Tf using Forward Euler
# with exact solution y(t) = y0*exp(-3t)
using Plots


function my_forward_euler(t0, Tf, Δt, y0, f, λ)
    N = Integer(Tf/Δt)  # N+1 total temporal nodes

    t = Vector{Float64}(undef, N+1)
    y = Vector{Float64}(undef, N+1)

    # fill in the initial condition:
    t[1] = t0
    y[1] = y0

    for n = 1:N # take N time steps
        y[n+1] = y[n] + Δt*f(t[n],y[n])
        t[n+1] = t[n] + Δt
    end
    
    return (t, y)
end


function test_equationRHS(t, y)
    return λ*y
end

function my_new_RHS(t, y)
    return cos(t)
end

Tf = 10
Δt = 0.1*2/3   #2/abs(λ) is limit of stability (λ = -3)
y0 = 17
λ = -3

(T, Y) = my_forward_euler(0, Tf, Δt, y0, my_new_RHS)

scatter(T, Y, label = "approx", shape = :circle, color = :green)

tfine = 0:Δt/20:Tf
#yexact = y0*exp.(λ*tfine)
yexact = sin.(tfine) .+ y0

plot!(tfine, yexact, label = "exact", color = :red)