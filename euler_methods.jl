using Plots

# write forward Euler for the IVP system y' = f(t, y)
# where y is a vector in R^n (i.e. has n components)


function my_forward_euler_linear(t0, Tf, Δt, y0, A, b)
    
    # y0 has N components 
    N = size(y0)

    M = Integer(Tf/Δt)  # M+1 total temporal nodes

    t = Vector{Float64}(undef, M+1)
    y = Matrix{Float64}(undef, N, M+1)

    # fill in the initial condition:
    t[1] = t0
    y[:, 1] = y0

    for n = 1:M # take M time steps
        y[:, n+1] = y[:, n] + Δt*(A*y[:, n] + b)
        t[n+1] = t[n] + Δt
    end
    
    return (t, y)
end

function my_forward_euler(t0, Tf, Δt, y0, f)
    
    # y0 has N components 
    N = size(y0)

    M = Integer(Tf/Δt)  # M+1 total temporal nodes

    t = Vector{Float64}(undef, M+1)
    y = Matrix{Float64}(undef, N, M+1)

    # fill in the initial condition:
    t[1] = t0
    y[:, 1] = y0

    for n = 1:M # take N time steps
        y[:, n+1] = y[:, n] + Δt*f(t[n],y[:, n])
        t[n+1] = t[n] + Δt
    end
    
    return (t, y)
end


#=
# Write forward Euler to solve the linear system IVP:
# y' = Ay + b on 0 ≤ t ≤ Tf
# with initial y0

# Do this on particular example, where A = [-3 13; -5 -1]
# y0 = [3; -10]

t0 = 0
y0 = [3; -10]

A = [-3 13;-5 -1]

Tf = 10
Δt = 1

N = Integer(Tf/Δt) # N+1 total temporal nodes


y = Matrix{Float64}(undef, 2, N+1)
t = Vector{Float64}(undef, N+1)
t[1] = 0
# fill in initial condition:
y[:, 1] = y0


for n = 1:N  # take N time steps

    y[:, n+1] = y[:, n] + Δt*A*y[:, n]  # Forward Euler
    t[n+1] = t[n] + Δt
end


#plot(t, y[1, :])  # plot first component of solution vector 
#plot!(t, y[2, :])  # plot second component 

# plot initial condition:
p = plot([y[1, 1]], [y[2, 1]], marker=(:circle,5), color = :blue, legend = false)

for n = 2:N+1
    p = plot!([y[1, n]], [y[2, n]], marker=(:circle,5), color = :blue, legend = false)
    display(p) 
    sleep(1)
end
=#

