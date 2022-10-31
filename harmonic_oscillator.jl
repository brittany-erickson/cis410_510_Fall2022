include("euler_methods.jl")

# Solve m x_tt = -kx - cx + F(t) 
# where F(t) = sin(t) 



# Write forward Euler to solve the linear system IVP:
# y' = Ay + b on 0 ≤ t ≤ Tf
# with initial y0

t0 = 0
y0 = [1; 0]  # initial displacement, initial velocity 

k = 1 # spring stiffness
m = 1 # block mass
c = 0 # damping 

A = [0 1;-k/m -c/m]

Tf = 200
Δt = .001

N = Integer(Tf/Δt) # N+1 total temporal nodes


y = Matrix{Float64}(undef, 2, N+1)
t = Vector{Float64}(undef, N+1)
t[1] = 0
# fill in initial condition:
y[:, 1] = y0

function F(t)
    return 0*sin(t)
end



for n = 1:N  # take N time steps

    b = [0; F(t[n])]

    y[:, n+1] = y[:, n] + Δt*(A*y[:, n] + b)  # Forward Euler
    t[n+1] = t[n] + Δt
end


plot(t, y[1, :])  # plot first component of solution vector 
plot!(t, y[2, :])  # plot second component 

# plot initial condition:
#=
p = plot([y[1, 1]], [y[2, 1]], marker=(:circle,5), color = :blue, legend = false)

for n = 2:N+1
    p = plot!([y[1, n]], [y[2, n]], marker=(:circle,5), color = :blue, legend = false)
    display(p) 
    sleep(.1)
end
=#
