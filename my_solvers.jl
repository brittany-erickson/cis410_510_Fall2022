using Plots # add Plots.jl from the package manager if you have not already done so.

# HW 1 starting script (if you want): contains function computeLU() - to compute an LU-factorization of square
# matrix A, namely, A = LU, where L and U are lower and upper triangular matrices. 


"""
    computeLU(A)
Compute and return LU factorization `LU = A` of square matrix `A`.  
Might not work on all matrices, since no pivoting is done!

# Examples (don't need examples, but fine to include)
'''
julia> A = [6 -2 2;12 -8 6;3 -13 3]
3×3 Array{Int64,2}:
  6   -2  2
 12   -8  6
  3  -13  3
julia> (L, U) = computeLU(A)
([1.0 0.0 0.0; 2.0 1.0 0.0; 0.5 3.0 1.0], [6.0 -2.0 2.0; 0.0 -4.0 2.0; 0.0 0.0 -4.0])
julia> norm(A - L*U)
0.0
'''
"""
function computeLU(A)

    N = size(A)[1]

    Id = Matrix{Float64}(I, N, N) # N x N identity matrix

    L = copy(Id)   # initialize
    U = copy(Id)   # initialize
    Ã  = copy(A) # initialize. Ã corresponds to A as it goes under elimination stages

    for k = 1:N-1 # march across columns

        (Lk, Lk_inv) = compute_Lk(Ã, k)

        Ã .= Lk * Ã
        L .= L * Lk_inv

    end

    U .= Ã

    return (L, U)

end


"""
    compute_Lk(A, k)
Compute Lk and its inverse from A, assuming first k-1 columns have undergone elimination.

"""
function compute_Lk(A, k)

    
    N = size(A)[1]

    Lk = Matrix{Float64}(I, N, N)       # initialize as identity matrix
    Lk_inv = Matrix{Float64}(I, N, N)   # initialize as identity matrix

    # now modify column k, strictly below diagonal (i = k+1:N)
    for i = k+1:N
        Lk[i,k] = ???       # fill me in (compute elimination factors)
        Lk_inv[i,k] = ???   # fill me in (compute elimination factors)
    end
 
    return (Lk, Lk_inv)

end

A = Matrix{Float64}(undef, 3, 3)
A .= [6 -2 2;12 -8 6;3 -13 3]
b = rand(3, 1)

(L, U) = computeLU(A)
@assert L*U ≈ A