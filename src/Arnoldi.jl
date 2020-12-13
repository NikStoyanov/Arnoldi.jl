module Arnoldi

import LinearAlgebra: norm

export arnoldi_iteration

function arnoldi_iteration(A::AbstractMatrix, b::AbstractVector, n::Int; nullvec=1e-5)
    m = size(A, 1)
    h = Matrix{eltype(A)}(undef, n, n - 1)
    Q = Matrix{eltype(A)}(undef, m, n)

    # First Krylov vector.
    q = b ./ norm(b)
    Q[:, 1] = q

    for k in 1:n - 1
        v = A'q
        # Gram–Schmidt process.
        for j in 1:k+1
            h[j, k] = conj(Q[:, j])'v
            v = v - h[j, k] * Q[:, j]
        end

        h[k + 1, k] = norm(v)
        # Check for zero vector.
        if h[k + 1, k] > nullvec
            q = v / h[k + 1, k]
            Q[:, k + 1] = q
        else
            return Q
        end
    end

    return Q
end

end # module
