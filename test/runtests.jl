using Arnoldi
using Test

@testset "Arnoldi iteration works" begin
    A = [1.0 -1.0 0.0;
         -1.0 2.0 0.0;
         0.0 0.0 3.0]
    b = [1.0; 2.0; 3.0]

    Q = [0.26726124 -0.77754191 0.56920998;
         0.53452248 -0.37186787 -0.75894664;
         0.80178373 0.50709255 0.31622777]

    @test arnoldi_iteration(A, b, 3) ≈ Q
end
