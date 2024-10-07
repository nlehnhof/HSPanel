using Test
include("HS_Panel.jl")


x = [1.0, 0.5, 0.0, 0.5, 1.0]
y = [0.0, -0.25, 0.0, 0.25, 0.0] 

@testset "find midpoints" begin
    x_mid, y_mid = find_midpoints(x, y)
    @test x_mid == [0.75, 0.25, 0.25, 0.75] 
    @test y_mid == [-0.125, -0.125, 0.125, 0.125]
end;


@testset "sin_theta and cos_theta" begin
    sin_theta, cos_theta = find_length(x, y)
    @test sin_theta == [-0.4472135954999579, 0.4472135954999579, 0.4472135954999579, -0.4472135954999579]
    @test cos_theta == [-0.8944271909999159, -0.8944271909999159, 0.8944271909999159, 0.8944271909999159]
end

@testset "r_ijs" begin
    r_ij = find_rijs(x, y, x_mid, y_mid)
    @test size(r_ij) == (4, 5)
    @test isapprox(r_ij[1, 1], 0.279508, rtol=1e-5)
    @test isapprox(r_ij[1, end], 0.279508, rtol=1e-5)
    @test isapprox(r_ij[2, 1], 0.760345, rtol=1e-5) 
    @test isapprox(r_ij[2, end], 0.760345, rtol=1e-5) 
    @test isapprox(r_ij[3, 1], 0.760345, rtol=1e-5) 
    @test isapprox(r_ij[3, end], 0.760345, rtol=1e-5) 
    @test isapprox(r_ij[4, 1], 0.279508, rtol=1e-5) 
    @test isapprox(r_ij[4, end], 0.279508, rtol=1e-5) 
end

@testset "sin_theta_ij cos_theta_ij" begin
    sin_theta_ij, cos_theta_ij = find_thetas(sin_theta, cos_theta)
    @test size(sin_theta_ij) == (4, 4)
    @test size(cos_theta_ij) == (4, 4)
    @test isapprox(sin_theta_ij[1, 1], -0.6, rtol=1e-5)
    @test isapprox(cos_theta_ij[1 , 1], 0.999999, rtol=1e-5)
    @test isapprox(sin_theta_ij[1, end], 0.999999, rtol=1e-5)
    @test isapprox(cos_theta_ij[1 , end], -0.6, rtol=1e-5)
    @test isapprox(sin_theta_ij[end, end], -0.6, rtol=1e-5)
    @test isapprox(cos_theta_ij[end, end], 0.999999, rtol=1e-5)
    @test isapprox(sin_theta_ij[end, 1], 0.999999, rtol=1e-5)
    @test isapprox(cos_theta_ij[end, 1], -0.6, rtol=1e-5)
end

@testset "find beta" begin
    beta = find_beta(x, y, x_mid, y_mid)
    @test size(beta) == (4, 4)
    @test isapprox(beta[1, 1], π , rtol=1e-5)
    @test isapprox(beta[1, end], -1.6951513213416582, rtol=1e-5)
    @test isapprox(beta[end, end], π, rtol=1e-5)
    @test isapprox(beta[end, 1], -1.6951513213416582, rtol=1e-5)
end

@testset "find A" begin
    A = find_A(r_ij, sin_theta_ij, cos_theta_ij, beta)
    @test size(A) == (5, 5)
    @test isapprox(A[1, 1], π, rtol=1e-5)
    @test isapprox(A[1, end-1], 0.5393350702912767, rtol=1e-10)
    @test isapprox(A[1, end], 4.851979365742334, rtol=1e-5)
    @test isapprox(A[end, 1], -3.2934534799873028, rtol=1e-5)
    @test isapprox(A[end, end-1], -3.8667603470037646, rtol=1e-5)
    @test isapprox(A[end, end], 9.19810144075646, rtol=1e-5)
end

@testset "find b" begin
    b = find_b(sin_theta, cos_theta, V_inf, alpha)
    @test length(b) == 5
    @test isapprox(b[1], 56.19851784832581, rtol=1e-5)
    @test isapprox(b[end], 0.0, rtol=1e-5)
end

@testset "find b" begin
    @test length(q_gamma) == 5
    @test isapprox(q_gamma[1], 60.779075233301334, rtol=1e-5)
    @test isapprox(q_gamma[end-1], -54.3044382858514, rtol=1e-5)
    @test isapprox(q_gamma[end], -1.1032507986492309e-15, rtol=1e-25)
end

@testset "find Vt" begin
    Vti = find_vt(r_ij, sin_theta_ij, cos_theta_ij, beta, q_gamma, V_inf, alpha)
    @test length(Vti) == 4
    @test isapprox(Vti[1], -2.269190481726291e-15, rtol=1e-10)
    @test isapprox(Vti[end], 8.944271909999149, rtol=1e-10)
end

@testset "find CP" begin
    CP = cpressure(Vti)
    @test length(CP) == 4
    @test isapprox(CP[1], 1.0, rtol=1e-10)
    @test isapprox(CP[end], 0.2, rtol=1e-10 )
end