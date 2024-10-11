using Plots

x = [1.0, 0.5, 0.0]
y = [0.0, -0.25, 0.0]

"""
find_midpoints(x::Vector, y::Vector)

    Find the midpoint between each coordinate. 
    Returns the x-coordinate of the midpoint and the y-coordinate of the midpoint.

# Examples
```julia-repl
julia> find_midpoints([1.0, 0.5], [0.0, -0.25])
([0.75],[-0.125])
```
"""
function find_midpoints(x, y)
    # Initialize new variables for midpoints
    n = length(x) - 1
    x_mid = zeros(Float64, n)
    y_mid = zeros(Float64, n)
    # Find midpoints
    for i in eachindex(x_mid)
        x_mid[i] = 0.5 * (x[i] + x[i+1])
        y_mid[i] = 0.5 * (y[i] + y[i+1])
    end
    return x_mid, y_mid
end

# Midpoints : vectors : Control points for each panel are located at the center of each panel -- length(n)
# x_mid, y_mid = find_midpoints(x, y)

"""
    find_sin_cos(x, y)

Computes the sin_theta and cos_theta of each panel with respect to the x-axis.

# Examples
```julia-repl
julia> find_sin_cos([1.0, 0.5], [0.0, -0.25])
([-0.4472], [-0.8944])
```
"""
function find_sin_cos(x, y)
    distance = similar(x_mid, length(x_mid))
    sin_theta = similar(x_mid, length(x_mid))
    cos_theta = similar(x_mid, length(x_mid))
    for i in eachindex(distance)
        distance[i] = sqrt((x[i+1] - x[i])^2 + (y[i+1] - y[i])^2)
        sin_theta[i] = (y[i+1] - y[i]) / distance[i]
        cos_theta[i] = (x[i+1] - x[i]) / distance[i]
        if distance[i] == 0.0
            error("Zero distance detected between points $i and $(i+1).")
        end
    end
    return sin_theta, cos_theta
end

# sin_theta, cos_theta : vectors : The sin and cos of each panel with respect to the x-axis -- length(n)
# sin_theta, cos_theta = find_sin_cos(x, y)

"""
    find_rijs(x, y, x_mid, y_mid)

Computes the distance from each midpoint to each (x, y) coordinate

# Examples
```julia-repl
julia> find_rijs([1.0, 0.5], [0.0, -0.25], [0.75], [-0.125])
[0.2795 0.2795]
```
"""
function find_rijs(x, y, x_mid, y_mid)  # lengths: 131, 131, 130, 130
    r_ij = similar(x_mid, length(x_mid), length(x_mid)+1)
    for i in eachindex(x_mid)
        for j in eachindex(x)
            r_ij[i, j] = sqrt((x_mid[i] - x[j])^2 + (y_mid[i] - y[j])^2)
        end
    end
    return r_ij   # 130 x 131
end

# r_ij : matrix : Find the distance between the control point (mid_point) of panel i with respect to source points j and j+1 -- size(n, n+1)
# r_ij = find_rijs(x, y, x_mid, y_mid)

"""
    find_thetas(sin_theta, cos_theta)

Computes sin(theta1 - theta2) and cos(theta1 - theta2)

# Examples
```julia-repl
julia> find_thetas([0.5, 0.6, 0.7], [0.4, 0.3, 0.2])
([0.0 -0.09 -0.1799; 0.09 0.0 -0.09; 0.1799 0.09 0.0], [0.4100 0.42 0.43; 0.42 0.4499 0.48; 0.43 0.48 0.5299])
```
"""
function find_thetas(sin_theta, cos_theta)
    cos_theta_ij = similar(cos_theta, length(cos_theta), length(cos_theta))
    sin_theta_ij = similar(sin_theta, length(sin_theta), length(sin_theta)) 
    for i in eachindex(sin_theta)
        for j in eachindex(cos_theta)
            sin_theta_ij[i, j] = sin_theta[i] * cos_theta[j] - cos_theta[i] * sin_theta[j]
            cos_theta_ij[i, j] = cos_theta[i] * cos_theta[j] + sin_theta[i] * sin_theta[j]
        end
    end
    return sin_theta_ij, cos_theta_ij   # 130 x 130 matrices
end

# sin_theta_ij, cos_theta_ij : using sin_theta and cos_theta to find sin(theta[i] - theta[j]) and cos(theta[i] - theta[j])
# sin_theta_ij, cos_theta_ij = find_thetas(sin_theta, cos_theta)

"""
    find_beta(x, y, x_mid, y_mid)

Computes the angle from each midpoint to each panel

# Examples
```julia-repl
julia> find_beta([1.0, 0.5, 0.0], [0.0, -0.25, 0.0], [0.75, 0.25], [-0.125, -0.125])
2×2 Matrix{Float64}:
  3.14159   -0.628796
 -0.628796   3.14159
```
"""
function find_beta(x, y, x_mid, y_mid)
    beta = similar(x_mid, length(x_mid),length(x_mid))
    for i in eachindex(x_mid)
        for j in eachindex(x_mid)
            if j == i
                beta[i, j] = π 
            else 
                numerator = (x[j] - x_mid[i]) * (y[j+1] - y_mid[i]) - (y[j] - y_mid[i]) * (x[j+1] - x_mid[i])
                denominator = ((x[j] - x_mid[i]) * (x[j+1] - x_mid[i]) + (y[j] - y_mid[i]) * (y[j+1] - y_mid[i]))
                beta[i, j] = atan(numerator, denominator)
            end
        end
    end
    return beta
end

# beta : matrix : use coordinates and midpoints to find the angle between control point of panel i and source coordinates j and j+1 -- size(n, n)
# beta = find_beta(x, y, x_mid, y_mid)

"""
    find_A(r_ij, sin_theta_ij, cos_theta_ij, beta)

Computes the influence of each panel on every other panel.
See Dr. Ning's textbook: eq. 2.223 and eq. 2.233

# Examples
```julia-repl
julia> find_A([0.279508  0.279508  0.760345; 0.760345  0.279508  0.279508], [0.0 0.7999; -0.7999 0.0], [0.9999 0.6; 0.6 0.9999], [3.14159   -0.628796; -0.628796   3.14159])
3×3 Matrix{Float64}:
 3.14159    0.423314   1.10348
 0.423314   3.14159   -1.10348
 1.10348   -1.10348    3.56491
```
"""
function find_A(r_ij, sin_theta_ij, cos_theta_ij, beta)
    A = similar(r_ij, size(r_ij,1)+1, size(r_ij,2))
    # These for loops give us the values for the [n x n-1] matrix
    for i in eachindex(sin_theta_ij[:, 1])
        for j in eachindex(cos_theta_ij[1, :])
            A[i, j] = log(ℯ, r_ij[i, j+1] / r_ij[i, j]) * sin_theta_ij[i, j] + beta[i, j] * cos_theta_ij[i, j] 
            # check to make sure the value is a real number
            if any(isnan, A[i, j]) || any(isinf, A[i, j])
                println("$i , $j")
                error("A[i ,end] contains NaN or Inf values.")
            end               
        end
    end

    # Gives the values for the nth column 
    for i in eachindex(sin_theta_ij[:, 1])
        for j in eachindex(cos_theta_ij[1, :])
            A[i, end] += log(ℯ, r_ij[i, j+1] / r_ij[i, j]) * cos_theta_ij[i, j] - beta[i, j] * sin_theta_ij[i, j]
            # check to make sure the value is a real number
            if any(isnan, A[i, end]) || any(isinf, A[i, end])
                println("$i , $j")
                println(A[i, end])
                error("A[i ,end] contains NaN or Inf values.")
            end
        end
    end

    # This gives us the n+1 row
    for j in eachindex(A[1,1:end-1])
        sin_theta_k1 = sin_theta_ij[1, j]
        cos_theta_k1 = cos_theta_ij[1, j]
        sin_theta_kn = sin_theta_ij[end, j]
        cos_theta_kn = cos_theta_ij[end, j]
        betak1 = beta[1, j]
        betakn = beta[end, j]
        r1j = r_ij[1, j]
        rnj = r_ij[end, j]
        r1j1 = r_ij[1, j+1]
        rnj1 = r_ij[end, j+1]
    
        k1 = betak1 * sin_theta_k1 - log(ℯ, r1j1 / r1j) * cos_theta_k1
        kn = betakn * sin_theta_kn - log(ℯ, rnj1 / rnj) * cos_theta_kn
        A[end, j] = k1 + kn
        if any(isnan, A[end, j]) || any(isinf, A[end, j])
            error("A[end, j] contains NaN or Inf values.")
        end
    end

    # Gives us the [n+1,n+1] value
    for j in eachindex(A[1, 1:end-1])
        sin_theta_k1 = sin_theta_ij[1, j]
        sin_theta_kn = sin_theta_ij[end, j]
        cos_theta_k1 = cos_theta_ij[1, j]
        cos_theta_kn = cos_theta_ij[end, j]
        betak1 = beta[1, j]
        betakn = beta[end, j]
        r1j = r_ij[1, j]
        rnj = r_ij[end, j]
        r1j1 = r_ij[1, j+1]
        rnj1 = r_ij[end, j+1]
    
        A[end, end] = betakn * cos_theta_kn + log(rnj1 / rnj) * sin_theta_kn + betak1 * cos_theta_k1 + log(r1j1 / r1j) * sin_theta_k1
        
        # check to make sure the value is a real number
        if any(isnan, A[end, end]) || any(isinf, A[end, end])
            error("A[end, end] contains NaN or Inf values.")
        end
    end
    return A 
end

# A : the matrix that is the influence of every panel on every other panel -- size(n+1, n+1)
# A = find_A(r_ij, sin_theta_ij, cos_theta_ij, beta)

"""
    find_b(sin_theta, cos_theta, V_inf, alpha)

Assembles the no-flow-through conditions and the Kutta condition in one matrix.

# Examples
```julia-repl
julia> find_b([0.0 0.7999; -0.7999 0.0], [0.9999 0.6; 0.6 0.9999], 10, 0.0)
3-element Vector{Float64}:
 -28.099258924162903
  28.099258924162903
 112.39703569665161
```
"""
function find_b(sin_theta, cos_theta, V_inf, alpha)
    b = similar(sin_theta, length(sin_theta)+1)
    for i in eachindex(sin_theta)
        b[i] = 2 * π * V_inf * (sin_theta[i] * cos(alpha) - cos_theta[i] * sin(alpha))
    end
    b[end] = -2 * π * V_inf * ((cos_theta[1] * cos(alpha) + sin_theta[1] * sin(alpha)) + (cos_theta[end] * cos(alpha) + sin_theta[end] * sin(alpha)))
    return b
end

# b : vector that provides the boundary conditions -- length(n+1)
# b = find_b(sin_theta, cos_theta, 10, 0.0) 

"""
    find_q_gamma(A, b)

Computes the source strength for each point and the circulation strength for the body.

# Examples
```julia-repl
julia> find_q_gamma([3.14159, 0.423314, 1.10348; 0.423314, 3.14159, -1.10348; 1.10348, -1.10348, 3.56491], [-28.099258924162903, 28.099258924162903, 112.39703569665161])
3-element Vector{Float64}:
 -30.90242799357707
  30.902427993577074
  50.6598146255635
```
"""
function find_q_gamma(A, b) 
    q_gamma = A \ b
    return q_gamma
end

# q_gamma = find_q_gamma(A, b)

"""
    find_vt(r_ij, sin_theta_ij, cos_theta_ij, beta, q_gamma, V_inf, alpha)

Computes the tangential velocity of each panel.

# Examples
```julia-repl
julia> find_vt([0.279508  0.279508  0.760345; 0.760345  0.279508  0.279508], [0.0 0.7999; -0.7999 0.0], [0.9999 0.6; 0.6 0.9999], [3.14159   -0.628796; -0.628796   3.14159],[-30.90242799357707, 30.902427993577074, 50.6598146255635], 10, 0.0)
2-element Vector{Float64}:
 14.37149415156244
 14.371494151562441
```
"""
function find_vt(r_ij, sin_theta_ij, cos_theta_ij, beta, q_gamma, V_inf, alpha)
    Vti = similar(q_gamma, length(q_gamma)-1)
    for i in eachindex(q_gamma[1:end-1])
        set1 = 0.0
        set2 = 0.0
        for j in eachindex(q_gamma[1:end-1])
            set1 += q_gamma[j] * (beta[i, j] * sin_theta_ij[i, j] - log(ℯ, r_ij[i, j+1] / r_ij[i, j]) * cos_theta_ij[i, j])
            set2 += beta[i, j] * cos_theta_ij[i, j] + log(ℯ, r_ij[i, j+1] / r_ij[i, j]) * sin_theta_ij[i, j] 
        end
        Vti[i] = V_inf * cos_theta[i] * cos(alpha) + (set1 / (2*π)) + (q_gamma[end]/(2*π)) * set2
    end
    return Vti
end

# Vti : vector : finds the tangential velocity at each panel since we have a no-flow through condition so there is no normal velocity -- length(n)
# Vti = find_vt(r_ij, sin_theta_ij, cos_theta_ij, beta, q_gamma, 10, 0.0)

"""
    cpressure(Vti, V_inf)

Computes the coefficient of pressure of each panel.

# Examples
```julia-repl
julia> cpressure([14.37149415156244, 14.371494151562441], 10)
2-element Vector{Float64}:
 -1.065398441483934
 -1.0653984414839344
```
"""
function cpressure(Vti, V_inf)
    CP = 1 .- (Vti ./ V_inf) .^2
    return CP 
end

# CP : vector : Coefficient of Pressure at each control point (midpoint of each panel) -- length(n)
# CP = cpressure(Vti, V_inf)

"""
    HS_Panel_CP(x, y, V_inf, alpha)

Computes the coefficient of pressure of each panel using the above equations.
Returns x_mid and CP

# Examples
```julia-repl
julia> HS_Panel_CP([1.0, 0.5, 0.0], [0.0, -0.25, 0.0], 10, 0.0)
([0.75, 0.25], [-1.065398441483934, -1.0653984414839344])
```
"""
function HS_Panel_CP(x, y, V_inf, alpha)
    x_mid, y_mid = find_midpoints(x, y)
    sin_theta, cos_theta = find_sin_cos(x, y)
    r_ij = find_rijs(x, y, x_mid, y_mid)
    sin_theta_ij, cos_theta_ij = find_thetas(sin_theta, cos_theta)
    beta = find_beta(x, y, x_mid, y_mid)
    A = find_A(r_ij, sin_theta_ij, cos_theta_ij, beta)
    b = find_b(sin_theta, cos_theta, V_inf, alpha)
    q_gamma = A \ b
    Vti = find_vt(r_ij, sin_theta_ij, cos_theta_ij, beta, q_gamma, V_inf, alpha)
    CP = cpressure(Vti, V_inf)

    return x_mid, CP
end


########## Validate with Joukowsky #################

import FLOWFoil.AirfoilTools as at

# - Parameters - #
center = [-0.1; 0.1]
radius = 1.0
alpha = 4.0
Vinf = 1.0

# - Joukowsky Geometry - #
x, y = at.joukowsky(center, radius)

# - Surface Values - #
surface_velocity, surface_pressure_coefficient, cl = at.joukowsky_flow(
    center, radius, alpha, Vinf
)

# - Your Stuff - #

alpha = deg2rad(alpha)

x_mid, CP = HS_Panel_CP(x, y, Vinf, alpha)

# - Plot Stuff - #
pl = plot(; xlabel="x", ylabel="cp", yflip=true)
plot!(
    pl,
    x[1:342],
    surface_pressure_coefficient[1:342];
    linestyle=:dash,
    linewidth=2,
    label="Analytic Solution",
)

plot!(pl, x[1:360], CP[1:360], label="Hess-Smith")

display(pl)
# savefig(pl, "Hess_Smith_vs_Analytic_Solution.png")
