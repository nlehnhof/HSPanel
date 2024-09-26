#=

Hess-Smith panel method is a numberical technique used in computational
fluid dynamics to analyze fluid flow about a body. In this example, we use an airfoil (NACA0012).

SET UP
1. Input Airfoil coordinates
2. Plot airfoil coordinates
3. Plot airfoil panels using the airfoil coordinates
4. Find midpoint-coordinates for each panel

MATH
5. Calculate influences - find how each panel vortex influences every other panel vortex
6. Input Boundary Conditions
7. Set up systems of equations (matrix) -- should have N + 1 equations where N = number of panels (one equation for every panel + 1 for Kutta condition)
8. Solve systems of equations to find the source strength of each point and our circulation strenght (gamma)
9. Find the tangential velocity at each panel
10. Find the Coefficient of Pressure at each point (x_mid, y_mid)

=#

using Plots, LinearAlgebra

# Freestream Velocity and AoA
alpha = 3 * pi/180
V_inf = 10.0

function get_coordinates(file)
    x, y = open(file, "r") do f
        x = Float64[]
        y = Float64[]
        for line in eachline(f)
            entries = split(chomp(line))
            push!(x, parse(Float64, entries[1]))
            push!(y, parse(Float64, entries[2]))
        end
        x, y
    end
end

x, y = get_coordinates("naca0012.txt")

function find_midpoints(x, y)
    # Initialize new variables for midpoints
    n = length(x)
    x_mid = zeros(n-1)
    y_mid = zeros(n-1)

    # Find midpoints
    for i in 1:n-1
        x_mid[i] = 0.5 * (x[i] + x[i+1])
        y_mid[i] = 0.5 * (y[i] + y[i+1])
    end
    x_mid, y_mid
end

x_mid, y_mid = find_midpoints(x, y)

function plot_coordinates()
    plot(x, y, marker = (:circle, 2), line = :solid, aspect_ratio=:equal, label = "coordinates", title = "NACA 0012")
    scatter!(x_mid, y_mid, marker = (:square, 1), label = "midpoints")
end

function find_length(x, y)
    n = length(x) - 1
    l = Vector{}(undef, n)
    distance = zeros(n)
    sin_theta = zeros(n)
    cos_theta = zeros(n)
    for i in 1:n
        l[i] = [x[i+1] - x[i], y[i+1] - y[i]]
        distance[i] = sqrt((x[i+1] - x[i])^2 + (y[i+1] - y[i])^2)
        sin_theta[i] = (y[i+1] - y[i]) / distance[i]
        cos_theta[i] = (x[i+1] - x[i]) / distance[i]
    end
    l, distance, sin_theta, cos_theta
end

l, distance, sin_theta, cos_theta = find_length(x, y)

# Find the distance from each midpoint to every other midpoint
# Because we chose to model control points at the center of the panels rather than at the center of the surface they are easily computed. The methodology doesn't actually need \theta but rather sin \theta and cos \theta.
# All works up to this point and has been confirmed

# We can now write our boundary conditions in equation form.
# Flow tangency condition is V dot n_hat = 0

# Calculate velocities while inputting boundary conditions such as Kutta condition and no-flow condition...

# Find rijs...

function find_rijs(x, y, x_mid, y_mid)  # lengths: 131, 131, 130, 130
    n = length(x)
    r_ij = zeros(Float64, n-1, n)

    for i in 1:n-1
        for j in 1:n
            r_ij[i, j] = sqrt((x_mid[i] - x[j])^2 + (y_mid[i] - y[j])^2)
        end
    end
    r_ij   # 130 x 131
end

r_ij = find_rijs(x, y, x_mid, y_mid)

# # Find sin(theta i - theta j) etc.

function find_thetas(x, y, sin_theta, cos_theta)
    n = length(x)
    sin_theta_ij = zeros(n-1, n-1) 
    cos_theta_ij = zeros(n-1, n-1)
    for i in 1:n-1
        for j in 1:n-1
            sin_theta_ij[i, j] = sin_theta[i] * sin_theta[j] - cos_theta[i] * cos_theta[j]
            cos_theta_ij[i, j] = cos_theta[i] * cos_theta[j] + sin_theta[i] * sin_theta[j]
        end
    end
    sin_theta_ij, cos_theta_ij   # 130 x 130 matrices
end

sin_theta_ij, cos_theta_ij = find_thetas(x, y, sin_theta, cos_theta)

# Find Beta

function find_beta(x, y, x_mid, y_mid)
    n = length(x)
    beta = zeros(n-1, n-1)
    for i in 1:n-1
        for j in 1:n-1
            if (j == i) beta[i, j] = π else beta[i, j] = (atan(((x[j] - x_mid[i]) * (y[j+1] - y_mid[i]) - (y[j] - y_mid[i]) * (x[j+1] - x_mid[i])) , ((x[j] - x_mid[i]) * (x[j+1] - x_mid[i]) + (y[j] - y_mid[i]) * (y[j+1] - y_mid[i])))) end
        end
    end
    beta   # 130 x 130
end

beta = find_beta(x, y, x_mid, y_mid)  # CORRECT

# Find Aij

function find_A(x, y, x_mid, y_mid, r_ij, sin_theta_ij, cos_theta_ij, beta)
    n = length(x) - 1  # 130
    A = zeros(n+1, n+1)   # 131 x 131

    # These for loops give us the values for the 130 x 129 matrix
    val = zeros(n, 130)
    for i in 1:n
        for j in 1:n
            A[i, j] = log(r_ij[i, j+1]/r_ij[i, j]) * sin_theta_ij[i, j] + beta[i, j] * cos_theta_ij[i, j]
            val[i, j] = log(r_ij[i, j+1]/r_ij[i, j]) * cos_theta_ij[i, j] - beta[i, j] * sin_theta_ij[i, j]
        end
        A[i, n+1] = sum(val[i, :]) 
    end

    # This supposedly gives us the 131st row  (A_N+1,j)
    for j in 1:n
        sin_theta_k1 = sin_theta_ij[1, j]
        cos_theta_k1 = cos_theta_ij[1, j]
        sin_theta_kn = sin_theta_ij[n, j]
        cos_theta_kn = cos_theta_ij[n, j]

        betak1 = beta[1, j]
        betakn = beta[n, j]

        r1j = r_ij[1, j]
        rnj = r_ij[n, j]
        r1j1 = r_ij[1, j+1]
        rnj1 = r_ij[1, j+1]
    
        k1 = betak1 * sin_theta_k1 - log(r1j1 / r1j) * cos_theta_k1
        kn = betakn * sin_theta_kn - log(rnj1 / rnj) * cos_theta_kn

        A[n+1, j] = k1 + kn
    end

    # Gives us (A_n+1,n+1)
    next = zeros(2, 130)
    for j in 1:n
        sin_theta_k1 = sin_theta_ij[1, j]
        cos_theta_k1 = cos_theta_ij[1, j]
        sin_theta_kn = sin_theta_ij[n, j]
        cos_theta_kn = cos_theta_ij[n, j]

        betak1 = beta[1, j]
        betakn = beta[n, j]

        r1j = r_ij[1, j]
        rnj = r_ij[n, j]
        r1j1 = r_ij[1, j+1]
        rnj1 = r_ij[1, j+1]
    
        next[1, j] = betak1 * cos_theta_k1 + log(r1j1 / r1j) * sin_theta_k1
        next[2, j] = betakn * cos_theta_kn + log(rnj1 / rnj) * sin_theta_kn
    end

    k1 = sum(next[1, :])
    kn = sum(next[2, :])

    A[n+1, n+1] = k1 + kn

    return A
end

A = find_A(x, y, x_mid, y_mid, r_ij, sin_theta_ij, cos_theta_ij, beta)

# Find b vector (boundary conditions)

function find_b(sin_theta, cos_theta, V_inf, alpha)
    n = length(x) - 1  # 131
    b = zeros(n + 1)

    for i in 1:n
        b[i] = 2 * π * V_inf * (sin_theta[i] * sin(alpha) - cos_theta[i] * cos(alpha))
    end

    b[n+1] = -2 * π * V_inf * ((cos_theta[1] * cos(alpha) + sin_theta[1] * sin(alpha)) + (cos_theta[n] * cos(alpha) + sin_theta[n] * sin(alpha)))

    return b
end

b = find_b(sin_theta, cos_theta, V_inf, alpha) # 131

q_gamma = A \ b

# Find tangential velocity at each panel

function find_vt(x, r_ij, sin_theta_ij, cos_theta_ij, beta, q_gamma, V_inf, alpha)
    n = length(x) - 1
    Vti = zeros(n)
    set1 = zeros(n+1, n+1)
    set2 = zeros(n+1, n+1)

    for i in 1:n
        for j in 1:n
            set1[i, j] = q_gamma[j] * (beta[i, j] * sin_theta_ij[i, j] - log(r_ij[i, j+1] / r_ij[i, j]) * cos_theta_ij[i, j])
            set2[i, j] = beta[i, j] * cos_theta_ij[i, j] + log(r_ij[i, j+1] / r_ij[i, j]) * sin_theta_ij[i, j] 
        end
        Vti[i] = V_inf * (cos_theta[i] * cos(alpha) + sin_theta[i] * sin(alpha)) + (1/(2π)) * sum(set1[i,:]) + (q_gamma[end]/(2π)) * sum(set2[i,:])
    end
    Vti
end

Vti = find_vt(x, r_ij, sin_theta_ij, cos_theta_ij, beta, q_gamma, V_inf, alpha)

# print(Vti)

# Find CP for each point

function cpressure(Vti)
    CP = 1 .- (Vti ./ V_inf) .^2
end

CP = cpressure(Vti)


# We have a problem. CP is mostly negative, but should be half negative and half positive