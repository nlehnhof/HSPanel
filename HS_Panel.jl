#=

SET UP
1. Input Airfoil coordinates
2. Plot airfoil coordinates
3. Plot airfoil panels using the airfoil coordinates
4. Find midpoint-coordinates for each panel
5. Find distance from each panel to every other panel

MATH
6. Calculate influences - find how each panel vortex influences every other panel vortex
7. Input Boundary Conditions
8. Set up systems of equations (matrix) -- should have N + 1 equations where N = number of panels (one equation for every panel + 1 for Kutta condition)
9. Solve systems of equations
10. Plot results

=#

using Plots, LinearAlgebra

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

function plot_coordinates(x, y)
    plot(x, y, marker = (:circle, 2), line = :solid, aspect_ratio=:equal)
end


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

# find_midpoints("naca0012.txt")

x_mid, y_mid = find_midpoints(x, y)

function plot_midpoints(x_mid, y_mid)
    scatter!(x_mid, y_mid, marker = (:square, 2))
end

# Find the distance from each midpoint to every other midpoint
# distnace is an Array. Each vector contains the distance from a specific midpoint/panel to every other panel. 

function distance_midpoints(x_mid, y_mid)
    n = length(x_mid)
    distance = zeros(Float64, n-1, n)
    for i in 1:n-1
        x_new = x_mid[i]
        y_new = y_mid[i]
        for j in 1:n-1
            distance[i, j] = sqrt((x_new + x_mid[j])^2 + (y_new + y_mid[j])^2) 
        end
    end
    distance
end


# Make each panel a vector -- simply point 2 - point 1
# Find the length of each panel -- will use this in our integrals

function panel_length(x, y)
    n = length(x)
    len = zeros(n-1) 
    for i in 1:n-1
        panel_x = x[i + 1] - x[i]
        panel_y = y[i + 1] - y[i]
        len[i] = sqrt((panel_x^2 + panel_y^2))
    end
    len
end

# Because we chose to model control points at the center of the panels rather than at the center of the surface they are easily computed. The methodology doesn't actually need \theta but rather sin \theta and cos \theta.

function sin_cos_theta(x, y)
    n = length(x)
    sin_theta = zeros(n - 1)
    cos_theta = zeros(n - 1)

    for i in 1:n-1
        sin_theta[i] = (y[i + 1] - y[i]) / sqrt((x[i + 1] + x[i])^2 + (y[i + 1] + y[i])^2)
        cos_theta[i] = (x[i + 1] - x[i]) / sqrt((x[i + 1] + x[i])^2 + (y[i + 1] + y[i])^2)
    end
    sin_theta, cos_theta
end

sin_theta, cos_theta = sin_cos_theta(x, y)


# We can now write our boundary conditions in equation form.
# Flow tangency condition is V dot n_hat = 0

# Calculate velocities while inputting boundary conditions such as Kutta condition and no-flow condition...
# Page 66 --- equations 2.195 and 2.196

# Find rijs...
# change from x_star to x_mid


function find_rijs(x, y, x_mid, y_mid)  # lengths: 131, 131, 130, 130
    n = length(x_mid)
    r_ij = zeros(Float64, n, n)

    for i in 1:n
        x_star = x_mid[i]
        y_star = y_mid[i]
        for j in 1:n
            r_ij[i, j] = sqrt((x_star - x[j])^2 + (y_star - y[j])^2)
        end
    end
    r_ij   # 130 x 130
end

r_ij = find_rijs(x, y, x_mid, y_mid)

# Find sin(theta i - theta j) etc.

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

# sin_theta_ij, cos_theta_ij = find_thetas(x, y, sin_theta, cos_theta)

# Find Beta

function find_beta(x, y, x_mid, y_mid)
    n = length(x)
    beta = zeros(n-1, n-1)
    for i in 1:n-1
        x_bar = x_mid[i]
        y_bar = y_mid[i]
        for j in 1:n-1
            if (j == i) beta[i, j] = π else beta[i, j] = (atan(((x[j] - x_bar) * (y[j+1] - y_bar) - (y[j] - y_bar) * (x[j+1] - x_bar)) , ((x[j] - x_bar) * (x[j+1] - x_bar) + (y[j] - y_bar) * (y[j+1] - y_bar)))) end
        end
    end
    beta   # length 130 x 130
end

# beta = find_beta(x, y, x_mid, y_mid)

# Find Aij
# rijs, sin_cos_theta, beta

function find_A(x, y, r_ij, sin_theta_ij, cos_theta_ij, beta)
    n = length(x) - 1  # 130
    A = zeros(n+1, n+1)   # 131 x 131

    # These for loops give us the values for the 130 x 129 matrix
    for i in 1:n
        val = zeros(130)
        for j in 1:n-1
            A[i, j] = log(r_ij[i, j+1]/r_ij[i, j]) * sin_theta_ij[i] + beta[i] * cos_theta_ij[i]
            val[j] = log(r_ij[i, j+1]/r_ij[i, n-2]) * cos_theta_ij[i] - beta[i] * sin_theta_ij[i]
        end
        A[i, n+1] = sum(val)  # HELP
    end

    # This supposedly gives us the 131st row  # HELP
    for j in 1:n
        A[n+1, j] = beta[j] * sin_theta_ij[j] - log(r_ij[j+1]/r_ij[j]) * cos_theta_ij[j]
    end

    # Gives us (131, 131) value  # HELP
        A[n+1, n+1] = beta[n] * cos_theta_ij[n] + log(r_ij[n+1]/r_ij[n]) * sin_theta_ij[n]

    return A
end

find_A(x, y, r_ij, sin_theta_ij, cos_theta_ij, beta)

# Help: A[:, n]
# I need j to go to n = 130, but the sin_cos_theta vectors restrict it to n-1 = 129


# Freestream Velocity and AoA
V_inf = 10.00
alpha = 5.0 * (180/pi)

function find_b(sin_theta, cos_theta, V_inf, alpha)
    n = length(x)  # 131
    b = zeros(n)

    for i in 1:n-1
        b[i] = 2 * π * V_inf * (sin_theta[i] * sin(alpha) - cos_theta[i] * cos(alpha))
    end
    
    b[n] = -2 * π * V_inf * ((cos_theta[1] * cos(alpha) + sin_theta[1] * sin(alpha)) + (cos_theta[n-1] * cos(alpha) - sin_theta[n-1] * sin(alpha)))

    return b
end

b = find_b(sin_theta, cos_theta, V_inf, alpha)