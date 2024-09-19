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

using Plots, LinearAlgebra, QuadGK

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

function plot_coordinates(file)
    x, y = get_coordinates(file)
    plot(x, y, marker = (:circle, 2), line = :solid, aspect_ratio=:equal)
end


function find_midpoints(file)
    # Get current coordinates
    x, y = get_coordinates(file)
    
    # Initialize new variables for midpoints
    x_mid = Float64[]
    y_mid = Float64[]

    # Find midpoints
    for i in range(1, length(x) - 1)
        push!(x_mid, 0.5 * (x[i] + x[i+1]))
        push!(y_mid, 0.5 * (y[i] + y[i+1]))
    end
    x_mid, y_mid
end

# find_midpoints("naca0012.txt")

function plot_midpoints(file)
    x_mid, y_mid = find_midpoints(file)
    scatter!(x_mid, y_mid, marker = (:square, 2))
end

# Find the distance from each midpoint to every other midpoint
# distnace is an Array. Each vector contains the distance from a specific midpoint/panel to every other panel. 

function distance_midpoints(file)
    x_mid, y_mid = find_midpoints(file)
    distance = []
    n = length(x_mid)
    for j in range(1, length(x_mid))
        x_new = x_mid[j]
        y_new = y_mid[j]
        length = []
        for i in 1:n
            pyth = sqrt((x_new + x_mid[i])^2 + (y_new + y_mid[i])^2)
            push!(length, pyth)
        end
        push!(distance, length)
    end
    distance
end

# Make each panel a vector -- simply point 2 - point 1
# Find the length of each panel -- will use this in our integrals

function make_sources(file)
    x, y = get_coordinates(file)
    n = length(x)
    sources = []
    panel_length = []
    panel_unit
    for each in eachindex(x)
        push!(sources, (x[each], y[each]))
    end

    for each in 1:n-1
        panel_x = x[each + 1] - x[each]
        panel_y = y[each + 1] - y[each]
        push!(panel_length, [panel_x, panel_y])
    end

    l = []
    for each in eachindex(panel_length)
        push!(l, sqrt((each[1]^2 + each[2]^2)))
    end
    sources, panel_length, l
end


# Because we chose to model control points at the center of the panels rather than at the center of the surface they are easily computed. The methodology doesn't actually need \theta but rather sin \theta and cos \theta.

function sin_theta(file)
    x, y = get_coordinates(file)
    n = length(n)
    sin_full = []
    cos_full = []
    tan_full = []
    norm_full = []

    for i in 1:n-1
        sin_theta = (y[i + 1] - y[i]) / sqrt((x[i + 1] + x[i])^2 + (y[i + 1] + y[i])^2)
        cos_theta = (x[i + 1] - x[i]) / sqrt((x[i + 1] + x[i])^2 + (y[i + 1] + y[i])^2)
        tan_i = [cos_theta, sin_theta]
        norm_i = [-sin_theta + cos_theta]

        push!(sin_full, sin_theta)
        push!(cos_full, cos_theta)
        push!(tan_full, tan_i)
        push!(norm_full, norm_i)
    end
    sin_full, cos_full, tan_full, norm_full
end

# We can now write our boundary conditions in equation form.
# Flow tangency condition is V dot n_hat = 0


# Freestream Velocity and AoA
V_inf = 00.00
AoA = 0.0 * (180/pi)

# Calculate velocities while inputting boundary conditions such as Kutta condition and no-flow condition...
# Page 66 --- equations 2.195 and 2.196

function velocities()
    
    x_mid, y_mid = find_midpoints(file)
    sources, panel_length, l = make_sources(file)
    
    u_velocity = []
    v_velocity = []

end


#=

\Beta
\thtea k - \theta j

r_kj+1
r_kj


=#


# _______________________________________________________________
# function integral_u(x_int, y_int, l)
#     integral_u = (x_int - x) / ((x_int - x)^2 + y_int^2)
#     result, error = quadgk(integral_u, 0, l)
#     return result
# end

# function integral_v(x_int, y_int, l)
#     integral_v = (y_int) / ((x_int - x)^2 + y_int^2)
#     result, error = quadgk(integral_v, 0, l)
#     return result
# end

# # Calculate velocities!!!

# function velocity(file)
#     x_mid, y_mid = find_midpoints(file)
#     panel_length = make_sources(file)
#     l = []
#     for each in eachindex(panel_length)
#         push!(l, sqrt((each[1]^2 + each[2]^2)))
#     end

#     u = []
#     for each in eachindex(x_mid)
#         int_u = (1/ 2pi) * q(x) * integral_u(x_mid[each], y_mid[each], l[each])
#         push!(u, int_u)
#     end

#     v = []
#     for each in eachindex(x_mid)
#         int_v = (1 / 2pi) * q(x) * integral_v(x_mid[each], y_mid[each], l[each])
#         push!(v, int_v)
#     end
#     u, v
# end

# WHAT IS q(x)????