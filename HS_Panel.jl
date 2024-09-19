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
    n = length(x)
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

# Find rijs...
# change from x_star to x_mid
function find_rijs(file)
    x, y = get_coordinates(file)
    n = length(x)
    r_ij_lst = []
    r_ij1_lst = []
    for i in 1:n-1
        x_star = x[i]
        y_star = y[i]
        for j in 1:n-1
            r_ij = [x_star - x[j], y_star - y[j]]
            r_ij1 = [x_star - x[j+1], y_star - y[j+1]]
            push!(r_ij_lst, r_ij)
            push!(r_ij1_lst, r_ij1)
        end
    end
    r_ij_lst, r_ij1_lst
end

# Find sin(theta i - theta j) etc.

function find_thetas(file)
    sin_full, cos_full, tan_full, norm_full = sin_theta(file)
    n = length(sin_full)
    sin_thetai_j = [] 
    cos_thetai_j = []
    for i in 1:n-1
        for j in 1:n-1
            sin_theta = sin_full[i] * sin_full[j] - cos_full[i] * cos_full[j]
            cos_theta = cos_full[i] * cos_full[j] + sin_full[i] * sin_full[j]
            push!(sin_thetai_j, sin_theta)
            push!(cos_thetai_j, cos_theta)
        end
    end
    sin_thetai_j, cos_thetai_j
end


# Find Beta

function find_beta(file)
    beta_list = []
    x, y = get_coordinates(file)
    x_mid, y_mid = find_midpoints(file)
    n = length(x)
    for i in 1:n
        x_bar = x_mid[1]
        y_bar = y_mid[1]

        for j in 1:n-1
            if j == i
                beta = π
                push!(beta_list, beta)
            else
                beta = (atan(((x[j] - x_bar) * (y[j+1] - y_bar) - (y[j] - y_bar) * (x[j+1] - x_bar)) , ((x[j] - x_bar) * (x[j+1] - x_bar) + (y[j] - y_bar) * (y[j+1] - y_bar))))
                push!(beta_list, beta)
            end
        end
    end
    beta_list
end

print(find_beta("naca0012.txt"))

#=

\Beta -- help!!
\thtea k - \theta j -- help-ish

r_kj+1 -- yes
r_kj -- yes


Remove all push!()
Don't keep reading in the file
Initialize matrix rather than do it twice (i and j)
Trianary Operator

=#