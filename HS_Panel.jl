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

using Plots, Printf

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


# distance_midpoints("naca0012.txt")