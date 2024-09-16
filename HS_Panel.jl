#=

1. Input Airfoil coordinates
2. Plot airfoil coordinates
3. Plot airfoil panels using the airfoil coordinates
4. Find midpoint-coordinates for each panel

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

function distance_midpoints(x_mid, y_mid)
    distance = Float64[]
    for i in range(1, length(x_mid))
        push!(distance, x[i])


# plot_midpoints("naca0012.txt")