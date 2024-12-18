using Plots

# Read the NACA0008 data from a file
function read_airfoil_data(filename)
    x_points = Float64[]
    y_points = Float64[]
    
    # Open the file and read line by line
    open(filename, "r") do file
        for line in eachline(file)
            # Split each line into x and y values, convert to Float64, and store them
            values = split(line)
            push!(x_points, parse(Float64, values[1]))
            push!(y_points, parse(Float64, values[2]))
        end
    end
    
    return x_points, y_points
end

# Read the data from the file 'naca0008.txt'
x_points, y_points = read_airfoil_data("C:\\Users\\nlehn\\HSPanel\\2D\\NACA0008.txt")

# Plot the airfoil
plot(x_points, y_points, label="NACA2414", xlabel="Chord", ylabel="Thickness", legend=:topright, aspect_ratio=:equal, linewidth=2)
title!("NACA2414 Airfoil")

savefig("NACA2414.png")
