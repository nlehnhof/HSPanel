using Xfoil, Plots

# Function to read airfoil coordinates from a text file
function read_airfoil_coordinates(file_path)
    x, y = open(file_path, "r") do f
        x = Float64[]  # Initialize empty arrays for x and y
        y = Float64[]  
        for line in eachline(f)
            entries = split(chomp(line))  # Split the line into coordinates
            push!(x, parse(Float64, entries[1]))  # Parse x-coordinate
            push!(y, parse(Float64, entries[2]))  # Parse y-coordinate
        end
        x, y  # Return x and y coordinates
    end
    return x, y
end