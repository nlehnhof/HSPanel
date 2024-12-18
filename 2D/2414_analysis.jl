using Plots

include("Hess_Smith_Panel_Method.jl")

# Read airfoil coordinates from a file
function extract_x_y(file)
    x = Float64[]  # Initialize empty arrays for x and y
    y = Float64[]

    open(file, "r") do f
        for line in eachline(f)
            entries = split(chomp(line))
            push!(x, parse(Float64, entries[1]))
            push!(y, parse(Float64, entries[2]))
        end
    end
    return x, y
end

# Extract coordinates from a file
x, y = extract_x_y("2D\\joukowsky_coordinates.txt")

# pl1 = scatter(x, y, aspect_ratio=1)
# display(pl1)

# Create the geometric object using the HS method
geo = analyze(x, y, 1.0, 4.0)

# Check if geo is nothing or if it contains valid data
if geo === nothing
    println("HS function returned nothing!")
else
    # println("geo object:", geo)
    # println("geo.x_mid:", geo.x_mid)
    # println("geo.CP:", geo.CP)
    pl = plot(geo.x_mid, geo.CP, yflip=true)
    display(pl)
end