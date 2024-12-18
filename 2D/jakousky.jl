using FLOWFoil
using Plots
using .AirfoilTools

# - Parameters - #
center = [-0.1; 0.1]
radius = 1.0

# - Joukowsky Geometry - #
x, y = FLOWFoil.AirfoilTools.joukowsky(center, radius)

pl = plot(x, y, aspect_ratio=1, label="Joukousky Airfoil")
xlabel!(pl, "Chord")
ylabel!(pl, "Thickness")
title!(pl, "Joukousky Geometry")

savefig(pl, "joukousky_geometry.png")

# # Define the output file name
# output_file = "joukowsky_airfoil.txt"

# # Open the file for writing
# open(output_file, "w") do file    
#     # Write each pair of x and y coordinates as columns
#     for i in eachindex(x)
#         println(file, "$(x[i]) $(y[i])")
#     end
# end

# println("Joukowsky airfoil coordinates have been written to $output_file")