import FLOWFoil.AirfoilTools as at
using Plots

# - Parameters - #
center = [-0.1; 0.1]
radius = 1.0
alpha = 4.0
Vinf = 1.0

# - Joukowsky Geometry - #
x, y = at.joukowsky(center, radius)

# # Open a file to write the coordinates
# open("joukowsky_coordinates.txt", "w") do file
#     for i in eachindex(x)
#         # Write each pair of coordinates to the file
#         println(file, "$(x[i]) $(y[i])")
#     end
# end

# plot(x, y, aspect_ratio=1)