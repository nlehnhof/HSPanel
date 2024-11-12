using LinearAlgebra

# Function to generate points on a semicircle between two points
function inscribe_semicircle(p1::Tuple{Float64, Float64}, p2::Tuple{Float64, Float64}, num_points::Int)
    # Calculate the midpoint and radius
    midpoint = [(p1[1] + p2[1]) / 2, (p1[2] + p2[2]) / 2]
    radius = norm([p2[1] - p1[1], p2[2] - p1[2]]) / 2
    
    # Calculate the angle between the points and the x-axis
    angle_offset = atan(p2[2] - p1[2], p2[1] - p1[1])

    # Generate points on the semicircle
    semicircle_points = []
    for i in 1:num_points
        theta = π * (i - 1) / (num_points - 1)  # Angle in radians for each point along the semicircle
        x = midpoint[1] + radius * cos(theta + angle_offset)
        y = midpoint[2] + radius * sin(theta + angle_offset)
        push!(semicircle_points, (x, y))
    end

    return semicircle_points
end

# Example usage
p1 = (1.0, 2.0)
p2 = (4.0, -1.0)
num_points = 10
points = inscribe_semicircle(p1, p2, num_points)

# # Output points
# println("Semicircle Points: ")
# for (i, point) in enumerate(points)
#     println("Point $i: $point")
# end

plot(first.(points), last.(points))