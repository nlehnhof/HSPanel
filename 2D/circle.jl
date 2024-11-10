using Plots


# Function to generate an arc between two points using the hinge point
function arc_coordinates(x_position_of_rotation, y_position_of_rotation, x_intersection, y_intersection, radius, num=10)
    height = abs(y_position_of_rotation - y_intersection)
    width = abs(x_position_of_rotation - x_intersection)
   
    radius = sqrt(height^2 + width^2)

    if y_intersection < 0
        start_angle = -pi/2
        stop_angle = -atan(height, width)
    else
        start_angle = pi/2
        stop_angle = atan(height, width)
    end
    
    angles = range(start_angle, stop_angle, length=num)

    x_coords = x_position_of_rotation .+ radius * cos.(angles)
    y_coords = y_position_of_rotation .+ radius * sin.(angles)
   
    arc_points = [(x, y) for (x, y) in zip(x_coords, y_coords)]

    return arc_points
end

arc = arc_coordinates(0.75, 0.0, 0.85, 0.1, 10)

plot(first.(arc), last.(arc))