using Plots
using Interpolations
# import Interpolations.linear_interpolation as lin_int

x = [1.0, 0.75, 0.5, 0.25, 0.0, 0.25, 0.5, 0.75, 1.0]
y = [0.0, -0.125, -0.25, -0.125, 0.0, 0.125, 0.25, 0.125, 0.0]

# Use 2D rotation matrix
function change_deflection_angle(x, y, angle_of_control_surface, percent_of_chord, percent_of_thickness; V_inf=[1.0], alpha=[0.0])

    @assert length(x) == length(y) "length(x) does not equal lenght(y)."    
    x_control_surface = similar(x, length(x)) * 0.0
    y_control_surface = similar(y, length(y)) * 0.0
    angle_of_control_surface = deg2rad(angle_of_control_surface)

    coordinates = [ [xi, yi] for (xi, yi) in zip(x, y) ]

    # Sort points
    sorted_indices = sortperm(x[1:round(Int, length(x)/2)])
    x_sorted = x[sorted_indices]
    y_sorted = y[sorted_indices]

    # Use interpolation to find the x and y coordinates about which to rotate the control surface
    x_position_of_rotation = percent_of_chord * x[1]
    interpolation = LinearInterpolation(x_sorted, y_sorted, extrapolation_bc=NaN)  # NaN for values outside the range
    y_position_of_rotation = abs(interpolation(x_position_of_rotation)) * (percent_of_thickness)
    y_position_of_hinge = abs(interpolation(x_position_of_rotation))

    first = (x[end] * percent_of_chord)
    start_index = findmin(abs.(x .- first))[2]
    end_index = length(x)+2 - start_index

    if x_position_of_rotation > (x[end] * percent_of_chord)
        insert!(coordinates, start_index, [x_position_of_rotation, -y_position_of_hinge])
        insert!(coordinates, end_index, [x_position_of_rotation, y_position_of_hinge])
    else
        insert!(coordinates, start_index+1, [x_position_of_rotation, -y_position_of_hinge])
        insert!(coordinates, end_index, [x_position_of_rotation, y_position_of_hinge])
    end

    # Rotation matrix
    matrix_rotation = [cos(angle_of_control_surface) -sin(angle_of_control_surface); sin(angle_of_control_surface) cos(angle_of_control_surface)]

    for i in eachindex(coordinates)
        if coordinates[i][1] > x_position_of_rotation
            coordinate = [coordinates[i][1] - x_position_of_rotation, coordinates[i][2] - y_position_of_rotation]
            new_coordinate = matrix_rotation * coordinate
            coordinates[i] = [new_coordinate[1] + x_position_of_rotation, new_coordinate[2] + y_position_of_rotation]
        end
    end

    x_control_surface = similar(x, length(coordinates))
    y_control_surface = similar(y, length(coordinates))

    for i in eachindex(coordinates)
        x_control_surface[i] = coordinates[i][1]
        y_control_surface[i] = coordinates[i][2]
    end

    # return x_position_of_rotation, y_position_of_rotation
    return x, y, x_control_surface, y_control_surface, coordinates, x_position_of_rotation, y_position_of_rotation
end

# change_deflection_angle(x, y, 10, 0.72, 0.0)
x, y, x_control_surface, y_control_surface, coordinates, x_position_of_rotation, y_position_of_rotation = change_deflection_angle(x, y, 30, 0.72, -0.5)

function plot_geometry()
    pl = plot(x, y, label = "Before", markers=true, aspect_ratio=1, color=:blue)
    title!("Lower Surface")
    xlabel!(pl, "Normalized Chord")
    ylabel!(pl, "Thickness")
    plot!(pl, x_control_surface, y_control_surface, label="After", markers=true, color=:green)
    scatter!(pl, (x_position_of_rotation, y_position_of_rotation), label="Coordinate of Rotation", markers=true, color=:orange)
    display(pl)
end

plot_geometry()
