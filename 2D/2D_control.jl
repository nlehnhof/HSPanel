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

    # Sort points
    sorted_indices = sortperm(x[1:round(Int, length(x)/2+1)])
    x_sorted = x[sorted_indices]
    y_sorted = y[sorted_indices]

    # Use interpolation to find the x and y coordinates about which to rotate the control surface
    x_position_of_rotation = percent_of_chord * x[1]
    interpolation = LinearInterpolation(x_sorted, y_sorted, extrapolation_bc=NaN)  # NaN for values outside the range
    y_position_of_rotation = abs(interpolation(x_position_of_rotation)) * (percent_of_thickness)

    first = (x[end] * percent_of_chord)
    start_index = findmin(abs.(x .- first))[2]
    
    # Rotation matrix
    matrix_rotation = [cos(angle_of_control_surface) -sin(angle_of_control_surface); sin(angle_of_control_surface) cos(angle_of_control_surface)]

    for i in eachindex(x_control_surface)
        if x[i] >= x_position_of_rotation
            coordinate = [x[i] - x_position_of_rotation, y[i] - y_position_of_rotation]
            new_coordinate = matrix_rotation * coordinate
            x_control_surface[i] = new_coordinate[1] + x_position_of_rotation
            y_control_surface[i] = new_coordinate[2] + y_position_of_rotation
        else
            x_control_surface[i] = x[i]
            y_control_surface[i] = y[i]
        end
    end

    println(x_control_surface)
    println(y_control_surface)

    # if x_position_of_rotation > (x[end] * percent_of_chord)
    #     splice!(x_control_surface, start_index:start_index+2, x_position_of_rotation)
    #     splice!(y_control_surface, start_index:start_index+2, y_position_of_rotation)
    # elseif x_position_of_rotation == (x[end] * percent_of_chord)
    #     splice!(x_control_surface, start_index-1:start_index+1, x_position_of_rotation)
    #     splice!(y_control_surface, start_index-1:start_index+1, y_position_of_rotation)
    # else
    #     splice!(x_control_surface, start_index-2:start_index, x_position_of_rotation)
    #     splice!(y_control_surface, start_index-2:start_index, y_position_of_rotation)
    # end

    return x, y, x_control_surface, y_control_surface, x_position_of_rotation, y_position_of_rotation
end

# change_deflection_angle(x, y, 10, 0.72, -0.5)
x, y, x_control_surface, y_control_surface, x_position_of_rotation, y_position_of_rotation = change_deflection_angle(x, y, 45, 0.72, 0.0)

# println(x_control_surface)
# println(y_control_surface)

function plot_geometry()
    pl = plot(x, y, label = "Before", markers=true, aspect_ratio=1)
    title!("Lower Surface")
    scatter!(pl, (x_position_of_rotation, y_position_of_rotation), label="Coordinate of Rotation")
    plot!(pl, x_control_surface, y_control_surface, label="After", markers=true)
    display(pl)
end

plot_geometry()
