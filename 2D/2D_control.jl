using Plots
using Interpolations

x = [1.0, 0.75, 0.5, 0.25, 0.0, 0.25, 0.5, 0.75, 1.0]
y = [0.0, -0.155, -0.25, -0.155, 0.0, 0.155, 0.25, 0.155, 0.0]

# Function to find intersections after divergence
function find_upper_intersection(new_coords, old_coords, index)
    # Extract x and y coordinates from the input pairs (using only the last few indices)
    x_coords1 = first.(new_coords)[length(new_coords) - index - 2:end]
    y_coords1 = last.(new_coords)[length(new_coords) - index - 2:end]
    x_coords2 = first.(old_coords)[length(old_coords) - index:end]
    y_coords2 = last.(old_coords)[length(old_coords) - index:end]

    sorted_indices = sortperm(x_coords1)
    x_sorted = x_coords1[sorted_indices]
    y_sorted = y_coords1[sorted_indices]

    sorted_indices2 = sortperm(x_coords2)
    x_sorted2 = x_coords2[sorted_indices2]
    y_sorted2 = y_coords2[sorted_indices2]

    # Create interpolation function for coords1
    interp1 = LinearInterpolation(x_sorted, y_sorted, extrapolation_bc=NaN)
    interp2 = LinearInterpolation(x_sorted2, y_sorted2, extrapolation_bc=NaN)

    tol = 1e-4

    upper_intersection = []

    # Loop through x_coords2 range to find intersections
    for i in x_coords1[1]:0.001:x_coords1[end]
        y1 = interp1(i)
        y2 = interp2(i)
        if !isnan(y1) && !isnan(y2) && abs(y2 - y1) < tol
            # Use push! to store the intersection point as a tuple
            push!(upper_intersection, (i, y1))
        end
    end

    # if length(upper_intersection) > 1
    #     x_point = sum(first.(upper_intersection)) / length(upper_intersection)
    #     y_point = sum(last.(upper_intersection)) / length(upper_intersection)
    #     push!(upper_intersection, (x_point, y_point))
    # end

    # upper = upper_intersection[end]
    return upper_intersection
end

function find_lower_intersection(new_coords, old_coords, index)
    # Extract x and y coordinates from the input pairs (using only the last few indices)
    x_coords1 = first.(new_coords)[1:index + 2]
    y_coords1 = last.(new_coords)[1:index + 2]
    x_coords2 = first.(old_coords)[1:index+1]
    y_coords2 = last.(old_coords)[1:index+1]

    sorted_indices = sortperm(x_coords1)
    x_sorted = x_coords1[sorted_indices]
    y_sorted = y_coords1[sorted_indices]

    sorted_indices2 = sortperm(x_coords2)
    x_sorted2 = x_coords2[sorted_indices2]
    y_sorted2 = y_coords2[sorted_indices2]

    # Create interpolation function for coords1
    interp1 = LinearInterpolation(x_sorted, y_sorted, extrapolation_bc=NaN)
    interp2 = LinearInterpolation(x_sorted2, y_sorted2, extrapolation_bc=NaN)

    tol = 1e-4

    lower_intersection = []

    # Loop through x_coords2 range to find intersections
    for i in x_coords1[end]:0.0001:x_coords1[1]
        y1 = interp1(i)
        y2 = interp2(i)
        if !isnan(y1) && !isnan(y2) && abs(y2 - y1) < tol
            # Use push! to store the intersection point as a tuple
            push!(lower_intersection, (i, y2))
        end
    end

    if length(lower_intersection) > 1
        x_point = sum(first.(lower_intersection)) / length(lower_intersection)
        y_point = sum(last.(lower_intersection)) / length(lower_intersection)
        push!(lower_intersection, (x_point, y_point))
    end

    lower = lower_intersection[end]
    return lower
end


# Use 2D rotation matrix
function change_deflection_angle(x, y, angle_of_control_surface, percent_of_chord, percent_of_thickness=[0.0], V_inf=[1.0], alpha=[0.0])

    angle_of_control_surface = deg2rad(angle_of_control_surface)
    original_coordinates = [(xi, yi) for (xi, yi) in zip(x, y)]
    coordinates = [(xi, yi) for (xi, yi) in zip(x, y)]

    # Sort points
    sorted_indices = sortperm(x[1:round(Int, length(x) / 2)])
    x_sorted = x[sorted_indices]
    y_sorted = y[sorted_indices]

    # Use interpolation to find the x and y coordinates about which to rotate the control surface
    x_position_of_rotation = percent_of_chord * x[1]
    interpolation = LinearInterpolation(x_sorted, y_sorted, extrapolation_bc=NaN)  # NaN for values outside the range
    y_position_of_rotation = abs(interpolation(x_position_of_rotation)) * (percent_of_thickness)
    y_position_of_hinge = abs(interpolation(x_position_of_rotation))

    @show start_index = findmin(abs.(x .- x_position_of_rotation))[2]

    if x_position_of_rotation < x[start_index]
        insert!(coordinates, start_index + 1, (x_position_of_rotation, -y_position_of_hinge))
        insert!(coordinates, length(coordinates) - start_index + 1, (x_position_of_rotation, y_position_of_hinge))
    else
        insert!(coordinates, start_index, (x_position_of_rotation, -y_position_of_hinge))
        insert!(coordinates, length(coordinates) - start_index + 2, (x_position_of_rotation, y_position_of_hinge))
    end

    # Rotation matrix
    matrix_rotation = [cos(angle_of_control_surface) -sin(angle_of_control_surface); sin(angle_of_control_surface) cos(angle_of_control_surface)]

    for i in eachindex(coordinates)
        if coordinates[i][1] >= x_position_of_rotation
            coordinate = [coordinates[i][1] - x_position_of_rotation, coordinates[i][2] - y_position_of_rotation]
            new_coordinate = matrix_rotation * coordinate
            coordinates[i] = (new_coordinate[1] + x_position_of_rotation, new_coordinate[2] + y_position_of_rotation)
        end
    end

    upper_intersection = find_upper_intersection(original_coordinates, coordinates, start_index)
    lower_intersection = find_lower_intersection(original_coordinates, coordinates, start_index)

    insert!(coordinates, start_index +1, lower_intersection)
    insert!(coordinates, length(coordinates) - start_index, upper_intersection)    

    # filter!(coordinates -> coordinates[1] >= x_position_of_rotation, coordinates)

    geometry_rotation = ((x=x, y=y, coordinates=coordinates, x_position_of_rotation=x_position_of_rotation, y_position_of_rotation=y_position_of_rotation, original_coordinates=original_coordinates))

    return geometry_rotation
end

# I have found the original geometry, used interpolation to find the hinge points on the upper and lower surfaces, inserted a point at each hingepoint, and have rotated the geometry (including the hinge points) about a point of rotation.
geometry_rotation = change_deflection_angle(x, y, -30, 0.72, 0.0)

function find_intersection(old_coords, rotated_coords, index)
end

function plot_geometry()
    pl = plot(geometry_rotation.x, geometry_rotation.y, label="Before Rotation", markers=true, aspect_ratio=1, color=:blue, legend=:topright)
    plot!(pl, first.(geometry_rotation.coordinates), last.(geometry_rotation.coordinates), label="After Rotation", markers=true, color=:green)
    # scatter!(pl, (geometry_rotation.intersection[1], geometry_rotation.intersection[2]), label="Point of Intersection")
    scatter!(pl, (geometry_rotation.x_position_of_rotation, geometry_rotation.y_position_of_rotation), label="Point of Rotation", markers=true)
    title!("Rotated Geometry")
    xlabel!(pl, "Normalized Chord")
    ylabel!(pl, "Thickness")
    display(pl)
end

plot_geometry()