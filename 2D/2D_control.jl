using Plots
using Interpolations

include("intersection2.jl")

x = [1.0, 0.75, 0.5, 0.25, 0.0, 0.25, 0.5, 0.75, 1.0]
y = [0.0, -0.125, -0.25, -0.125, 0.0, 0.125, 0.25, 0.125, 0.0]

# Function to find intersections after divergence
function find_intersections(new_coords, old_coords, index)
    # Extract x and y coordinates from the input pairs (using only the last few indices)
    x_coords1 = first.(new_coords)[index+1:end]
    y_coords1 = last.(new_coords)[index+1:end]
    x_coords2 = first.(old_coords)[index-2:end]
    y_coords2 = last.(old_coords)[index-2:end]

    # Create interpolation function for coords1
    interp1 = LinearInterpolation(x_coords1, y_coords1, extrapolation_bc=NaN)
    interp2 = LinearInterpolation(x_coords2, y_coords2, extrapolation_bc=NaN)

    # Initialize as an empty array to store intersection points

    tol = 1e-5
    intersections = []
    # Loop through x_coords2 range to find intersections
    for x in x_coords1[1]:0.0001:x_coords1[end]
        y1 = interp1(x)
        y2 = interp2(x)
        if !isnan(y1) && !isnan(y2) && abs(y2 - y1) < tol
            # Use push! to store the intersection point as a tuple
            push!(intersections, x)
            push!(intersections, y1)
        end
    end

    return intersections
end

# Use 2D rotation matrix
function change_deflection_angle(x, y, angle_of_control_surface, percent_of_chord, percent_of_thickness=[0.0], V_inf=[1.0], alpha=[0.0])

    @assert length(x) == length(y) "length(x) does not equal lenght(y)."    
    angle_of_control_surface = deg2rad(angle_of_control_surface)
    original_coordinates = [ (xi, yi) for (xi, yi) in zip(x, y) ]
    coordinates = [ (xi, yi) for (xi, yi) in zip(x, y) ]
    
    # Sort points
    sorted_indices = sortperm(x[1:round(Int, length(x)/2)])
    x_sorted = x[sorted_indices]
    y_sorted = y[sorted_indices]

    # Use interpolation to find the x and y coordinates about which to rotate the control surface
    x_position_of_rotation = percent_of_chord * x[1]
    interpolation = LinearInterpolation(x_sorted, y_sorted, extrapolation_bc=NaN)  # NaN for values outside the range
    y_position_of_rotation = abs(interpolation(x_position_of_rotation)) * (percent_of_thickness)
    y_position_of_hinge = abs(interpolation(x_position_of_rotation))

    first1 = (x[end] * percent_of_chord)
    start_index = findmin(abs.(x .- first1))[2]
    end_index = length(x)+2 - start_index

    if x_position_of_rotation > (x[end] * percent_of_chord)
        insert!(coordinates, start_index, (x_position_of_rotation, -y_position_of_hinge))
        insert!(coordinates, end_index, (x_position_of_rotation, y_position_of_hinge))
    else
        insert!(coordinates, start_index+1, (x_position_of_rotation, -y_position_of_hinge))
        insert!(coordinates, end_index, (x_position_of_rotation, y_position_of_hinge))
    end

    # Rotation matrix
    matrix_rotation = [cos(angle_of_control_surface) -sin(angle_of_control_surface); sin(angle_of_control_surface) cos(angle_of_control_surface)]

    for i in eachindex(coordinates)
        if coordinates[i][1] > x_position_of_rotation
            coordinate = [coordinates[i][1] - x_position_of_rotation, coordinates[i][2] - y_position_of_rotation]
            new_coordinate = matrix_rotation * coordinate
            coordinates[i] = (new_coordinate[1] + x_position_of_rotation, new_coordinate[2] + y_position_of_rotation)
        end
    end

    index = findlast(==(x_position_of_rotation), first.(coordinates))
    intersection = find_intersections(coordinates, original_coordinates, index)
    geometry_rotation = ((x=x, y=y, coordinates=coordinates, x_position_of_rotation=x_position_of_rotation, y_position_of_rotation=y_position_of_rotation, original_coordinates=original_coordinates, index=index, intersection=intersection))

    return geometry_rotation
end

geometry_rotation = change_deflection_angle(x, y, 30, 0.72, 0.5)

println(geometry_rotation.intersection)

function plot_geometry()
    pl = plot(geometry_rotation.x[(geometry_rotation.index -2):end], geometry_rotation.y[(geometry_rotation.index - 2):end], label = "Before", markers=true, aspect_ratio=1, color=:blue, legend=:topleft)
    title!("Lower Surface")
    xlabel!(pl, "Normalized Chord")
    ylabel!(pl, "Thickness")
    plot!(pl, first.(geometry_rotation.coordinates)[(geometry_rotation.index):end], last.(geometry_rotation.coordinates)[(geometry_rotation.index):end], label="After", markers=true, color=:green)
    scatter!(pl, (geometry_rotation.intersection[1], geometry_rotation.intersection[2]), label="Point of Intersection")
    scatter!(pl, (geometry_rotation.x_position_of_rotation, geometry_rotation.y_position_of_rotation), label="Point of Rotation")
    plot!(pl, geometry_rotation.x, geometry_rotation.y, label="Original Geometry")
    plot!(pl, first.(geometry_rotation.coordinates), last.(geometry_rotation.coordinates), label = "Full Rotated Geometry")
    
    display(pl)
end

plot_geometry()