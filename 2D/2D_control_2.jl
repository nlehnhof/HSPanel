using Plots, Interpolations
import LinearAlgebra.norm as norm

"""
    find_upper_intersection(new_coords::Vector{Tuple{Float64, Float64}}, old_coords::Vector{Tuple{Float64, Float64}}, index::Int)

Finds the point of intersection on the upper surface between the original airfoil geometry and the control surface.
Uses interpolation and returns the point of intersection.

# Arguments:
- `new_coords::Vector{Tuple{Float64, Float64}}` : the coordinates of the airfoil geometry with the deflected control surface
- `old_coords::Vector{Tuple{Float64, Float64}}` : the coordinates of the original airfoil geometry
- `index::Int` : the index of the first hinge point in the original geometry after hinge points have been inserted

# Returns:
- `upper_intersection::Vector{Tuple{Float64, Float64}}` : the coordinate where the original geometry crosses the new geoemtry 
"""
function find_upper_intersection(new_coords, old_coords, index)
    # Extract x and y coordinates from the input pairs (using only the last few indices)
    x_coords1 = first.(new_coords)[length(new_coords) - index - 2:end]
    y_coords1 = last.(new_coords)[length(new_coords) - index - 2:end]
    x_coords2 = first.(old_coords)[length(old_coords) - index:end]
    y_coords2 = last.(old_coords)[length(old_coords) - index:end]

    # Process the new coordinates for interpolation
    sorted_indices = sortperm(x_coords1)
    x_sorted = x_coords1[sorted_indices]
    y_sorted = y_coords1[sorted_indices]

    # Process the old coordinates for interpolation
    sorted_indices2 = sortperm(x_coords2)
    x_sorted2 = x_coords2[sorted_indices2]
    y_sorted2 = y_coords2[sorted_indices2]

    # Create interpolation function for coords1
    interp1 = LinearInterpolation(x_sorted, y_sorted, extrapolation_bc=NaN)
    interp2 = LinearInterpolation(x_sorted2, y_sorted2, extrapolation_bc=NaN)

    # Sets tolerance value
    tol = 1e-4
    upper_intersection = []

    # Loop through x_coords2 range to find intersections
    for i in x_coords1[1]:0.0001:x_coords1[end]
        y1 = interp1(i)
        y2 = interp2(i)
        if !isnan(y1) && !isnan(y2) && abs(y2 - y1) < tol
            # Use push! to store the intersection point as a tuple
            push!(upper_intersection, (i, y1))
        end
    end

    # Find the average of possible intersection points since they will be right next to each other
    if length(upper_intersection) > 1
        x_point = sum(first.(upper_intersection)) / length(upper_intersection)
        y_point = sum(last.(upper_intersection)) / length(upper_intersection)
        empty!(upper_intersection)
        push!(upper_intersection, (x_point, y_point))
    end

    return upper_intersection
end

"""
    find_lower_intersection(new_coords::Vector{Tuple{Float64, Float64}}, old_coords::Vector{Tuple{Float64, Float64}}, index::Int)

Finds the point of intersection on the lower surface between the original airfoil geometry and the control surface.
Uses interpolation and returns the point of intersection.

# Arguments:
- `new_coords::Vector{Tuple{Float64, Float64}}` : the coordinates of the airfoil geometry with the deflected control surface
- `old_coords::Vector{Tuple{Float64, Float64}}` : the coordinates of the original airfoil geometry
- `index::Int` : the index of the first hinge point in the original geometry after hinge points have been inserted

# Returns:
- `upper_intersection::Vector{Tuple{Float64, Float64}}` : the coordinate where the original geometry crosses the new geoemtry 
"""
function find_lower_intersection(new_coords, old_coords, index)
    # Extract x and y coordinates from the input pairs (using only the last few indices)
    x_coords1 = first.(new_coords)[1:index + 2]
    y_coords1 = last.(new_coords)[1:index + 2]
    x_coords2 = first.(old_coords)[1:index+1]
    y_coords2 = last.(old_coords)[1:index+1]

    # Process the x and y coordinates for interpolation
    sorted_indices = sortperm(x_coords1)
    x_sorted = x_coords1[sorted_indices]
    y_sorted = y_coords1[sorted_indices]

    sorted_indices2 = sortperm(x_coords2)
    x_sorted2 = x_coords2[sorted_indices2]
    y_sorted2 = y_coords2[sorted_indices2]

    # Create interpolation function for coords1
    interp1 = LinearInterpolation(x_sorted, y_sorted, extrapolation_bc=NaN)
    interp2 = LinearInterpolation(x_sorted2, y_sorted2, extrapolation_bc=NaN)

    # Set a tolerance value
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
        empty!(lower_intersection)
        push!(lower_intersection, (x_point, y_point))
    end

    return lower_intersection
end

"""
    inscribe_semicircle(p1::Tuple{Float64, Float64}, p2::Tuple{Float64, Float64}, num_points::Int=10)

Function to generate points on a semicircle between two points. In this case, used to find the points along the semicircle between the hinge points.

# Arguments:
- `p1::Tuple{Float64, Float64}` : the first endpoint of the semicircle; the first hinge point of the airfoil
- `p2::Tuple{Float64, Float64}` : the second endpoint of the semicircle; the second hinge point of the airfoil

# Key Word Arguments:
- `num_points::Int=10` : The number of points desired to form the semicircle

# Returns:
- `semicircle_points::Vector{Tuple{Float64, Float64}}` : the coordinates that form the semicircle 
"""
function inscribe_semicircle(p1, p2, num_points=10)
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


"""
    change_deflection_angle(x::Vector, y::Vector, angle_of_control_surface::Float64, percent_of_chord::Float64, percent_of_thickness::Float64=0.0, V_inf::Float64=1.0, alpha::Float64=0.0)

This function takes the original geometry of the airfoil and rotates it about the point of rotation as defined by the user, before updating the geometry.

# Arguments:
- `x::Vector` : the x coordinates of the original geometry
- `y::Vector` : the y coordinates of the original geometry
- `angle_of_control_surface::Float64` : user-defined angle of the control surface
- `percent_of_chord::Float64` : user-defined x-position of the point of rotation given as a decimal between 0.0 and 1.0 where 0.0 is the leading edge and 1.0 the trailing edge. It should be between 0.5 and 1.0 for trailing edge control surfaces.
- `percent_of_thickness::Float64` : user-defined y-position of the point of rotation given as a decimal between -1.0 and 1.0, where 0.0 is at the camber line

# Returns:
- `geometry_rotation::NamedTuple` : Contains relevant points and geometry
"""
function change_deflection_angle(x, y, angle_of_control_surface, percent_of_chord, percent_of_thickness)

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
    
    # y position of the upper and lower surface hinge points
    y_position_of_hinge = abs(interpolation(x_position_of_rotation))
    y_lower_hingeower_hinge = -1 * y_position_of_hinge

    # Find the index of the value that is closest to the position of rotation
    start_index = findmin(abs.(x .- x_position_of_rotation))[2]

    # Use the start_index to insert the hinge points into the geomtry of the airfoil
    if x_position_of_rotation < x[start_index]
        insert!(coordinates, start_index + 1, (x_position_of_rotation, -y_position_of_hinge))
        insert!(coordinates, length(coordinates) - start_index + 1, (x_position_of_rotation, y_position_of_hinge))
        index_low = start_index + 1
        index_high = length(coordinates) -start_index
    else
        insert!(coordinates, start_index, (x_position_of_rotation, -y_position_of_hinge))
        insert!(coordinates, length(coordinates) - start_index + 2, (x_position_of_rotation, y_position_of_hinge))
        index_low = start_index
        index_high = length(coordinates) - start_index + 2
    end

    # Rotation matrix
    matrix_rotation = [cos(angle_of_control_surface) -sin(angle_of_control_surface); sin(angle_of_control_surface) cos(angle_of_control_surface)]

    # Rotated the geometry of the control surface about the point of rotation
    for i in eachindex(coordinates)
        if coordinates[i][1] >= x_position_of_rotation
            coordinate = [coordinates[i][1] - x_position_of_rotation, coordinates[i][2] - y_position_of_rotation]
            new_coordinate = matrix_rotation * coordinate
            coordinates[i] = (new_coordinate[1] + x_position_of_rotation, new_coordinate[2] + y_position_of_rotation)
        end
    end

    # Find the points that make up the semicircle between the two hinge points
    points = inscribe_semicircle(coordinates[index_low], coordinates[index_high], 20)

    # Find the indexes of the hinge points in the new, rotated geoemtry
    x_hinge_rotated = coordinates[start_index][1]
    x_top_hinge_rotated = coordinates[length(coordinates)-start_index][1]
    index_lower_hinge = findfirst(t -> t == x_hinge_rotated, first.(coordinates))
    index_upper_hinge = findfirst(t -> t == x_top_hinge_rotated, first.(coordinates))

    # Find where the new, rotated geometry (coordinates), intersects the old geometry (original_coordinates)
    lower_intersection = find_lower_intersection(original_coordinates, coordinates, start_index)
    upper_intersection = find_upper_intersection(original_coordinates, coordinates, start_index)

    # Insert the hinge points and the intersection points into the rotated geoemtry (coordinates)
    insert!(coordinates, start_index, lower_intersection[1])
    insert!(coordinates, length(coordinates)-start_index+2, upper_intersection[1])
    insert!(coordinates, index_lower_hinge+3, (x_position_of_rotation, y_lower_hinge))
    insert!(coordinates, index_upper_hinge+2, (x_position_of_rotation, y_position_of_hinge))
    
    # Remove values of the rotated geomtry that are inside the airfoil geometry
    # Insert those values of the semicircle that bridge the gap between the rotated hinge point and the original hinge point
    if angle_of_control_surface > 0.0
        points = reverse(points)
        for i in eachindex(points)
            if first.(points)[i] > x_position_of_rotation
                insert!(coordinates, start_index + 2 + i, points[i])
            end
        end
        deleteat!(coordinates, length(coordinates)-start_index)
        deleteat!(coordinates, length(coordinates)-start_index)
    end

    # Remove values of the rotated geomtry that are inside the airfoil geometry
    # Insert those values of the semicircle that bridge the gap between the rotated hinge point and the original hinge point
    if angle_of_control_surface < 0.0
        deleteat!(coordinates, start_index+1)
        deleteat!(coordinates, start_index+1)

        points = reverse(points)
        for i in eachindex(points)
            if first.(points)[i] > x_position_of_rotation
                insert!(coordinates, length(coordinates)-start_index, points[i])
            end
        end
    end

    position_of_rotation = (x_position_of_rotation, y_position_of_rotation)
    upper_hinge = (x_position_of_rotation, y_position_of_hinge)
    lower_hinge = (x_position_of_rotation, y_lower_hinge)

    # Return the NamedTuple geometry_rotation that consists of the original geometry, updated geometry,
    # and relevant points for graphing and understanding, including the point of rotation, hinge points, and the intersection points.
    geometry_rotation = (
        original_coordinates=original_coordinates,      # Original geometry
        coordinates=coordinates,                        # Updated geometry
        position_of_rotation=position_of_rotation,      # Coordinate of the point of rotation
        upper_hinge=upper_hinge,                        # Coordinate of the upper hinge point
        lower_hinge=lower_hinge,                        # Coordinate of the lower hinge point 
        lower_intersection=lower_intersection,          # Intersection point at the lower surface
        upper_intersection=upper_intersection,          # Intersection point at the upper surface  
        points=points                                   # Points of the semicircle with the two hinge points as endpoints 
        )                                  

    return geometry_rotation
end

################# PLOTTING ############################
#=
function plot_geometry(geometry_rotation)
    pl = plot(; aspect_ratio=:equal, color=:blue, legend=:topleft)
    # plot!(geometry_rotation.x, geometry_rotation.y, markers=true, label="Before Rotation")
    plot!(pl, first.(geometry_rotation.coordinates), last.(geometry_rotation.coordinates), markers=true, color=:green, label="Final Rotation")
    # plot!(pl, first.(geometry_rotation.rotated_coordinates), last.(geometry_rotation.rotated_coordinates), markers=true, label="Rotated Geometry")
    # plot!(pl, first.(geometry_rotation.points), last.(geometry_rotation.points), color=:red, label="circle")
    # plot!(pl, first.(geometry_rotation.new_control), last.(geometry_rotation.new_control), color=:purple, label="Sort of Full-Geo")
    # plot!(pl, first.(geometry_rotation.control_surface), last.(geometry_rotation.control_surface), markers=true, color=:orange, label="Only Control Surface")
    # plot!(pl, first.(geometry_rotation.full_geo), last.(geometry_rotation.full_geo), color=:purple, markers=true, label="Final Geometry")
    # scatter!(pl, (geometry_rotation.x_position_of_rotation, geometry_rotation.y_position_of_rotation), label="Point of Rotation", markers=true)
    # scatter!(pl, (geometry_rotation.x_position_of_rotation, geometry_rotation.y_position_of_hinge), color=:red, label = "Point of Hinge")
    # scatter!(pl, (geometry_rotation.x_position_of_rotation, geometry_rotation.lower_hinge), label="lower hinge")
    # scatter!(pl, geometry_rotation.lower_intersection[1], label="lower intersection")
    # scatter!(pl, geometry_rotation.upper_intersection[1], label="upper intersection")
    title!("Full Airfoil Geometry with Control Surface")
    xlabel!(pl, "Normalized Chord")
    ylabel!(pl, "Thickness")
    display(pl)
    savefig(pl, "Full_Geometry_at_negative_20_degrees_Deflection.png")
end

plot_geometry(geometry_rotation)
############################################
=#