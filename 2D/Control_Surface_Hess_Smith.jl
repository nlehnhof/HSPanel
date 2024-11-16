using Plots, Interpolations
import LinearAlgebra.norm as norm

"""
    find_intersection(new_coords::Vector{Tuple{Float64, Float64}}, old_coords::Vector{Tuple{Float64, Float64}}, index::Int, surface::String="upper" || "lower")

Finds the point of intersection on the lower or upper surface between the original airfoil geometry and the control surface.
Uses interpolation and returns the point of intersection.

# Arguments:
- `new_coords::Vector{Tuple{Float64, Float64}}` : the coordinates of the airfoil geometry with the deflected control surface
- `old_coords::Vector{Tuple{Float64, Float64}}` : the coordinates of the original airfoil geometry
- `index::Int` : the index of the first hinge point in the original geometry after hinge points have been inserted
- `surface::String="upper" || "lower"` : specifices whether to find the intersection point of the upper or lower surface

# Returns:
- `intersection_points::Vector{Tuple{Float64, Float64}}` : the coordinate where the original geometry crosses the new geoemtry 
"""
function find_intersection(new_coords, old_coords, index, surface::String)
    if surface == "upper"
        # Extract x and y coordinates for the upper surface
        x_coords1 = first.(new_coords)[length(new_coords) - index - 2:end]
        y_coords1 = last.(new_coords)[length(new_coords) - index - 2:end]
        x_coords2 = first.(old_coords)[length(old_coords) - index:end]
        y_coords2 = last.(old_coords)[length(old_coords) - index:end]
    elseif surface == "lower"
        # Extract x and y coordinates for the lower surface
        x_coords1 = first.(new_coords)[1:index + 2]
        y_coords1 = last.(new_coords)[1:index + 2]
        x_coords2 = first.(old_coords)[1:index + 1]
        y_coords2 = last.(old_coords)[1:index + 1]
    else
        error("Invalid surface specified. Use 'upper' or 'lower'.")
    end

    # Process the coordinates for interpolation
    sorted_indices1 = sortperm(x_coords1)
    x_sorted1 = x_coords1[sorted_indices1]
    y_sorted1 = y_coords1[sorted_indices1]

    sorted_indices2 = sortperm(x_coords2)
    x_sorted2 = x_coords2[sorted_indices2]
    y_sorted2 = y_coords2[sorted_indices2]

    # Create interpolation functions
    interp1 = LinearInterpolation(x_sorted1, y_sorted1, extrapolation_bc=NaN)
    interp2 = LinearInterpolation(x_sorted2, y_sorted2, extrapolation_bc=NaN)

    # Set a tolerance value
    tol = 1e-4
    intersection_points = []

    # Determine loop range based on surface
    if surface == "upper"
        loop_range = x_coords1[1]:0.0001:x_coords1[end]
    else
        loop_range = x_coords1[end]:-0.0001:x_coords1[1]
    end

    # Find intersections
    for i in loop_range
        y1 = interp1(i)
        y2 = interp2(i)
        if !isnan(y1) && !isnan(y2) && abs(y2 - y1) < tol
            push!(intersection_points, (i, y1))
        end
    end

    # Average intersection points if there are multiple close ones
    if length(intersection_points) > 1
        x_point = sum(first.(intersection_points)) / length(intersection_points)
        y_point = sum(last.(intersection_points)) / length(intersection_points)
        empty!(intersection_points)
        push!(intersection_points, (x_point, y_point))
    end

    return intersection_points
end

"""
    inscribe_semicircle(p1::Tuple{Float64, Float64}, p2::Tuple{Float64, Float64}, num_points::Int=10)

Function to generate points on a semicircle between two points. In this case, used to find the points along the semicircle between the hinge points.

# Arguments:
- `p1::Tuple{Float64, Float64}` : the first endpoint of the semicircle; the first hinge point of the airfoil
- `p2::Tuple{Float64, Float64}` : the second endpoint of the semicircle; the second hinge point of the airfoil

# Key Word Arguments:
- `num_points::Int=10` : The number of points that form the semicircle

# Returns:
- `semicircle_points::Vector{Tuple{Float64, Float64}}` : the coordinates that form the semicircle with p1 and p2 as the endpoints
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
    generate_panel_geometry(x::Vector, y::Vector, angle_of_control_surface::Float64, percent_of_chord::Float64, percent_of_thickness::Float64=0.0, num_points::Int=10)

This function takes the original geometry of the airfoil and rotates it about the point of rotation as defined by the user, before updating the geometry.

# Arguments:
- `x::Vector` : the x coordinates of the original geometry
- `y::Vector` : the y coordinates of the original geometry
- `angle_of_control_surface::Float64` : user-defined angle of the control surface
- `percent_of_chord::Float64` : user-defined x-position of the point of rotation given as a decimal between 0.0 and 1.0 where 0.0 is the leading edge and 1.0 the trailing edge. It should be between 0.5 and 1.0 for trailing edge control surfaces.
- `percent_of_thickness::Float64` : user-defined y-position of the point of rotation given as a decimal between -1.0 and 1.0, where 0.0 is at the camber line

# Key Word Arugments:
- `num_points::Int=10` : number of points to form the semicircle

# Returns:
- `geometry_rotation::NamedTuple` : Contains relevant points and geometry
    - `original_coordinates::Vector{Tupe{Float64, Float64}}` : Original geometry 
    - `rotated_coordinates::Vector{Tupe{Float64, Float64}}` : Rotated geometry
    - `points::Vector{Tuple{Float64, Float64}}` : Points of the semicircle with the two hinge points as endpoints
    - `index::Int` : index of the hinge point on the lower surface of the original geometry
    - `point_of_rotation::Tuple{Float64, Float64}` : Coordinate of the point of rotation
    - `upper_hinge::Tuple{Float64, Float64}` : Coordinate of the upper hinge point
    - `lower_hinge::Tuple{Float64, Float64}` : Coordinate of the lower hinge point 
"""
function generate_panel_geometry(x, y, angle_of_control_surface, percent_of_chord, percent_of_thickness, num_points)
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
    y_lower_hinge = -1 * y_position_of_hinge

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
    points = inscribe_semicircle(coordinates[index_low], coordinates[index_high], num_points)

    point_of_rotation = (x_position_of_rotation, y_position_of_rotation)
    upper_hinge = (x_position_of_rotation, y_position_of_hinge)
    lower_hinge = (x_position_of_rotation, y_lower_hinge)

    geometry_rotation = (
        original_coordinates=original_coordinates,
        rotated_coordinates=coordinates,
        points=points,
        index=start_index,
        point_of_rotation=point_of_rotation,
        upper_hinge=upper_hinge,
        lower_hinge=lower_hinge
        )

    return geometry_rotation
end

"""
    update_geometry(geometry_rotation::NamedTuple)

This function takes the original geometry of the airfoil and rotates it about the point of rotation as defined by the user, before updating the geometry.

# Arguments:
- `x::Vector` : the x coordinates of the original geometry
- `y::Vector` : the y coordinates of the original geometry
- `angle_of_control_surface::Float64` : user-defined angle of the control surface
- `percent_of_chord::Float64` : user-defined x-position of the point of rotation given as a decimal between 0.0 and 1.0 where 0.0 is the leading edge and 1.0 the trailing edge. It should be between 0.5 and 1.0 for trailing edge control surfaces.
- `percent_of_thickness::Float64` : user-defined y-position of the point of rotation given as a decimal between -1.0 and 1.0, where 0.0 is at the camber line

# Returns:
- `system_geometry::NamedTuple` : Contains relevant points and geometry
    - `original_coordinates::Vector{Tupe{Float64, Float64}}` : Original geometry 
    - `updated_coordinates::Vector{Tupe{Float64, Float64}}` : updated geometry
    - `point_of_rotation::Tuple{Float64, Float64}` : Coordinate of the point of rotation
    - `lower_hinge::Tuple{Float64, Float64}` : Coordinate of the lower hinge point 
    - `upper_hinge::Tuple{Float64, Float64}` : Coordinate of the upper hinge point
    - `lower_intersection::Tuple{Float64, Float64}` : Coordinate of the intersection point of the lower surface
    - `upper_intersection::Tuple{Float64, Float64}` : Coordinate of the intersection point of the upper surface
    - `points::Vector{Tuple{Float64, Float64}}` : Coordinates of the semicircle formed between the two hinge points
"""
function update_geometry(geometry_rotation)
    new_coords = geometry_rotation.rotated_coordinates
    old_coords = geometry_rotation.original_coordinates
    index = geometry_rotation.index
    lower_hinge = geometry_rotation.lower_hinge
    upper_hinge = geometry_rotation.upper_hinge 
    point_of_rotation = geometry_rotation.point_of_rotation
    points = geometry_rotation.points

    # Find the indexes of the hinge points in the new, rotated geoemtry
    x_hinge_rotated = new_coords[index][1]
    x_top_hinge_rotated = new_coords[length(new_coords)-index][1]
    index_lower_hinge = findfirst(t -> t == x_hinge_rotated, first.(new_coords))
    index_upper_hinge = findfirst(t -> t == x_top_hinge_rotated, first.(new_coords))

    # Find where the new, rotated geometry (new_coords), intersects the old geometry (original_coordinates)
    lower_intersection = find_lower_intersection(old_coords, new_coords, index, "lower")
    upper_intersection = find_upper_intersection(old_coords, new_coords, index, "upper")

    # Insert the hinge points and the intersection points into the rotated geoemtry (new_coords)
    insert!(new_coords, index, lower_intersection[1])
    insert!(new_coords, length(new_coords)-index+2, upper_intersection[1])
    insert!(new_coords, index_lower_hinge+3, lower_hinge)
    insert!(new_coords, index_upper_hinge+2, upper_hinge)
    
    # Remove values of the rotated geomtry that are inside the airfoil geometry
    # Insert those values of the semicircle that bridge the gap between the rotated hinge point and the original hinge point
    if angle_of_control_surface > 0.0
        points = reverse(points)
        for i in eachindex(points)
            if first.(points)[i] > x_position_of_rotation
                insert!(new_coords, index + 2 + i, points[i])
            end
        end
        deleteat!(new_coords, length(new_coords)-index)
        deleteat!(new_coords, length(new_coords)-index)
    end

    # Remove values of the rotated geomtry that are inside the airfoil geometry
    # Insert those values of the semicircle that bridge the gap between the rotated hinge point and the original hinge point
    if angle_of_control_surface < 0.0
        deleteat!(new_coords, index+1)
        deleteat!(new_coords, index+1)

        points = reverse(points)
        for i in eachindex(points)
            if first.(points)[i] > x_position_of_rotation
                insert!(new_coords, length(new_coords)-index, points[i])
            end
        end
    end

    system_geometry = (
        original_coordinates=original_coordinates,
        updated_geometry=new_coords,
        point_of_rotation=point_of_rotation,
        lower_hinge=lower_hinge,
        upper_hinge=upper_hinge,
        lower_intersection=lower_intersection,
        upper_intersection=upper_intersection,
        points=points
        )

    return system_geometry
end


"""
    analyze(x::Vector, y::Vector, angle_of_control_surface::Float64, percent_of_chord::Float64, percent_of_thickness::Float64, num_points::Int=10)

Convenience function that updates the airfoil geometry to include a control surface.

# Arguments
- `x::Vector` : the x coordinates of the original geometry
- `y::Vector` : the y coordinates of the original geometry
- `angle_of_control_surface::Float64` : user-defined angle of the control surface
- `percent_of_chord::Float64` : user-defined x-position of the point of rotation given as a decimal between 0.0 and 1.0 where 0.0 is the leading edge and 1.0 the trailing edge. It should be between 0.5 and 1.0 for trailing edge control surfaces.
- `percent_of_thickness::Float64` : user-defined y-position of the point of rotation given as a decimal between -1.0 and 1.0, where 0.0 is at the camber line

# Key Word Arguments
- `num_points::Int=10` : the number of points that will form the semicircle between the hinge points

# Returns
- `system_geometry::NamedTuple` : NamedTuple that contains the original geometry, updated geoemtry with the control surface, and relevant points (point of rotation, lower hinge, upper hinge, lower intersection, upper intersection, and points (of the semicircle))
    - `original_coordinates::Vector{Tupe{Float64, Float64}}` : Original geometry 
    - `updated_coordinates::Vector{Tupe{Float64, Float64}}` : updated geometry
    - `point_of_rotation::Tuple{Float64, Float64}` : Coordinate of the point of rotation
    - `lower_hinge::Tuple{Float64, Float64}` : Coordinate of the lower hinge point 
    - `upper_hinge::Tuple{Float64, Float64}` : Coordinate of the upper hinge point
    - `lower_intersection::Tuple{Float64, Float64}` : Coordinate of the intersection point of the lower surface
    - `upper_intersection::Tuple{Float64, Float64}` : Coordinate of the intersection point of the upper surface
    - `points::Vector{Tuple{Float64, Float64}}` : Coordinates of the semicircle formed between the two hinge points
"""
function analyze(
    x, y, angle_of_control_surface, percent_of_chord, percent_of_thickness, num_points=10
)

    # Generate Rotated Panel Geometry
    geometry_rotation = generate_panel_geometry(x, y, angle_of_control_surface, percent_of_chord, percent_of_thickness, num_points)

    # Update geometry with intersection points
    system_geometry = update_geometry(geometry_rotation)

    return system_geometry
end