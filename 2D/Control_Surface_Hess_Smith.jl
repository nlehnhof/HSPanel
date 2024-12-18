#=

Convenience functions wrapping problem, system, solution, and post processing
steps into single functions for user convenience. See analyze().

Author: Nate Lehnhof

=#

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
        x_coords1 = first.(new_coords)[length(new_coords) - index - 4:end]
        y_coords1 = last.(new_coords)[length(new_coords) - index - 4:end]
        x_coords2 = first.(old_coords)[length(old_coords) - index - 4:end]
        y_coords2 = last.(old_coords)[length(old_coords) - index - 4:end]
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
    intersection_points = []

    # Determine loop range based on surface
    if surface == "upper"
        tol = 1e-4
        loop_range = x_coords1[1]:0.0001:x_coords1[end]
    else
        tol = 1e-4
        loop_range = x_coords1[length(x_coords1)]:0.0001:x_coords1[1]
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
    
    # Split and sort upper and lower surfaces
    mid_index = argmin(x)
    x_upper, y_upper = x[1:mid_index], y[1:mid_index]
    x_lower, y_lower = x[mid_index+1:end], y[mid_index+1:end]

    x_sorted_upper, y_sorted_upper = sort(x_upper), y_upper[sortperm(x_upper)]
    x_sorted_lower, y_sorted_lower = sort(x_lower), y_lower[sortperm(x_lower)]

    # Interpolation for hinge points
    interpolation_upper = LinearInterpolation(x_sorted_upper, y_sorted_upper, extrapolation_bc=NaN)
    interpolation_lower = LinearInterpolation(x_sorted_lower, y_sorted_lower, extrapolation_bc=NaN)
    x_position_of_rotation = percent_of_chord * x[1]

    y_lower_hinge = interpolation_upper(x_position_of_rotation)
    y_upper_hinge = interpolation_lower(x_position_of_rotation)
    y_position_of_rotation = y_lower_hinge + (y_upper_hinge - y_lower_hinge) * percent_of_thickness

    # Find the index of the value that is closest to the position of rotation
    start_index = findmin(abs.(x[1:round(Int, length(x)/2)] .- x_position_of_rotation))[2]

    # Use the start_index to insert the hinge points into the geomtry of the airfoil
    if x_position_of_rotation < x[start_index]
        insert!(coordinates, start_index+1, (x_position_of_rotation, y_lower_hinge))
        insert!(coordinates, length(coordinates) - start_index-1, (x_position_of_rotation, y_upper_hinge))
    else
        insert!(coordinates, start_index, (x_position_of_rotation, y_lower_hinge))
        insert!(coordinates, length(coordinates) - start_index -3, (x_position_of_rotation, y_upper_hinge))
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

    point_of_rotation = (x_position_of_rotation, y_position_of_rotation)
    upper_hinge = (x_position_of_rotation, y_upper_hinge)
    lower_hinge = (x_position_of_rotation, y_lower_hinge)

    points = inscribe_semicircle(lower_hinge, upper_hinge, num_points)

    for i in eachindex(points)
        coordinate = [points[i][1] - x_position_of_rotation, points[i][2] - y_position_of_rotation]
        new_coordinate = matrix_rotation * coordinate
        points[i] = (new_coordinate[1] + x_position_of_rotation, new_coordinate[2] + y_position_of_rotation)
    end

    geometry_rotation = (
        original_coordinates=original_coordinates,
        rotated_coordinates=coordinates,
        index=start_index,
        point_of_rotation=point_of_rotation,
        points=points,
        upper_hinge=upper_hinge,
        lower_hinge=lower_hinge,
        angle_of_control_surface=angle_of_control_surface
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
    - `angle_of_control_surface::Float64` : user-defined angle of the control surface 
"""
function update_geometry(geometry_rotation)
    new_coords = geometry_rotation.rotated_coordinates
    old_coords = geometry_rotation.original_coordinates
    index = geometry_rotation.index
    lower_hinge = geometry_rotation.lower_hinge
    upper_hinge = geometry_rotation.upper_hinge 
    point_of_rotation = geometry_rotation.point_of_rotation
    angle_of_control_surface = geometry_rotation.angle_of_control_surface
    points=geometry_rotation.points

    # Find where the new, rotated geometry (new_coords), intersects the old geometry (original_coordinates)
    if angle_of_control_surface < 0
        lower_intersection = find_intersection(old_coords, new_coords, index, "lower")
        upper_intersection = 1.0
        low_inter_index = findmin(abs.(first.(new_coords[1:round(Int, length(new_coords)/2)]) .- lower_intersection[1][1]))[2]
        if lower_intersection[1][1] < first.(new_coords)[low_inter_index]
            insert!(new_coords, low_inter_index + 1, lower_intersection[1])
            low_index = low_inter_index + 1
        elseif lower_intersection[1][1] > first.(new_coords)[low_inter_index]
            insert!(new_coords, low_inter_index, lower_intersection[1])
            low_index = low_inter_index
        else
            low_index = low_inter_index
        end
        n = 0
        for i in low_index-5:low_index
            j = i - n
            if last.(new_coords)[j] > lower_intersection[1][2]
                deleteat!(new_coords, j)
                n += 1
            end
        end
    elseif angle_of_control_surface == 0.0
        lower_intersection = 1.0
        upper_intersection = 1.0
    else
        upper_intersection = find_intersection(old_coords, new_coords, index, "upper")
        lower_intersection = 1.0
        upper_inter_index = findmin(abs.(first.(new_coords)[round(Int, length(new_coords)/2):end] .- upper_intersection[1][1]))[2] + round(Int, length(new_coords)/2) - 1
        if upper_intersection[1][1] < first.(new_coords)[upper_inter_index]
            insert!(new_coords, upper_inter_index, upper_intersection[1])
            up_index = upper_inter_index + 1
        elseif upper_intersection[1][1] > first.(new_coords)[upper_inter_index]
            insert!(new_coords, upper_inter_index + 1, upper_intersection[1])
            up_index = upper_inter_index + 1
        else
            up_index = upper_inter_index
        end
        n = 0
        for i in up_index-5:up_index
            j = i - n
            if last.(new_coords)[j] < upper_intersection[1][2]
                deleteat!(new_coords, j)
                n += 1
            end
        end
    end

    angle_of_control_surface = rad2deg(angle_of_control_surface)
    
    system_geometry = (
        original_coordinates=old_coords,
        updated_geometry=new_coords,
        point_of_rotation=point_of_rotation,
        lower_hinge=lower_hinge,
        upper_hinge=upper_hinge,
        lower_intersection=lower_intersection,
        upper_intersection=upper_intersection,
        points=points,
        index=index,
        angle_of_control_surface=angle_of_control_surface
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
    - `angle_of_control_surface::Float64` : Angle of control surface
"""
function analyze_cs(
    x, y, angle_of_control_surface, percent_of_chord, percent_of_thickness, num_points=10
)

    # Generate Rotated Panel Geometry
    geometry_rotation = generate_panel_geometry(x, y, angle_of_control_surface, percent_of_chord, percent_of_thickness, num_points)

    # Update geometry with intersection points
    system_geometry = update_geometry(geometry_rotation)

    return system_geometry
end

########## Validate with Joukowsky #################
include("Hess_Smith_Panel_Method.jl")

system_geometry = analyze_cs(x, y, -20.0, 0.8, 0.5, 10)

geo = analyze(first.(system_geometry.updated_geometry), last.(system_geometry.updated_geometry), Vinf, alpha)

# - Plot Stuff - #
pl2 = plot(geo.x[10:350], geo.CP[10:350], markers=false, label=false, yflip=true, title="Coefficient of Pressure at $(system_geometry.angle_of_control_surface) Degrees", xlabel="Chord", ylabel="Thickness")
display(pl2)
savefig(pl2, "CP_flap_$(system_geometry.angle_of_control_surface)_degrees.png")

# function plot_geometry()
#     pl = plot(; aspect_ratio=1, color=:blue, legend=:topleft)
#     plot!(pl, first.(system_geometry.updated_geometry), last.(system_geometry.updated_geometry), markers=true, markersize=1, label=false)
#     title!(pl, "Jakousky with Flap at $(system_geometry.angle_of_control_surface) Degree Deflection")
#     xaxis!(pl, "Normalized Chord")
#     yaxis!(pl, "Airfoil Thickness")
#     display(pl)
#     # savefig(pl, "j_flap_-10.png")
# end

# plot_geometry()