using Interpolations
using Plots

coords1 = [(0.9274871130596428, 0.1306217782649107), (0.7734807621135331, -0.10263139720814413), (0.72, -0.14), (0.5, -0.25), (0.25, -0.125), (0.0, 0.0), (0.25, 0.125), (0.5, 0.25), (0.72, 0.14), (0.6484807621135331, 0.11387495373796555), (0.9274871130596428, 0.1306217782649107)]    
coords2 = [(1.0, 0.0), (0.75, -0.125), (0.5, -0.25), (0.25, -0.125), (0.0, 0.0), (0.25, 0.125), (0.5, 0.25), (0.75, 0.125), (1.0, 0.0)]

function find_divergent_intersection(coords1, coords2)
    # Sort and remove duplicates manually
    coords1 = sort(coords1, by=x->x[1])
    coords1 = unique(coords1)  # Removes duplicates based on exact (x, y) pairs
    
    coords2 = sort(coords2, by=x->x[1])
    coords2 = unique(coords2)  # Removes duplicates based on exact (x, y) pairs

    x_coords1 = first.(coords1)
    y_coords1 = last.(coords1)
    x_coords2 = first.(coords2)
    y_coords2 = last.(coords2)

    # Interpolations with extrapolation set to NaN
    interp1 = LinearInterpolation(x_coords1, y_coords1, extrapolation_bc=NaN)
    interp2 = LinearInterpolation(x_coords2, y_coords2, extrapolation_bc=NaN)
    
    tolerance = 1e-3
    divergent = false
    intersection_after_divergence = []

    # Loop through x values within the range
    for x in range(minimum(x_coords1), stop=maximum(x_coords1), length=500)
        y1 = interp1(x)
        y2 = interp2(x)

        if !isnan(y1) && !isnan(y2)
            # Check for divergence
            if !divergent && abs(y1 - y2) > tolerance
                divergent = true
            end

            # If they have diverged, look for an intersection
            if divergent && abs(y1 - y2) < tolerance
                push!(intersection_after_divergence, (x, y1))
            end
        end
    end

    return intersection_after_divergence
end

# Find intersection points after divergence
intersections = find_divergent_intersection(coords1, coords2)

# Plot the original curves
plot(first.(coords1), last.(coords1), label="Coords1", lw=2)
plot!(first.(coords2), last.(coords2), label="Coords2", lw=2)

# Plot intersection points after divergence, if found
if !isempty(intersections)
    scatter!(first.(intersections), last.(intersections), color=:red, label="Intersections after Divergence", ms=4)
end

# println("Intersections after divergence: ", intersections)