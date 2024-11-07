using Interpolations
using Plots

# Function to find intersections after divergence
function find_intersections_after_divergence(new_coords, old_coords, index2)
    # Extract x and y coordinates from the input pairs (using only the last few indices)
    x_coords1 = first.(new_coords)[index2+1:end]
    y_coords1 = last.(new_coords)[index2+1:end]
    x_coords2 = first.(old_coords)[index2-1:end]
    y_coords2 = last.(old_coords)[index2-1:end]

    # Create interpolation function for coords1
    interp1 = LinearInterpolation(x_coords1, y_coords1, extrapolation_bc=NaN)
    interp2 = LinearInterpolation(x_coords2, y_coords2, extrapolation_bc=NaN)

    tol = 1e-3

    intersections = []

    for x in 0.0:0.01:1.0
        y1 = interp1(x)
        y2 = interp2(x)
        if abs(y2 - y1) < tol
            push!(intersections, (x, interp1(x)))
        end
    end

    return intersections 
end