using Interpolations

function find_intersection_line_interpolation(coords1, coords2)
    # coords1 and coords2 are the coordinates of the two lines
    # Each is a list of tuples like [(x1, y1), (x2, y2), ...]
    
    # We will use linear interpolation to find where the two lines intersect.
    
    function interpolate_line(coords, t)
        # Interpolate the coordinates for a given t value (from 0 to 1)
        n = length(coords)
        idx = floor(Int, t * (n - 1)) + 1
        t1, t2 = coords[idx], coords[idx + 1]
        x1, y1 = t1
        x2, y2 = t2
        
        # Linear interpolation formula for x and y
        x_interpolated = x1 + (x2 - x1) * (t - (idx - 1) / (n - 1))
        y_interpolated = y1 + (y2 - y1) * (t - (idx - 1) / (n - 1))
        
        return x_interpolated, y_interpolated
    end
    
    # Check each t value for intersection
    tolerance = 1e-5  # A small tolerance to consider two points equal
    for t1 in 0:0.01:1
        x1_interpolated, y1_interpolated = interpolate_line(coords1, t1)
        
        for t2 in 0:0.01:1
            x2_interpolated, y2_interpolated = interpolate_line(coords2, t2)
            
            # If the two points are close enough, return the intersection point
            if abs(x1_interpolated - x2_interpolated) < tolerance && abs(y1_interpolated - y2_interpolated) < tolerance
                return (x1_interpolated, y1_interpolated)
            end
        end
    end
    
    # If no intersection found
    return nothing
end
