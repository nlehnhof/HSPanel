mutable struct Panel3D
    # Source points (x, y of airfoil)
    Panel_le_point::Array{Vector{Float64}, 2}

    # Midpoint of each chordwise panel connecting airfoils
    Panel_chord_mid::Array{Vector{Float64}, 2}

    # Midpoint of each spanwise panel (within airfoil)
    Panel_span_mid::Array{Vector{Float64}, 2}

    # Control point of each panel
    Panel_control_point::Array{Vector{Float64}, 2}

    # Theta of each panel (pitch)
    Panel_theta_y::Vector{Float64}
end

function panel_setup_3D(x, y, z; graph = false)
    n = length(x) - 1  # number of panels in x direction (chordwise)
    m = length(z)      # number of layers in z direction (spanwise)
    
    # Initialize arrays of size n × m (chordwise × spanwise)
    le_point = Array{Vector{Float64}, 2}(undef, n, m)
    span_mid = Array{Vector{Float64}, 2}(undef, n, m)
    chord_mid = Array{Vector{Float64}, 2}(undef, n, m - 1)
    theta_y = zeros(Float64, n)
    control_point = Array{Vector{Float64}, 2}(undef, n, m - 1)

    # Loop over spanwise and chordwise points to assign values
    for j in 1:n  # chordwise panels
        for i in 1:m  # spanwise layers
            le_point[j, i] = [x[j], y[j], z[i]]  # 3D leading edge point
        end
    end

    for j in 1:n
        theta_y[j] = atan((y[j + 1] - y[j]) / (x[j + 1] - x[j]))  # pitch angle
        
        for i in 1:(m - 1)  # as you go along the span
            chord_mid[j, i] = [(x[j + 1] - x[j]), (y[j + 1] - y[j]), z[i]]
            span_mid[j, i] = [x[j], y[j], (z[i + 1] - z[i])]
            control_point[j, i] = [0.75 * (x[j + 1] - x[j]), 0.75 * (y[j + 1] - y[j]), 0.75 * (z[i + 1] - z[i])]
        end
    end

    # Create the panel3D struct with all computed points
    panel_data_3D = Panel3D(le_point, chord_mid, span_mid, control_point, theta_y)

    # Optional: Plot the 3D geometry if graph=true
    if graph
        plot1 = plot(x, y, z)
        plot!(legend = false)
        xlims!(-2.0, 2.0)
        ylims!(-2.0, 2.0)
        zlims!(-2.0, 2.0)
        display(plot1)
    end

    return panel_data_3D
end

# Example inputs
x = [1.0, 0.5, 0.0, 0.5, 1.0]
y = [0.0, -0.25, 0.0, 0.25, 0.0]
z = [0.0, 0.25, 0.5, 0.75]

panel_data = panel_setup_3D(x, y, z)
