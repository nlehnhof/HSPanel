#=

Input airfoil Geometry
User input control surface (% of chord, deflection angle)
Update Geometry
Run Hess-Smith with new surface

=#

using Plots

x = [1.0, 0.75, 0.5, 0.25, 0.0, 0.25, 0.5, 0.75, 1.0]
y = [0.0, -0.125, -0.25, -0.125, 0.0, 0.125, 0.25, 0.125, 0.0]

# """
#     get_geometry(x, y)

#     get_geometry takes in the x and y coordinates of a body (airfoil) and returns the x-midpoints and theta of each panel

#     # Arguments
#     x::Vector y::Vector
#     x and y are vectors where (x[i], y[i]) corresponds to the ith panel of the body

#     # Returns
#     x_mid::Vector
#     y_mid::Vector
#     angle_of_panel::vector
#     Returns the vectors x_mid, y_mid, and angle_of_panel 
#     where (x_mid[i], y_mid[i]) are the coordinates of 
#     the midpoint of the ith panel and angle_of_panel[i] is 
#     the angle of the ith panel with respect to the x-axis

# Example:
# julia> x = [1.0, 0.75, 0.5, 0.25, 0.0, 0.25, 0.5, 0.75, 1.0]
# julia> y = [0.0, -0.125, -0.25, -0.125, 0.0, 0.125, 0.25, 0.125, 0.0]
# julia> x_mid, y_mid, angle_of_panel = get_geometry(x, y)
# ([0.875, 0.625, 0.375, 0.125, 0.125, 0.375, 0.625, 0.875], [-0.0625, -0.1875, -0.1875, -0.0625, 0.0625, 0.1875, 0.1875, 0.0625], [-0.44838647994305775, -0.44838647994305775, 0.44838647994305775, 0.44838647994305775, 0.44838647994305775, 0.44838647994305775, -0.44838647994305775, -0.44838647994305775])

# """
# function get_geometry(x, y)
#     @assert length(x) == length(y) "length(x) does not equal length(y)."
#     x_mid = similar(x, length(x)-1) .* zeros(length(x)-1)
#     y_mid = similar(y, length(y)-1) .* zeros(length(y)-1)
#     # Find midpoints
#     for i in eachindex(x_mid)
#         x_mid[i] = 0.5 * (x[i] + x[i+1])
#         y_mid[i] = 0.5 * (y[i] + y[i+1])
#     end
        
#     angle_of_panel = similar(x_mid, length(x_mid)) .* zeros(length(x_mid))
    
#     for i in eachindex(x_mid)
#         distance = sqrt((x[i+1] - x[i])^2 + (y[i+1] - y[i])^2)
#         @assert distance > 0.0 "Distance must be positive and greater than 0.0. Points $i and $(i+1)."
#         angle_of_panel[i] = asin(y[i+1] - y[i]) / distance
#     end
#     return x_mid, y_mid, angle_of_panel
# end

# x_mid, y_mid, angle_of_panel = get_geometry(x, y)

function change_geometry(x, y, angle_of_control_surface, percent_of_chord)

    @assert length(x) == length(y) "length(x) does not equal lenght(y)."    
    x_control_surface = similar(x, length(x)) * 0.0
    y_control_surface = similar(y, length(y)) * 0.0

    percent_of_chord = 1.0 - percent_of_chord
    angle_of_control_surface = deg2rad(angle_of_control_surface)
    @assert percent_of_chord < 1.0 && percent_of_chord > 0.0 "percent_of_chord is a float and must be within 0.0 and 1.0"
    start_index = Int(round(length(y) * percent_of_chord, digits=0))

    x_control_surface[start_index:length(x) - start_index] = x[start_index:length(x) - start_index]
    y_control_surface[start_index:length(y) - start_index] = y[start_index:length(y) - start_index]
    x_mid_control = similar(x, length(x) - 1) * 0.0
    y_mid_control = similar(y, length(y) - 1) * 0.0
    angle_of_panel = similar(x_mid_control, length(x_mid)) * 0.0

    plot(x[1:start_index], y[1:start_index], label="before")

    for i in 1:start_index
        distance = sqrt((x[i+1] - x[i])^2 + (y[i+1] - y[i])^2)
        @assert distance > 0.0 "Distance must be positive and greater than 0.0. Points $i and $(i+1)."
        angle1 = asin((y[i+1] - y[i]) / distance) + (pi/2)
        angle2 = angle_of_control_surface
        angle3 = 2pi - angle1 - angle2
        height = distance * sin(angle2) / sin(angle3)
        y_control_surface[i] = y[i] + height
        x_control_surface[i] = x[i]
        length = sqrt((y_control_surface[i+1] - y_control_surface[i])^2 + (x_control_surface[i+1] - x_control_surface[i])^2)
        angle_of_panel[i] = asin(y_control_surface[i] - y_control_surface[i+1]) / length
    end

    plot!(x_control_surface[1:start_index], y_control_surface[1:start_index], label="After")
    xlims!(0.0, 1.0)
    ylims!(-0.5, 0.5)

    # for i in (length(x) - start_index):length(x)
    #     radius = sqrt((x[i] - x_control_surface[start_index])^2 + (y[i])^2)
    #     x_control_surface[i] = radius * cos(angle_of_control_surface)
    #     y_control_surface[i] = radius * sin(angle_of_control_surface)
    # end

    # for i in eachindex(angle_of_panel)
    #     distance = sqrt((x[i+1] - x[i])^2 + (y[i+1] - y[i])^2)
    #     @assert distance > 0.0 "Distance must be positive and greater than 0.0. Points $i and $(i+1)."
    #     angle_of_panel[i] = asin(y_control_surface[i+1] - y_control_surface[i]) / distance
    # end

    # return x_control_surface, y_control_surface, angle_of_panel

end

change_geometry(x, y, 5, 0.75)
# x_control_surface, y_control_surface, angle_of_panel = change_geometry(x, y, 5, 0.9) 