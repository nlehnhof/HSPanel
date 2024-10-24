#=

Input airfoil Geometry
User input control surface (% of chord, deflection angle)
Update Geometry
Run Hess-Smith with new surface

=#

using Plots
include("C:\\Users\\nlehn\\HSPanel\\revised_HS_Panel.jl")


x = [1.0, 0.75, 0.5, 0.25, 0.0, 0.25, 0.5, 0.75, 1.0]
y = [0.0, -0.125, -0.25, -0.125, 0.0, 0.125, 0.25, 0.125, 0.0]

"""
get_geometry(x_control_surface, y_control_surface)

    get_geometry takes in the x and y coordinates of a body 
    (airfoil) at a specific angle of attack and deflection angle
    and returns the x-midpoints and theta of each panel

# Arguments
    x_control_surface::Vector y_control_surface::Vector
    x_control_surface and y_control_surface are vectors 
    where (x_control_surface[i], y_control_surface[i]) 
    corresponds to the ith panel of the body

# Returns
    x_mid::Vector
    y_mid::Vector
    angle_of_panel::vector
    Returns the vectors x_mid, y_mid, and angle_of_panel 
    where (x_mid[i], y_mid[i]) are the coordinates of 
    the midpoint of the ith panel and angle_of_panel[i] is 
    the angle of the ith panel with respect to the x-axis

# Example:
```
julia> x_control_surface = [1.0, 0.75, 0.5, 0.25, 0.0, 0.25, 0.5, 0.75, 1.0]
julia> y_control_surface = [0.0, -0.125, -0.25, -0.125, 0.0, 0.125, 0.25, 0.125, 0.0]
julia> x_mid, y_mid, angle_of_panel = get_geometry(x_control_surface, y_control_surface)
([0.8690771033397395, 0.625, 0.375, 0.125, 0.125, 0.375, 0.625, 0.8690771033397395], [-0.051843363525808636, -0.1875, -0.1875, -0.0625, 0.0625, 0.1875, 0.1875, 0.07315663647419136], [-0.5253521943066973, -0.44838647994305775, 0.44838647994305775, 0.44838647994305775, 0.44838647994305775, 0.44838647994305775, -0.44838647994305775, -0.39990270559192903])
```
"""
function get_geometry(x_control_surface, y_control_surface)
    @assert length(x_control_surface) == length(y_control_surface) "length(x) does not equal length(y)."
    x_mid = similar(x_control_surface, length(x_control_surface)-1) .* zeros(length(x_control_surface)-1)
    y_mid = similar(y_control_surface, length(y_control_surface)-1) .* zeros(length(y_control_surface)-1)
    # Find midpoints
    for i in eachindex(x_mid)
        x_mid[i] = 0.5 * (x_control_surface[i] + x_control_surface[i+1])
        y_mid[i] = 0.5 * (y_control_surface[i] + y_control_surface[i+1])
    end
        
    angle_of_panel = similar(x_mid, length(x_mid)) .* zeros(length(x_mid))
    
    for i in eachindex(x_mid)
        distance = sqrt((x_control_surface[i+1] - x_control_surface[i])^2 + (y_control_surface[i+1] - y_control_surface[i])^2)
        @assert distance > 0.0 "Distance must be positive and greater than 0.0. Points $i and $(i+1)."
        angle_of_panel[i] = asin(y_control_surface[i+1] - y_control_surface[i]) / distance
    end
    return x_mid, y_mid, angle_of_panel
end

"""
change_deflection_angle(x, y, angle_of_control_surface, percent_of_chord, V_inf, alpha)
    
    Takes the x and y coordinates of a body, the user input deflection angle of the control surface, the 
    percentage of the chord at which the user wants the control surface to begin,
    and the freestream velocity and angle of attack of the body.

# Arguments
    x::Vector (vector of the x-coordinates for the body)
    y::Vector (vector of the y-coordinates for the body)
    angle_of_control_surface::Float64 (user-defined angle of the control surface)
    percent_of_chord::Float64 (user-defined point of the chord that is the start of the control surface)
    V_inf::Float64 (user-defined freestream velocity)
    alpha::Float64 (user-defined angle-of-attack of the body)

# Returns
    x_mid_CP, y_mid_CP, CP (all Vectors) --> these values are found by using the Hess-Smith Panel Method

"""
function change_deflection_angle(x, y, angle_of_control_surface, percent_of_chord, V_inf, alpha)

    @assert length(x) == length(y) "length(x) does not equal lenght(y)."    
    x_control_surface = similar(x, length(x)) * 0.0
    y_control_surface = similar(y, length(y)) * 0.0

    @assert percent_of_chord < 1.0 "percent_of_chord is a float and must be within 0.5 and 1.0"
    @assert percent_of_chord >= 0.5 "percent_of_chord must be between 0.5 and 1.0."

    angle_of_control_surface = deg2rad(angle_of_control_surface)
    first = (x[end] * percent_of_chord)
    start_index = findmin(abs.(x .- first))[2]

    percent_of_chord = 1.0 - percent_of_chord

    for i in eachindex(x)
        if i < start_index
            radius = sqrt((x[i] - x[start_index])^2 + (y[i] - y[start_index])^2)
            x_diff = x[i] - x[start_index]
            y_diff = y[i] - y[start_index]
            initial_theta = atan(y_diff/x_diff)
            final_theta = initial_theta + angle_of_control_surface
            x_cs = radius * cos(final_theta)
            y_cs = radius * sin(final_theta)
            x_control_surface[i] = x[start_index] + x_cs
            y_control_surface[i] = y[start_index] + y_cs                
        elseif i > length(y)-start_index+1
            radius = sqrt((x[i] - x[start_index])^2 + (y[i] - y[start_index])^2)
            x_diff = x[i] - x[start_index]
            y_diff = y[i] - y[start_index]
            initial_theta = atan(y_diff/x_diff)
            final_theta = initial_theta + angle_of_control_surface
            x_cs = radius * cos(final_theta)
            y_cs = radius * sin(final_theta)
            x_control_surface[i] = x[start_index] + x_cs
            y_control_surface[i] = y[start_index] + y_cs
        else
            x_control_surface[i] = x[i]
            y_control_surface[i] = y[i]
        end
    end

    # x_mid, y_mid, angle_of_panel = get_geometry(x_control_surface, y_control_surface)

    x_mid_CP, y_mid_CP, CP = HS_Panel_CP(x_control_surface, y_control_surface, V_inf, alpha)
    return x_mid_CP, y_mid_CP, CP
end


angle_of_control_surface = 0.0
percent_of_chord = 0.75 
V_inf = 1.0
alpha = 0.0

x_mid_CP, y_mid_CP, CP = change_deflection_angle(x, y, angle_of_control_surface, percent_of_chord, V_inf, alpha)

"""
plot_geometry()

    plots the original geometry of the body with 
    the new geometry of the body with its control surface

# Inputs
    (x, y) # Original body
    (x_control_surface, y_control_surface) # New body
    (x_mid_CP, y_mid_CP) # Panel midpoints on New body

# Outputs
    Plot of the new body versus the old body
"""
function plot_geometry()
    pl = plot(x, y, label = "Before", markers=true)
    # xlims!(pl, 0.0, 1.0)
    # ylims!(pl, -0.5, 0.5)
    plot!(pl, x_control_surface, y_control_surface, label="After", markers=true)
    scatter!(pl, x_mid_CP, y_mid_CP, label = "Midpoints")
    display(pl)
end

"""
plot_CP()

    Plots the Coefficient of Pressure that corresponds to each panel.

# Inputs
    CP::Vector of the coefficients of pressures of each panel
    x_mid_CP::Vector of the control points (midpoints) of each panel

# Outputs
    Plot of the Coefficient of Pressure of the body with 
    the control surface deflection
"""
function plot_CP()
    pl1 = plot(x_mid_CP[1:4], CP[1:4], markers=true, yflip=true, label="lower")
    plot!(pl1, x_mid_CP[5:end], CP[5:end], markers=true, label="upper")
    # xlims!(pl1, 0.0, 1.0)
    # ylims!(pl1, -1.0, 1.0)
    display(pl1)
    savefig("zero_degree_AoA_and_deflection_test.png")
end

plot_CP()