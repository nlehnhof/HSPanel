#=

Input airfoil Geometry
User input control surface (% of chord, deflection angle)
Update Geometry
Run Hess-Smith with new surface

=#

using Plots
include("C:\\Users\\nlehn\\HSPanel\\HS_Panel_2.jl")


x = [1.0, 0.75, 0.5, 0.25, 0.0, 0.25, 0.5, 0.75, 1.0]
y = [0.0, -0.125, -0.25, -0.125, 0.0, 0.125, 0.25, 0.125, 0.0]

"""
change_deflection_angle(x, y, angle_of_control_surface, percent_of_chord, V_inf, alpha)
    
    Takes the x and y coordinates of a body, the user input deflection angle of the control surface, the 
    percentage of the chord at which the user wants the control surface to begin,
    and the freestream velocity and angle of attack of the body.

# Arguments
    x::Vector(n+1) (vector of the x-coordinates for the body)
    y::Vector(n+1) (vector of the y-coordinates for the body)
    angle_of_control_surface::Float64 (user-defined angle of the control surface)
    percent_of_chord::Float64 (user-defined point of the chord that is the start of the control surface)
    V_inf::Float64 (user-defined freestream velocity)
    alpha::Float64 (user-defined angle-of-attack of the body)

# Returns
    x_mid_CP, y_mid_CP, tangential_velocity, CP (all Vectors of size n) --> these values are found by using the Hess-Smith Panel Method
    x_control_surface::Vector(n)
    y_control_surface::Vector(n)
        Where (x_control_surface[i], y_control_surface[i]) corresponds to the ith panel of the body

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

    x, y, x_mid_CP, y_mid_CP, tangential_velocity, CP = HS_Panel_CP(x_control_surface, y_control_surface, V_inf, alpha)
    
    return x_control_surface, y_control_surface, x_mid_CP, y_mid_CP, tangential_velocity, CP
end


angle_of_control_surface = 5.0
percent_of_chord = 0.75 
V_inf = 1.0
alpha = 0.0

x_control_surface, y_control_surface, x_mid_CP, y_mid_CP, tangential_velocity, CP = change_deflection_angle(x, y, angle_of_control_surface, percent_of_chord, V_inf, alpha)

"""
plot_geometry()

    plots the original geometry of the body with 
    the new geometry of the body with its control surface

# Arguments
    (x, y) # Original body
    (x_control_surface, y_control_surface) # New body
    (x_mid_CP, y_mid_CP) # Panel midpoints on New body

# Returns
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

# plot_geometry()

"""
plot_CP()

    Plots the Coefficient of Pressure that corresponds to each panel.

# Arguments
    CP::Vector(n) of the coefficients of pressures of each panel
    x_mid_CP::Vector(n) of the control points (midpoints) of each panel

# Returns
    Plot of the Coefficient of Pressure of the body with 
    the control surface deflection
"""
function plot_CP()
    pl1 = plot(x_mid_CP[1:4], CP[1:4], markers=true, yflip=true, label="lower")
    plot!(pl1, x_mid_CP[5:end], CP[5:end], markers=true, label="upper")
    # xlims!(pl1, 0.0, 1.0)
    # ylims!(pl1, -1.0, 1.0)
    display(pl1)
    # savefig("zero_degree_AoA_and_deflection_test.png")
end