using Plots

# x = [1.0, 0.75, 0.5, 0.25, 0.0, 0.25, 0.5, 0.75, 1.0]
# y = [0.0, -0.125, -0.25, -0.125, 0.0, 0.125, 0.25, 0.125, 0.0]

"""
find_midpoints(x::Vector, y::Vector)

    Find the midpoint between each coordinate. 
    Returns the x-coordinate of the midpoint and the y-coordinate of the midpoint.

# Arguments
    x::Vector(n+1)
    y::Vector(n+1)
    The x and y coordinates of each point on the body.

# Returns
    x_mid::Vector(n) the x coordinate of the control point (midpoint) of each panel
    y_mid::Vector(n) the y coordinate of the control point (midpoint) of each panel
"""
function find_midpoints(x, y)
    # Initialize new variables for midpoints
    x_mid = zeros(Float64, length(x) - 1)
    y_mid = zeros(Float64, length(x) - 1)
    # Find midpoints
    for i in eachindex(x_mid)
        x_mid[i] = 0.5 * (x[i] + x[i+1])
        y_mid[i] = 0.5 * (y[i] + y[i+1])
    end
    return x_mid, y_mid
end

"""
    find_sin_cos_of_panel(x, y)

Computes the sin_panel and cos_panel of each panel with respect to the x-axis.
See Fig. 2.28 on pg. 67 in Computational Aerodynamics by Dr. Ning

# Arguments
    x::Vector(n+1)
    y::Vector(n+1)
    The x and y coordinates of each point on the body.

# Returns
    sin_panel::Vector(n)
    cos_panel::Vector(n)
    The sin and cos value of each panel with respect to the x-axis.
"""
function find_sin_cos_of_panel(x, y, x_mid)
    sin_panel = similar(x_mid, length(x_mid))
    cos_panel = similar(x_mid, length(x_mid))
    for i in eachindex(sin_panel)
        distance = sqrt((x[i+1] - x[i])^2 + (y[i+1] - y[i])^2)
        sin_panel[i] = (y[i+1] - y[i]) / distance
        cos_panel[i] = (x[i+1] - x[i]) / distance
    end
    return sin_panel, cos_panel
end

"""
    find_distance_between_panels(x, y, x_mid, y_mid)

Computes the distance from each midpoint to each (x, y) coordinate
See Fig. 2.29 on pg. 69 in Computational Aerodynamics by Dr. Ning

# Arguments
    x::Vector and y::Vector are the x and y coordinates of each point on the body.
    x_mid::Vector and y_mid::Vector are the (x_mid, y_mid) coordinates for the control point on each panel.

# Returns
    r_panel::Array(n, n+1) Gives the distance between each control point (x_mid, y_mid) and each coordinate (x, y).
"""
function find_distance_between_panels(x, y, x_mid, y_mid)
    r_panel = similar(x_mid, length(x_mid), length(x_mid) + 1)
    for i in eachindex(x_mid)
        for j in eachindex(x)
            r_panel[i, j] = sqrt((x_mid[i] - x[j])^2 + (y_mid[i] - y[j])^2)
        end
    end
    return r_panel
end

"""
    find_angle_between_panels(sin_panel, cos_panel)

Computes sin(theta1 - theta2) and cos(theta1 - theta2)
See Fig. 2.29 on pg. 69 in Computational Aerodynamics by Dr. Ning
See Eq. 2.224 and 2.225 on pg. 72
Note that the identities are not correct. 
Eq. 224 should be sin i * cos j - cos i sin j

# Arguments
    sin_panel::Vector(n)
    cos_panel::Vector(n)
    The sin and cos values of each panel with respect to the x-axis 

# Returns
    sin_angle_panels::Array(n, n)
    cos_angle_panels::Array(n, n)
    The sin and cos values for the angle between the control point of a panel and the end points of every other panel
"""
function find_angle_between_panels(sin_panel, cos_panel)
    cos_angle_panels = similar(cos_panel, length(cos_panel), length(cos_panel))
    sin_angle_panels = similar(sin_panel, length(sin_panel), length(sin_panel))
    for i in eachindex(sin_panel)
        for j in eachindex(cos_panel)
            sin_angle_panels[i, j] = sin_panel[i] * cos_panel[j] - cos_panel[i] * sin_panel[j]
            cos_angle_panels[i, j] = cos_panel[i] * cos_panel[j] + sin_panel[i] * sin_panel[j]
        end
    end
    return sin_angle_panels, cos_angle_panels
end

"""
    find_beta(x, y, x_mid, y_mid)

Computes the angle from each midpoint to each panel
See Fig. 2.29 and Equation 2.212 on pg. 69 in Computational Aerodynamics by Dr. Ning

# Arguments
    x::Vector(n+1)
    y::Vector(n+1) 
    The x and y coordinates of each point on the body.
    x_mid::Vector(n)
    y_mid::Vector(n) 
    The (x_mid, y_mid) coordinates for the control point on each panel.

# Returns
    beta::Array(n, n) that gives us the beta angle 
    between the control point of a panel (x_mid, y_mid) and another panel. See figure ___ in ____ textbook.
"""
function find_beta(x, y, x_mid, y_mid)
    beta = similar(x_mid, length(x_mid), length(x_mid))
    for i in eachindex(x_mid)
        for j in eachindex(x_mid)
            if j == i
                beta[i, j] = π
            else
                numerator = (x[j] - x_mid[i]) * (y[j+1] - y_mid[i]) - (y[j] - y_mid[i]) * (x[j+1] - x_mid[i])
                denominator = ((x[j] - x_mid[i]) * (x[j+1] - x_mid[i]) + (y[j] - y_mid[i]) * (y[j+1] - y_mid[i]))
                beta[i, j] = atan(numerator, denominator)
            end
        end
    end
    return beta
end

"""
    find_A(r_panel, sin_angle_panels, cos_angle_panels, beta)

Computes the influence of each panel on every other panel.
See Eq. 2.223 and 2.233 and Fig. 2.234 on pg. 72 and 76 in Computational Aerodynamics by Dr. Ning

# Arguments
    r_panel::Array(n, n+1) that is the distance between the control point of each panel to the source points of every other panel. 
    sin_angle_panels::Array(n, n)
    cos_angle_panels::Array(n, n)
        The sin and cos values for the angle between the control point of a panel and the end points of every other panel
    beta::Array(n, n) that gives us the beta angle between the control point of a panel (x_mid, y_mid) and another panel. See figure ___ in ____ textbook.

# Returns
    A::Matrix(n+1, n+1) that gives us the matrix of how each panel influences every other panel.
"""
function find_A(r_panel, sin_angle_panels, cos_angle_panels, beta)
    A = similar(r_panel, size(r_panel, 1) + 1, size(r_panel, 2)) * 0.0

    # These for loops give us the values for the [n x n] matrix
    for i in eachindex(sin_angle_panels[:, 1])
        k1 = beta[1, i] * sin_angle_panels[1, i] - log(r_panel[1, i+1] / r_panel[1, i]) * cos_angle_panels[1, i]
        kn = beta[end, i] * sin_angle_panels[end, i] - log(r_panel[end, i+1] / r_panel[end, i]) * cos_angle_panels[end, i]
        A[end, i] = k1 + kn
        A[end, end] = beta[end, i] * cos_angle_panels[end, i] + log(r_panel[end, i+1] / r_panel[end, i]) * sin_angle_panels[end, i] + beta[1, i] * cos_angle_panels[1, i] + log(r_panel[1, i+1] / r_panel[1, i]) * sin_angle_panels[1, i]

        for j in eachindex(cos_angle_panels[1, :])
            A[i, j] = log(r_panel[i, j+1] / r_panel[i, j]) * sin_angle_panels[i, j] + beta[i, j] * cos_angle_panels[i, j]
            A[i, end] += log(r_panel[i, j+1] / r_panel[i, j]) * cos_angle_panels[i, j] - beta[i, j] * sin_angle_panels[i, j]
        end
    end
    return A
end

"""
    find_b(sin_theta, cos_panel, V_inf, alpha)

Assembles the no-flow-through conditions and the Kutta condition in one matrix.
See Eq. 2.223 and 2.233 and Fig. 2.234 on pg. 72 and 76 in Computational Aerodynamics by Dr. Ning

# Arguments
    sin_panel::Vector(n)
    cos_panel::Vector(n)
    The sin and cos values of each panel with respect to the x-axis 
    V_inf::Float64 user-defined freestream velocity
    alpha::Float64 user-defined angle of attack in degrees

# Returns
    b::Vector(n+1) the boundary conditions including
    the no-flow through condition and the Kutta-Conidtion
"""
function find_b(sin_panel, cos_panel, V_inf, alpha)
    b = similar(sin_panel, length(sin_panel) + 1)
    for i in eachindex(sin_panel)
        b[i] = 2 * π * V_inf * (sin_panel[i] * cos(alpha) - cos_panel[i] * sin(alpha))
    end
    b[end] = -2 * π * V_inf * ((cos_panel[1] * cos(alpha) + sin_panel[1] * sin(alpha)) + (cos_panel[end] * cos(alpha) + sin_panel[end] * sin(alpha)))
    return b
end

"""
    find_q_gamma(A, b)

Computes the source strength for each point and the circulation strength for the body.

# Arguments
    A::Array(n+1, n+1) influences of every panel on every other panel
    b::Vector(n+1) boundary conditions: no-flow through and Kutta-Condition
    
# Returns
    q_gamma::Vector(n+1)
        q_gamma[1:n] gives us the source strength at each (x, y) coordinate
        q_gamma[n+1] gives us the circulation strength
```
"""
function find_q_gamma(A, b)
    q_gamma = A \ b
    return q_gamma
end

"""
    find_vt(r_panel, sin_angle_panels, cos_angle_panels, beta, q_gamma, V_inf, alpha)

Computes the tangential velocity of each panel.
See Eq. 2.237 on pg. 78 in Computational Aerodynamics by Dr. Ning

# Arguments
    sin_panel::Vector(n)
    cos_panel::Vector(n)
        The sin and cos values of each panel with respect to the x-axis 
    sin_angle_panels::Array(n, n)
    cos_angle_panels::Array(n, n)
        The sin and cos values for the angle between the control point of a panel and the end points of every other panel
    r_panel::Array(n, n+1) that is the distance between the control point of each panel to the source points of every other panel.
    sin_angle_panels::Array(n, n)
    cos_angle_panels::Array(n, n)
        The sin and cos values for the angle between the control point of a panel and the end points of every other panel
    beta::Array(n, n) that gives us the beta angle 
        between the control point of a panel (x_mid, y_mid) and another panel. See figure ___ in ____ textbook.
    q_gamma::Vector(n+1)
        q_gamma[1:n] gives us the source strength at each (x, y) coordinate
        q_gamma[n+1] gives us the circulation strength 
    V_inf::Float64 user-defined freestream velocity
    alpha::Float64 user-defined angle of attack in degrees

# Returns
    tangential_velocity::Vector(n) gives the tangential velocity at each panel
"""
function find_vt(cos_panel, sin_panel, r_panel, sin_angle_panels, cos_angle_panels, beta, q_gamma, V_inf, alpha)
    tangential_velocity = similar(q_gamma, length(q_gamma) - 1)
    for i in eachindex(q_gamma[1:end-1])
        set1 = 0.0
        set2 = 0.0
        for j in eachindex(q_gamma[1:end-1])
            set1 += q_gamma[j] * (beta[i, j] * sin_angle_panels[i, j] - log(r_panel[i, j+1] / r_panel[i, j]) * cos_angle_panels[i, j])
            set2 += beta[i, j] * cos_angle_panels[i, j] + log(r_panel[i, j+1] / r_panel[i, j]) * sin_angle_panels[i, j]
        end
        tangential_velocity[i] = V_inf * (cos_panel[i] * cos(alpha) + sin_panel[i] * sin(alpha)) + (set1 / (2 * π)) + (q_gamma[end] / (2 * π)) * set2
    end
    return tangential_velocity
end

"""
    cpressure(tangential_velocity, V_inf)

Computes the coefficient of pressure of each panel.
See Eq. 2.238 on pg. 78 in Computational Aerodynamics by Dr. Ning

# Arguments
    tangential_velocity::Vector(n) gives the tangential velocity at each panel
    V_inf::Float64 user-defined freestream velocity

# Returns
    CP::Vector(n) gives the coefficient of pressure at each panel
"""
function cpressure(tangential_velocity, V_inf)
    CP = 1 .- (tangential_velocity ./ V_inf) .^ 2
    return CP
end

"""
    HS_Panel_CP(x, y, V_inf, alpha)

Computes the coefficient of pressure of each panel.

# Arguments
    x::Vector(n+1)
    y::Vector(n+1)
    The x and y coordinates of each point on the body.
    V_inf::Float64 user-defined freestream velocity
    alpha::Float64 user-defined angle of attack in degrees

# Returns
    x::Vector(n+1)
    y::Vector(n+1) 
        The x and y coordinates of each point on the body.
    x_mid::Vector(n)
    y_mid::Vector(n) 
        The (x_mid, y_mid) coordinates for the control point on each panel.
    tangential_velocity::Vector(n) gives the tangential velocity at each panel
    CP::Vector(n) gives the coefficient of pressure at each panel
"""
function HS_Panel_CP(x, y, V_inf, alpha)
    x_mid, y_mid = find_midpoints(x, y)

    sin_panel, cos_panel = find_sin_cos(x, y, x_mid)

    r_panel = find_distance_between_panels(x, y, x_mid, y_mid)

    sin_angle_panels, cos_angle_panels = find_thetas(sin_panel, cos_panel)

    beta = find_beta(x, y, x_mid, y_mid)

    A = find_A(r_panel, sin_angle_panels, cos_angle_panels, beta)

    b = find_b(sin_panel, cos_panel, V_inf, alpha)

    q_gamma = A \ b

    tangential_velocity = find_vt(cos_panel, sin_panel, r_panel, sin_angle_panels, cos_angle_panels, beta, q_gamma, V_inf, alpha)

    CP = cpressure(tangential_velocity, V_inf)

    return x, y, x_mid, y_mid, tangential_velocity, q_gamma, CP
end


########## Validate with Joukowsky #################

import FLOWFoil.AirfoilTools as at

# - Parameters - #
center = [-0.1; 0.1]
radius = 1.0
alpha = 4.0
Vinf = 1.0

# - Joukowsky Geometry - #
x, y = at.joukowsky(center, radius)

# - Surface Values - #
surface_velocity, surface_pressure_coefficient, cl = at.joukowsky_flow(
    center, radius, alpha, Vinf
)

# - Your Stuff - #

alpha = deg2rad(alpha)

x, y, x_mid, y_mid, tangential_velocity, q_gamma, CP = HS_Panel_CP(x, y, Vinf, alpha)

# - Plot Stuff - #
pl = plot(; xlabel="x", ylabel="cp", yflip=true)
plot!(
    pl,
    x[1:342],
    surface_pressure_coefficient[1:342];
    linestyle=:dash,
    linewidth=2,
    label="Analytic Solution",
)

plot!(pl, x[1:360], CP[1:360], label="Hess-Smith")

display(pl)
# savefig(pl, "Hess_Smith_vs_Analytic_Solution.png")