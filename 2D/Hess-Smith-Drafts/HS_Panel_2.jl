using Plots

# x = [1.0, 0.75, 0.5, 0.25, 0.0, 0.25, 0.5, 0.75, 1.0]
# y = [0.0, -0.125, -0.25, -0.125, 0.0, 0.125, 0.25, 0.125, 0.0]

"""
    find_midpoints(x::Vector, y::Vector)

Finds the midpoint between each consecutive coordinate pair.

Returns the x and y coordinates of the midpoint for each panel along the body.

# Arguments:
- `x::Vector(n+1)` : x-coordinates of each point on the body
- `y::Vector(n+1)` : y-coordinates of each point on the body

# Returns:
- `x_mid::Vector(n)` : x-coordinate of the midpoint for each panel
- `y_mid::Vector(n)` : y-coordinate of the midpoint for each panel
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
    find_sin_cos_of_panel(x::Vector, y::Vector, x_mid::Vector)

Computes the sine and cosine of the angle of each panel relative to the x-axis.

This function determines the orientation of each panel, which is useful for various aerodynamic calculations, referencing Fig. 2.28 on pg. 67 of *Computational Aerodynamics* by Dr. Ning.

# Arguments:
- `x::Vector(n+1)` : x-coordinates of each point on the body
- `y::Vector(n+1)` : y-coordinates of each point on the body
- `x_mid::Vector(n)` : x-coordinate midpoints of each panel

# Returns:
- `sin_panel::Vector(n)` : sine values of the angle of each panel relative to the x-axis
- `cos_panel::Vector(n)` : cosine values of the angle of each panel relative to the x-axis
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
    find_distance_between_panels(x::Vector, y::Vector, x_mid::Vector, y_mid::Vector)

Computes the distance from each midpoint (control point) to each boundary coordinate.

Useful for determining relative positioning between control points and boundary points, based on Fig. 2.29 in *Computational Aerodynamics* by Dr. Ning.

# Arguments:
- `x::Vector(n+1)` : x-coordinates of each point on the body
- `y::Vector(n+1)` : y-coordinates of each point on the body
- `x_mid::Vector(n)` : x-coordinates of control points
- `y_mid::Vector(n)` : y-coordinates of control points

# Returns:
- `r_panel::Array(n, n+1)` : distances between each control point `(x_mid, y_mid)` and each boundary point `(x, y)`
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
    find_angle_between_panels(sin_panel::Vector, cos_panel::Vector)

Computes the sine and cosine of the angle differences between each pair of panels.

These values, based on equations 2.224 and 2.225 in *Computational Aerodynamics* by Dr. Ning, represent the angles between control points and endpoints, with a correction on Eq. 224 for accuracy.

# Arguments:
- `sin_panel::Vector(n)` : sine values of each panel's angle relative to the x-axis
- `cos_panel::Vector(n)` : cosine values of each panel's angle relative to the x-axis

# Returns:
- `sin_angle_panels::Array(n, n)` : sine values of angle differences between panel control points and endpoints
- `cos_angle_panels::Array(n, n)` : cosine values of angle differences between panel control points and endpoints
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
    find_beta(x::Vector, y::Vector, x_mid::Vector, y_mid::Vector)

Computes the angle from each control point (midpoint) to each panel.

This function calculates the `beta` angle, as referenced in Fig. 2.29 and Equation 2.212 on pg. 69 of *Computational Aerodynamics* by Dr. Ning.

# Arguments:
- `x::Vector(n+1)` : x-coordinates of each point on the body
- `y::Vector(n+1)` : y-coordinates of each point on the body
- `x_mid::Vector(n)` : x-coordinates of control points
- `y_mid::Vector(n)` : y-coordinates of control points

# Returns:
- `beta::Array(n, n)` : angle between each control point `(x_mid, y_mid)` and each panel, useful for aerodynamic calculations
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
    find_A(r_panel::Array, sin_angle_panels::Array, cos_angle_panels::Array, beta::Array)

Computes the influence of each panel on every other panel.

This function determines the interaction matrix `A`, where each entry represents the influence of one panel on another. Refer to Eq. 2.223 and 2.233, as well as Fig. 2.234 on pages 72 and 76 in *Computational Aerodynamics* by Dr. Ning.

# Arguments:
- `r_panel::Array(n, n+1)` : distances between each control point and the source points of all other panels
- `sin_angle_panels::Array(n, n)` : sine values for the angle between each panel's control point and endpoints of all other panels
- `cos_angle_panels::Array(n, n)` : cosine values for the angle between each panel's control point and endpoints of all other panels
- `beta::Array(n, n)` : angle `beta` between each panel's control point `(x_mid, y_mid)` and every other panel

# Returns:
- `A::Matrix(n+1, n+1)` : matrix representing the influence of each panel on every other panel
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
    find_b(sin_panel::Vector, cos_panel::Vector, V_inf::Float64, alpha::Float64)

Assembles the no-flow-through conditions and the Kutta condition in one matrix.

This function constructs vector `b`, which incorporates both the no-flow-through boundary conditions and the Kutta condition. See Eq. 2.223 and 2.233, and Fig. 2.234 on pages 72 and 76 in *Computational Aerodynamics* by Dr. Ning.

# Arguments:
- `sin_panel::Vector(n)` : sine values of each panel's angle relative to the x-axis
- `cos_panel::Vector(n)` : cosine values of each panel's angle relative to the x-axis
- `V_inf::Float64` : freestream velocity defined by the user
- `alpha::Float64` : angle of attack in degrees, defined by the user

# Returns:
- `b::Vector(n+1)` : boundary conditions vector, including the no-flow-through and Kutta conditions
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
    find_q_gamma(A::Array, b::Vector)

Computes the source strength for each boundary point and the circulation strength for the body.

This function calculates `q_gamma`, where each entry represents the source strength at boundary points and the final entry represents the circulation strength of the body.

# Arguments:
- `A::Array(n+1, n+1)` : influence matrix representing how each panel affects every other panel
- `b::Vector(n+1)` : boundary conditions vector, including no-flow-through and Kutta conditions

# Returns:
- `q_gamma::Vector(n+1)` : vector of source and circulation strengths
    - `q_gamma[1:n]` : source strengths at each boundary point `(x, y)`
    - `q_gamma[n+1]` : circulation strength around the body
"""
function find_q_gamma(A, b)
    q_gamma = A \ b
    return q_gamma
end

"""
    find_vt(r_panel::Array, sin_panel::Vector, cos_panel::Vector, sin_angle_panels::Array, cos_angle_panels::Array, beta::Array, q_gamma::Vector, V_inf::Float64, alpha::Float64)

Computes the tangential velocity of each panel.

This function calculates the tangential velocities based on Eq. 2.237 on pg. 78 in *Computational Aerodynamics* by Dr. Ning.

# Arguments:
- `r_panel::Array(n, n+1)` : distances between each control point and source points of all other panels
- `sin_panel::Vector(n)` : sine values of each panel's angle relative to the x-axis
- `cos_panel::Vector(n)` : cosine values of each panel's angle relative to the x-axis
- `sin_angle_panels::Array(n, n)` : sine values of angles between each panel's control point and endpoints of all other panels
- `cos_angle_panels::Array(n, n)` : cosine values of angles between each panel's control point and endpoints of all other panels
- `beta::Array(n, n)` : angle `beta` between each panel's control point `(x_mid, y_mid)` and every other panel
- `q_gamma::Vector(n+1)` : source and circulation strengths, where:
    - `q_gamma[1:n]` : source strengths at each boundary point `(x, y)`
    - `q_gamma[n+1]` : circulation strength around the body
- `V_inf::Float64` : freestream velocity defined by the user
- `alpha::Float64` : angle of attack in degrees, defined by the user

# Returns:
- `tangential_velocity::Vector(n)` : tangential velocity at each panel
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
    cpressure(tangential_velocity::Vector, V_inf::Float64)

Computes the coefficient of pressure for each panel.

This function calculates the pressure coefficient `CP` based on Eq. 2.238 on pg. 78 in *Computational Aerodynamics* by Dr. Ning.

# Arguments:
- `tangential_velocity::Vector(n)` : tangential velocity at each panel
- `V_inf::Float64` : freestream velocity defined by the user

# Returns:
- `CP::Vector(n)` : coefficient of pressure at each panel
"""
function cpressure(tangential_velocity, V_inf)
    CP = 1 .- (tangential_velocity ./ V_inf) .^ 2
    return CP
end

"""
    HS_Panel_CP(x::Vector, y::Vector, V_inf::Float64, alpha::Float64)

Computes the coefficient of pressure for each panel.

This function calculates the coefficient of pressure `CP` along with the tangential velocities and midpoint coordinates for each panel on the body.

# Arguments:
- `x::Vector(n+1)` : x-coordinates of each point on the body
- `y::Vector(n+1)` : y-coordinates of each point on the body
- `V_inf::Float64` : freestream velocity defined by the user
- `alpha::Float64` : angle of attack in degrees, defined by the user

# Returns:
- `x::Vector(n+1)` : x-coordinates of each point on the body
- `y::Vector(n+1)` : y-coordinates of each point on the body
- `x_mid::Vector(n)` : x-coordinates for the control points of each panel
- `y_mid::Vector(n)` : y-coordinates for the control points of each panel
- `tangential_velocity::Vector(n)` : tangential velocity at each panel
- `CP::Vector(n)` : coefficient of pressure at each panel
"""
function HS_Panel_CP(x, y, V_inf, alpha)
    x_mid, y_mid = find_midpoints(x, y)

    sin_panel, cos_panel = find_sin_cos_of_panel(x, y, x_mid)

    r_panel = find_distance_between_panels(x, y, x_mid, y_mid)

    sin_angle_panels, cos_angle_panels = find_angle_between_panels(sin_panel, cos_panel)

    beta = find_beta(x, y, x_mid, y_mid)

    A = find_A(r_panel, sin_angle_panels, cos_angle_panels, beta)

    b = find_b(sin_panel, cos_panel, V_inf, alpha)

    q_gamma = A \ b

    tangential_velocity = find_vt(cos_panel, sin_panel, r_panel, sin_angle_panels, cos_angle_panels, beta, q_gamma, V_inf, alpha)

    CP = cpressure(tangential_velocity, V_inf)

    return x, y, x_mid, y_mid, tangential_velocity, q_gamma, CP
end

HS_Panel_CP(x, y, 1.0, 0.0)


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
