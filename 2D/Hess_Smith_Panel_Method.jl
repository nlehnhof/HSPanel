#=

Convenience functions wrapping problem, system, solution, and post processing
steps into single functions for user convenience. See analyze().

Author: Nate Lehnhof

=#

using Plots
using LinearAlgebra

"""
    find_midpoints(x::Vector, y::Vector)

Finds the midpoint between each consecutive coordinate pair.
Returns the x and y coordinates of the midpoint for each panel along the body.

# Arguments:
- `x::Vector` : x-coordinates of each point on the body
- `y::Vector` : y-coordinates of each point on the body

# Returns:
- `x_mid::Vector` : x-coordinate of the midpoint for each panel
- `y_mid::Vector` : y-coordinate of the midpoint for each panel
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
This function determines the orientation of each panel, which is useful for various aerodynamic calculations, 
    referencing Fig. 2.28 on pg. 67 of *Computational Aerodynamics* by Dr. Ning.

# Arguments:
- `x::Vector` : x-coordinates of each point on the body
- `y::Vector` : y-coordinates of each point on the body
- `x_mid::Vector` : x-coordinate midpoints of each panel

# Returns:
- `sin_panel::Vector` : sine values of the angle of each panel relative to the x-axis
- `cos_panel::Vector` : cosine values of the angle of each panel relative to the x-axis
"""
function find_sin_cos_of_panel(x, y, x_mid)
    # Initialize new variables for the sin and cos of the panels
    sin_panel = similar(x_mid, length(x_mid))
    cos_panel = similar(x_mid, length(x_mid))
    # Find sin and cos values for each panel
    for i in eachindex(sin_panel)
        distance = sqrt((x[i+1] - x[i])^2 + (y[i+1] - y[i])^2)
        sin_panel[i] = (y[i+1] - y[i]) / distance
        cos_panel[i] = (x[i+1] - x[i]) / distance
    end
    return sin_panel, cos_panel
end

"""
    geneterate_panel_geometry(x::Vector, y::Vector)

Finds the values necessary for the panel geometry.

# Arguments
- `x::Vector` : x-coordinates of each point on the body
- `y::Vector` : y-coordinates of each point on the body

# Returns
- `panel_geometry::NamedTuple` : NamedTuple including the geometry of each panel (x, y, x_mid, y_mid, sin_panel, cos_panel) 

"""
function generate_panel_geometry(x, y)
    x_mid, y_mid = find_midpoints(x, y)
    sin_panel, cos_panel = find_sin_cos_of_panel(x, y, x_mid)
    # Create a NamedTuple with the panel geoemtry
    panel_geometry = (x_m=x_mid, y_m=y_mid, s_p=sin_panel, c_p=cos_panel)
    return panel_geometry
end

"""
    find_distance_between_panels(x::Vector, y::Vector, x_mid::Vector, y_mid::Vector)

Computes the distance from each midpoint (control point) to each boundary coordinate.
Useful for determining relative positioning between control points and boundary points, based on Fig. 2.29 in *Computational Aerodynamics* by Dr. Ning.

# Arguments:
- `x::Vector` : x-coordinates of each point on the body
- `y::Vector` : y-coordinates of each point on the body
- `x_mid::Vector` : x-coordinates of control points
- `y_mid::Vector` : y-coordinates of control points

# Returns:
- `r_panel::Array` : distances between each control point `(x_mid, y_mid)` and each boundary point `(x, y)`
"""
function find_distance_between_panels(x, y, x_mid, y_mid)
    r_panel = similar(x_mid, length(x_mid), length(x_mid) + 1)
    # Use distance formula to find the distance from each midpoint to each source point 
    for i in eachindex(x_mid)
        for j in eachindex(x)
            r_panel[i, j] = sqrt((x_mid[i] - x[j])^2 + (y_mid[i] - y[j])^2)
        end
    end
    return r_panel
end

"""
    find_angle_between_panels(sin_panel::Vector, cos_panel::Vector)

Computes the sine and cosine of the angle differences between each pair of panels
These values, based on equations 2.224 and 2.225 in *Computational Aerodynamics* by Dr. Ning, 
represent the angles between control points and endpoints, with a correction on Eq. 224 for accuracy.

# Arguments:
- `sin_panel::Vector` : sine values of each panel's angle relative to the x-axis
- `cos_panel::Vector` : cosine values of each panel's angle relative to the x-axis

# Returns:
- `sin_angle_panels::Array` : sine values of angle differences between panel control points and endpoints
- `cos_angle_panels::Array` : cosine values of angle differences between panel control points and endpoints
"""
function find_angle_between_panels(sin_panel, cos_panel)
    cos_angle_panels = similar(cos_panel, length(cos_panel), length(cos_panel))
    sin_angle_panels = similar(sin_panel, length(sin_panel), length(sin_panel))
    # Use the sin(theta1 - theta2) and cos(theta1 theta2) identities to find the sin and cos values of a panel with respect to the other panels.
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
This function calculates the `beta` angle, as referenced in Fig. 2.29 
and Equation 2.212 on pg. 69 of *Computational Aerodynamics* by Dr. Ning.

# Arguments:
- `x::Vector` : x-coordinates of each point on the body
- `y::Vector` : y-coordinates of each point on the body
- `x_mid::Vector` : x-coordinates of control points
- `y_mid::Vector` : y-coordinates of control points

# Returns:
- `beta::Array` : angle between each control point `(x_mid, y_mid)` and each panel, useful for aerodynamic calculations
"""
function find_beta(x, y, x_mid, y_mid)
    beta = similar(x_mid, length(x_mid), length(x_mid))
    # The `beta` angle as referenced in Fig. 2.29 
    for i in eachindex(x_mid)
        for j in eachindex(x_mid)
            if j == i
                beta[i, j] = π
            else
                # Equation 2.212 on pg. 69 of *Computational Aerodynamics* by Dr. Ning.
                numerator = (x[j] - x_mid[i]) * (y[j+1] - y_mid[i]) - (y[j] - y_mid[i]) * (x[j+1] - x_mid[i])
                denominator = ((x[j] - x_mid[i]) * (x[j+1] - x_mid[i]) + (y[j] - y_mid[i]) * (y[j+1] - y_mid[i]))
                beta[i, j] = atan(numerator, denominator)
            end
        end
    end
    return beta
end

"""
    generate_system_geometry(x::Vector, y::Vector, generate_panel_geometry::NamedTuple)

Compute the system geometry.

# Arguments
- `x::Vector` : x-coordinates of each point on the body
- `y::Vector` : y-coordinates of each point on the body
- `panel_geometry::NamedTuple` : NamedTuple including the geometry of each panel (x, y, x_mid, y_mid, sin_panel, cos_panel) 

# Returns
- `system_geometry::NamedTuple` : NamedTuple with the system geoemtry of each panel and every other panel (sin_angle_panels, cos_angle_panels, r_panel, beta)
"""
function generate_system_geometry(x, y, panel_geometry)
    x_mid = panel_geometry.x_m
    y_mid = panel_geometry.y_m
    sin_panel = panel_geometry.s_p
    cos_panel = panel_geometry.c_p

    r_panel = find_distance_between_panels(x, y, x_mid, y_mid)
    sin_angle_panels, cos_angle_panels = find_angle_between_panels(sin_panel, cos_panel)
    beta = find_beta(x, y, x_mid, y_mid)

    # Create a NamedTuple with the system geometry
    system_geometry = ((r_p=r_panel, s_a_p=sin_angle_panels, c_a_p=cos_angle_panels, beta=beta))

    return system_geometry
end

"""
    find_A(r_panel::Array, sin_angle_panels::Array, cos_angle_panels::Array, beta::Array)

Computes the influence of each panel on every other panel.
This function determines the interaction matrix `A`, where each entry represents the influence of one panel on another. 
Refer to Eq. 2.223 and 2.233, as well as Fig. 2.234 on pages 72 and 76 in *Computational Aerodynamics* by Dr. Ning.

# Arguments:
- `r_panel::Array` : distances between each control point and the source points of all other panels
- `sin_angle_panels::Array` : sine values for the angle between each panel's control point and endpoints of all other panels
- `cos_angle_panels::Array` : cosine values for the angle between each panel's control point and endpoints of all other panels
- `beta::Array` : angle `beta` between each panel's control point `(x_mid, y_mid)` and every other panel

# Returns:
- `A::Matrix(n+1, n+1)` : matrix representing the influence of each panel on every other panel
"""
function find_A(r_panel, sin_angle_panels, cos_angle_panels, beta)
    A = similar(r_panel, size(r_panel, 1) + 1, size(r_panel, 2)) * 0.0

    # This function determines the interaction matrix `A`, where each entry represents the influence of one panel on another. 
    # Refer to Eq. 2.223 and 2.233, as well as Fig. 2.234 on pages 72 and 76 in *Computational Aerodynamics* by Dr. Ning.
    for i in eachindex(sin_angle_panels[:, 1])
        # Refer to Eq. 2.233 in Dr. Ning's textbook
        k1 = beta[1, i] * sin_angle_panels[1, i] - log(r_panel[1, i+1] / r_panel[1, i]) * cos_angle_panels[1, i]
        kn = beta[end, i] * sin_angle_panels[end, i] - log(r_panel[end, i+1] / r_panel[end, i]) * cos_angle_panels[end, i]
        A[end, i] = k1 + kn
        A[end, end] = beta[end, i] * cos_angle_panels[end, i] + log(r_panel[end, i+1] / r_panel[end, i]) * sin_angle_panels[end, i] + beta[1, i] * cos_angle_panels[1, i] + log(r_panel[1, i+1] / r_panel[1, i]) * sin_angle_panels[1, i]

        # Refer to Eq. 2.223 in Dr. Ning's textbook
        for j in eachindex(cos_angle_panels[1, :])
            A[i, j] = log(r_panel[i, j+1] / r_panel[i, j]) * sin_angle_panels[i, j] + beta[i, j] * cos_angle_panels[i, j]
            A[i, end] += log(r_panel[i, j+1] / r_panel[i, j]) * cos_angle_panels[i, j] - beta[i, j] * sin_angle_panels[i, j]
        end
    end
    
    return A
end

"""
    find_b(sin_panel::Vector, cos_panel::Vector, V_inf::Float64, AoA::Float64)

Assembles the no-flow-through conditions and the Kutta condition in one Vector.
This function constructs vector `b`, which incorporates both the no-flow-through boundary conditions and the Kutta condition. See Eq. 2.223 and 2.233, and Fig. 2.234 on pages 72 and 76 in *Computational Aerodynamics* by Dr. Ning.

# Arguments:
- `sin_panel::Vector` : sine values of each panel's angle relative to the x-axis
- `cos_panel::Vector` : cosine values of each panel's angle relative to the x-axis

# Keyword Arguments
- `V_inf::Float64=1.0` : freestream velocity defined by the user
- `AoA::Float64=0.0` : angle of attack in degrees, defined by the user

# Returns:
- `b::Vector` : boundary conditions vector, including the no-flow-through and Kutta conditions
"""
function find_b(sin_panel, cos_panel, V_inf, AoA)
    b = similar(sin_panel, length(sin_panel) + 1)
    # Refer to Eq. 2.223 in Dr. Ning's textbook
    for i in eachindex(sin_panel)
        b[i] = 2 * π * V_inf * (sin_panel[i] * cos(AoA) - cos_panel[i] * sin(AoA))
    end
    # Refer to Eq. 2.233 in Dr. Ning's textbook
    b[end] = -2 * π * V_inf * ((cos_panel[1] * cos(AoA) + sin_panel[1] * sin(AoA)) + (cos_panel[end] * cos(AoA) + sin_panel[end] * sin(AoA)))
    return b
end

"""
    generate_system_matrices(panel_geometry::NamedTuple, system_geometry::NamedTuple)

Assembles the Linear System (matrix A and matrix b)

# Arguments
- `panel_geometry::NamedTuple` : NamedTuple including the geometry of each panel (x, y, x_mid, y_mid, sin_panel, cos_panel) 
- `system_geometry::NamedTuple` : NamedTuple with the system geoemtry of each panel and every other panel (sin_angle_panels, cos_angle_panels, r_panel, beta)

# Returns
- `system_matrices::NamedTuple` : NamedTuple of matrix A (incluence on each panel) and matrix b (boundary conditions).
"""
function generate_system_matrices(panel_geometry, system_geometry, V_inf, AoA)
    sin_panel = panel_geometry.s_p
    cos_panel = panel_geometry.c_p

    r_panel = system_geometry.r_p
    sin_angle_panels = system_geometry.s_a_p
    cos_angle_panels = system_geometry.c_a_p
    beta = system_geometry.beta

    @assert all(isfinite, sin_panel) "sin_panel"
    @assert all(isfinite, cos_panel) "cos_panel"
    @assert all(isfinite, r_panel) "r_panel"
    @assert all(isfinite, sin_angle_panels) "sin_angle_panels"
    @assert all(isfinite, cos_angle_panels) "cos_angle_panels"
    @assert all(isfinite, beta) "beta"

    A = find_A(r_panel, sin_angle_panels, cos_angle_panels, beta)
    b = find_b(sin_panel, cos_panel, V_inf, AoA)

    # Create NamedTuple for the system matrices
    system_matrices = ((a_matrix=A, b_matrix=b))

    return system_matrices
end

"""
    solve(system_matrices::NamedTuple)

Computes the source strength for each boundary point and the circulation strength for the body.
This function calculates `strengths`, where each entry represents the source strength at boundary points and the final entry represents the circulation strength of the body.

# Arguments:
- `system_matrices::NamedTuple` : NamedTuple of matrix A (incluence on each panel) and matrix b (boundary conditions).

# Returns:
- `strengths::Vector` : vector of source and circulation strengths
"""
function solve(system_matrices)
    # Solve the systems of equations
    A_reg = system_matrices.a_matrix + 1e-5 * I
    strengths = A_reg \ system_matrices.b_matrix
    return strengths
end

"""
    find_vt(r_panel::Array, sin_panel::Vector, cos_panel::Vector, sin_angle_panels::Array, cos_angle_panels::Array, beta::Array, q_gamma::Vector, V_inf::Float64, AoA::Float64)

Computes the tangential velocity of each panel.

This function calculates the tangential velocities based on Eq. 2.237 on pg. 78 in *Computational Aerodynamics* by Dr. Ning.

# Arguments:
- `cos_panel::Vector` : cosine values of each panel's angle relative to the x-axis
- `sin_panel::Vector` : sine values of each panel's angle relative to the x-axis
- `r_panel::Array` : distances between each control point and source points of all other panels
- `sin_angle_panels::Array` : sine values of angles between each panel's control point and endpoints of all other panels
- `cos_angle_panels::Array` : cosine values of angles between each panel's control point and endpoints of all other panels
- `beta::Array` : angle `beta` between each panel's control point `(x_mid, y_mid)` and every other panel
- `q_gamma::Vector` : source and circulation strengths

# Keyword Arguments
- `V_inf::Float64` : freestream velocity defined by the user
- `AoA::Float64` : angle of attack in degrees, defined by the user

# Returns:
- `tangential_velocity::Vector` : tangential velocity at each panel
"""
function find_vt(cos_panel, sin_panel, r_panel, sin_angle_panels, cos_angle_panels, beta, strengths, V_inf, AoA)
    tangential_velocity = similar(strengths, length(strengths) - 1)
    # See eq. 2.237 in Dr. Ning's textbook
    
    for i in eachindex(strengths[1:end-1])
        set1 = 0.0
        set2 = 0.0

        for j in eachindex(strengths[1:end-1])    
            set1 += strengths[j] * (beta[i, j] * sin_angle_panels[i, j] - log(r_panel[i, j+1] / r_panel[i, j]) * cos_angle_panels[i, j])
            set2 += beta[i, j] * cos_angle_panels[i, j] + log(r_panel[i, j+1] / r_panel[i, j]) * sin_angle_panels[i, j]
        end
        tangential_velocity[i] = V_inf * (cos_panel[i] * cos(AoA) + sin_panel[i] * sin(AoA)) + (set1 / (2 * π)) + (strengths[end] / (2 * π)) * set2
    end
    return tangential_velocity
end

"""
    cpressure(tangential_velocity::Vector, V_inf::Float64)

Calculates the pressure coefficient `CP` at each panel based on Eq. 2.238 on pg. 78 in *Computational Aerodynamics* by Dr. Ning.

# Arguments:
- `tangential_velocity::Vector` : tangential velocity at each panel

# Keyword Argument
- `V_inf::Float64=1.0` : freestream velocity defined by the user

# Returns:
- `CP::Vector` : coefficient of pressure at each panel
"""
function cpressure(tangential_velocity, V_inf)
    # See eq. 238 in Dr. Ning's textbook
    CP = 1 .- (tangential_velocity ./ V_inf) .^ 2
    return CP
end

"""
    post_process(panel_geometry::NamedTuple, system_geometry::NamedTuple, strengths::NamedTuple, V_inf::Float64, AoA::Float64)

# Arguments
- `panel_geometry::NamedTuple` : NamedTuple including the geometry of each panel (x, y, x_mid, y_mid, sin_panel, cos_panel) 
- `system_geometry::NamedTuple` : NamedTuple with the system geoemtry of each panel and every other panel (sin_angle_panels, cos_angle_panels, r_panel, beta)
- `strengths::Vector` : The solved system that gives the strength of every source and the strength of the circulation about the airfoils

# Keyword Arguments
- `V_inf::Float64=1.0` : User-defined freestream velocity of the airfoil
- `AoA::Float64=0.0` : User-defined angle of attack of the airfoil in degrees

# Returns
- `x::Vector` : x-coordinates of each point on the body
- `y::Vector` : y-coordinates of each point on the body
- `x_mid::Vector` : x-coordinate of the midpoint for each panel
- `y_mid::Vector` : y-coordinate of the midpoint for each panel
- `V_inf::Float64=1.0` : User-defined freestream velocity.
- `AoA::Float64=0.0` : Angle of attack in degrees 
- `tangential_velocity::Vector` : tangential velocity at each panel
- `CP::Vector` : coefficient of pressure at each panel
"""
function post_process(panel_geometry, system_geometry, strengths, V_inf, AoA)
    x_mid = panel_geometry.x_m
    y_mid = panel_geometry.y_m
    cos_panel = panel_geometry.c_p
    sin_panel = panel_geometry.s_p
    r_panel = system_geometry.r_p
    sin_angle_panels = system_geometry.s_a_p
    cos_angle_panels = system_geometry.c_a_p
    beta = system_geometry.beta
    
    tangenetial_velocity = find_vt(cos_panel, sin_panel, r_panel, sin_angle_panels, cos_angle_panels, beta, strengths, V_inf, AoA)
    CP = cpressure(tangenetial_velocity, V_inf)

    return x, y, x_mid, y_mid, V_inf, AoA, tangenetial_velocity, CP
end

"""
    analyze(x::Vector, y::Vector, V_inf::Float64, AoA::Float64)

Convenience function for setting up, solving, and post-processing airfoils and airfoil systems.

# Arguments
- `x::Vector(n+1)` : x-coordinates of each point on the body
- `y::Vector(n+1)` : y-coordinates of each point on the body

# Keyword Arguments
- `V_inf::Float64=1.0` : User-defined freestream velocity.
- `AoA::Float64=0.0` : Angle of attack in degrees 

# Returns
- `geo::NamedTuple` : NamedTuple including the following
    - `x::Vector` : x-coordinates of each point on the body
    - `y::Vector` : y-coordinates of each point on the body
    - `x_mid::Vector` : x-coordinate of the midpoint for each panel
    - `y_mid::Vector` : y-coordinate of the midpoint for each panel
    - `V_inf::Float64=1.0` : User-defined freestream velocity.
    - `AoA::Float64=0.0` : Angle of attack in degrees 
    - `tangential_velocity::Vector` : tangential velocity at each panel
    - `CP::Vector` : coefficient of pressure at each panel
"""
function analyze(
    x, y, V_inf=1.0, AoA=0.0
)

    # Generate Panel Geometry
    panel_geometry = generate_panel_geometry(x, y)

    # Generate Influence Mesh
    system_geometry = generate_system_geometry(x, y, panel_geometry)

    # Assemble Linear System
    system_matrices = generate_system_matrices(panel_geometry, system_geometry, V_inf, AoA)

    # Solve System
    strengths = solve(system_matrices)

    # Post Process Solution
    x, y, x_mid, y_mid, V_inf, AoA, tangenetial_velocity, CP = post_process(panel_geometry, system_geometry, strengths, V_inf, AoA)

    geo = (x=x, y=y, x_mid=x_mid, y_mid=y_mid, V_inf=V_inf, AoA=AoA, tangenetial_velocity=tangenetial_velocity, CP=CP)

    return geo
end

# ########## Validate with Joukowsky #################
# import FLOWFoil.AirfoilTools as at

# # - Parameters - #
# center = [-0.1; 0.1]
# radius = 1.0
# alpha = 4.0
# Vinf = 1.0

# # - Joukowsky Geometry - #
# x, y = at.joukowsky(center, radius)

# # - Surface Values - #
# surface_velocity, surface_pressure_coefficient, cl = at.joukowsky_flow(
#     center, radius, alpha, Vinf
# )

# # - Your Stuff - #
# alpha = deg2rad(alpha)

# # - Plot Stuff - #
# pl = plot(; xlabel="X", ylabel="CP", yflip=true)
# plot!(
#     pl,
#     x[7:360],
#     surface_pressure_coefficient[7:360];
#     linestyle=:dash,
#     linewidth=2,
#     label="Analytic Solution",
#     title = "Coefficient of Pressure"
# )

# x, y, x_mid, y_mid, V_inf, alpha, tangenetial_velocity, CP = analyze(x, y, Vinf, alpha)

# plot!(pl, x[10:350], CP[10:350], label="Hess-Smith")

# display(pl)
# # savefig(pl, "Hess_Smith_vs_Analytic_Solution.png")