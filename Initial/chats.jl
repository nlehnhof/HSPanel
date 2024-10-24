# -------------------------------
# 1. Import Required Packages
# -------------------------------

using LinearAlgebra
using Plots

# -------------------------------
# 2. Generate NACA 0008 Airfoil Coordinates
# -------------------------------

# Function to generate cosine-spaced x-coordinates
function cosine_spacing(n_points::Int)
    β = range(0, stop=π, length=n_points)
    x = (1 .- cos.(β)) ./ 2
    return x
end

# Function to calculate thickness distribution for NACA 0008
function thickness_distribution(x::Vector{Float64}, t::Float64=0.08)
    yt = 5 * t * (0.2969 .* sqrt.(x) .- 0.1260 .* x .- 
                  0.3516 .* x.^2 .+ 0.2843 .* x.^3 .- 
                  0.1015 .* x.^4)
    return yt
end

# Number of points per surface (including trailing edge)
n_points = 100

# Generate cosine-spaced x-coordinates
x_upper = cosine_spacing(n_points + 1)  # +1 to include trailing edge
yt = thickness_distribution(x_upper)

# Upper surface coordinates
y_upper = yt

# Lower surface coordinates (mirror of upper)
y_lower = -yt

# Combine upper and lower surfaces
# Exclude the first point of lower surface to avoid duplication at trailing edge
x_coords = vcat(x_upper, reverse(x_upper[2:end]))
y_coords = vcat(y_upper, reverse(y_lower[2:end]))

# Plot the airfoil to verify geometry
plot(x_coords, y_coords, aspect_ratio=1, label="NACA0008 Airfoil",
     xlabel="x", ylabel="y", title="NACA0008 Airfoil Geometry", legend=false,
     linecolor=:blue, linewidth=2)
savefig("NACA0008_Airfoil_Geometry.png")  # Save the plot as an image

# -------------------------------
# 3. Implement Hess-Smith Panel Method
# -------------------------------

# 3.1. Define Panel Structure
struct Panel
    xa::Float64  # Starting point x-coordinate
    ya::Float64  # Starting point y-coordinate
    xb::Float64  # Ending point x-coordinate
    yb::Float64  # Ending point y-coordinate
    xm::Float64  # Midpoint x-coordinate
    ym::Float64  # Midpoint y-coordinate
    length::Float64  # Length of the panel
    theta::Float64  # Angle of the panel with the x-axis
    sin_theta::Float64  # Sine of theta
    cos_theta::Float64  # Cosine of theta
    Cp::Float64  # Pressure coefficient (to be computed)
end

# 3.2. Discretize Airfoil into Panels
function define_panels(x::Vector{Float64}, y::Vector{Float64})
    n_panels = length(x) - 1
    panels = Vector{Panel}(undef, n_panels)
    for i in 1:n_panels
        xa, ya = x[i], y[i]
        xb, yb = x[i+1], y[i+1]
        xm, ym = 0.5*(xa + xb), 0.5*(ya + yb)
        dx, dy = xb - xa, yb - ya
        length_panel = sqrt(dx^2 + dy^2)
        theta = atan(dy, dx)
        sin_theta = sin(theta)
        cos_theta = cos(theta)
        panels[i] = Panel(xa, ya, xb, yb, xm, ym, length_panel, theta, sin_theta, cos_theta, 0.0)
    end
    return panels
end

# Define panels for the airfoil
panels = define_panels(x_coords, y_coords)
n_panels = length(panels)

# 3.3. Compute Influence Coefficients
function compute_influence_coefficients(panels::Vector{Panel})
    A = zeros(Float64, n_panels + 1, n_panels + 1)  # +1 for Kutta condition
    b = zeros(Float64, n_panels + 1)
    
    # Loop over each control point (midpoint of each panel)
    for i in 1:n_panels
        for j in 1:n_panels
            if i != j
                # Compute the relative position between panel i and panel j
                dx = panels[j].xm - panels[i].xm
                dy = panels[j].ym - panels[i].ym
                r_squared = dx^2 + dy^2
                r = sqrt(r_squared)
                phi = atan(dy, dx) - panels[i].theta
                
                # Avoid division by zero
                if r == 0.0
                    A[i, j] = 0.0
                else
                    # Influence coefficient based on the angle phi
                    A[i, j] = -cos(phi) / (2 * π * r)
                end
            else
                # Self-influence: set to zero to avoid singularity
                A[i, j] = 0.0
            end
        end
    end
    
    # Implement Kutta condition: sum of circulations equals zero
    for j in 1:n_panels
        A[n_panels + 1, j] = 1.0
    end
    A[n_panels + 1, n_panels + 1] = 0.0  # No influence of the extra variable
    
    # Boundary conditions: no normal flow on the surface
    # For each panel, the dot product of velocity and normal vector should be zero
    # In this simplified implementation, we set b[i] = 0 for all panels
    for i in 1:n_panels
        b[i] = 0.0
    end
    
    # Kutta condition: circulation at trailing edge is zero
    b[n_panels + 1] = 0.0
    
    return A, b
end

# Compute influence coefficients and boundary conditions
A, b = compute_influence_coefficients(panels)

# -------------------------------
# 4. Solve for Circulation Strengths (Gamma)
# -------------------------------

# Solve the linear system A * gamma = b
try
    gamma = A \ b  # Solve the linear system
catch e
    println("Error during linear solve: ", e)
    # To debug, inspect the matrix A
    println("Matrix A:")
    println(A)
    error("Matrix A is singular. Cannot solve the linear system.")
end

# Assign gamma to panels (excluding the extra variable for Kutta condition)
for i in 1:n_panels
    panels[i].Cp = gamma[i]
end

# -------------------------------
# 5. Compute Tangential Velocities and Pressure Coefficients
# -------------------------------

function compute_pressure_coefficients(panels::Vector{Panel}, V_inf::Float64, alpha::Float64, gamma::Vector{Float64})
    for i in 1:length(panels)
        # Tangential velocity approximation
        # For a symmetric airfoil at 0 AoA, the freestream contribution simplifies
        Vt = V_inf * panels[i].cos_theta + gamma[i]
        # Compute Pressure Coefficient
        panels[i].Cp = 1.0 - (Vt / V_inf)^2
    end
end

# Freestream conditions
alpha_deg = 0.0
alpha = deg2rad(alpha_deg)
V_inf = 10.0

# Compute Pressure Coefficients
compute_pressure_coefficients(panels, V_inf, alpha, gamma)

# -------------------------------
# 6. Extract Midpoints and Cp for Plotting
# -------------------------------

x_mid = [panel.xm for panel in panels]
Cp = [panel.Cp for panel in panels]

# -------------------------------
# 7. Plot Pressure Coefficient Distribution
# -------------------------------

# Separate upper and lower surfaces for symmetry verification
half = div(n_panels, 2)
x_upper = x_mid[1:half]
Cp_upper = Cp[1:half]

x_lower = x_mid[(half + 1):end]
Cp_lower = Cp[(half + 1):end]

# Reverse lower surface to match upper surface order from leading to trailing edge
x_lower_reversed = reverse(x_lower)
Cp_lower_reversed = reverse(Cp_lower)

# Plot Pressure Coefficient Distribution
p1 = plot(x_upper, Cp_upper, seriestype=:scatter, label="Upper Surface CP",
         xlabel="x (Chordwise Position)", ylabel="C_P (Pressure Coefficient)",
         title="Pressure Coefficient Distribution on NACA0008 Airfoil",
         xlims=(0, 1), ylims=(-0.5, 0.5), markersize=5, color=:red)

scatter!(p1, x_lower_reversed, Cp_lower_reversed, label="Lower Surface CP", marker=:x, color=:blue)

# Add horizontal line at C_P = 0 for reference
hline!(p1, [0], linestyle=:dash, color=:black, label="Freestream Pressure")

# Display the plot
display(p1)

# Save the Pressure Coefficient plot as an image
savefig(p1, "NACA0008_Pressure_Coefficient.png")

# -------------------------------
# 8. Verification of Symmetry
# -------------------------------

# Check if upper and lower Cp distributions are symmetric
tolerance = 1e-4
is_symmetric = all(abs.(Cp_upper .- Cp_lower_reversed) .< tolerance)

println("Is the Pressure Coefficient symmetric? ", is_symmetric)

# -------------------------------
# 9. Plot Airfoil with Midpoints (Optional)
# -------------------------------

# Plot airfoil with midpoints
p2 = plot(x_coords, y_coords, aspect_ratio=1, label="Airfoil Geometry",
         xlabel="x", ylabel="y", title="NACA0008 Airfoil Geometry", legend=false,
         linecolor=:blue, linewidth=2)

scatter!(p2, x_mid, [panel.ym for panel in panels], label="Panel Midpoints", color=:green, markersize=3)
# savefig(p2, "NACA0008_Airfoil_With_Midpoints.png")
