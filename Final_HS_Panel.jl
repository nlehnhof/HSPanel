using Plots
using LinearAlgebra

################ Get Coordinates from File ###########

# function get_coordinates(file)
#     x, y = open(file, "r") do f
#         x = Float64[]
#         y = Float64[]
#         for line in eachline(f)
#             entries = split(chomp(line))
#             push!(x, parse(Float64, entries[1]))
#             push!(y, parse(Float64, entries[2]))
#         end
#         x, y
#     end
# end

# # x, y : vectors : coordinates of the body -- size(n+1)
# x, y = get_coordinates("naca0008.txt")

########## Panel Method ###################
#=
    Takes the coordinates of a body and returns the Coefficient of Pressure
=#

function HS_Panel_CP(x, y, V_inf, alpha)

    function find_midpoints(x, y)
        # Initialize new variables for midpoints
        n = length(x) - 1
        x_mid = zeros(Float64, n)
        y_mid = zeros(Float64, n)
        # Find midpoints
        for i in eachindex(x_mid)
            x_mid[i] = 0.5 * (x[i] + x[i+1])
            y_mid[i] = 0.5 * (y[i] + y[i+1])
        end
        return x_mid, y_mid
    end

    # Midpoints : vectors : Control points for each panel are located at the center of each panel -- length(n)
    x_mid, y_mid = find_midpoints(x, y)

    function find_sin_cos(x, y)
        distance = similar(x_mid, length(x_mid))
        sin_theta = similar(x_mid, length(x_mid))
        cos_theta = similar(x_mid, length(x_mid))
        for i in eachindex(distance)
            distance[i] = sqrt((x[i+1] - x[i])^2 + (y[i+1] - y[i])^2)
            sin_theta[i] = (y[i+1] - y[i]) / distance[i]
            cos_theta[i] = (x[i+1] - x[i]) / distance[i]
            if distance[i] == 0.0
                error("Zero distance detected between points $i and $(i+1).")
            end
        end
        return sin_theta, cos_theta
    end

    # sin_theta, cos_theta : vectors : The sin and cos of each panel with respect to the x-axis -- length(n)
    sin_theta, cos_theta = find_sin_cos(x, y)

    # Find rijs...

    function find_rijs(x, y, x_mid, y_mid)  # lengths: 131, 131, 130, 130
        r_ij = similar(x_mid, length(x_mid), length(x_mid)+1)
        for i in eachindex(x_mid)
            for j in eachindex(x)
                r_ij[i, j] = sqrt((x_mid[i] - x[j])^2 + (y_mid[i] - y[j])^2)
            end
        end
        return r_ij   # 130 x 131
    end

    # r_ij : matrix : Find the distance between the control point (mid_point) of panel i with respect to source points j and j+1 -- size(n, n+1)
    r_ij = find_rijs(x, y, x_mid, y_mid)

    # Find sin(theta i - theta j) etc.
    function find_thetas(sin_theta, cos_theta)
        cos_theta_ij = similar(cos_theta, length(cos_theta), length(cos_theta))
        sin_theta_ij = similar(sin_theta, length(sin_theta), length(sin_theta)) 

        for i in eachindex(sin_theta)
            for j in eachindex(cos_theta)
                sin_theta_ij[i, j] = sin_theta[i] * cos_theta[j] - cos_theta[i] * sin_theta[j]
                cos_theta_ij[i, j] = cos_theta[i] * cos_theta[j] + sin_theta[i] * sin_theta[j]
            end
        end
        return sin_theta_ij, cos_theta_ij   # 130 x 130 matrices
    end

    # sin_theta_ij, cos_theta_ij : using sin_theta and cos_theta to find sin(theta[i] - theta[j]) and cos(theta[i] - theta[j])
    sin_theta_ij, cos_theta_ij = find_thetas(sin_theta, cos_theta)

    # Find Beta
    function find_beta(x, y, x_mid, y_mid)
        beta = similar(x_mid, length(x_mid),length(x_mid))
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

    # beta : matrix : use coordinates and midpoints to find the angle between control point of panel i and source coordinates j and j+1 -- size(n, n)
    beta = find_beta(x, y, x_mid, y_mid)

    # Find Aij
    function find_A(r_ij, sin_theta_ij, cos_theta_ij, beta)
        A = similar(r_ij, size(r_ij,1)+1, size(r_ij,2))

        # These for loops give us the values for the [n x n-1] matrix
        for i in eachindex(sin_theta_ij[:, 1])
            for j in eachindex(cos_theta_ij[1, :])
                rs = r_ij[i, j+1] / r_ij[i, j]
                A[i, j] = log(rs) * sin_theta_ij[i, j] + beta[i, j] * cos_theta_ij[i, j] 

                # check to make sure the value is a real number
                if any(isnan, A[i, j]) || any(isinf, A[i, j])
                    println("$i , $j")
                    error("A[i ,end] contains NaN or Inf values.")
                end               
            end
        end

        # Gives the values for the nth column 
        for i in eachindex(sin_theta_ij[:, 1])
            A[i, end] = 0.0
            for j in eachindex(cos_theta_ij[1, :])
                rs = r_ij[i, j+1] / r_ij[i, j]
                A[i,end] += log(rs) * cos_theta_ij[i, j] - beta[i, j] * sin_theta_ij[i, j]
                
                # check to make sure the value is a real number
                if any(isnan, A[i, end]) || any(isinf, A[i, end])
                    println("$i , $j")
                    println(A[i, end])
                    error("A[i ,end] contains NaN or Inf values.")
                end
            end
        end

        # This gives us the n+1 row
        for j in eachindex(A[1,1:end-1])
            sin_theta_k1 = sin_theta_ij[1, j]
            cos_theta_k1 = cos_theta_ij[1, j]
            sin_theta_kn = sin_theta_ij[end, j]
            cos_theta_kn = cos_theta_ij[end, j]

            betak1 = beta[1, j]
            betakn = beta[end, j]

            r1j = r_ij[1, j]
            rnj = r_ij[end, j]
            r1j1 = r_ij[1, j+1]
            rnj1 = r_ij[end, j+1]
        
            k1 = betak1 * sin_theta_k1 - log(r1j1 / r1j) * cos_theta_k1
            kn = betakn * sin_theta_kn - log(rnj1 / rnj) * cos_theta_kn

            A[end, j] = k1 + kn

            if any(isnan, A[end, j]) || any(isinf, A[end, j])
                error("A[end, j] contains NaN or Inf values.")
            end

        end

        # Gives us the [n+1,n+1] value
        A[end, end] = 0.0
        for j in eachindex(A[1, 1:end-1])
            sin_theta_k1 = sin_theta_ij[1, j]
            cos_theta_k1 = cos_theta_ij[1, j]
            sin_theta_kn = sin_theta_ij[end, j]
            cos_theta_kn = cos_theta_ij[end, j]

            betak1 = beta[1, j]
            betakn = beta[end, j]

            r1j = r_ij[1, j]
            rnj = r_ij[end, j]
            r1j1 = r_ij[1, j+1]
            rnj1 = r_ij[end, j+1]
        
            A[end, end] += betak1 * cos_theta_k1 + log(r1j1 / r1j) * sin_theta_k1
            A[end, end] += betakn * cos_theta_kn + log(rnj1 / rnj) * sin_theta_kn

            # check to make sure the value is a real number
            if any(isnan, A[end, end]) || any(isinf, A[end, end])
                error("A[end, end] contains NaN or Inf values.")
            end
        end

        return A 
    end

    # A : the matrix that is the influence of every panel on every other panel -- size(n+1, n+1)
    A = find_A(r_ij, sin_theta_ij, cos_theta_ij, beta)

    # Find b vector (boundary conditions -- using no-flow condition and Kutta condition)
    function find_b(sin_theta, cos_theta, V_inf, alpha)
        b = similar(sin_theta, length(sin_theta)+1)

        for i in eachindex(sin_theta)
            b[i] = 2 * π * V_inf * (sin_theta[i] * cos(alpha) - cos_theta[i] * sin(alpha))
        end

        b[end] = -2 * π * V_inf * ((cos_theta[1] * cos(alpha) + sin_theta[1] * sin(alpha)) + (cos_theta[end] * cos(alpha) + sin_theta[end] * sin(alpha)))

        return b
    end

    # b : vector that provides the boundary conditions -- length(n+1)
    b = find_b(sin_theta, cos_theta, V_inf, alpha) 

    # q_gamma : vector : solves for the source strengths of each point and the circulation strength of the body : length(n+1) where (n+1) = gamma
    q_gamma = A \ b

    # Find tangential velocity at each panel
    function find_vt(r_ij, sin_theta_ij, cos_theta_ij, beta, q_gamma, V_inf, alpha)
        Vti = similar(q_gamma, length(q_gamma)-1)
        for i in eachindex(q_gamma[1:end-1])
            set1 = 0.0
            set2 = 0.0
            for j in eachindex(q_gamma[1:end-1])
                set1 += q_gamma[j] * (beta[i, j] * sin_theta_ij[i, j] - log(r_ij[i, j+1] / r_ij[i, j]) * cos_theta_ij[i, j])
                set2 += beta[i, j] * cos_theta_ij[i, j] + log(r_ij[i, j+1] / r_ij[i, j]) * sin_theta_ij[i, j] 
            end
            Vti[i] = V_inf * cos_theta[i] * cos(alpha) + (set1 / (2π)) + (q_gamma[end]/(2π)) * set2
        end
        return Vti
    end

    # Vti : vector : finds the tangential velocity at each panel since we have a no-flow through condition so there is no normal velocity -- length(n)
    Vti = find_vt(r_ij, sin_theta_ij, cos_theta_ij, beta, q_gamma, V_inf, alpha)

    # Find CP for each point
    function cpressure(Vti, V_inf)
        CP = 1 .- (Vti ./ V_inf) .^2
        return CP 
    end

    # CP : vector : Coefficient of Pressure at each control point (midpoint of each panel) -- length(n)
    CP = cpressure(Vti, V_inf)

    #= 
    plot each coefficient of pressure at its corresponding panel control point (x_mid); 
    flip the y_axis so the top curve is the upper surface and the bottom curve
    is the lower surface for a positively cambered airfoil
    =#
    # plcp = plot(x_mid, CP, yflip = true, markers = true)

    return x_mid, CP
end

########### Validate with Joukowsky #################

using FLOWFoil
using Plots
using .AirfoilTools

# - Parameters - #
center = [-0.1; 0.1]
radius = 1.0
alpha = 4.0
Vinf = 1.0

# - Joukowsky Geometry - #
x, y = FLOWFoil.AirfoilTools.joukowsky(center, radius, N = 720)

# - Surface Values - #
surface_velocity, surface_pressure_coefficient, cl = FLOWFoil.AirfoilTools.joukowsky_flow(
    center, radius, alpha, Vinf
)

# - Your Stuff - #
x_mid, CP = HS_Panel_CP(x, y, Vinf, deg2rad(alpha)) 

# - Plot Stuff - #
pl = plot(; xlabel="x", ylabel="cp", yflip=true)
# plot!(
#     pl,
#     x,
#     surface_pressure_coefficient;
#     linestyle=:dash,
#     linewidth=2,
#     label="Analytic Solution",
# )
plot!(pl, x_mid[2:end-1], CP[2:end-1], label="Hess-Smith")

################################################