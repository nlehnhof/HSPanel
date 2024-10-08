using Plots
using LinearAlgebra

# Freestream Velocity and AoA
alpha = 0.0 * π/180
V_inf = 10.0

################ Get Coordinates from File ############

function get_coordinates(file)
    x, y = open(file, "r") do f
        x = Float64[]
        y = Float64[]
        for line in eachline(f)
            entries = split(chomp(line))
            push!(x, parse(Float64, entries[1]))
            push!(y, parse(Float64, entries[2]))
        end
        x, y
    end
end

x, y = get_coordinates("naca0008.txt")

# ###### Generate NACA 0012 coordinates ######

# # Function to generate cosine-spaced x-coordinates
# function cosine_spacing(n_points)
#     β = range(0, stop=π, length=n_points)
#     x = (1 .- cos.(β)) ./ 2
#     return x
# end

# # Function to calculate thickness distribution
# function thickness_distribution(x, t=0.12)
#     yt = 5 * t * (0.2969 .* sqrt.(x) .- 0.1260 .* x .- 
#                   0.3516 .* x.^2 .+ 0.2843 .* x.^3 .- 
#                   0.1015 .* x.^4)
#     return yt
# end

# # Number of points per surface (excluding the trailing point to avoid duplication)
# n_points = 100

# # Generate cosine-spaced x-coordinates
# x = cosine_spacing(n_points + 1)  # +1 to include the trailing edge

# # Calculate thickness distribution
# yt = thickness_distribution(x)

# # Upper surface coordinates
# x_upper = x
# y_upper = yt

# # Lower surface coordinates
# x_lower = x
# y_lower = -yt

# # Combine upper and lower surfaces
# # Exclude the first and last points of the lower surface to avoid duplication at leading and trailing edges
# x_coords = vcat(x_upper, reverse(x_lower[2:end-1]))
# y_coords = vcat(y_upper, reverse(y_lower[2:end-1]))

# # Plot the airfoil to verify
# plot(x_coords, y_coords, aspect_ratio=1, label="NACA0012 Airfoil", xlabel="x", ylabel="y", title="NACA0012 Airfoil Coordinates", legend=false)


########## Panel Method ###################
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

x_mid, y_mid = find_midpoints(x, y)

function find_length(x, y)
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

sin_theta, cos_theta = find_length(x, y)

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

r_ij = find_rijs(x, y, x_mid, y_mid)

# Find sin(theta i - theta j) etc.

function find_thetas(sin_theta, cos_theta)
    cos_theta_ij = similar(cos_theta, length(cos_theta), length(cos_theta))
    sin_theta_ij = similar(sin_theta, length(sin_theta), length(sin_theta)) 
    
    for i in eachindex(sin_theta)
        for j in eachindex(cos_theta)
            sin_theta_ij[i, j] = sin_theta[i] * sin_theta[j] - cos_theta[i] * cos_theta[j]
            cos_theta_ij[i, j] = cos_theta[i] * cos_theta[j] + sin_theta[i] * sin_theta[j]
        end
    end
    return sin_theta_ij, cos_theta_ij   # 130 x 130 matrices
end

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
    return beta   # 130 x 130
end

beta = find_beta(x, y, x_mid, y_mid)  # CORRECT

# Find Aij

function find_A(r_ij, sin_theta_ij, cos_theta_ij, beta)
    A = similar(r_ij, size(r_ij,1)+1, size(r_ij,2))   # 131 x 131

    # These for loops give us the values for the 130 x 129 matrix
    for i in eachindex(sin_theta_ij[:, 1])
        for j in eachindex(cos_theta_ij[1, :])
            rs = r_ij[i, j+1] / r_ij[i, j]
            A[i, j] = log(rs) * sin_theta_ij[i, j] + beta[i, j] * cos_theta_ij[i, j] 
            if any(isnan, A[i, j]) || any(isinf, A[i, j])
                println("$i , $j")
                error("A[i ,end] contains NaN or Inf values.")
            end               
        end
    end

    for i in eachindex(sin_theta_ij[:, 1])
        A[i, end] = 0.0
        for j in eachindex(cos_theta_ij[1, :])
            rs = r_ij[i, j+1] / r_ij[i, j]
            A[i,end] += log(rs) * cos_theta_ij[i, j] - beta[i, j] * sin_theta_ij[i, j]
            if any(isnan, A[i, end]) || any(isinf, A[i, end])
                println("$i , $j")
                println(A[i, end])
                error("A[i ,end] contains NaN or Inf values.")
            end
        end
    end

    # This gives us the 131st row  (A_N+1,j)
 
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

    # Gives us (A_n+1,n+1)

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
        if any(isnan, A[end, end]) || any(isinf, A[end, end])
            error("A[end, end] contains NaN or Inf values.")
        end
    end

    return A 
end

A = find_A(r_ij, sin_theta_ij, cos_theta_ij, beta)

# Find b vector (boundary conditions)

function find_b(sin_theta, cos_theta, V_inf, alpha)
    b = similar(sin_theta, length(sin_theta)+1)

    for i in eachindex(sin_theta)
        b[i] = 2 * π * V_inf * (sin_theta[i] * sin(alpha) - cos_theta[i] * cos(alpha))
    end

    b[end] = -2 * π * V_inf * ((cos_theta[1] * cos(alpha) + sin_theta[1] * sin(alpha)) + (cos_theta[end] * cos(alpha) + sin_theta[end] * sin(alpha)))

    return b
end

b = find_b(sin_theta, cos_theta, V_inf, alpha) # 131

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

Vti = find_vt(r_ij, sin_theta_ij, cos_theta_ij, beta, q_gamma, V_inf, alpha)


# Find CP for each point

function cpressure(Vti, V_inf)
    CP1 = 1 .- (Vti ./ V_inf) .^2
    CP = [round(CP1[i], digits = 4) for i in 1:length(Vti)]
    return CP 
end

CP = cpressure(Vti, V_inf)
println(CP)

plot(x_mid, CP, markers = true)

