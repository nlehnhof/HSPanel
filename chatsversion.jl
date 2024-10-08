using Plots

# -------------------------------
# 1. Input NACA0012 Coordinates
# -------------------------------

# NACA0012 Coordinates
coordinates = [
    (1.000000,  0.001260),
    (0.999416,  0.001342),
    (0.997666,  0.001587),
    (0.994753,  0.001994),
    (0.990685,  0.002560),
    (0.985471,  0.003280),
    (0.979123,  0.004152),
    (0.971656,  0.005169),
    (0.963087,  0.006324),
    (0.953437,  0.007611),
    (0.942728,  0.009022),
    (0.930985,  0.010549),
    (0.918235,  0.012182),
    (0.904509,  0.013914),
    (0.889837,  0.015735),
    (0.874255,  0.017635),
    (0.857800,  0.019605),
    (0.840508,  0.021635),
    (0.822421,  0.023714),
    (0.803581,  0.025834),
    (0.784032,  0.027983),
    (0.763820,  0.030152),
    (0.742992,  0.032329),
    (0.721596,  0.034506),
    (0.699682,  0.036670),
    (0.677303,  0.038811),
    (0.654509,  0.040917),
    (0.631354,  0.042978),
    (0.607892,  0.044980),
    (0.584179,  0.046912),
    (0.560268,  0.048762),
    (0.536217,  0.050516),
    (0.512082,  0.052162),
    (0.487918,  0.053687),
    (0.463783,  0.055077),
    (0.439732,  0.056320),
    (0.415822,  0.057403),
    (0.392108,  0.058314),
    (0.368646,  0.059042),
    (0.345492,  0.059575),
    (0.322698,  0.059903),
    (0.300318,  0.060017),
    (0.278404,  0.059910),
    (0.257008,  0.059576),
    (0.236180,  0.059008),
    (0.215968,  0.058205),
    (0.196419,  0.057164),
    (0.177579,  0.055886),
    (0.159492,  0.054372),
    (0.142201,  0.052625),
    (0.125745,  0.050651),
    (0.110163,  0.048457),
    (0.095492,  0.046049),
    (0.081765,  0.043437),
    (0.069015,  0.040631),
    (0.057272,  0.037641),
    (0.046563,  0.034479),
    (0.036913,  0.031156),
    (0.028344,  0.027683),
    (0.020877,  0.024071),
    (0.014529,  0.020330),
    (0.009315,  0.016471),
    (0.005247,  0.012501),
    (0.002334,  0.008429),
    (0.000584,  0.004260),
    (0.000000,  0.000000),
    (0.000584, -0.004260),
    (0.002334, -0.008429),
    (0.005247, -0.012501),
    (0.009315, -0.016471),
    (0.014529, -0.020330),
    (0.020877, -0.024071),
    (0.028344, -0.027683),
    (0.036913, -0.031156),
    (0.046563, -0.034479),
    (0.057272, -0.037641),
    (0.069015, -0.040631),
    (0.081765, -0.043437),
    (0.095492, -0.046049),
    (0.110163, -0.048457),
    (0.125745, -0.050651),
    (0.142201, -0.052625),
    (0.159492, -0.054372),
    (0.177579, -0.055886),
    (0.196419, -0.057164),
    (0.215968, -0.058205),
    (0.236180, -0.059008),
    (0.257008, -0.059576),
    (0.278404, -0.059910),
    (0.300318, -0.060017),
    (0.322698, -0.059903),
    (0.345492, -0.059575),
    (0.368646, -0.059042),
    (0.392108, -0.058314),
    (0.415822, -0.057403),
    (0.439732, -0.056320),
    (0.463783, -0.055077),
    (0.487918, -0.053687),
    (0.512082, -0.052162),
    (0.536217, -0.050516),
    (0.560268, -0.048762),
    (0.584179, -0.046912),
    (0.607892, -0.044980),
    (0.631354, -0.042978),
    (0.654509, -0.040917),
    (0.677303, -0.038811),
    (0.699682, -0.036670),
    (0.721596, -0.034506),
    (0.742992, -0.032329),
    (0.763820, -0.030152),
    (0.784032, -0.027983),
    (0.803581, -0.025834),
    (0.822421, -0.023714),
    (0.840508, -0.021635),
    (0.857800, -0.019605),
    (0.874255, -0.017635),
    (0.889837, -0.015735),
    (0.904509, -0.013914),
    (0.918235, -0.012182),
    (0.930985, -0.010549),
    (0.942728, -0.009022),
    (0.953437, -0.007611),
    (0.963087, -0.006324),
    (0.971656, -0.005169),
    (0.979123, -0.004152),
    (0.985471, -0.003280),
    (0.990685, -0.002560),
    (0.994753, -0.001994),
    (0.997666, -0.001587),
    (0.999416, -0.001342),
    (1.000000, -0.001260)
]

# Extract x and y coordinates
x = [coord[1] for coord in coordinates]
y = [coord[2] for coord in coordinates]

# -------------------------------
# 2. Implement the Panel Method
# -------------------------------

# 2.1. Find Midpoints
function find_midpoints(x, y)
    n = length(x) - 1
    x_mid = zeros(Float64, n)
    y_mid = zeros(Float64, n)
    for i in 1:n
        x_mid[i] = 0.5 * (x[i] + x[i+1])
        y_mid[i] = 0.5 * (y[i] + y[i+1])
    end
    return x_mid, y_mid
end

x_mid, y_mid = find_midpoints(x, y)

# 2.2. Calculate Panel Lengths and Angles
function find_length(x, y)
    n = length(x) - 1
    distance = zeros(Float64, n)
    sin_theta = zeros(Float64, n)
    cos_theta = zeros(Float64, n)
    for i in 1:n
        dx = x[i+1] - x[i]
        dy = y[i+1] - y[i]
        distance[i] = sqrt(dx^2 + dy^2)
        sin_theta[i] = dy / distance[i]
        cos_theta[i] = dx / distance[i]
    end
    return distance, sin_theta, cos_theta
end

distance, sin_theta, cos_theta = find_length(x, y)

# 2.3. Compute Influence Coefficients (A Matrix)
function find_rijs(x_mid, y_mid, x, y)
    n = length(x_mid)
    m = length(x)
    r_ij = zeros(Float64, n, m)
    for i in 1:n
        for j in 1:m
            r_ij[i, j] = sqrt((x_mid[i] - x[j])^2 + (y_mid[i] - y[j])^2)
        end
    end
    return r_ij
end

r_ij = find_rijs(x_mid, y_mid, x, y)

# 2.4. Find Theta Angles for Influence Coefficients
function find_thetas(sin_theta, cos_theta)
    n = length(sin_theta)
    sin_theta_ij = zeros(Float64, n, n)
    cos_theta_ij = zeros(Float64, n, n)
    for i in 1:n
        for j in 1:n
            sin_theta_ij[i, j] = sin_theta[i] * sin_theta[j] - cos_theta[i] * cos_theta[j]
            cos_theta_ij[i, j] = cos_theta[i] * cos_theta[j] + sin_theta[i] * sin_theta[j]
        end
    end
    return sin_theta_ij, cos_theta_ij
end

sin_theta_ij, cos_theta_ij = find_thetas(sin_theta, cos_theta)

# 2.5. Compute Beta Angles
function find_beta(x, y, x_mid, y_mid)
    n = length(x_mid)
    beta = zeros(Float64, n, n)
    for i in 1:n
        for j in 1:n
            if j == i
                beta[i, j] = π
            else
                numerator = (x[j] - x_mid[i]) * (y[j + 1] - y_mid[i]) - (y[j] - y_mid[i]) * (x[j + 1] - x_mid[i])
                denominator = (x[j] - x_mid[i]) * (x[j + 1] - x_mid[i]) + (y[j] - y_mid[i]) * (y[j + 1] - y_mid[i])
                beta[i, j] = atan(numerator, denominator)
            end
        end
    end
    return beta
end

beta = find_beta(x, y, x_mid, y_mid)

# 2.6. Assemble the A Matrix
function find_A(r_ij, sin_theta_ij, cos_theta_ij, beta)
    n = size(r_ij, 1)
    m = size(r_ij, 2)
    A = zeros(Float64, n + 1, m)  # Extra row for Kutta condition

    # Fill A matrix for each panel
    for i in 1:n
        for j in 1:(m - 1)
            A[i, j] = log(r_ij[i, j + 1] / r_ij[i, j]) * sin_theta_ij[i, j] + beta[i, j] * cos_theta_ij[i, j]
            A[i, end] += log(r_ij[i, j + 1] / r_ij[i, j]) * cos_theta_ij[i, j] - beta[i, j] * sin_theta_ij[i, j]
        end
    end

    # Apply Kutta condition (ensuring smooth flow at trailing edge)
    for j in 1:(m - 1)
        sin_theta_k1 = sin_theta_ij[1, j]
        cos_theta_k1 = cos_theta_ij[1, j]
        sin_theta_kn = sin_theta_ij[end, j]
        cos_theta_kn = cos_theta_ij[end, j]

        betak1 = beta[1, j]
        betakn = beta[end, j]

        r1j = r_ij[1, j]
        rnj = r_ij[end, j]
        r1j1 = r_ij[1, j + 1]
        rnj1 = r_ij[end, j + 1]

        k1 = betak1 * sin_theta_k1 - log(r1j1 / r1j) * cos_theta_k1
        kn = betakn * sin_theta_kn - log(rnj1 / rnj) * cos_theta_kn

        A[end, j] = k1 + kn
    end

    # Compute A[end, end]
    for j in 1:(m - 1)
        sin_theta_k1 = sin_theta_ij[1, j]
        cos_theta_k1 = cos_theta_ij[1, j]
        sin_theta_kn = sin_theta_ij[end, j]
        cos_theta_kn = cos_theta_ij[end, j]

        betak1 = beta[1, j]
        betakn = beta[end, j]

        r1j = r_ij[1, j]
        rnj = r_ij[end, j]
        r1j1 = r_ij[1, j + 1]
        rnj1 = r_ij[end, j + 1]

        A[end, end] += betak1 * cos_theta_k1 + log(r1j1 / r1j) * sin_theta_k1
        A[end, end] += betakn * cos_theta_kn + log(rnj1 / rnj) * sin_theta_kn
    end

    return A
end

A = find_A(r_ij, sin_theta_ij, cos_theta_ij, beta)

# 2.7. Assemble the Boundary Conditions Vector (b)
function find_b(sin_theta, cos_theta, V_inf, alpha)
    n = length(sin_theta)
    b = zeros(Float64, n + 1)
    for i in 1:n
        b[i] = 2 * π * V_inf * (sin_theta[i] * sin(alpha) - cos_theta[i] * cos(alpha))
    end
    b[end] = -2 * π * V_inf * ((cos_theta[1] * cos(alpha) + sin_theta[1] * sin(alpha)) +
                                (cos_theta[end] * cos(alpha) + sin_theta[end] * sin(alpha)))
    return b
end

# Freestream conditions
alpha = 0.0 * π / 180  # Angle of attack in radians
V_inf = 10.0           # Freestream velocity

b = find_b(sin_theta, cos_theta, V_inf, alpha)

# 2.8. Solve for Circulation (gamma)
q_gamma = A \ b

# -------------------------------
# 3. Compute Pressure Coefficient (C_P)
# -------------------------------

# 3.1. Find Tangential Velocity at Each Panel
function find_vt(r_ij, sin_theta_ij, cos_theta_ij, beta, q_gamma, V_inf, alpha)
    n = length(q_gamma) - 1
    Vti = zeros(Float64, n)
    for i in 1:n
        set1 = 0.0  # Reset for each panel
        set2 = 0.0  # Reset for each panel
        for j in 1:n
            set1 += q_gamma[j] * (beta[i, j] * sin_theta_ij[i, j] - log(r_ij[i, j + 1] / r_ij[i, j]) * cos_theta_ij[i, j])
            set2 += beta[i, j] * cos_theta_ij[i, j] + log(r_ij[i, j + 1] / r_ij[i, j]) * sin_theta_ij[i, j]
        end
        Vti[i] = V_inf * cos_theta[i] * cos(alpha) + (set1 / (2π)) + (q_gamma[end] / (2π)) * set2
    end
    return Vti
end

Vti = find_vt(r_ij, sin_theta_ij, cos_theta_ij, beta, q_gamma, V_inf, alpha)

# 3.2. Calculate Pressure Coefficient (C_P)
function cpressure(Vti, V_inf)
    CP = 1 .- (Vti ./ V_inf) .^ 2
    return CP
end

CP = cpressure(Vti, V_inf)

plot(x_mid, CP, yflip = true, markers = true)