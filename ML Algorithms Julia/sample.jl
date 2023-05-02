using Random, Distributions
using DelimitedFiles, Formatting
"For a given privacy level, sample Gaussian noise."
function Gaussian_noise(epsilon::Float64, delta::Float64, sensitivity::Float64)
    if epsilon <= 0 || delta <= 0 || sensitivity < 0
        error("Revise the privacy parameters.")
    end
    mu_Gauss = 0.0
    var_gauss = (2*log(1.25 / delta)*(sensitivity^2)) / (epsilon^2)
    return rand(Normal(mu_Gauss, sqrt(var_gauss)), 1)[1] #sample from normal dist
end

"For a given privacy level, sample Truncated Laplace noise."
function TLAP_noise(epsilon::Float64, delta::Float64, sensitivity::Float64)
    if epsilon <= 0 || delta <= 0 || sensitivity < 0
        error("Revise the privacy parameters.")
    end
    lam = sensitivity/epsilon
    A = lam * log(1 + (exp(epsilon) - 1)/(2 * delta))
    B = 1/(2 * lam * (1 - exp(- A / lam)))
    u = rand(Uniform(0,1),1)[1]
    return sign(u - 0.5)*lam*log(min(u, 1-u)/(B*lam) + exp(-A/lam))
end

"For a given privacy level, call the corresponding distribution."
function call_optimal_dist(epsilon::Float64, delta::Float64, sensitivity::Float64)
    if epsilon <= 0 || delta <= 0 || sensitivity < 0
        error("Revise the privacy parameters.")
    end
    #read the noise break-points
    X = vec(readdlm("Distributions/eps"*string(round(epsilon,digits=1))*"delta"*string(round(delta,digits=1))*
    "/X_use"*sprintf1("%.2f", sensitivity )*".csv", ',', Float64))
    beta = round(X[2] - X[1],digits = 6) #take the beta
    #read the optimized mixture weights
    p = vec(readdlm("Distributions/eps"*string(round(epsilon,digits=1))*"delta"*string(round(delta,digits=1))*
    "/p_use"*sprintf1("%.2f", sensitivity )*".csv", ',', Float64))
    return X, p, beta
end

"For a given privacy level, call the corresponding multiple distribution."
function call_optimal_dist_mult(epsilon::Float64, delta::Float64, sensitivity::Float64, F_lower::Float64, F_upper::Float64, truth::Float64)
    #read the noise break-points
    X = vec(readdlm("Distributions/new/X_sens"*sprintf1("%.2f", sensitivity)*"lower"*sprintf1("%.2f", F_lower)*
        "upper"*sprintf1("%.2f", F_upper)*".csv", ',', Float64))
    beta = round(X[2] - X[1],digits = 6) #take the beta
    #read range break-points
    ranges_breaks = vec(readdlm("Distributions/new/ranges_sens"*sprintf1("%.2f", sensitivity)*"lower"*sprintf1("%.2f", F_lower)*
        "upper"*sprintf1("%.2f", F_upper)*".csv", ',', Float64))
    dist_index = findfirst(isequal(1), truth - 10^(-6) .<= ranges_breaks) #first
    if typeof(dist_index) == Nothing
        dist_index = length(ranges_breaks)
        # error("Sensitivity "*string(sensitivity)*" truth "*string(truth)*" F_lower "*string(F_lower)*" F_upper "*string(F_upper)*
        # "range does not fit the truth.") #later see the source of error message and fix
    end
    #read the optimized mixture weights
    p = readdlm("Distributions/new/p_sens"*sprintf1("%.2f", sensitivity)*"lower"*sprintf1("%.2f", F_lower)*
        "upper"*sprintf1("%.2f", F_upper)*".csv", ',', Float64)
    p = vec(p[dist_index, :])
    return X, p, beta
end

"For a given privacy level, call the corresponding multiple distribution."
function call_optimal_dist_mult_efficient(epsilon::Float64, delta::Float64, sensitivity::Float64, F_lower::Float64, F_upper::Float64, truth::Float64, Xs,ranges,ps,betas)
    lookup = sprintf1("%.2f", sensitivity)*"lower"*sprintf1("%.2f", F_lower)*"upper"*sprintf1("%.2f", F_upper)*".csv"
    #read the noise break-points
    X,beta,p,ranges_breaks = Xs[lookup], betas[lookup], ps[lookup], ranges[lookup]
    dist_index = findfirst(isequal(1), truth - 10^(-6) .<= ranges_breaks) #first
    if typeof(dist_index) == Nothing
        dist_index = length(ranges_breaks)
    end
    #read the optimized mixture weights
    p = vec(p[dist_index, :])
    return X, p, beta
end

"For a given privacy level, sample from the optimized noise."
function optimal_noise(X::Vector{Float64}, p::Vector{Float64}, beta::Float64)
    u = rand(Uniform(0,1),1)[1]
    #find first index where p gets larger than u, means we are in that interval
    cumul_sum = 0.0
    index_where = 0
    for (index, value_p) in enumerate(p)
        cumul_sum += value_p
        if cumul_sum >= u - 10^-6
            index_where = index
            break
        end
    end
    #index_where = findfirst(cumsum(p) .>= u) #slower version of the above code
    sample = X[index_where] + (beta * rand(Uniform(0,1),1)[1]) #go to the interval that starts with X[index_where] and move a bit
    return sample
end
