"Compute sensitivities of the numerical columns -- to be used later. sensitivities_mean[1] will give the sensitivities of col 1 under different labels."
function widths(df_cont::DataFrame, prior_counts::Dict{Int64, Int64})
    n,m = size(df_cont) #n is the number of rows, m is the number of numerical columns
    maxs = [maximum(df_cont[:,col_index]) for col_index in 1:m] #column-wise max values
    mins = [minimum(df_cont[:,col_index]) for col_index in 1:m]
    max_diff = maxs-mins #note: we can call C++ by using these
    #above was all about the width of the numerical columns. But the sensitivities will depend on the labels. Below we do that.
    prior_vector = [prior_counts[prior] for prior in 1:length(prior_counts)] #prior counts to use -- make them in order

    sensitivities_mean = [ceil.(max_diff[col_index]./prior_vector) for col_index = 1:m] #TODO change later to 2 digit precision
    sensitivities_sd = [ceil.(max_diff[col_index]./sqrt.(prior_vector)) for col_index = 1:m]#TODO change later to 2 digit precision

    max_sds_possible = [std([repeat([mins[col_index]], floor(Int64, n/2)); repeat([maxs[col_index]], ceil(Int64, n/2)) ] , corrected = false) for col_index = 1:m]

    return mins, maxs, max_diff, sensitivities_mean, sensitivities_sd, max_sds_possible #TODO power way more than amplitude -- report.
end

"Return noisy variants of the counts."
function noisy_counts(epsilon, delta, df_counts_given, prior_counts_given, method)
    ## Case 1: no noise
    if method == "null"
        return df_counts, prior_counts
    end
    df_counts = copy(df_counts_given)
    prior_counts = copy(prior_counts_given)
    ## compute noisy df counts
    X, p, beta = call_optimal_dist(epsilon, delta, 1.0) #optimal dist
    df_counts = convert(Vector{Vector{Dict{Int64,Float64}}}, df_counts) #convert
    for (col_index, col_dicts) in enumerate(df_counts)
        for (label_index, label_dict) in enumerate(col_dicts)
            for (key, value) in label_dict
                if method == "optimal"
                    df_counts[col_index][label_index][key] = max(value + optimal_noise(X, p, beta), 10^-6)
                elseif method == "truncated"
                    df_counts[col_index][label_index][key] = max(value +  TLAP_noise(epsilon, delta, 1.0),  10^-6)
                elseif method == "gaussian"
                    df_counts[col_index][label_index][key] = max(value + Gaussian_noise(epsilon, delta, 1.0),  10^-6)
                elseif method == "analytic" #new
                    df_counts[col_index][label_index][key] = max(value + rand(Normal(0.0, 1.0858777651918035), 1)[1],  10^-6)
                else
                    error("the method is not defined.")
                end
            end
        end
    end
    #add noise to prior_counts
    prior_counts = convert(Dict{Int64, Float64}, prior_counts) #convert
    for (key, value) in prior_counts
        if method == "optimal"
            prior_counts[key] = max(value + optimal_noise(X, p, beta), 10^-6)
        elseif method == "truncated"
            prior_counts[key] = max(value + TLAP_noise(epsilon, delta, 1.0), 10^-6)
        elseif method == "gaussian"
            prior_counts[key] = max(value + Gaussian_noise(epsilon, delta, 1.0), 10^-6)
        elseif method == "analytic" #new
            prior_counts[key]  = max(value + rand(Normal(0.0, 1.0858777651918035), 1)[1],  10^-6)
        else
            error("the method is not defined.")
        end
    end
    return df_counts, prior_counts
end


function noisy_means_sds(epsilon, delta, sensitivities_mean, sensitivities_sd, df_means_given, df_sds_given, method)
    ## Case 1: no noise
    if method == "null"
        return df_means, df_sds
    end
    df_means = copy.(df_means_given) #otherwise not immutable
    df_sds = copy.(df_sds_given)
    ## compute noisy mean/sds
    for (col_index, col_vec) in enumerate(df_means)
        for (label_index, label_val) in enumerate(col_vec)
            if method == "truncated"
                df_means[col_index][label_index] = df_means[col_index][label_index] + TLAP_noise(epsilon, delta, sensitivities_mean[col_index][label_index])
                df_sds[col_index][label_index] = max(df_sds[col_index][label_index] + TLAP_noise(epsilon, delta, sensitivities_sd[col_index][label_index]),10^0)
            elseif method == "gaussian"
                df_means[col_index][label_index] = df_means[col_index][label_index] + Gaussian_noise(epsilon, delta, sensitivities_mean[col_index][label_index])
                df_sds[col_index][label_index] = max(df_sds[col_index][label_index] + Gaussian_noise(epsilon, delta, sensitivities_sd[col_index][label_index]),10^0)
            elseif method == "optimal"
                X, p, beta = call_optimal_dist(epsilon, delta, sensitivities_mean[col_index][label_index]) #optimal dist
                df_means[col_index][label_index] = df_means[col_index][label_index] + optimal_noise(X, p, beta)
                X, p, beta = call_optimal_dist(epsilon, delta, sensitivities_sd[col_index][label_index]) #optimal dist
                df_sds[col_index][label_index] = max(df_sds[col_index][label_index] + optimal_noise(X, p, beta), 10^0)
            elseif method == "analytic"
                df_means[col_index][label_index] = df_means[col_index][label_index] + rand(Normal(0.0, calibrateAnalyticGaussianMechanism(epsilon, delta, sensitivities_mean[col_index][label_index])), 1)[1]
                df_sds[col_index][label_index] = max(df_sds[col_index][label_index] + rand(Normal(0.0, calibrateAnalyticGaussianMechanism(epsilon, delta, sensitivities_sd[col_index][label_index])), 1)[1],10^0)
            end
        end
    end
    return df_means, df_sds
end

"Return noisy variants of the counts."
function noisy_counts_mult(epsilon, delta, df_counts_given, prior_counts_given, F_uppers,Xs, ranges, ps, betas)
    ## Case 1: no noise
    df_counts = copy(df_counts_given)
    prior_counts = copy(prior_counts_given)
    ## compute noisy df counts
    df_counts = convert(Vector{Vector{Dict{Int64,Float64}}}, df_counts) #convert
    for (col_index, col_dicts) in enumerate(df_counts)
        for (label_index, label_dict) in enumerate(col_dicts)
            for (key, value) in label_dict
                #X, p, beta = call_optimal_dist_mult(epsilon, delta, 1.0, 0.0, 1.0*F_uppers[1], 1.0*value)
                X, p, beta = call_optimal_dist_mult_efficient(epsilon, delta, 1.0, 0.0, 1.0*F_uppers[1], 1.0*value,Xs, ranges, ps, betas)
                df_counts[col_index][label_index][key] = max(value + optimal_noise(X, p, beta), 10^-6)
            end
        end
    end
    #add noise to prior_counts
    prior_counts = convert(Dict{Int64, Float64}, prior_counts) #convert
    for (key, value) in prior_counts
        #X, p, beta = call_optimal_dist_mult(epsilon, delta, 1.0, 0.0, 1.0*F_uppers[1], 1.0*value)
        X, p, beta = call_optimal_dist_mult_efficient(epsilon, delta, 1.0, 0.0, 1.0*F_uppers[1], 1.0*value,Xs, ranges, ps, betas)
        prior_counts[key] = max(value + optimal_noise(X, p, beta), 10^-6)
    end
    return df_counts, prior_counts
end

function noisy_means_sds_mult(epsilon, delta, sensitivities_mean, sensitivities_sd, df_means_given, df_sds_given, mins, maxs, max_sds_possible,Xs, ranges, ps, betas)
    df_means = copy.(df_means_given) #otherwise not immutable
    df_sds = copy.(df_sds_given)
    ## compute noisy mean/sds
    for (col_index, col_vec) in enumerate(df_means)
        for (label_index, label_val) in enumerate(col_vec)
            #noisy means below
            # X, p, beta = call_optimal_dist_mult(epsilon, delta, sensitivities_mean[col_index][label_index],
            #     mins[col_index], round(maxs[col_index], digits = 3), 1.0*df_means[col_index][label_index])
            X, p, beta = call_optimal_dist_mult_efficient(epsilon, delta, sensitivities_mean[col_index][label_index],
                mins[col_index], round(maxs[col_index], digits = 3), 1.0*df_means[col_index][label_index], Xs, ranges, ps, betas)

            df_means[col_index][label_index] = df_means[col_index][label_index] + optimal_noise(X, p, beta)
            #noisy sds below
            # X, p, beta = call_optimal_dist_mult(epsilon, delta, sensitivities_sd[col_index][label_index],
            #     0.0, round(max_sds_possible[col_index], digits = 3), 1.0*df_sds[col_index][label_index])
            X, p, beta = call_optimal_dist_mult_efficient(epsilon, delta, sensitivities_sd[col_index][label_index],
                0.0, round(max_sds_possible[col_index], digits = 3), 1.0*df_sds[col_index][label_index], Xs, ranges, ps, betas)
            df_sds[col_index][label_index] = max(df_sds[col_index][label_index] + optimal_noise(X, p, beta), 10^0)
        end
    end
    return df_means, df_sds
end

function noisy_classify_training_dataset(df, df_cont, labels, noisy_df_counts, noisy_prior_counts, noisy_df_means,noisy_df_sds)
    ######################## First, the "training" stage
    #noise these
    n = size(labels)[1] #number of test instances
    nr_categorical = size(df)[2] #number of cateorical variables
    nr_cont = size(df_cont)[2] #number of continuous variables
    ######################## Second, classify the points
    predictions = zeros(Int64, n) #predictions to fill
    for i = 1:n
        if nr_categorical > 0
            test_categorical = df[i, :]
        else
            test_categorical = DataFrame([])
        end
        if nr_cont > 0
            test_continuous = df_cont[i,:]
        else
            test_continuous = DataFrame([])
        end
        _, predictions[i] = test_error(test_categorical, test_continuous, n, noisy_df_counts, noisy_prior_counts, noisy_df_means, noisy_df_sds)
    end
    msc = sum(predictions .!= labels[:,1])
    mscrate = msc/n
    return predictions, msc, mscrate
end
