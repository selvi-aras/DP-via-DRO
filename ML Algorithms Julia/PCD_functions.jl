function raw_data(dataset_index)
    datasets = ["cylinder-bands", "abalone", "ecoli", "dermatology", "absent","spect","annealing","breast-cancer","spambase", "post-operative","contraceptive","bank","adult","colon-cancer", "adult-v2"] #dataset names
    d_name = datasets[dataset_index] #pick which one to take
    df = DataFrame(CSV.File("./SGDData/"*d_name*"-cooked.csv", header=false, types = Float64)) #read the specified dataset
    insertcols!(df, 1, :constant => ones(size(df)[1])) #add the intercept
    df = Matrix(df)
    labels = vec(df[:,end:end]) #last column of the categorcial dataset is labels
    df = df[:, 1:end-1] #first columns of the categorical dataset are categorical values
    df[:,2:end] = mapslices(x -> x .- mean(x) , df[:,2:end], dims=1) #except for intercept -- demean
    df[:,2:end] = mapslices(x -> x / maximum(abs.(x)), df[:,2:end], dims=1)#except for intercept
    return df, labels
end

#prox operator
function prox(val,lr; reg_const = 10^-2) #ell_1 prox with learning rate
    compare_with = reg_const * lr
    if val >= compare_with
        return val - compare_with
    elseif val <= -compare_with
        return val + compare_with
    end
    return 0
end

function PCD_train(df_train, labels_train, sim, X, p, beta; epsilon = 1.0, delta = 10^(-2), lambda=10^-2, T = 100, K = 10, method = "naive")
    Random.seed!(sim)
    n_train = length(labels_train) #number of rows in the training set
    d = size(df_train)[2] #number of variables
    #new addition -- save the whole trajectory
    xi_collection = zeros(d, T)
    #initalize the ξ first
    xi_bar = (rand(d) .- 0.5).*2 #starting ξ is in [-1,1]^n
    lr = 1.0 #learning rate, normally η_j, but for now we keep it constant because 1 / M_j and M_j is Lips.const. which is 1.
    # error_collection = []
    for t in 1:T
        xi_matrix = zeros(d, K+1) #each column will be alternatively, # xi_matrix = Array{Float64}(undef, d, T+1)
        xi_matrix[:,1] = xi_bar#initialize
        inner_products = df_train * xi_matrix[:,1] #for efficiency -- no need to re-compute the inner product of ξ^k^\top x^i for all i. Instead update gradually.
        J = randperm(d) #shuffled elements to update -- add a random seed
        #figure out which noise to sample from
        sensitivity = 2.0
        if method == "naive"
            noises = zeros(K,1)
        elseif method == "truncated"
            noises = [TLAP_noise(epsilon, delta, sensitivity) for k in 1:K]#noises to be used
        elseif method == "gaussian"
            noises = [Gaussian_noise(epsilon, delta, sensitivity) for k in 1:K]#noises to be used
        elseif method == "optimal"
            noises = [optimal_noise(X,p,beta) for k in 1:K]#noises to be used
        elseif method == "analytic"
            noises = rand(Normal(0.0, sigma_gauss), K)
        end
        #now iterate
        for k in 1:K #because all the coordinates will be updated
            xi_matrix[:, k+1] = xi_matrix[:, k] #copy the previous
            j = J[k] #take the next coordinate to flip
            # quantity_to_prox = xi_matrix[j, k] - (lr/n_train)*
            #     sum([ (1.0/(1.0+exp(labels_train[i] * (df_train[i,:]'xi_matrix[:, k]))))*(-labels_train[i]*df_train[i,j]) for i in 1:n_train])
            #append `error_collection` with the number 3
            # push!(error_collection, sum([ (1.0/(1.0+exp(labels_train[i] * (df_train[i,:]'xi_matrix[:, k]))))*(-labels_train[i]*df_train[i,j]) for i in 1:n_train])) #uncomment later for histogram
            sum_errors = sum([ (1.0/(1.0+exp(labels_train[i] * inner_products[i])))*(-labels_train[i]*df_train[i,j]) for i in 1:n_train])
            quantity_to_prox = xi_matrix[j, k] - (lr/n_train)*(noises[k] + sum_errors)
            xi_matrix[j, k+1] = prox(quantity_to_prox, lr; reg_const = lambda)
            inner_products = inner_products + (xi_matrix[j, k+1] - xi_matrix[j, k]).*df_train[:,j] #update inner products
        end
        xi_bar = (1/K) .* sum(xi_matrix[:,2:end], dims = 2) #average but skip the first column
        #make the t-th column of xi_collection the current xi_bar
        xi_collection[:,t] = xi_bar
    end
    return xi_collection #, error_collection
end

function PCD_train_return_errors(df_train, labels_train, sim, X, p, beta; epsilon = 1.0, delta = 10^(-2), lambda=10^-2, T = 100, K = 10, method = "naive")
    #same but errors are reported for the histogram supplied in GitHub
    Random.seed!(sim)
    n_train = length(labels_train) #number of rows in the training set
    d = size(df_train)[2] #number of variables
    #new addition -- save the whole trajectory
    xi_collection = zeros(d, T)
    #initalize the ξ first
    xi_bar = (rand(d) .- 0.5).*2 #starting ξ is in [-1,1]^n
    lr = 1.0 #learning rate, normally η_j, but for now we keep it constant because 1 / M_j and M_j is Lips.const. which is 1.
    error_collection = []
    for t in 1:T
        xi_matrix = zeros(d, K+1) #each column will be alternatively, # xi_matrix = Array{Float64}(undef, d, T+1)
        xi_matrix[:,1] = xi_bar#initialize
        inner_products = df_train * xi_matrix[:,1] #for efficiency -- no need to re-compute the inner product of ξ^k^\top x^i for all i. Instead update gradually.
        J = randperm(d) #shuffled elements to update -- add a random seed
        #figure out which noise to sample from
        sensitivity = 2.0
        if method == "naive"
            noises = zeros(K,1)
        elseif method == "truncated"
            noises = [TLAP_noise(epsilon, delta, sensitivity) for k in 1:K]#noises to be used
        elseif method == "gaussian"
            noises = [Gaussian_noise(epsilon, delta, sensitivity) for k in 1:K]#noises to be used
        elseif method == "optimal"
            noises = [optimal_noise(X,p,beta) for k in 1:K]#noises to be used
        elseif method == "analytic"
            noises = rand(Normal(0.0, sigma_gauss), K)
        end
        #now iterate
        for k in 1:K #because all the coordinates will be updated
            xi_matrix[:, k+1] = xi_matrix[:, k] #copy the previous
            j = J[k] #take the next coordinate to flip
            # quantity_to_prox = xi_matrix[j, k] - (lr/n_train)*
            #     sum([ (1.0/(1.0+exp(labels_train[i] * (df_train[i,:]'xi_matrix[:, k]))))*(-labels_train[i]*df_train[i,j]) for i in 1:n_train])
            #append `error_collection` with the number 3
            push!(error_collection, sum([ (1.0/(1.0+exp(labels_train[i] * (df_train[i,:]'xi_matrix[:, k]))))*(-labels_train[i]*df_train[i,j]) for i in 1:n_train])) #uncomment later for histogram
            sum_errors = sum([ (1.0/(1.0+exp(labels_train[i] * inner_products[i])))*(-labels_train[i]*df_train[i,j]) for i in 1:n_train])
            quantity_to_prox = xi_matrix[j, k] - (lr/n_train)*(noises[k] + sum_errors)
            xi_matrix[j, k+1] = prox(quantity_to_prox, lr; reg_const = lambda)
            inner_products = inner_products + (xi_matrix[j, k+1] - xi_matrix[j, k]).*df_train[:,j] #update inner products
        end
        xi_bar = (1/K) .* sum(xi_matrix[:,2:end], dims = 2) #average but skip the first column
        #make the t-th column of xi_collection the current xi_bar
        xi_collection[:,t] = xi_bar
    end
    return xi_collection, error_collection
end

function PCD_train_mult(df_train, labels_train, sim, X_mult, p_mult,beta_mult; epsilon = 1.0, delta = 10^(-2), lambda=10^-2, T = 100, K = 10, method = "mult")
    #for data dependent
    Random.seed!(sim)
    n_train = length(labels_train) #number of rows in the training set
    d = size(df_train)[2] #number of variables
    #new addition -- save the whole trajectory
    xi_collection = zeros(d, T)
    #initalize the ξ first
    xi_bar = (rand(d) .- 0.5).*2 #starting ξ is in [-1,1]^n
    lr = 1.0 #learning rate, normally η_j, but for now we keep it constant because 1 / M_j and M_j is Lips.const. which is 1.
    # error_collection = []
    for t in 1:T
        xi_matrix = zeros(d, K+1) #each column will be alternatively, # xi_matrix = Array{Float64}(undef, d, T+1)
        xi_matrix[:,1] = xi_bar#initialize
        inner_products = df_train * xi_matrix[:,1] #for efficiency -- no need to re-compute the inner product of ξ^k^\top x^i for all i. Instead update gradually.
        J = randperm(d) #shuffled elements to update -- add a random seed
        #figure out which noise to sample from
        sensitivity = 2.0
        #create a vector of vectors
        noises = Array{Array{Float64,1},1}(undef, size(p_mult)[1])
        for i in 1:size(p_mult)[1]
            #append `noises`
            noises[i] = [optimal_noise(X_mult,vec(p_mult[i,:]),beta_mult) for k in 1:K]
        end
        # noises[1] = noises[21]
        # noises[42] = noises[21]

        #now iterate
        for k in 1:K #because all the coordinates will be updated
            xi_matrix[:, k+1] = xi_matrix[:, k] #copy the previous
            j = J[k] #take the next coordinate to flip
            # quantity_to_prox = xi_matrix[j, k] - (lr/n_train)*
            #     sum([ (1.0/(1.0+exp(labels_train[i] * (df_train[i,:]'xi_matrix[:, k]))))*(-labels_train[i]*df_train[i,j]) for i in 1:n_train])
            #append `error_collection` with the number 3
            # push!(error_collection, sum([ (1.0/(1.0+exp(labels_train[i] * (df_train[i,:]'xi_matrix[:, k]))))*(-labels_train[i]*df_train[i,j]) for i in 1:n_train])) #uncomment later for histogram
            sum_errors = sum([ (1.0/(1.0+exp(labels_train[i] * inner_products[i])))*(-labels_train[i]*df_train[i,j]) for i in 1:n_train])
            # println("Sumerr: ", sum_errors)
            if method == "mult" #this will alwasy be true
                distnr = 0
                if sum_errors < -10 #as anything less than -10 was dist#1
                    # println(sum_errors)
                    distnr = 1
                elseif sum_errors > 10 #as anything greater than 10 was dist#42
                    distnr = 42 #size(p_mult)[1]
                else
                    distnr = Int(div(sum_errors + 10, 0.5) + 2)  #Int(2 + trunc((sum_errors + 10) * 2)) #this is only because the noise F_k intervals are of 0.5 length btw -10, 10 and -10< or >10 are 1 and 42
                end

                quantity_to_prox = xi_matrix[j, k] - (lr/n_train)*(noises[distnr][k] + sum_errors)
            end
            # println("distribution: ", distnr)
            xi_matrix[j, k+1] = prox(quantity_to_prox, lr; reg_const = lambda)
            inner_products = inner_products + (xi_matrix[j, k+1] - xi_matrix[j, k]).*df_train[:,j] #update inner products
        end
        xi_bar = (1/K) .* sum(xi_matrix[:,2:end], dims = 2) #average but skip the first column
        #make the t-th column of xi_collection the current xi_bar
        xi_collection[:,t] = xi_bar
    end
    return xi_collection #, error_collection
end

function mscs(xi, df_train, df_test, labels_train, labels_test)
    n_train = length(labels_train)
    predictions_train = round.((df_train*xi .>= 0.0).*2 .- 1)
    error_in_sample = sum(predictions_train .!= labels_train)/n_train

    n_test = length(labels_test) #number of rows in the training set
    predictions_test = round.((df_test*xi .>= 0.0).*2 .- 1)
    error_out_sample = sum(predictions_test .!= labels_test)/n_test
    return error_in_sample, error_out_sample
end

function generalized_mscs(xi_collection, df_train, df_test, labels_train, labels_test)
    n_train = length(labels_train)
    n_test = length(labels_test) #number of rows in the training set
    nr_cols = size(xi_collection)[2]
    errors_in_sample = zeros(nr_cols)
    errors_out_sample = zeros(nr_cols)
    #iterate over the columns of xi_collection
    for xi_index in 1:nr_cols
        xi = xi_collection[:,xi_index] #take the relevant column as your ξ
        
        predictions_train = round.((df_train*xi .>= 0.0).*2 .- 1)
        errors_in_sample[xi_index] = sum(predictions_train .!= labels_train)/n_train
        
        predictions_test = round.((df_test*xi .>= 0.0).*2 .- 1)
        errors_out_sample[xi_index] = sum(predictions_test .!= labels_test)/n_test
    end
    return errors_in_sample, errors_out_sample
end

function summarize(results_in_sample, results_out_sample)
    println("IN SAMPLE ERRORS")
    for method in keys(results_in_sample)
        println(method*": "*sprintf1("%.2f", 100*mean(results_in_sample[method]))*"%")
    end
    println("OUT SAMPLE ERRORS")
    for method in keys(results_in_sample)
        println(method*": "*sprintf1("%.2f", 100*mean(results_out_sample[method]))*"%")
    end
end

#below has all the functions for trajectory analysis -- in case you want to see what happened in each iteration
function average_simulated_errors(results_in_sample, results_out_sample)
    in_sample_average,out_sample_average = Dict([]),Dict([])
    for method in keys(results_in_sample)
        in_sample_average[method] = mean(results_in_sample[method])
        out_sample_average[method] = mean(results_out_sample[method])
    end
    return in_sample_average, out_sample_average
end

function summarize_last_iterate(in_sample_average, out_sample_average) #summarize the final results
    println("IN SAMPLE ERRORS")
    for method in keys(in_sample_average)
        println(method*": "*sprintf1("%.2f", 100*in_sample_average[method][end])*"%")
    end
    println("OUT SAMPLE ERRORS")
    for method in keys(out_sample_average)
        println(method*": "*sprintf1("%.2f", 100*out_sample_average[method][end])*"%")
    end

    #second best in sample
    second_in = in_sample_average["truncated"][end]
    our_in = in_sample_average["optimal"][end]
    base_in = in_sample_average["naive"][end]
    println("In Sample improvement: ", round(1 - (our_in - base_in)/(second_in - base_in),digits = 2)*100, "%")
    in_sample_imp = 1 - (our_in - base_in)/(second_in - base_in)
    #second best out sample
    second_out = out_sample_average["truncated"][end]
    our_out = out_sample_average["optimal"][end]
    base_out = out_sample_average["naive"][end]
    println("Out Sample improvement: ", round(1 - (our_out - base_out)/(second_out - base_out),digits = 2)*100, "%")
    out_sample_imp = 1 - (our_out - base_out)/(second_out - base_out)

    #
    println("Copy this to Latex:")
    latex = sprintf1("%.2f", 100*in_sample_average["gaussian"][end])*"\\%"*"("*sprintf1("%.2f", 100*in_sample_average["analytic"][end])*"\\%) & "* 
        sprintf1("%.2f", 100*in_sample_average["truncated"][end])*"\\% & \\textbf{"*sprintf1("%.2f", 100*in_sample_average["optimal"][end])*"\\%} {\\color{red} (TBA)} & "*
        sprintf1("%.2f", 100*in_sample_average["naive"][end])*"\\% & "*"\\textit{"*sprintf1("%.2f", 100*in_sample_imp)*"\\%} & "*
        sprintf1("%.2f", 100*out_sample_average["gaussian"][end])*"\\%"*"("*sprintf1("%.2f", 100*out_sample_average["analytic"][end])*"\\%) & "* 
        sprintf1("%.2f", 100*out_sample_average["truncated"][end])*"\\% & \\textbf{"*sprintf1("%.2f", 100*out_sample_average["optimal"][end])*"\\%} {\\color{red} (TBA)} & "*
        sprintf1("%.2f", 100*out_sample_average["naive"][end])*"\\% & "*"\\textit{"*sprintf1("%.2f", 100*out_sample_imp)*"\\%}"
    print(latex)
end