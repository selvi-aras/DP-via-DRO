## Import packages and scripts
include("./NB.jl")
include("./sample.jl")
include("./train_test.jl")
include("./SGD_functions.jl")
include("./analytic_gaussian.jl")
using JLD2, LinearAlgebra
## Read a dataset and split tr:te
job_nr = 1 #parse(Int64, ENV["PBS_ARRAY_INDEX"])
nr_sim = 100 #100 tr:te split
dataset_index = ceil(Int, job_nr / nr_sim) #for each dataset we'll have 100 tr:te splits
split_nr = mod(job_nr, nr_sim)
# dataset_index, split_nr = 1, 51 #parallelize this via HPC
datasets = ["cylinder-bands", "abalone", "ecoli", "dermatology", "absent","spect","annealing",
    "breast-cancer","spambase","post-operative","contraceptive","bank","adult","colon-cancer", "adult-v2"] #dataset names #missing: breast_cancer, spambase
dataset_name = datasets[dataset_index]
#
df, labels = raw_data(dataset_index) #dataset_index
df_train,df_test, labels_train, labels_test = train_test_split_simple(df, labels, split_nr) #take a tr:te split
d = size(df_test)[2]
##improve this later #optimized already. Sensitivity always 2.0. I manually optimized epsilon/d, delta/d
epsilon = 1.0
delta = 0.1
X = vec(readdlm("Distributions/SGD/X_use.csv", ',', Float64)) 
beta = round(X[2] - X[1],digits = 6) #take the beta
p = vec(readdlm("Distributions/SGD/p_use.csv", ',', Float64))

#Analytic Gaussian parameters
sigma_gauss = calibrateAnalyticGaussianMechanism(epsilon, delta, 2.0) #we need this as analytic gaussian method will always use this sigma
## Optimize & Predict
# noise_methods = ["naive", "gaussian", "analytic", "truncated", "optimal", "mult"]
noise_methods = ["naive"]
# noise_methods = ["naive"]
nr_sim = 10^2
results_in_sample, results_out_sample = Dict([]), Dict([])
errors = []
for method in noise_methods
    mscs_out_sample, mscs_in_sample = [zeros(100) for _ in 1:nr_sim],[zeros(100) for _ in 1:nr_sim] #100 is from T, change otherwise
    for sim in 1:nr_sim
        println("method: ", method, " simul: ", sim)
        xi_collection, error_collection = PCD_train_return_errors(df_train, labels_train, sim, X, p, beta; epsilon = epsilon, delta = delta, lambda=10^-8, K = Int(ceil(d/4)), T = 100, method = method)
        #append `errors` with the elements of `error_collection`
        append!(errors, error_collection)
        # mscs_in_sample[sim], mscs_out_sample[sim] = mscs(xi_opt, df_train, df_test, labels_train, labels_test)
        mscs_in_sample[sim], mscs_out_sample[sim] = generalized_mscs(xi_collection, df_train, df_test, labels_train, labels_test)
    end
    results_in_sample[method], results_out_sample[method]= mscs_in_sample, mscs_out_sample
end

#plot a histogram of errors
using Plots
#use bins from -20 to 20 via 0.5  length intervals
b = Plots.histogram(errors, bins = range(-50, 50, length = 101), xlabel = "components of the sum of gradients",ylabel = "frequency throughout PCD iterations", yticks = 0:0.02:1,normalize=:pdf, legend = :none, color=:darkgrey, xlims = (-50,50))
#save figure as pdf
Plots.savefig(b, "./Images/histogram.pdf")

#cout how many times "errors" is between -10 and 10
sum(abs.(errors).<=10)/length(errors)

# summarize(results_in_sample, results_out_sample)
in_sample_average, out_sample_average = average_simulated_errors(results_in_sample, results_out_sample);
summarize_last_iterate(in_sample_average, out_sample_average)
# sum(abs.(errors).<=10)/length(errors)

# jldsave("./SGDResults/dataset_"*string(dataset_name)*"_split_"*string(split_nr); results_in_sample, results_out_sample, epsilon, delta,
                                                                                                # nr_sim, dataset_index, split_nr)


# jldsave("./SGDResults/dataset_"*string(dataset_name)*"_split_"*string(split_nr); in_sample_average, out_sample_average, epsilon, delta,
#                                                                                                 nr_sim, dataset_index, split_nr)


