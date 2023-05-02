## Import packages and scripts
include("./NB.jl")
include("./sample.jl")
include("./train_test.jl")
include("./SGD_functions.jl")
include("./analytic_gaussian.jl")
using JLD2, LinearAlgebra
## Read a dataset and split tr:te
job_nr = 1101 #parse(Int64, ENV["PBS_ARRAY_INDEX"])
nr_sim = 100 #100 tr:te split
dataset_index = ceil(Int, job_nr / nr_sim) #for each dataset we'll have 100 tr:te splits
split_nr = mod(job_nr, nr_sim)
# dataset_index, split_nr = 1, 51 #parallelize this via HPC
datasets = ["cylinder-bands", "abalone", "ecoli", "dermatology", "absent","spect","annealing",
    "breast-cancer","spambase","post-operative","contraceptive","bank","adult","colon-cancer"]
dataset_name = datasets[dataset_index]
#
df, labels = raw_data(dataset_index) #dataset_index
df_train,df_test, labels_train, labels_test = train_test_split_simple(df, labels, split_nr) #take a tr:te split
d = size(df_test)[2]
##improve this later #optimized already. Sensitivity always 2.0. I manually optimized epsilon/d, delta/d
epsilon = 1.0
delta = 0.1
#read multi parameters
X_mult = vec(readdlm("Distributions/SGD/dep/X_mult.csv", ',', Float64)) 
beta_mult = round(X_mult[2] - X_mult[1],digits = 6) #take the beta
p_mult =  readdlm("Distributions/SGD/dep/p_mult.csv", ',', Float64) #now a matrix

nr_sim = 10^3
results_in_sample, results_out_sample = Dict([]), Dict([])
errors = []
#before this was a for loop over methods
T = 100
method = "mult"
mscs_out_sample, mscs_in_sample = [zeros(T) for _ in 1:nr_sim],[zeros(T) for _ in 1:nr_sim] #100 is from T, change otherwise
for sim in 1:nr_sim
    println("method: ", method, " simul: ", sim)
    xi_collection = PCD_train_mult(df_train, labels_train, sim, X_mult, p_mult,beta_mult; epsilon = epsilon, delta = delta, lambda=10^-8, K = Int(ceil(d/4)), T = 100, method = method)
    #append `errors` with the elements of `error_collection`
    # append!(errors, error_collection)
    # mscs_in_sample[sim], mscs_out_sample[sim] = mscs(xi_opt, df_train, df_test, labels_train, labels_test)
    mscs_in_sample[sim], mscs_out_sample[sim] = generalized_mscs(xi_collection, df_train, df_test, labels_train, labels_test)
end
results_in_sample[method], results_out_sample[method]= mscs_in_sample, mscs_out_sample

# summarize(results_in_sample, results_out_sample)
in_sample_average, out_sample_average = average_simulated_errors(results_in_sample, results_out_sample);

println("IN SAMPLE ERRORS")
println(method*": "*sprintf1("%.2f", 100*in_sample_average[method][end])*"%")

println("OUT SAMPLE ERRORS")
println(method*": "*sprintf1("%.2f", 100*out_sample_average[method][end])*"%")


# sum(abs.(errors).<=10)/length(errors)

# jldsave("./SGDResults/multi/dataset_"*string(dataset_name)*"_split_"*string(split_nr); results_in_sample, results_out_sample, epsilon, delta,
                                                                                                # nr_sim, dataset_index, split_nr)


# jldsave("./SGDResults/multi/dataset_"*string(dataset_name)*"_split_"*string(split_nr); in_sample_average, out_sample_average, epsilon, delta,
#                                                                                                 nr_sim, dataset_index, split_nr)


