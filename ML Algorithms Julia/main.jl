## Import packages and scripts
include("./NB.jl")
include("./sample.jl")
include("./sensitivities.jl")
include("./train_test.jl")
using JLD2 #to save results
## Figure out the job number
job_nr = parse(Int64, ENV["PBS_ARRAY_INDEX"])
nr_sim = 10^3
dataset_index = ceil(Int, job_nr / nr_sim)
split_nr = mod(job_nr, nr_sim)
## Specify datasets and take whichever you like.
datasets = ["post-operative", "adult", "breast-cancer", "contraceptive", "dermatology", "cylinder-bands",
                "annealing", "spect","bank", "abalone", "spambase", "ecoli", "absent"]#dataset names
d = datasets[dataset_index] #pick which one to take
df = DataFrame(CSV.File("./Datasets/"*d*"-cooked.csv", header=false, types = Int64)) #read the specified dataset
labels = df[:,end:end] #last column of the categorcial dataset is labels
df = df[:, 1:end-1] #first columns of the categorical dataset are categorical values
df_cont = DataFrame(CSV.File("./Datasets/"*d*"-cooked-cont.csv", header=false, types = Float64)) #the numerical features of this dataset
## Train/test-set split
df_train,df_cont_train, labels_train, df_test, df_cont_test, labels_test = train_test_split(df,df_cont, labels, split_nr)
## STEP 0: Count the classes (categorical), and compute the mean/sds (numerical). This step will be used by all the private or non-private mechanisms.
df_counts, prior_counts = counts_of_df(df_train, labels_train) #count the relevant conditional values
df_means, df_sds = means_sds_of_df(df_cont_train, labels_train)
## STEP 1: Non-private Naive Bayes
predictions_naive, msc_naive, mscrate_naive = classify_training_dataset(df_test, df_cont_test, labels_test, df_counts, prior_counts, df_means, df_sds)
predictions_naive_in_sample, msc_naive_in_sample, mscrate_naive_in_sample = classify_training_dataset(df_train, df_cont_train, labels_train, df_counts, prior_counts, df_means, df_sds)
## STEP 2: Private Naive Bayes methods
epsilon = 1.0 #desired epsilon guarantee
delta = 0.1 #desirec delta guarantee
#sensitivities
mins, maxs, max_diff, sensitivities_mean, sensitivities_sd, max_sds_possible =  widths(df_cont, countmap(labels[:,1])) #compute sensitivities of numerical cols
#note: above, we actually need widths(df_cont_train, prior_counts), but this will change the sensitivities all the same so we keep it constant now.
#loop over results
noise_methods = ["gaussian", "truncated", "optimal"]
nr_sim = 10^3
results = Dict([])
results_in_sample = Dict([])
for method in noise_methods
    mscs = zeros(nr_sim)
    mscs_in_sample = zeros(nr_sim)
    for sim in 1:nr_sim
        println("method: ", method, " simul: ", sim)
        Random.seed!(sim) #give the job nr if you want
        #noisy counts
        noisy_df_counts, noisy_prior_counts = noisy_counts(epsilon, delta, df_counts, prior_counts, method)
        noisy_df_means, noisy_df_sds = noisy_means_sds(epsilon, delta, sensitivities_mean, sensitivities_sd, df_means, df_sds, method)
        #now give the noisy counts to classify
        _, _, mscrate = noisy_classify_training_dataset(df_test, df_cont_test,labels_test, noisy_df_counts, noisy_prior_counts,noisy_df_means,noisy_df_sds)
        mscs[sim] = mscrate
        #below is the in-sample
        _, _, mscrate_in_sample = noisy_classify_training_dataset(df_train, df_cont_train, labels_train, noisy_df_counts,
                                    noisy_prior_counts,noisy_df_means,noisy_df_sds)
        mscs_in_sample[sim] = mscrate_in_sample
    end
    results[method] = mscs
    results_in_sample[method] = mscs_in_sample
end
jldsave("./Results/dataset_"*string(d)*"_split_"*string(split_nr); predictions_naive, msc_naive, mscrate_naive,
                                        predictions_naive_in_sample, msc_naive_in_sample, mscrate_naive_in_sample,
                                        epsilon, delta, nr_sim,
                                        results, results_in_sample, split_nr)

println("Gaussian error: "*string(round(mean(results["gaussian"])*100, digits = 4))*"%")
println("Truncated error: "*string(round(mean(results["truncated"])*100, digits = 4))*"%")
println("Optimal error: "*string(round(mean(results["optimal"])*100, digits = 4))*"%")
println("Non-Private error: "*string(round(mscrate_naive*100, digits = 4))*"%")
