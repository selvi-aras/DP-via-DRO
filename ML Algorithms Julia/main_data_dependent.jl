## Import packages and scripts
using JLD2 #to save results
include("./NB.jl")
include("./sample.jl")
include("./sensitivities.jl")
include("./train_test.jl")
include("./read_distributions.jl")
## Figure out the job number
job_nr = parse(Int64, ENV["PBS_ARRAY_INDEX"])
nr_sim = 100
dataset_index = ceil(Int, job_nr / nr_sim)
split_nr = mod(job_nr, nr_sim)
## Specify datasets and take whichever you like.
datasets = ["post-operative", "adult", "breast-cancer", "contraceptive", "dermatology", "cylinder-bands", "credit",
                "annealing", "hepatitis", "spect", "glass", "bank", "banknote", "abalone", "mushroom", "spambase", "ecoli","page-blocks", "absent", "parkinsons", "wholesale"] #dataset names
d = datasets[dataset_index] #pick which one to take
df = DataFrame(CSV.File("./Datasets/"*d*"-cooked.csv", header=false, types = Int64)) #read the specified dataset
labels = df[:,end:end] #last column of the categorcial dataset is labels
df = df[:, 1:end-1] #first columns of the categorical dataset are categorical values
df_cont = DataFrame(CSV.File("./Datasets/"*d*"-cooked-cont.csv", header=false, types = Float64)) #the numerical features of this dataset
## Train/test-set split
#split_nr = 67
df_train,df_cont_train, labels_train, df_test, df_cont_test, labels_test = train_test_split(df,df_cont, labels, split_nr)
## STEP 0: Count the classes (categorical), and compute the mean/sds (numerical). This step will be used by all the private or non-private mechanisms.
df_counts, prior_counts = counts_of_df(df_train, labels_train) #count the relevant conditional values
#now quickly add the min/max/sensitivity of this noise for the data_dependent case -- the ones below are common for all count queries
F_lowers = vec([0.0])
F_uppers = vec([size(df)[1]])
f_overlines = vec([1.0])
#mean and sensitivities
df_means, df_sds = means_sds_of_df(df_cont_train, labels_train)
## STEP 1: Non-private Naive Bayes
predictions_naive, msc_naive, mscrate_naive = classify_training_dataset(df_test, df_cont_test, labels_test, df_counts, prior_counts, df_means, df_sds)
predictions_naive_in_sample, msc_naive_in_sample, mscrate_naive_in_sample = classify_training_dataset(df_train, df_cont_train, labels_train, df_counts, prior_counts, df_means, df_sds)
## STEP 2: Private Naive Bayes methods
epsilon = 1.0 #desired epsilon guarantee
delta = 0.1 #desirec delta guarantee
#sensitivities
mins, maxs, max_diff, sensitivities_mean, sensitivities_sd, max_sds_possible =  widths(df_cont, countmap(labels[:,1])) #compute sensitivities of numerical cols
#now append the lower and upper bounds wrt the mean
for (attribute, row) in enumerate(sensitivities_mean)
    for (lab, specific_value) in enumerate(row)
        global F_lowers = [F_lowers; mins[attribute]]
        global F_uppers = [F_uppers; maxs[attribute]]
        global f_overlines = [f_overlines; sensitivities_mean[attribute][lab]]
    end
end
#now append the lower and upper bounds wrt the st.deviation
for (attribute, row) in enumerate(sensitivities_sd)
    for (lab, specific_value) in enumerate(row)
        global F_lowers = [F_lowers; 0.0] #lower bound on SD is always 0
        global F_uppers = [F_uppers; max_sds_possible[attribute]] #upper bound on SD is what we have computed
        global f_overlines = [f_overlines; sensitivities_sd[attribute][lab]] #
    end
end
#get rid of the duplicates
matrix_together = round.(hcat(F_lowers, F_uppers, f_overlines), digits = 3) #matrix that combines
matrix_together = reduce(hcat, unique(eachrow(matrix_together)))' #reduce to unique rows
F_lowers = vec(matrix_together[:,1])
F_uppers = vec(matrix_together[:,2])
f_overlines = vec(matrix_together[:,3])
Xs, ranges, ps, betas = read_once(epsilon, delta, f_overlines, F_lowers, F_uppers)
# uncomment below to figure CSV-out parameters
st_F_lowers = join(sprintf1.("%.2f", F_lowers), ',')
st_F_uppers = join(sprintf1.("%.2f", F_uppers), ',')
st_f_overlines = join(sprintf1.("%.2f", f_overlines), ',')
open("./Datasets/mins_read_params.csv", "w") do file
    write(file, st_F_lowers)
end
open("Datasets/maxs_read_params.csv", "w") do file
    write(file, st_F_uppers)
end
open("Datasets/sens_read_params.csv", "w") do file
    write(file, st_f_overlines)
end

#loop over results
noise_methods = ["multiple_optimal"] #'multiple' stands for the multi-noise setting
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
        noisy_df_counts, noisy_prior_counts = noisy_counts_mult(epsilon, delta, df_counts, prior_counts, F_uppers,Xs, ranges, ps, betas)
        noisy_df_means, noisy_df_sds =  noisy_means_sds_mult(epsilon, delta, sensitivities_mean, sensitivities_sd,
                                                                        df_means, df_sds, mins, maxs, max_sds_possible,Xs, ranges, ps, betas)
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
jldsave("./Results/new/dataset_"*string(d)*"_split_"*string(split_nr); predictions_naive, msc_naive, mscrate_naive,
                                        predictions_naive_in_sample, msc_naive_in_sample, mscrate_naive_in_sample,
                                        epsilon, delta, nr_sim,
                                        results, results_in_sample, split_nr)

println("Multiple_optimal OOS error: "*string(round(mean(results["multiple_optimal"])*100, digits = 4))*"%")
println("Multiple_optimal IS error: "*string(round(mean(results_in_sample["multiple_optimal"])*100, digits = 4))*"%")
