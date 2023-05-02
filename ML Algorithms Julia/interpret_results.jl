datasets = ["post-operative", "adult", "breast-cancer", "contraceptive", "dermatology", "cylinder-bands",
                "annealing", "spect","bank", "abalone", "spambase", "ecoli", "absent"]#dataset names
dataset_index = 5 
d = datasets[dataset_index] #pick which one to take
#read
println("Dataset: "*string(d))
df = DataFrame(CSV.File("./Datasets/"*d*"-cooked.csv", header=false, types = Int64)) #read the specified dataset
labels = df[:,end:end] #last column of the categorcial dataset is labels
df = df[:, 1:end-1] #first columns of the categorical dataset are categorical values
df_cont = DataFrame(CSV.File("./Datasets/"*d*"-cooked-cont.csv", header=false, types = Float64)) #the numerical features of
println("n: ", size(df)[1])
println("V_cont: ", size(df_cont)[2])
println("V_cat: ", size(df)[2])
println("C: ", maximum(labels[:,1]))
#
nr_sim = 100 #number of tr:te splits. #10^2 # DO NOT CHANGe
#means to fill -- out-of-sample
means_truncated = zeros(nr_sim)
means_gaussian = zeros(nr_sim)
means_optimal = zeros(nr_sim)
means_naive = zeros(nr_sim)
#means to fill -- in-sample
means_truncated_in_sample = zeros(nr_sim)
means_gaussian_in_sample = zeros(nr_sim)
means_optimal_in_sample = zeros(nr_sim)
means_naive_in_sample = zeros(nr_sim)

for sim in 1:nr_sim
    obj = load("./Results/dataset_"*string(d)*"_split_"*string(mod(sim, nr_sim))) #load the file
    #read out-of-sample
    means_naive[sim] = obj["mscrate_naive"]
    means_gaussian[sim] = mean(obj["results"]["gaussian"])
    means_truncated[sim] = mean(obj["results"]["truncated"])
    means_optimal[sim] = mean(obj["results"]["optimal"])
    #read in-sample
    means_naive_in_sample[sim] = obj["mscrate_naive_in_sample"]
    means_gaussian_in_sample[sim] = mean(obj["results_in_sample"]["gaussian"])
    means_truncated_in_sample[sim] = mean(obj["results_in_sample"]["truncated"])
    means_optimal_in_sample[sim] = mean(obj["results_in_sample"]["optimal"])
end
#means of out-of-sample erros
mean_gaussian = mean(means_gaussian)
mean_truncated = mean(means_truncated)
mean_optimal = mean(means_optimal)
mean_naive = mean(means_naive)
#means of in-sample erros
mean_gaussian_in_sample = mean(means_gaussian_in_sample)
mean_truncated_in_sample = mean(means_truncated_in_sample)
mean_optimal_in_sample = mean(means_optimal_in_sample)
mean_naive_in_sample = mean(means_naive_in_sample)
#hypothesis testing -- out of sample
diff_vec =  100*(means_optimal - means_truncated) #reject that this is >= 0
mu = mean(diff_vec)
sdev = sqrt(sum((diff_vec .- mu).^2)/(100 - 1))
t_stat = (mu - (0.0))/(sdev/sqrt(100))
t_dist = TDist(100-1)
p_val_out = cdf(t_dist, t_stat)
println("In sample p_val "*string(p_val_out))
#hypothesis testing -- out of sample
diff_vec =  100*(means_optimal_in_sample - means_truncated_in_sample) #reject that this is >= 0
mu = mean(diff_vec)
sdev = sqrt(sum((diff_vec .- mu).^2)/(100 - 1))
t_stat = (mu - 0.0)/(sdev/sqrt(100))
t_dist = TDist(100-1)
p_val_in = cdf(t_dist, t_stat)
println("Out-of-sample p_val "*string(p_val_in))
#compute P()
#print
println("**********"*"Out of Sample Errors"*"**********")
println("Gaussian error: "*string(round(mean_gaussian*100, digits = 3))*"%")
println("Truncated error: "*string(round(mean_truncated*100, digits = 3))*"%")
println("Optimal error: "*string(round(mean_optimal*100, digits = 3))*"%")
println("Non-Private error: "*string(round(mean_naive*100, digits = 3))*"%")
println("**********"*"In Sample Errors"*"**********")
println("Gaussian error: "*string(round(mean_gaussian_in_sample*100, digits = 3))*"%")
println("Truncated error: "*string(round(mean_truncated_in_sample*100, digits = 3))*"%")
println("Optimal error: "*string(round(mean_optimal_in_sample*100, digits = 3))*"%")
println("Non-Private error: "*string(round(mean_naive_in_sample*100, digits = 3))*"%")
# difference closed
in_sample_diff = mean_truncated_in_sample - mean_naive_in_sample
in_sample_opt_diff = mean_optimal_in_sample - mean_naive_in_sample
in_perentage_closed = 100*((in_sample_diff - in_sample_opt_diff)/in_sample_diff)
# difference closed
out_sample_diff = mean_truncated - mean_naive
out_sample_opt_diff = mean_optimal - mean_naive
out_perentage_closed = 100*((out_sample_diff - out_sample_opt_diff)/out_sample_diff)
#for latex table
println(string(size(df)[1])*" & "*string(size(df_cont)[2])*" & "*string(size(df)[2])*" & "*string( maximum(labels[:,1]))*" & "*string(round(mean_gaussian_in_sample*100, digits = 2))*"\\% & "*
string(round(mean_truncated_in_sample*100, digits = 2))*"\\% & "*string(round(mean_optimal_in_sample*100, digits = 2))*"\\% & "*
string(round(mean_naive_in_sample*100, digits = 2))*"\\% & "*string(round(mean_gaussian*100, digits = 2))*"\\% & "*string(round(mean_truncated*100, digits = 3))*"\\% & "
*string(round(mean_optimal*100, digits = 2))*"\\% & "*string(round(mean_naive*100, digits = 2))*"\\%")
# truncated_rel_err = (mean_truncated - mean_naive)/mean_naive
# optimal_rel_err = (mean_optimal - mean_naive)/mean_naive
# println("**************************************************")
# println("Truncated relative error: "*string(round(truncated_rel_err*100, digits=4))*"%")
# println("Optimal relative error: "*string(round(optimal_rel_err*100, digi ts=4))*"%")

println(round(in_perentage_closed,digits=  2))
println(round(out_perentage_closed,digits= 2))
