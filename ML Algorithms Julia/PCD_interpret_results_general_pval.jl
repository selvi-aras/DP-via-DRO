using Plots, Measures, GR
## Dataset to read
datasets = ["cylinder-bands","abalone", "ecoli", "dermatology", "absent","spect","annealing","breast-cancer","spambase", "post-operative","contraceptive","bank","adult","colon-cancer"] #dataset names
dataset_index = 11 #waiting for 9 (is bad),12,13, 14
dataset_name = datasets[dataset_index]

sim_nr = 10^2
T  = 100 #T is number of xi's in the trajectory of the training
df, labels = raw_data(dataset_index) #dataset_index
size(df)
## Obtain a new dictionary that has the average behavior of each method
in_sample_list = Dict("truncated" => zeros(sim_nr), "optimal" => zeros(sim_nr))
out_sample_list = Dict("truncated" => zeros(sim_nr), "optimal" => zeros(sim_nr))

for split_nr in 1:sim_nr
    if isfile("./SGDResults/dataset_"*dataset_name*"_split_"*string(mod(split_nr, sim_nr)))
        result_sim = load("./SGDResults/dataset_"*dataset_name*"_split_"*string(mod(split_nr, sim_nr)))
        in_sample, out_sample = result_sim["in_sample_average"], result_sim["out_sample_average"]
        #next take average of all the methods
        for method in ["truncated", "optimal"]
            in_sample_list[method][split_nr] = in_sample[method][end] #I think this needs to be divided by sim_nr
            out_sample_list[method][split_nr] = out_sample[method][end]
        end
    else
        println("File not found")
    end
end

#hypothesis testing -- out of sample
diff_vec =  100*(out_sample_list["optimal"] - out_sample_list["truncated"]) #reject that this is >= 0
mu = mean(diff_vec)
sdev = sqrt(sum((diff_vec .- mu).^2)/(100 - 1))
t_stat = (mu - 0.0)/(sdev/sqrt(100))
t_dist = TDist(100-1)
p_val_in = cdf(t_dist, t_stat)
println("Out-of-sample p_val "*string(p_val_in))
#in 
diff_vec =  100*(in_sample_list["optimal"] - in_sample_list["truncated"]) #reject that this is >= 0
mu = mean(diff_vec)
sdev = sqrt(sum((diff_vec .- mu).^2)/(100 - 1))
t_stat = (mu - 0.0)/(sdev/sqrt(100))
t_dist = TDist(100-1)
p_val_in = cdf(t_dist, t_stat)
println("In-sample p_val "*string(p_val_in))
