## Dataset to read
datasets = ["cylinder-bands","abalone", "ecoli", "dermatology", "absent","spect","annealing","breast-cancer","spambase", "post-operative","contraceptive","bank","adult","colon-cancer", "adult-v2"] #dataset names
dataset_index = 1
dataset_name = datasets[dataset_index]
sim_nr = 100
## Initialize list of average errors with each 100 elements
in_sample_collection = Dict("naive" => zeros(sim_nr), "gaussian" => zeros(sim_nr), "analytic" => zeros(sim_nr),
                                                            "truncated" => zeros(sim_nr), "optimal" => zeros(sim_nr))
out_sample_collection = Dict("naive" => zeros(sim_nr), "gaussian" => zeros(sim_nr), "analytic" => zeros(sim_nr),
                                                            "truncated" => zeros(sim_nr), "optimal" => zeros(sim_nr))
## Go over every tr:te split, average the 10^3 simulations, save the results
for split_nr in 1:sim_nr
    result_sim = load("./SGDResults/dataset_"*dataset_name*"_split_"*string(mod(split_nr, sim_nr)))
    in_sample, out_sample = result_sim["results_in_sample"], result_sim["results_out_sample"]
    #next take average of all the methods
    for method in keys(in_sample)
        in_sample_collection[method][split_nr] = mean(in_sample[method])
        out_sample_collection[method][split_nr] = mean(out_sample[method])
    end
end
## Average the results to have final performances
for method in keys(in_sample_collection)
    in_sample_collection[method] = [mean(in_sample_collection[method])]
    out_sample_collection[method] = [mean(out_sample_collection[method])]
end

summarize(in_sample_collection, out_sample_collection)

