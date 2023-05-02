using Plots, Measures, GR
## Dataset to read
datasets = ["cylinder-bands","abalone", "ecoli", "dermatology", "absent","spect","annealing","breast-cancer","spambase", "post-operative","contraceptive","bank","adult","colon-cancer"] #dataset names
dataset_index = 1 
dataset_name = datasets[dataset_index]
sim_nr = 10^2
T  = 100 #T is number of xi's in the trajectory of the training
df, labels = raw_data(dataset_index) #dataset_index
size(df)
## Obtain a new dictionary that has the average behavior of each method
in_sample_trajectory = Dict("mult" => zeros(T))
out_sample_trajectory = Dict("mult" => zeros(T))

for split_nr in 1:sim_nr
    if isfile("./SGDResults/multi/dataset_"*dataset_name*"_split_"*string(mod(split_nr, sim_nr)))
        result_sim = load("./SGDResults/multi/dataset_"*dataset_name*"_split_"*string(mod(split_nr, sim_nr)))
        in_sample, out_sample = result_sim["in_sample_average"], result_sim["out_sample_average"]
        #next take average of all the methods
        for method in keys(in_sample) #only "mult" :shhh:
            in_sample_trajectory[method] = in_sample_trajectory[method] + (in_sample[method]./sim_nr)
            out_sample_trajectory[method] = out_sample_trajectory[method] + (out_sample[method]./sim_nr)
        end
    else
        println("File not found")
    end
end

println("The final classifier's performances compara as the following:")

println("IN SAMPLE ERRORS")
for method in keys(in_sample_trajectory)
    println(method*": "*sprintf1("%.2f", 100*in_sample_trajectory[method][end])*"%")
end
println("OUT SAMPLE ERRORS")
for method in keys(out_sample_trajectory)
    println(method*": "*sprintf1("%.2f", 100*out_sample_trajectory[method][end])*"%")
end
