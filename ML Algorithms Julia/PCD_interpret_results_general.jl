using Plots, Measures, GR
## Dataset to read
datasets = ["cylinder-bands","abalone", "ecoli", "dermatology", "absent","spect","annealing","breast-cancer","spambase", "post-operative","contraceptive","bank","adult","colon-cancer"] #dataset names
dataset_index = 11
dataset_name = datasets[dataset_index]
sim_nr = 10^2
T  = 100 #T is number of xi's in the trajectory of the training
df, labels = raw_data(dataset_index) #dataset_index
size(df)
## Obtain a new dictionary that has the average behavior of each method
in_sample_trajectory = Dict("naive" => zeros(T), "gaussian" => zeros(T), "analytic" => zeros(T), "truncated" => zeros(T), "optimal" => zeros(T))
out_sample_trajectory = Dict("naive" => zeros(T), "gaussian" => zeros(T), "analytic" => zeros(T), "truncated" => zeros(T), "optimal" => zeros(T))

for split_nr in 1:sim_nr
    if isfile("./SGDResults/dataset_"*dataset_name*"_split_"*string(mod(split_nr, sim_nr)))
        result_sim = load("./SGDResults/dataset_"*dataset_name*"_split_"*string(mod(split_nr, sim_nr)))
        in_sample, out_sample = result_sim["in_sample_average"], result_sim["out_sample_average"]
        #next take average of all the methods
        for method in keys(in_sample)
            in_sample_trajectory[method] = in_sample_trajectory[method] + (in_sample[method]./T) #I think this needs to be divided by sim_nr
            out_sample_trajectory[method] = out_sample_trajectory[method] + (out_sample[method]./T)
        end
    else
        error("File not found")
    end
end

println("The final classifier's performances comparea as the following:")
summarize_last_iterate(in_sample_trajectory, out_sample_trajectory) 

#convert everything into a df
col_names = keys(in_sample_trajectory)
df_in_sample = DataFrame(v_1 = zeros(T), v_2 = zeros(T), v_3 = zeros(T), v_4 = zeros(T), v_5 = zeros(T))
#copy df_in_sample and nam it df_out_sample
df_out_sample = deepcopy(df_in_sample)
rename!(df_in_sample, Symbol.(col_names))
for method in col_names
    df_in_sample[!, method] = in_sample_trajectory[method]
    df_out_sample[!, method] = out_sample_trajectory[method]
end

# to_show = 1:T #10:10:100

# df_to_plot = df_in_sample

# var = df_to_plot.naive
# plot_fig= Plots.plot(to_show,100*var[to_show],xlims = (to_show[1], to_show[end]), x_ticks = 10:10:100, xaxis=("PCD iteration", font(12)),xguidefontsize=50, yaxis=("in-sample error (%)", font(12)), c = 1,label = "No noise (" * sprintf1("%.2f", 100*var[end])*"%)",margin = 3mm)

# var = df_to_plot.optimal
# plot_fig= plot!(to_show,100*var[to_show],c = 2,label = "Optimized noise (" * sprintf1("%.2f", 100*var[end])*"%)",margin = 3mm)

# var = df_to_plot.truncated
# plot_fig= plot!(to_show,100*var[to_show],c = 3,label = "Truncated Laplace (" * sprintf1("%.2f", 100*var[end])*"%)",margin = 3mm)

# var = df_to_plot.analytic
# plot_fig= plot!(to_show,100*var[to_show],c = 4,label = "Analytic Gaussian (" * sprintf1("%.2f", 100*var[end])*"%)",margin = 3mm)

# var = df_to_plot.gaussian
# plot_fig= plot!(to_show,100*var[to_show],c = 5,label = "Gaussian (" * sprintf1("%.2f", 100*var[end])*"%)",margin = 3mm)

#save plot as an eps file
# Plots.savefig(plot_fig, "./Figures/PCD_"*dataset_name*"_out_sample.pdf")

