using DelimitedFiles
using CSV, DataFrames
using StatsBase, Statistics

"Will retun 'df_counts' where df_counts[1][2] gives the value counts of the 1st column conditioned on the 2nd label."
function counts_of_df(df::DataFrame, labels::DataFrame)
    nr_labels = maximum(labels[:,1]) #how many possible labels there are
    #coolest single liner below
    df_counts = [[countmap(df[:, col_index][labels[:,1] .== label]) for label in 1:nr_labels] for col_index in 1:size(df)[2]]
    prior_counts = countmap(labels[:,1])
    return df_counts, prior_counts
end

"Will return df_means and df_sds where df_means[1][2] will give the mean of 1st column of df_cont when conditioned on the 2nd label."
function means_sds_of_df(df_cont::DataFrame, labels::DataFrame)
    nr_labels = maximum(labels[:,1]) #how many possible labels there are
    #coolest single liner below
    df_means = [[mean(df_cont[:, col_index][labels[:,1] .== label]) for label in 1:nr_labels] for col_index in 1:size(df_cont)[2]]
    df_sds = [[max(std(df_cont[:, col_index][labels[:,1] .== label], corrected= false), 10^0) for label in 1:nr_labels] for col_index in 1:size(df_cont)[2]]
    #we took max above as sometimes by luck we can have 0 sd as we split train/test
    return df_means, df_sds
end

"For a given single test instance, this method return the likelihoods of the test data under each class. In log scale."
function test_error(test_categorical, test_continuous, n, df_counts, prior_counts, df_means, df_sds)
    # likelihoods = vec([])
    log_likelihoods = vec([])
    for label in 1:length(prior_counts) #this length is the max label
        # val_to_push = (prior_counts[label]/n)*
        # (prod([(1/prior_counts[label])*df_counts[col_index][label][test_categorical[col_index]] for col_index in 1:length(df_counts)]))*
        # (prod([(1/(sqrt(2*pi)*df_sds[col_index][label]))*exp((test_continuous[col_index] - df_means[col_index][label]^2)/(2*df_sds[col_index][label]^2)) for col_index in 1:length(df_means)]))
        # push!(likelihoods, val_to_push)
        val_to_push = 0
        try #we are "try"ing because the conditional count may not exist hence this will give an error
            val_to_push = log(prior_counts[label]/n)
            if length(df_counts) > 0
                val_to_push += sum([log((1/prior_counts[label])*df_counts[col_index][label][test_categorical[col_index]]) for col_index in 1:length(df_counts)])
            end
            if length(df_means) > 0
                val_to_push += sum([log((1/(sqrt(2*pi)*df_sds[col_index][label]))*exp(-(test_continuous[col_index] - df_means[col_index][label])^2/(2*df_sds[col_index][label]^2))) for col_index in 1:length(df_means)])
            end
        catch e
            val_to_push = log(0) #if there is an error it means in the training set one of the conditions do not exist for a category hence -Inf.
        end
        push!(log_likelihoods, val_to_push)
    end
    return log_likelihoods, argmax(log_likelihoods) #likelihoods
end

"Inputs are: df (test categorical dataframe), df_cont (test numerical dataframe), labels (test labels), df_counts (categorical training df counts), df_means/sds (numerical training df means/sds). Return msc."
function classify_training_dataset(df::DataFrame, df_cont::DataFrame, labels::DataFrame, df_counts, prior_counts, df_means, df_sds)
    ######################## First, specify the test-set values
    n = size(labels)[1] #number of instances in the test-set
    nr_categorical = size(df)[2] #number of categorical features
    nr_cont = size(df_cont)[2] #number of numerical features
    ######################## Second, classify the points
    predictions = zeros(Int64, n) #this will be filled
    for i = 1:n #iterate over test points
        if nr_categorical > 0 #if there is at least one categorical column
            test_categorical = df[i, :] #then take the test row
        else
            test_categorical = DataFrame([]) #otherwise the test row is empty
        end
        if nr_cont > 0 #similarly, if there is at least one numerical column
            test_continuous = df_cont[i,:] #take the test row
        else
            test_continuous = DataFrame([]) #otherwise empty
        end
        _, predictions[i] = test_error(test_categorical, test_continuous, n, df_counts, prior_counts, df_means, df_sds) #send the test row to holiday (classification)
    end
    msc = sum(predictions .!= labels[:,1])
    mscrate = msc/n
    return predictions, msc, mscrate
end
